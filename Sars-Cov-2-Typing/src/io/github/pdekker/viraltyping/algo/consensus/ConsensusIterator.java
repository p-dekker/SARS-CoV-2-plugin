package io.github.pdekker.viraltyping.algo.consensus;

import java.util.Map;
import java.util.TreeMap;

import com.clcbio.api.base.math.FrequencyDistribution;
import com.clcbio.api.base.process.Activity;
import com.clcbio.api.base.process.NullActivity;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.symbol.SymbolSource;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.LocalCursor;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.MatchIntersection;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.MatchList;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.MatchListIntersections;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.MatchListSession;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.UnalignedEnds;
import com.clcbio.api.genomics.base.algo.variation.iterator.AbstractMatchListIterator;
import com.clcbio.api.genomics.base.algo.variation.iterator.MatchListIteratorParameters;
import com.clcbio.api.genomics.base.algo.variation.iterator.PairedMode;
import com.clcbio.api.genomics.base.algo.variation.iterator.PositionHandler;
import com.clcbio.api.genomics.base.algo.variation.iterator.PositionInfo;

import io.github.pdekker.viraltyping.algo.consensus.ConsensusBuilder.BuilderSession;

class ConsensusIterator extends AbstractMatchListIterator {

	private final static int MAX_UNALIGNED_END = 5;

	private static MatchListIteratorParameters getMatchListIteratorParameters(final ConsensusInterpreter par) {
		return new MatchListIteratorParameters() {
			@Override
			public boolean getSkipRepeats() {
				return par.ignoreNonSpecificMatches.get() != ConsensusInterpreter.IgnoreNonSpecificType.NONE;
			}

			@Override
			public boolean getSkipDuplicates() {
				return false;
			}

			@Override
			public int getSkipRepeatRegionMinimumLength() {
				return par.ignoreNonSpecificMatches.get() == ConsensusInterpreter.IgnoreNonSpecificType.REGION
						? par.minimumIgnoreReadLength.get()
						: -1;
			}

			@Override
			public PairedMode getPairedMode() {
				return par.ignoreBrokenPairs.get() ? PairedMode.IGNORE_BROKEN_PAIRS : PairedMode.USE_ALL;
			}

			@Override
			public boolean optimizeWrites() {
				return true;
			}
		};
	}

	private final ConsensusBuilder builder;
	private BuilderSession session;
	private final PositionFilter overlapFilter;
	private final PositionFilter qualityFilter;
	private final int mainLength;
	private final double minBreakPointRatio;

	private final Map<Integer, int[]> unalignedSequencesFor = new TreeMap<Integer, int[]>();
	private final Map<Integer, int[]> unalignedSequencesRev = new TreeMap<Integer, int[]>();

	ConsensusIterator(final ConsensusBuilder builder, final MatchList matchList, final SymbolSource symbolSource,
			final ConsensusInterpreter parms) {
		super(matchList, symbolSource, 1, false, getMatchListIteratorParameters(parms));
		this.builder = builder;
		this.matchList = matchList;
		this.mainLength = matchList.getMainSequenceLength();
		this.overlapFilter = PositionFilter.createOverlapFilter();
		if (parms.useQualityFilter.get()) {
			this.qualityFilter = PositionFilter.createQualityFilter(parms);
		} else {
			this.qualityFilter = PositionFilter.createDummyFilter(parms);
		}
		this.minBreakPointRatio = parms.inDelResolution.get() ? parms.minBreakpoint.get() : 2.0; // ratio of 2.0 is
																									// impossible so we
																									// won't call
																									// breakpoints...
	}

	int searchStart() {
		if (matchList.isCircular()) {
			return 0;
		}
		final MatchListSession session = MatchListSession.forMatchList(matchList);
		final MatchIntersection intersection = MatchListIntersections.mainRange(session, 0, 1);
		final FrequencyDistribution fd = new FrequencyDistribution();

		while (intersection.findChild(UnalignedEnds.INCLUDE)) {
			final int pos = intersection.currentChild().getMainStartPosition();
			fd.add(pos);
		}
		final Long l = fd.getMode();
		return l == null ? 0 : l.intValue();
	}

	int searchEnd() {
		if (matchList.isCircular()) {
			return matchList.getMainSequenceLength();
		}
		final FrequencyDistribution fd = new FrequencyDistribution();
		final MatchListSession session = MatchListSession.forMatchList(matchList);
		final int len = matchList.getMainSequenceLength() - 1;
		final MatchIntersection intersection = MatchListIntersections.mainRange(session, len, len + 1);
		while (intersection.findChild(UnalignedEnds.INCLUDE)) {
			final int pos = intersection.currentChild().getMainEndPosition();
			fd.add(pos);
		}
		final Long l = fd.getMode();
		return l == null ? matchList.getMainSequenceLength() : l.intValue();
	}

	int fixStart(int mainPos) {
		if (matchList.isCircular()) {
			return 0;
		}
		clear();
		final MatchListSession session = MatchListSession.forMatchList(matchList);
		final MatchIntersection intersection = MatchListIntersections.mainRange(session, mainPos, mainPos + 1);
		int lowestPos = 0;

		while (intersection.findChild(UnalignedEnds.INCLUDE)) {
			final LocalCursor cursor = intersection.currentChildLocalCursor(UnalignedEnds.INCLUDE);
			final Map<Integer, int[]> map = cursor.matchIsReversed() ? unalignedSequencesRev : unalignedSequencesFor;
			cursor.moveToMain(mainPos);
			int pos = mainPos;
			while (cursor.prev()) {
				pos--;
				byte sym = cursor.matchSymbol();
				if (!map.containsKey(pos)) {
					map.put(pos, new int[6]);
				}
				final int[] counts = map.get(pos);
				if (sym > 4) {
					sym = 4;
				}
				counts[sym + 1]++;
				if (lowestPos > pos) {
					lowestPos = pos;
				}
			}
		}
		int startPos = 0;
		for (int i = lowestPos; i < mainPos; i++) {
			int[] countsFor = unalignedSequencesFor.get(i);
			if (countsFor == null) {
				countsFor = new int[6];
			}
			int[] countsRev = unalignedSequencesRev.get(i);
			if (countsRev == null) {
				countsRev = new int[6];
			}
			if (getMaxCoverage(countsFor, countsRev) >= builder.getMinCoverageExtend()) {
				if (startPos == 0) {
					startPos = i;
				}
			} else {
				startPos = 0;
			}
		}

		if (startPos < mainPos) {
			final BuilderSession bs = builder.createSession(startPos, mainPos, NullActivity.INSTANCE);
			for (int i = startPos; i < mainPos; i++) {
				int[] countsFor = unalignedSequencesFor.get(i);
				if (countsFor == null) {
					countsFor = new int[6];
				}
				int[] countsRev = unalignedSequencesRev.get(i);
				if (countsRev == null) {
					countsRev = new int[6];
				}
				bs.add(i, countsFor, countsRev, false, false);
			}
			builder.add(bs);
			return bs.rows.size();
		}
		return 0;
	}

	private void clear() {
		unalignedSequencesFor.clear();
		unalignedSequencesRev.clear();

	}

	int fixEnd(int mainEnd) {
		if (matchList.isCircular()) {
			return 0;
		}
		clear();
		final MatchListSession session = MatchListSession.forMatchList(matchList);
		final int len = mainEnd - 1;
		final MatchIntersection intersection = MatchListIntersections.mainRange(session, len, len + 1);
		int highestPos = len;
		while (intersection.findChild(UnalignedEnds.INCLUDE)) {
			final LocalCursor cursor = intersection.currentChildLocalCursor(UnalignedEnds.INCLUDE);
			final Map<Integer, int[]> map = cursor.matchIsReversed() ? unalignedSequencesRev : unalignedSequencesFor;
			cursor.moveToMain(len);
			int pos = len;
			while (cursor.next()) {
				pos++;
				byte sym = cursor.matchSymbol();
				if (!map.containsKey(pos)) {
					map.put(pos, new int[6]);
				}
				final int[] counts = map.get(pos);
				if (sym > 4) {
					sym = 4;
				}
				counts[sym + 1]++;
				if (highestPos < pos) {
					highestPos = pos;
				}
			}
		}
		int endPos = len;
		for (int i = highestPos; i > len; i--) {
			int[] countsFor = unalignedSequencesFor.get(i);
			if (countsFor == null) {
				countsFor = new int[6];
			}
			int[] countsRev = unalignedSequencesRev.get(i);
			if (countsRev == null) {
				countsRev = new int[6];
			}
			if (getMaxCoverage(countsFor, countsRev) >= builder.getMinCoverageExtend()) {
				if (endPos == len) {
					endPos = i;
				}
			} else {
				endPos = len;
			}

		}

		if (endPos > len) {
			final BuilderSession bs = builder.createSession(mainEnd, endPos, NullActivity.INSTANCE);
			for (int i = len + 1; i < endPos; i++) {
				int[] countsFor = unalignedSequencesFor.get(i);
				if (countsFor == null) {
					countsFor = new int[6];
				}
				int[] countsRev = unalignedSequencesRev.get(i);
				if (countsRev == null) {
					countsRev = new int[6];
				}
				bs.add(i, countsFor, countsRev, false, false);
			}
			builder.add(bs);
			return bs.rows.size();
		}
		return 0;
	}

	private static int getMaxCoverage(int[] forCount, int[] revCount) {
		int max = 0;
		for (int i = 0; i < 6; i++) {
			max = Math.max(forCount[i] + revCount[i], max);
		}
		return max;
	}

	@Override
	public void iterate(final int start, final int end, final boolean skipLastInsert, final Activity activity) {
		session = builder.createSession(start, end, activity);
		// activity is not used by iterate function so the builder will
		// take care of it..
		super.iterate(start, end, skipLastInsert, activity);
		builder.add(session);
	}

	@Override
	protected PositionHandler createPositionHandler() {
		return new PositionHandler() {
			private int position;

			private int[] symbolForCounts;
			private int[] symbolRevCounts;
			private int rightBreakPoint;
			private int leftBreakPoint;

			@Override
			public void init(final int position, final int subPosition, final int maxInsert, final byte mainSymbol) {
				this.position = position;

				symbolForCounts = new int[6]; // 0 -> gap 1,2,3,4 -> symbol, 5
												// -> N
				symbolRevCounts = new int[6];
				rightBreakPoint = 0;
				leftBreakPoint = 0;
			}

			@Override
			public void handleMatch(final PositionInfo info, final int round) {
				overlapFilter.filterInfo(info);
				qualityFilter.filterInfo(info);

				if (info.getAlternativeSymbol1() == -2) {
					return;
				}

				if (info.isUnalignedEndAfterPosition()) {
					final int readPos = info.getReadPosition();
					final int readLen = info.getSegmentEnd();
					final int unalignedLength = readLen - readPos;
					if (unalignedLength > MAX_UNALIGNED_END
							&& position + unalignedLength < ConsensusIterator.this.mainLength - MAX_UNALIGNED_END) {
						rightBreakPoint++;
					}
				}

				if (info.isUnalignedEndBeforePosition()) {
					final int readPos = info.getReadPosition();
					final int unalignedLength = readPos - info.getSegmentStart();
					if (unalignedLength > MAX_UNALIGNED_END && position - unalignedLength > MAX_UNALIGNED_END) {
						leftBreakPoint++;
					}
				}

				final byte forwardSymbol = info.getAlternativeSymbol2();
				if (0 <= forwardSymbol && forwardSymbol < 4) {
					if (info.isReverse()) {
						symbolRevCounts[forwardSymbol + 1]++;
					} else {
						symbolForCounts[forwardSymbol + 1]++;
					}
				} else {
					if (info.isReverse()) {
						symbolRevCounts[5]++;
					} else {
						symbolForCounts[5]++;
					}
				}
			}

			@Override
			public void handleDelete(final PositionInfo info, final int round) {
				overlapFilter.filterInfo(info);
				qualityFilter.filterInfo(info);

				if (info.getAlternativeSymbol1() == -2) {
					return;
				}
				if (info.isReverse()) {
					symbolRevCounts[0]++;
				} else {
					symbolForCounts[0]++;
				}
			}

			@Override
			public void done(final int round) {
				final double ratioRight = (double) rightBreakPoint / getCoverage();
				final boolean right = ratioRight > minBreakPointRatio;
				final double ratioLeft = (double) leftBreakPoint / getCoverage();
				final boolean left = ratioLeft > minBreakPointRatio;
				session.add(position, symbolForCounts, symbolRevCounts, left, right);
			}

			private int getCoverage() {
				int sum = 0;
				for (final int n : symbolForCounts) {
					sum += n;
				}
				for (final int n : symbolRevCounts) {
					sum += n;
				}
				return sum;
			}

		};
	}
}
