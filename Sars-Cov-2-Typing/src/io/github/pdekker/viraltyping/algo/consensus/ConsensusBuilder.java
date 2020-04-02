package io.github.pdekker.viraltyping.algo.consensus;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.stream.Collectors;

import com.clcbio.api.base.algorithm.AlgoOutputNamingTools;
import com.clcbio.api.base.math.FrequencyDistribution;
import com.clcbio.api.base.misc.Cleanupable;
import com.clcbio.api.base.process.Activity;
import com.clcbio.api.base.session.FactoryManager;
import com.clcbio.api.base.util.Comprehension;
import com.clcbio.api.base.util.CreateList;
import com.clcbio.api.base.util.CreateMap;
import com.clcbio.api.base.util.IteratorTools;
import com.clcbio.api.base.util.StringTools;
import com.clcbio.api.clc.datatypes.bioinformatics.variant.LinkageGroupImpl;
import com.clcbio.api.clc.datatypes.bioinformatics.variant.SequenceAlterationFactory;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.BasicSequence;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.Sequence;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.SequenceBuilder;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.alphabet.Alphabet;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.alphabet.AlphabetTools;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.feature.Feature;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.feature.FeatureTypes;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.interval.Interval;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.interval.SimpleInterval;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.region.Region;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.region.RegionTools;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.symbol.SymbolSource;
import com.clcbio.api.free.datatypes.bioinformatics.variant.MappingVariant;
import com.clcbio.api.free.datatypes.bioinformatics.variant.SequenceAlteration;
import com.clcbio.api.free.datatypes.bioinformatics.variant.Variant;

import io.github.pdekker.viraltyping.algo.consensus.ConsensusInterpreter.ConflictResolution;

/**
 * @author pdekker
 * 
 */
public class ConsensusBuilder implements Cleanupable {
	private final static int BUFFER_SIZE = 10_000;
	public final static String LOW_COVERAGE = "Low Coverage";
	private final static String COVERAGE = "Coverage";
	private final static String VAR = "var";

	public final static String FAILURES = "Failures";
	public final static String UNSURE = "Unsure";

	private final String name;
	private final int minCoverage;
	private final int minCoverageExtend;
	private final boolean addConflicts;
	private final double minFrequency;
	private CoverageInformation coverageInformation;
//	private final double minProblematicVariantFreq;

	private final int IGNORE_FAILURES_CLOSE_TO_END = 25;

	private final ConflictResolution conflictResolution;

	private static final byte GAP = (byte) -1;
	private static final Alphabet DNA = AlphabetTools.getDnaAlphabet();
	private static final byte N = DNA.getMostAmbiguousSymbolIndex();

	public static ConsensusBuilder createBuilder(SymbolSource main, boolean extend, int minCoverage,
			int minCoverageExtend, double minFrequency, boolean addConflicts, ConflictResolution conflictResolution,
			double minProblematicVariantFreq) {
		return new ConsensusBuilder(main.getName(), minCoverage, minCoverageExtend, minFrequency, addConflicts,
				conflictResolution, minProblematicVariantFreq);
	}

	List<BuilderSession> sessions;

	private ConsensusBuilder(String name, int minCoverage, int minCoverageExtend, double minFrequency,
			boolean addConflicts, ConflictResolution conflictResolution, double minProblematicVariantFreq) {
		sessions = CreateList.of();
		this.name = name;
		this.minCoverage = minCoverage;
		this.minCoverageExtend = minCoverageExtend;
		this.minFrequency = minFrequency;
		this.addConflicts = addConflicts;
		this.conflictResolution = conflictResolution;
//		this.minProblematicVariantFreq = minProblematicVariantFreq;
	}

	int getMinCoverageExtend() {
		return minCoverageExtend;
	}

	public Iterator<DataPoint> iterator() {
		if (sessions.size() == 0) {
			return IteratorTools.getEmptyIterator();
		}
		if (sessions.size() == 1) {
			return sessions.get(0).iterator();
		}
		final List<Iterator<DataPoint>> iterators = CreateList.of();
		for (final BuilderSession session : sessions) {
			iterators.add(session.iterator());
		}
		return Comprehension.flattenIteratorOfIterators(iterators.iterator());
	}

	public int getNumberOfDataPoints() {
		return sessions.stream().mapToInt(s -> s.size()).sum();
	}

	public Iterator<DataPoint> gapFilteredIterator() {
		return IteratorTools.getFilteredIterator(iterator(), dp -> {
			final List<PreVariant> prevariant = dp.getAlternativeSymbols(minFrequency, minCoverage);
			return !DataPoint.willResultInGap(prevariant, conflictResolution);
		});
	}

	public static Region getBreakPointRegion(BasicSequence bs) {
		final Optional<Feature> of = bs.getFeatures().stream()
				.filter(f -> f.getType().equals(FAILURES) && f.getName().equals(UNSURE)).findFirst();
		return of.isPresent() ? of.get().getRegion() : null;
	}

	public Sequence getConsensus() {
		final Iterator<DataPoint> it = iterator();
		final SequenceBuilder seqBuilder = FactoryManager.getInstance().getSequenceFactory().createBuilder();
		seqBuilder.setAlphabet(AlphabetTools.getDnaAlphabet());
		seqBuilder.setName(AlgoOutputNamingTools.createRetaggedName(name, "cons"));
		final byte[] tmpBuffer = new byte[BUFFER_SIZE];

		int i = 0;
		int pos = 0;
		int lastPos = Integer.MIN_VALUE;
		int count = 0;
		final int maxCount = getNumberOfDataPoints() - IGNORE_FAILURES_CLOSE_TO_END;

		final Map<Region, List<List<PreVariant>>> insertions = CreateMap.of();

		final List<Interval> breakpoints = new ArrayList<>();

		while (it.hasNext()) {
			final DataPoint dp = it.next();
			count++;
			final List<PreVariant> prevariant = dp.getAlternativeSymbols(minFrequency, minCoverage);
			if (DataPoint.willResultInGap(prevariant, conflictResolution)) {
				if (prevariant.size() > 1) {
					final Region variantRegion = new Region(pos, pos);
					if (addConflicts) {
						if (!insertions.containsKey(variantRegion)) {
							insertions.put(variantRegion, CreateList.<List<PreVariant>>of());
						}
						insertions.get(variantRegion).add(prevariant);
					}
				}
				continue;
			}
			if (dp.position >= lastPos) {
				lastPos = dp.position;
			} else {
				throw new IllegalStateException("Data points not sorted!");
			}
			final byte consensus = DataPoint.getConsensusSymbol(prevariant, conflictResolution);
			tmpBuffer[i] = consensus;
			if (prevariant.size() > 1 && addConflicts) {
				final Region variantRegion = new Region(pos, pos + 1);
				addVariants(seqBuilder, variantRegion, prevariant, consensus);
			}
			if (count > IGNORE_FAILURES_CLOSE_TO_END && count < maxCount) {
				if (dp.isBreakpointPosition()) {
					breakpoints.add(new SimpleInterval(pos, pos + 1));
				}
			}
			i++;
			pos++;
			if (i == BUFFER_SIZE) {
				seqBuilder.addSequenceData(tmpBuffer, 0, BUFFER_SIZE);
				i = 0;
			}
		}
		if (i > 0) {
			seqBuilder.addSequenceData(tmpBuffer, 0, i);
		}

		// finish insertions...
		for (final Entry<Region, List<List<PreVariant>>> e : insertions.entrySet()) {
			addInsertionVariants(seqBuilder, e.getKey(), e.getValue());
		}

		final List<Feature> lowCoverageRegions = getCoverageInformation().lowCoverageRegions;
		for (final Feature f : lowCoverageRegions) {
			seqBuilder.addFeature(f);
		}
		if (breakpoints.size() > 1) {
			final Region r = RegionTools.getOrderedRegionWithoutOverlaps(new Region(breakpoints));
			seqBuilder.addFeature(new Feature(UNSURE, r, FAILURES));
		}
		return seqBuilder.finish();
	}

	private void addVariants(SequenceBuilder seqBuilder, Region variantRegion, List<PreVariant> prevariant,
			byte consensus) {
		seqBuilder.addFeature(asFeature(prevariant, variantRegion, consensus));
	}

	private void addInsertionVariants(SequenceBuilder seqBuilder, Region variantRegion,
			List<List<PreVariant>> prevariant) {
		if (prevariant.size() == 1) {
			seqBuilder.addFeature(asFeature(prevariant.get(0), variantRegion, GAP));
			return;
		}

		final Feature f = new Feature(FeatureTypes.CONFLICT, variantRegion, FeatureTypes.CONFLICT);
		f.addAnnotation("Conflict resolution",
				"" + conflictResolution.getShortName() + " called '" + getBaseName(GAP) + "'");
		f.addAnnotation("Reference position", variantRegion.toUserString());
		f.addAnnotation("Start position", variantRegion.getFirstPos().getMin() + 1);
		f.addAnnotation("Ref", getBaseName(GAP));

		final StringBuilder variants = new StringBuilder();

		int maxVariants = -1;
		int coverage = 0;
		for (final List<PreVariant> pv : prevariant) {
			Collections.sort(pv, (o1, o2) -> Double.compare(o2.freq, o1.freq));
			maxVariants = Math.max(maxVariants, pv.size());
			coverage += pv.get(0).coverage;
		}

		coverage /= prevariant.size();

		final StringBuilder sb = new StringBuilder();
		for (int i = 0; i < maxVariants; i++) {
			sb.setLength(0);
			int count = 0;
			int gaps = 0;
			int gapCount = 0;

			for (final List<PreVariant> pvs : prevariant) {
				if (i < pvs.size()) {
					final PreVariant pv = pvs.get(i);
					if (pv.isGap()) {
						gaps++;
						gapCount += pv.count;
					} else {
						sb.append(getBaseName(pv.symbol));
						count += pv.count;
					}
				}
			}

			if (sb.length() == 0) {
				// gap..
				if (gaps > 0) {
					gapCount /= gaps;

					final double freq = (double) gapCount / coverage;

					f.addAnnotation("gap count", count);
					f.addAnnotation("gap perc", freq * 100);
				}
			} else {
				final String insertion = sb.toString();
				if (variants.length() > 0) {
					variants.append(',');
				}
				variants.append(insertion);

				count /= insertion.length();
				coverage /= insertion.length();
				final double freq = (double) count / coverage;

				f.addAnnotation(insertion + " count", count);
				f.addAnnotation(insertion + " perc", freq * 100);
			}
		}
		f.addAnnotation(COVERAGE, coverage);
		f.addAnnotation(VAR, variants.toString());
		seqBuilder.addFeature(f);
	}

	public CoverageInformation getCoverageInformation() {
		if (coverageInformation != null) {
			return coverageInformation;
		}
		final Iterator<DataPoint> it = gapFilteredIterator();
		final FrequencyDistribution fd = new FrequencyDistribution();

		int pos = 0;

		final List<Feature> lowCoverageRegions = CreateList.of();
		int startLow = -1;
		int coverage = 0;

		while (it.hasNext()) {
			final DataPoint dp = it.next();
			final int cov = dp.getCoverage();
			final boolean lowCoverage = cov < minCoverage;
			fd.add(cov);
			if (lowCoverage) {
				if (startLow == -1) {
					startLow = pos;
					coverage = 0;
				}
				coverage += cov;
			} else {
				if (startLow != -1) {
					final Feature f = new Feature(LOW_COVERAGE, new Region(startLow, pos), FAILURES);
					coverage /= pos - startLow;
					f.addAnnotation("Coverage", coverage);
					lowCoverageRegions.add(f);
					startLow = -1;
				}
			}
			pos++;
		}
		if (startLow != -1) {
			final Feature f = new Feature(LOW_COVERAGE, new Region(startLow, pos), FAILURES);
			coverage /= pos - startLow;
			f.addAnnotation("Coverage", coverage);
			lowCoverageRegions.add(f);
		}
		coverageInformation = new CoverageInformation(lowCoverageRegions, fd);
		return coverageInformation;
	}

	public void finish() {
	}

	@Override
	public void cleanup() {
	}

	public BuilderSession createSession(int start, int end, Activity activity) {
		return new BuilderSession(start, end, activity);
	}

	public void add(BuilderSession session) {
		sessions.add(session);
	}

	final static class BuilderSession {
		final List<DataPoint> rows;
		int lastPositionFoundByDone = 0;
		final int start;
		final int end;
		final Activity activity;

		private BuilderSession(int start, int end, Activity activity) {
			rows = new ArrayList<DataPoint>((int) ((end - start) * 1.2));
			this.activity = activity;
			this.start = start;
			this.end = end;
		}

		public Iterator<DataPoint> iterator() {
			return rows.iterator();
		}

		public void add(int pos, int[] symbolForCounts, int[] symbolRevCounts, boolean left, boolean right) {
			rows.add(new DataPoint(pos, symbolForCounts, symbolRevCounts, left, right));
			if (pos > lastPositionFoundByDone) {
				lastPositionFoundByDone = pos + 100; // we add 100 so we don't
				activity.setProgress((double) (pos - start) / end - start);
			}
		}

		public int size() {
			return rows.size();
		}
	}

	private static class PreVariant {
		private final byte symbol;
		private final int count;
		private final double freq;
		private final int coverage;

		PreVariant(int symbol, int count, int coverage) {
			if (symbol == 5) {
				this.symbol = N;
			} else {
				this.symbol = (byte) (symbol - 1);
			}
			this.count = count;
			this.coverage = coverage;
			this.freq = (double) count / coverage;
		}

		public boolean isGap() {
			return symbol == -1;
		}
	}

	private static class DataPoint {

		private final int position;
		private final int[] symbolForCounts;
		private final int[] symbolRevCounts;
		private final boolean left;
		private final boolean right;

		public DataPoint(int position, int[] symbolForCounts, int[] symbolRevCounts, boolean left, boolean right) {
			this.position = position;
			this.symbolForCounts = symbolForCounts;
			this.symbolRevCounts = symbolRevCounts;
			this.left = left;
			this.right = right;
		}

		public static byte getConsensusSymbol(List<PreVariant> prevariant, ConflictResolution conflictResolution) {
			if (prevariant.size() == 0) {
				return N; // too low coverage
			}
			if (prevariant.size() == 1) {
				return prevariant.get(0).symbol;
			}
			switch (conflictResolution) {
			case IUPAC_CONSENSUS:
				final List<Byte> list = CreateList.of();
				for (final PreVariant pv : prevariant) {
					if (pv.symbol >= 0) { // ignore gaps
						list.add(pv.symbol);
					}
				}
				return DNA.getConsensus(list);
			case MOST_AMBIGUOUS_IF_AMBIGUOUS:
				return N;
			case VOTE_UNAMBIGUOUS:
				PreVariant best = null;
				for (final PreVariant pv : prevariant) {
					if (best == null) {
						best = pv;
					} else {
						if (best.count < pv.count) {
							best = pv;
						}
					}
				}
				return best.symbol;
			default:
				throw new AssertionError("Enum missing");
			}
		}

		public static boolean willResultInGap(List<PreVariant> prevariant, ConflictResolution conflictResolution) {

			if (prevariant.size() == 0) {
				return false; // too low coverage
			}
			if (prevariant.size() == 1) {
				return prevariant.get(0).isGap();
			}
			switch (conflictResolution) {
			case IUPAC_CONSENSUS:
			case MOST_AMBIGUOUS_IF_AMBIGUOUS:
				return false;
			case VOTE_UNAMBIGUOUS:
				PreVariant best = null;
				for (final PreVariant pv : prevariant) {
					if (best == null) {
						best = pv;
					} else {
						if (best.count < pv.count) {
							best = pv;
						}
					}
				}
				return best.isGap();
			default:
				throw new AssertionError("Enum missing");
			}
		}

		public List<PreVariant> getAlternativeSymbols(double freq, int minCoverage) {
			if (isLowCoverage(minCoverage)) {
				return Collections.emptyList();
			}
			int sum = 0;
			// int max = 0;
			for (int i = 0; i < 6; i++) {
				final int c = symbolForCounts[i] + symbolRevCounts[i];
				sum += c;
				// max = Math.max(c, max);
			}
			final List<PreVariant> variants = CreateList.of();
			final int minCount = (int) (sum * freq);
			for (int i = 0; i < 5; i++) {
				// note i == 5 is count for N but we don't want to report N.
				final int c = symbolForCounts[i] + symbolRevCounts[i];
				if (c > minCount) {
					variants.add(new PreVariant(i, c, sum));
				}
			}
			return variants;
		}

		private int getCoverage() {
			int sum = 0;
			for (int i = 0; i < 6; i++) {
				sum += symbolForCounts[i] + symbolRevCounts[i];
			}
			return sum;
		}

		public boolean isLowCoverage(int minCoverage) {
			return getCoverage() < minCoverage;
		}

		public boolean isBreakpointPosition() {
			return left || right;
		}

	}

	private Feature asFeature(List<PreVariant> pv, Region variantRegion, byte consensus) {
		final Feature f = new Feature(FeatureTypes.CONFLICT, variantRegion, FeatureTypes.CONFLICT);
		f.addAnnotation("Conflict resolution",
				"" + conflictResolution.getShortName() + " called '" + getBaseName(consensus) + "'");

		final StringBuilder sb = new StringBuilder();

		for (final PreVariant p : pv) {
			final String symb = getBaseName(p.symbol);
			f.addAnnotation(symb + " count", p.count);
			f.addAnnotation(symb + " perc", p.freq * 100);
			if (sb.length() > 0) {
				sb.append(',');
			}
			sb.append(symb);
		}
		if (pv.size() > 0) {
			final PreVariant p = pv.get(0);
			f.addAnnotation("Coverage", p.coverage);
			f.addAnnotation("Reference position", variantRegion.toUserString());
			f.addAnnotation("Start position", variantRegion.getFirstPos().getMin() + 1);

			f.addAnnotation("Ref", getBaseName(consensus));
		}
		f.addAnnotation("Var", sb.toString());
		return f;
	}

	private static String getBaseName(final byte b) {
		return DNA.getSymbol(b).getShortName();
	}

	public static List<Variant> asVariant(Feature f) {
		if (!f.getType().equals(FeatureTypes.CONFLICT)) {
			throw new IllegalArgumentException("Feature type should be 'conflict' but got " + f.getType());
		}
		final List<Variant> variants = CreateList.of();

		try {
			final int startPos = (Integer) f.getFirstAnnotation("Start position") - 1;
			final String ref = (String) f.getFirstAnnotation("Conflict resolution");
			final int from = ref.indexOf('\'');
			final int to = ref.lastIndexOf('\'');
			final String refBase = ref.substring(from + 1, to);

			final int coverage = (Integer) f.getFirstAnnotation("Coverage");

			final byte[] refs = refBase.equals("gap") ? null : AlphabetTools.convertToBytes(refBase, DNA);

			final String[] vars = StringTools.split((String) f.getFirstAnnotation("Var"), ',');
			for (final String var : vars) {
				if (f.getFirstAnnotation(var + " count") != null) {
					final int count = (Integer) f.getFirstAnnotation(var + " count");
					final double freq = (Double) f.getFirstAnnotation(var + " perc");

					final byte[] varb = var.equals("gap") ? null : AlphabetTools.convertToBytes(var, DNA);

					final SequenceAlteration sa = SequenceAlterationFactory.createHeuristicSequenceAlteration(startPos,
							varb, refs);
					final MappingVariant mv = new MappingVariant(sa, LinkageGroupImpl.createLinkageGroupNone());
					mv.setCoverage(coverage);
					mv.setCount(count);
					mv.setFrequency(freq);
					variants.add(mv);
				}
			}
		} catch (NullPointerException | NumberFormatException | ClassCastException e) {
			// ignore, we had invalid features, so we return empty list of
			// variants
		}
		return variants;
	}

	public static boolean hasVariant(BasicSequence seq) {
		return seq.getFeatureCount(FeatureTypes.CONFLICT) > 0;
	}

	public static boolean hasFailures(BasicSequence seq) {
		return seq.getFeatureCount(FAILURES) > 0;
	}

	public static class CoverageInformation {
		private final static int MAX_ENUM_REGIONS = 10;

		final List<Feature> lowCoverageRegions;
		final FrequencyDistribution fd;

		public CoverageInformation(List<Feature> lowCoverage, FrequencyDistribution fd) {
			this.lowCoverageRegions = lowCoverage;
			this.fd = fd;
		}

		public String lowCoverageRegions() {
			if (lowCoverageRegions.isEmpty()) {
				return "-";
			}
			final String s = lowCoverageRegions.stream().limit(MAX_ENUM_REGIONS).map(f -> f.getRegion().toUserString())
					.collect(Collectors.joining(","));
			if (lowCoverageRegions.size() > MAX_ENUM_REGIONS) {
				return s + "...";
			}
			return s;
		}
	}
}
