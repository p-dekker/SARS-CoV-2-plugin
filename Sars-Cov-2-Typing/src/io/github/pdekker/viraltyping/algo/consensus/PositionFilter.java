package io.github.pdekker.viraltyping.algo.consensus;

import com.clcbio.api.genomics.base.algo.variation.iterator.PositionInfo;

public abstract class PositionFilter {
	public abstract void filterInfo(PositionInfo info);

	public static PositionFilter createOverlapFilter() {
		return new PositionFilter() {

			@Override
			public void filterInfo(final PositionInfo info) {
				if (info.getRemainingAlignmentsAtPosition() != 1) {
					// Only report one site in overlapping regions
					info.setAlternativeSymbol12((byte) -2);
				} else if (info.isConflict()) {
					// Pretend the conflict position is an N
					info.setAlternativeSymbol12((byte) 4);
				}
			}
		};
	}

	public static PositionFilter createQualityFilter(final ConsensusInterpreter p) {
		final int minCentralQuality = p.qualityMinCentral.get();
		final int minRegionQuality = p.qualityMinRegion.get();
		final int regionRadius = p.qualityRadius.get();

		return new PositionFilter() {

			@Override
			public void filterInfo(final PositionInfo info) {
				final byte[] qualities = info.getQualities();
				int readPos = info.getReadOrientedPosition();
				final int segStart = info.getSegmentStart();
				final int segEnd = info.getSegmentEnd();
				int start;
				int end;

				if (info.isDelete()) {
					if (info.isReverse()) {
						readPos++;
					}

					start = readPos - regionRadius;
					end = readPos + regionRadius;
				} else {
					if (qualities != null && qualities[readPos] < minCentralQuality) {
						info.setAlternativeSymbol12((byte) -2);
						return;
					}
					start = readPos - regionRadius;
					end = readPos + regionRadius + 1;
				}
				// shift the region so it is included in the read
				if (start < segStart) {
					end -= start - segStart;
					start = segStart;
				}
				if (end > segEnd) {
					start -= end - segEnd;
					end = segEnd;
				}
				// this happens when the read is shorter than the region
				if (start < segStart) {
					start = segStart;
				}
				if (qualities != null) {
					int totalQual = 0;
					for (int i = start; i < end; i++) {
						totalQual += qualities[i];
					}
					if (totalQual / (end - start) < minRegionQuality) {
						info.setAlternativeSymbol12((byte) -2);
					}
				}
			}

		};
	}

	public static PositionFilter createDummyFilter(final ConsensusInterpreter parms) {
		return new PositionFilter() {

			@Override
			public void filterInfo(final PositionInfo info) {
				// do nothing
			}
		};
	}
}
