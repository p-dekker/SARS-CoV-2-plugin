package io.github.pdekker.viraltyping.algo.consensus;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import com.clcbio.api.base.algorithm.ParameterValidationHandler;
import com.clcbio.api.base.algorithm.parameter.AlgoInputSelection;
import com.clcbio.api.base.algorithm.parameter.AlgoParameters;
import com.clcbio.api.base.algorithm.parameter.AlgoParametersInterpreter;
import com.clcbio.api.base.algorithm.parameter.ParameterGroup;
import com.clcbio.api.base.algorithm.parameter.keys.BooleanKey;
import com.clcbio.api.base.algorithm.parameter.keys.ByteKey;
import com.clcbio.api.base.algorithm.parameter.keys.CachingKeyChecker;
import com.clcbio.api.base.algorithm.parameter.keys.ClcObjectKey;
import com.clcbio.api.base.algorithm.parameter.keys.DoubleKey;
import com.clcbio.api.base.algorithm.parameter.keys.EnumKey;
import com.clcbio.api.base.algorithm.parameter.keys.IntegerKey;
import com.clcbio.api.base.algorithm.parameter.keys.Key;
import com.clcbio.api.base.algorithm.parameter.keys.KeyChecker;
import com.clcbio.api.base.algorithm.parameter.keys.KeyContainer;
import com.clcbio.api.base.algorithm.parameter.keys.Keys;
import com.clcbio.api.base.process.Activity;
import com.clcbio.api.base.session.ApplicationContext;
import com.clcbio.api.base.util.CreateSet;
import com.clcbio.api.base.util.Described;
import com.clcbio.api.base.util.Named;
import com.clcbio.api.free.datatypes.bioinformatics.trim.TrimAdapterList;
import com.clcbio.api.free.datatypes.framework.history.HistoryEntry;

public class ConsensusInterpreter extends AlgoParametersInterpreter {
	protected final KeyContainer keys;

	public final ParameterGroup firstPageGroup = ParameterGroup.topLevel("input-tracks", "Input tracks");

	public final ParameterGroup secondPageGroup = ParameterGroup.topLevel("read-filters", "Read filters");

	public final ParameterGroup coverageSettingsGroup = ParameterGroup.childOf(secondPageGroup, "Read Coverage");

	public final ParameterGroup qualityFilterGroup = ParameterGroup.childOf(secondPageGroup, "Quality filters");

	public final ParameterGroup extensionSettingsGroup = ParameterGroup.childOf(firstPageGroup, "Extension");

	public final ParameterGroup conflictResolutionGroup = ParameterGroup.childOf(firstPageGroup, "Conflict resolution");

	public final ParameterGroup inDelResolutionGroup = ParameterGroup.childOf(firstPageGroup, "Large InDel resolution");

	public final IntegerKey minCoverage = Keys.newIntegerKey(this, "minCoverage").defaultsTo(50).minMax(1, null)
			.labelled("Minimum coverage").withOptionKey("min-coverage").mandatory().inGroup(conflictResolutionGroup)
			.done();

	public final DoubleKey minFrequency = Keys.newDoubleKey(this, "minFrequency").defaultsTo(0.05)
			.minMax(0.0, false, 1.0, true).labelled("Minimum frequency").withOptionKey("min-frequency").mandatory()
			.inGroup(conflictResolutionGroup).done();

	public final BooleanKey addConflictAnnotations = Keys.newBooleanKey(this, "add_conflict_annotations")
			.labelled("Add conflict annotations").describedAs("Add conflict annotations")
			.withOptionKey("add-conflict-annotations").defaultsTo(true).inGroup(conflictResolutionGroup).done();

	public final EnumKey<ConflictResolution> conflictResolution = Keys
			.newEnumKey(this, "conflict_resolution", ConflictResolution.class).labelled("Conflict resolution")
			.describedAs("Conflict resolution").defaultsTo(ConflictResolution.VOTE_UNAMBIGUOUS)
			.inGroup(conflictResolutionGroup).withOptionKey("conflict-resolution").mandatory().done();

	public final BooleanKey extendStartEnd = Keys.newBooleanKey(this, "extend_sequence").labelled("Extend sequence")
			.describedAs("Extend the consensus sequence based on the unaligned end of the reads")
			.withOptionKey("extend-sequence").defaultsTo(false).inGroup(extensionSettingsGroup).done();

	public final IntegerKey minCoverageExtend = Keys.newIntegerKey(this, "minCoverageExtend").defaultsTo(20)
			.minMax(1, null).labelled("Minimum coverage extend step").withOptionKey("min-coverage-extend").mandatory()
			.inGroup(extensionSettingsGroup).done();

	public final EnumKey<PrimerMode> trimPrimers = Keys.newEnumKey(this, "trim_primers", PrimerMode.class)
			.labelled("Primers").withOptionKey("extend-sequence").defaultsTo(PrimerMode.REMOVE)
			.inGroup(extensionSettingsGroup).done();

	public final ClcObjectKey trimLinkerList = Keys.newClcObjectKey(this, "TrimLinkerLibrary")
			.withOptionKey("trim-adapter-list").labelled("Trim adapter list")
			.describedAs(
					"List of adapters to be used for trimming. Can be created in the Workbench in File -> New -> Adapter List")
			.allowedTypesDescription("Trim Adapter List").allowedTypes(TrimAdapterList.class).withMultiplicity(0, 1)
			.inGroup(extensionSettingsGroup).done();

	public final BooleanKey inDelResolution = Keys.newBooleanKey(this, "indel_resolution").labelled("InDel resolution")
			.describedAs("Perform local de Novo assembly to resolve larger InDels (experimental)")
			.withOptionKey("indel-resolution").defaultsTo(false).inGroup(inDelResolutionGroup).done();

	public final DoubleKey minBreakpoint = Keys.newDoubleKey(this, "minBreakpoint").defaultsTo(0.25)
			.minMax(0.0, false, 1.0, true).labelled("Breakpoint frequency").withOptionKey("min-frequency-breakpoint")
			.inGroup(inDelResolutionGroup).done();

	public static enum ConflictResolution implements Named, Described {
		VOTE_UNAMBIGUOUS("Vote", "Vote (A, C, G, T)", "Insert the most common nucleotide when a conflict occurs"),
		MOST_AMBIGUOUS_IF_AMBIGUOUS("Unknown", "Unknown nucleotide (N)",
				"Always insert an unknown nucleotide symbol (N) when a conflict occurs"),
		IUPAC_CONSENSUS("Ambiguous", "Ambiguity nucleotides (R, Y, etc.)",
				"Insert a specific ambiguity nucleotide symbol when a conflict occurs");

		private final String name;
		private final String shortName;
		private final String desc;

		private ConflictResolution(String shortName, String name, String desc) {
			this.name = name;
			this.desc = desc;
			this.shortName = shortName;
		}

		public String getShortName() {
			return shortName;
		}

		@Override
		public String getDescription() {
			return desc;
		}

		@Override
		public String getName() {
			return name;
		}
	}

	public static enum PrimerMode implements Named, Described {
		REMOVE("Remove primers"), IGNORE("Ignore primers");

		private String desc;

		PrimerMode(String desc) {
			this.desc = desc;
		}

		@Override
		public String getDescription() {
			return desc;
		}

		@Override
		public String getName() {
			return desc;
		}

	}

	public final BooleanKey ignoreBrokenPairs = Keys.newBooleanKey(this, "ignoreBrokenPairs").defaultsTo(true)
			.withOptionKey("ignore-broken-pairs").labelled("Ignore broken pairs")
			.describedAs("When ticked, broken pair reads will be ignored.").inGroup(coverageSettingsGroup).done();

	public final EnumKey<IgnoreNonSpecificType> ignoreNonSpecificMatches = Keys
			.newEnumKey(this, "ignoreNonSpecificMatches", IgnoreNonSpecificType.class)
			.defaultsTo(IgnoreNonSpecificType.READ).withOptionKey("ignore-nonspecific-matches")
			.labelled("Ignore non-specific matches").describedAs("").inGroup(coverageSettingsGroup).done();

	public final IntegerKey minimumIgnoreReadLength = Keys.newIntegerKey(this, "minimumIgnoreReadLength").defaultsTo(20)
			.minMax(1, null).labelled("Minimum read length").withOptionKey("min-ignore-read-length")
			.describedAs(
					"When ignoring regions with non-specific matches, this is the minimum read length needed to ignore the region.")
			.mandatory().inGroup(coverageSettingsGroup).done();

	// Quality filters
	public final BooleanKey useQualityFilter = Keys.newBooleanKey(this, "useQualityFilter").defaultsTo(false)
			.withOptionKey("use-quality-filter").labelled("Base quality filter")
			.describedAs(
					"When ticked, reads with bases that do not meet the quality requirements specified below will be ignored.")
			.inGroup(qualityFilterGroup).done();

	public final IntegerKey qualityRadius = Keys.newIntegerKey(this, "qualityRadius").defaultsTo(5).minMax(1, null)
			.labelled("Neighborhood radius").withOptionKey("neighborhood-radius")
			.describedAs(
					"Region size for quality filter defined as the distance to both sides of the central nucleotide")
			.inGroup(qualityFilterGroup).done();

	public final ByteKey qualityMinCentral = Keys.newByteKey(this, "qualityMinCentral").defaultsTo((byte) 20)
			.minMax((byte) 0, (byte) 64).labelled("Minimum central quality").withOptionKey("min-central-quality")
			.describedAs("Minimum quality score for central nucleotide").inGroup(qualityFilterGroup).done();

	public final ByteKey qualityMinRegion = Keys.newByteKey(this, "qualityMinRegion").defaultsTo((byte) 15)
			.minMax((byte) 0, (byte) 64).labelled("Minimum neighborhood quality")
			.withOptionKey("min-neighborhood-quality")
			.describedAs("Minimum average quality score for the region around the central nucleotide")
			.inGroup(qualityFilterGroup).done();

	public final BooleanKey createReport = Keys.newBooleanKey(this, "createReport").defaultsTo(true)
			.withOptionKey("create-report").labelled("Create report").describedAs("Create a report")
			.inGroup(ParameterGroup.OUTPUT_OPTIONS).done();

	@Override
	protected void validateKeys(final ParameterValidationHandler validatorHandler,
			final ApplicationContext applicationContext) {
		final Set<Key<?>> validateKeys = CreateSet.from(getKeyObjects());
		for (final Key<?> key : validateKeys) {
			key.validateCurrent(validatorHandler, applicationContext);
		}
	}

	@Override
	public KeyChecker createKeyChecker(final ApplicationContext applicationContext) {
		return new CachingKeyChecker(getAlgoParameters()) {
			@Override
			public boolean isEnabled(final Key<?> key, final AlgoInputSelection input, final Activity activity)
					throws InterruptedException {
				if (key == minimumIgnoreReadLength) {
					return ignoreNonSpecificMatches.get() == IgnoreNonSpecificType.REGION;
				}
				if (key == qualityMinCentral || key == qualityMinRegion || key == qualityRadius) {
					return useQualityFilter.get();
				}
				if (key == trimLinkerList) {
					return !trimPrimers.get().equals(PrimerMode.IGNORE);
				}
				if (key == minCoverageExtend) {
					return extendStartEnd.get();
				}

				if (key == minBreakpoint) {
					return inDelResolution.get();
				}
				return true;
			}

			@Override
			public void validate(final Collection<? extends Key<?>> keys, final AlgoInputSelection input,
					final ParameterValidationHandler handler, final Activity activity) throws InterruptedException {

				ConsensusInterpreter.this.validate(handler, applicationContext, activity);

				if (trimPrimers.get() == PrimerMode.REMOVE && trimLinkerList.getClcObject(applicationContext) == null) {
					handler.postMissingParameter(trimLinkerList);
				}
			}
		};
	}

	public ConsensusInterpreter(final AlgoParameters parameters) {
		super(parameters);
		keys = new KeyContainer(minCoverage, addConflictAnnotations, conflictResolution, minFrequency, extendStartEnd,
				minCoverageExtend, trimPrimers, trimLinkerList, inDelResolution, minBreakpoint, ignoreBrokenPairs,
				ignoreNonSpecificMatches, minimumIgnoreReadLength, useQualityFilter, qualityMinCentral,
				qualityMinRegion, qualityRadius, createReport);
	}

	@Override
	public void setToDefault() {
		keys.setToDefault();
	}

	@Override
	public String getClassKey() {
		return "sars_cov2_censensus_algo_interpreter";
	}

	@Override
	public Collection<Key<?>> getKeyObjects() {
		return keys.getKeys();
	}

	@Override
	public Set<String> getKeys() {
		return keys.getKeySet();
	}

	public enum IgnoreNonSpecificType {
		READ("Reads"), REGION("Regions"), NONE("No");

		private String name;

		private IgnoreNonSpecificType(final String name) {
			this.name = name;
		}

		@Override
		public String toString() {
			return name;
		}
	}

	public HistoryEntry enrichHistoryEntry(final HistoryEntry entry) {
		final Set<Key<?>> ignored = getIgnoredHistoryKeys();
		for (final Key<?> key : getKeyObjects()) {
			if (ignored.contains(key)) {
				continue;
			}
			entry.addParameterEntry(key.getShortDescription(), key.getUserString());
		}
		return entry;
	}

	public Set<Key<?>> getIgnoredHistoryKeys() {
		final Set<Key<?>> ignored = new HashSet<Key<?>>();
		if (ignoreNonSpecificMatches.get() == IgnoreNonSpecificType.NONE) {
			ignored.add(minimumIgnoreReadLength);
		}
		if (!useQualityFilter.get()) {
			ignored.add(qualityMinCentral);
			ignored.add(qualityMinRegion);
			ignored.add(qualityRadius);
		}

		if (!inDelResolution.get()) {
			ignored.add(minBreakpoint);
		}
		return ignored;
	}
}
