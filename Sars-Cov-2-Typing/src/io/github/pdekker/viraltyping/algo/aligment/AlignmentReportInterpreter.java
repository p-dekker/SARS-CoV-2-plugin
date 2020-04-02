package io.github.pdekker.viraltyping.algo.aligment;

import java.util.Collection;
import java.util.Set;

import com.clcbio.api.base.algorithm.ParameterValidationHandler;
import com.clcbio.api.base.algorithm.parameter.AlgoInputSelection;
import com.clcbio.api.base.algorithm.parameter.AlgoParameters;
import com.clcbio.api.base.algorithm.parameter.AlgoParametersInterpreter;
import com.clcbio.api.base.algorithm.parameter.ParameterGroup;
import com.clcbio.api.base.algorithm.parameter.keys.BooleanKey;
import com.clcbio.api.base.algorithm.parameter.keys.CachingKeyChecker;
import com.clcbio.api.base.algorithm.parameter.keys.EnumKey;
import com.clcbio.api.base.algorithm.parameter.keys.Key;
import com.clcbio.api.base.algorithm.parameter.keys.KeyChecker;
import com.clcbio.api.base.algorithm.parameter.keys.KeyContainer;
import com.clcbio.api.base.algorithm.parameter.keys.Keys;
import com.clcbio.api.base.algorithm.parameter.keys.StringKey;
import com.clcbio.api.base.process.Activity;
import com.clcbio.api.base.session.ApplicationContext;
import com.clcbio.api.base.util.Described;
import com.clcbio.api.base.util.Named;

public class AlignmentReportInterpreter extends AlgoParametersInterpreter {
	protected final KeyContainer keys;

	public final ParameterGroup firstPageGroup = ParameterGroup.topLevel("Setting", "Settings");

	public final ParameterGroup referenceGroup = ParameterGroup.childOf(firstPageGroup, "Reference");
	public final ParameterGroup optionsGroup = ParameterGroup.childOf(firstPageGroup, "Options");
	public final ParameterGroup translateGroup = ParameterGroup.childOf(firstPageGroup, "Translate to proteins");

	public final StringKey sequenceName = Keys.newStringKey(this, "sequence-name").labelled("Reference")
			.describedAs("Sequence used as reference").inGroup(referenceGroup).defaultsTo("").optional().done();

	public final EnumKey<ReferenceDetermination> referenceSelection = Keys
			.newEnumKey(this, "reference-selection", ReferenceDetermination.class).labelled("Reference Selection")
			.describedAs("Method to Select a reference sequence").defaultsTo(ReferenceDetermination.FIRST)
			.inGroup(referenceGroup).mandatory().done();

	public final BooleanKey ignoreGapsAtEnd = Keys.newBooleanKey(this, "ïgnore-gaps-at-end")
			.labelled("Ignore gaps at begin/end").inGroup(optionsGroup).defaultsTo(true).done();

	public AlignmentReportInterpreter(final AlgoParameters parameters) {
		super(parameters);
		keys = new KeyContainer(referenceSelection, sequenceName, ignoreGapsAtEnd);
	}

	@Override
	public void setToDefault() {
		keys.setToDefault();
	}

	@Override
	public String getClassKey() {
		return "sars_cov2_Alignment_report_interpreter";
	}

	@Override
	public Collection<Key<?>> getKeyObjects() {
		return keys.getKeys();
	}

	@Override
	public Set<String> getKeys() {
		return keys.getKeySet();
	}

	public static enum ReferenceDetermination implements Described, Named {
		AUTOMATIC("Automatic"), FIRST("First sequence"), SELECT("Select sequence");

		private String desc;

		ReferenceDetermination(String desc) {
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

	@Override
	public KeyChecker createKeyChecker(final ApplicationContext applicationContext) {
		return new CachingKeyChecker(getAlgoParameters()) {
			@Override
			public boolean isEnabled(Key<?> key, AlgoInputSelection input, Activity activity)
					throws InterruptedException {
				if (key == sequenceName) {
					return referenceSelection.get() == ReferenceDetermination.SELECT;
				}
				return true;
			}

			@Override
			public void validate(Collection<? extends Key<?>> keys, AlgoInputSelection input,
					ParameterValidationHandler handler, Activity activity) throws InterruptedException {
				final boolean sequenceNameShouldBeSet = referenceSelection.get() == ReferenceDetermination.SELECT;
				final boolean sequenceNameIsSet = sequenceName.get() != null;
				if (sequenceNameShouldBeSet && !sequenceNameIsSet) {
					handler.postInconsistentParameters("Reference should be set", sequenceName, referenceSelection);
				}
			}
		};
	}
}
