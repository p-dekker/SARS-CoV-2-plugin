package io.github.pdekker.viraltyping.algo.transferannot;

import java.util.Collection;
import java.util.Set;

import com.clcbio.api.base.algorithm.ParameterValidationHandler;
import com.clcbio.api.base.algorithm.parameter.AlgoParameters;
import com.clcbio.api.base.algorithm.parameter.AlgoParametersInterpreter;
import com.clcbio.api.base.algorithm.parameter.ParameterGroup;
import com.clcbio.api.base.algorithm.parameter.keys.BooleanKey;
import com.clcbio.api.base.algorithm.parameter.keys.ClcObjectKey;
import com.clcbio.api.base.algorithm.parameter.keys.Key;
import com.clcbio.api.base.algorithm.parameter.keys.KeyContainer;
import com.clcbio.api.base.algorithm.parameter.keys.Keys;
import com.clcbio.api.base.session.ApplicationContext;
import com.clcbio.api.base.util.CreateSet;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.NucleotideSequence;

public class TransferAnnotationInterpreter extends AlgoParametersInterpreter {
	protected final KeyContainer keys;

	public final ParameterGroup firstPageGroup = ParameterGroup.topLevel("settings", "Settings");

	public final ParameterGroup inputGroup = ParameterGroup.childOf(firstPageGroup, "Input");
	public final ParameterGroup optionGroup = ParameterGroup.childOf(firstPageGroup, "Options");

	public final ClcObjectKey reference = Keys.newClcObjectKey(this, "reference").withMultiplicity(1, 1)
			.withOptionKey("reference").labelled("Annotated Sequence").inGroup(inputGroup)
			.allowedType(NucleotideSequence.class).done();

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

	public TransferAnnotationInterpreter(final AlgoParameters parameters) {
		super(parameters);
		keys = new KeyContainer(reference, createReport);
	}

	@Override
	public void setToDefault() {
		keys.setToDefault();
	}

	@Override
	public String getClassKey() {
		return "sars_cov2_transfer_annotation_parameters";
	}

	@Override
	public Collection<Key<?>> getKeyObjects() {
		return keys.getKeys();
	}

	@Override
	public Set<String> getKeys() {
		return keys.getKeySet();
	}

}
