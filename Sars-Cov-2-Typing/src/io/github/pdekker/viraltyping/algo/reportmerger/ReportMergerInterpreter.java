package io.github.pdekker.viraltyping.algo.reportmerger;

import java.util.Collection;
import java.util.Set;

import com.clcbio.api.base.algorithm.parameter.AlgoParameters;
import com.clcbio.api.base.algorithm.parameter.AlgoParametersInterpreter;
import com.clcbio.api.base.algorithm.parameter.ParameterGroup;
import com.clcbio.api.base.algorithm.parameter.keys.BooleanKey;
import com.clcbio.api.base.algorithm.parameter.keys.Key;
import com.clcbio.api.base.algorithm.parameter.keys.KeyContainer;
import com.clcbio.api.base.algorithm.parameter.keys.Keys;
import com.clcbio.api.base.algorithm.parameter.keys.StringKey;

public class ReportMergerInterpreter extends AlgoParametersInterpreter {
	protected final KeyContainer keys;

	public final ParameterGroup firstPageGroup = ParameterGroup.topLevel("settings", "Settings");
	public final ParameterGroup inputGroup = ParameterGroup.childOf(firstPageGroup, "Elements");
	public final ParameterGroup optionGroup = ParameterGroup.childOf(firstPageGroup, "Options");

	public final StringKey elements = Keys.newStringKey(this, "elements").withMultiplicity(1, Integer.MAX_VALUE)
			.withOptionKey("elements").labelled("Elements").inGroup(inputGroup).done();

	public final BooleanKey createFrontPage = Keys.newBooleanKey(this, "frontPage").withOptionKey("frontpage")
			.labelled("Create frontpage").inGroup(optionGroup).defaultsTo(false).done();

	public ReportMergerInterpreter(final AlgoParameters parameters) {
		super(parameters);
		keys = new KeyContainer(elements, createFrontPage);
	}

	@Override
	public void setToDefault() {
		keys.setToDefault();
	}

	@Override
	public String getClassKey() {
		return "sars_cov2_report_merger_interpreter";
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
