package io.github.pdekker.viraltyping.action.consensus;

import com.clcbio.api.base.algorithm.Algo;
import com.clcbio.api.base.algorithm.parameter.AlgoParameters;
import com.clcbio.api.base.algorithm.parameter.keys.Key;
import com.clcbio.api.base.algorithm.parameter.keys.KeyChecker;
import com.clcbio.api.base.util.CreateList;
import com.clcbio.api.base.util.CreateSet;
import com.clcbio.api.clc.gui.wizard.WizardBuilder;
import com.clcbio.api.clc.gui.wizard.WizardContentsFactory;
import com.clcbio.api.clc.gui.wizard.WizardFacade;
import com.clcbio.api.clc.gui.wizard.WizardGroupAppender;
import com.clcbio.api.clc.gui.wizard.WizardGroupBuilder;
import com.clcbio.api.clc.gui.wizard.WizardState;
import com.clcbio.api.clc.gui.wizard.WizardStepBuilder;
import com.clcbio.api.free.actions.framework.ActionGroup;
import com.clcbio.api.free.algorithm.AlgoAction;
import com.clcbio.api.free.algorithm.wizard.AlgoSaveWizardStepModel;
import com.clcbio.api.free.datatypes.ClcObject;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.ReadMappingObject;
import com.clcbio.api.free.gui.components.MultiSelectClassRestrictor;
import com.clcbio.api.free.gui.components.MultiSelectRestrictor;
import com.clcbio.api.free.gui.icon.ClcIcon;
import com.clcbio.api.free.gui.icon.DefaultClcIcon;
import com.clcbio.api.free.wizard.dynamic.ClcWizardStepModel;

import io.github.pdekker.viraltyping.action.actiongroup.ViralTypingActionGroup;
import io.github.pdekker.viraltyping.algo.consensus.ConsensusAlgo;
import io.github.pdekker.viraltyping.algo.consensus.ConsensusInterpreter;

public class ExtractConsensusAction extends AlgoAction {
	private static final long serialVersionUID = 1L;

	public static final String PLUGIN_GROUP = "genomics";

	@Override
	protected void addToActionGroup() {
		final ActionGroup ag = manager.getActionManager().findActionGroup(ViralTypingActionGroup.CLASS_KEY);
		if (ag != null) {
			ag.addAction(this);
		}
	}

	@Override
	public Algo createAlgo() {
		return new ConsensusAlgo(getManager());
	}

	@Override
	public String getClassKey() {
		return "com.clcbio.plugins.ps.action.consensus.ExtractConsensusAction";
	}

	@Override
	protected int getOutputObjectsCount(final AlgoParameters parameters, final ClcObject[] selectedObjects) {
		return 1;
	}

	@Override
	public ClcWizardStepModel getFirstStep(final AlgoParameters parameters, final ClcWizardStepModel nextStep) {

		final ConsensusInterpreter p = new ConsensusInterpreter(parameters);
		final KeyChecker keyChecker = p.createKeyChecker(getManager());
		final WizardFacade facade = WizardFacade.getInstance();

		final WizardContentsFactory wcf = facade.getWizardContentsFactory();
		final WizardBuilder builder = facade.createWizardBuilder(parameters, nextStep);

		// Second wizard step...
		// Read coverage filters group
		final WizardStepBuilder variantStep = builder.prependStep(p.secondPageGroup.getId(),
				p.secondPageGroup.getTitle());
		final WizardState state = variantStep.getState();

		final WizardGroupBuilder readFiltersBuilder = variantStep.appendGroup(p.coverageSettingsGroup.getTitle());
		readFiltersBuilder.appendWidget(wcf.widget(state, p.ignoreBrokenPairs));
		readFiltersBuilder.appendWidget(wcf.widget(state, p.ignoreNonSpecificMatches));
		final WizardGroupAppender readFilterAppender = facade.createDefaultGroupAppender(keyChecker,
				CreateSet.of(p.minimumIgnoreReadLength));
		readFiltersBuilder.append(readFilterAppender, state, wcf);

		// Read Qualtiy filters group
		final WizardGroupBuilder qualityFiltersBuilder = variantStep.appendGroup(p.qualityFilterGroup.getTitle());
		qualityFiltersBuilder.appendWidget(wcf.widget(state, p.useQualityFilter));
		qualityFiltersBuilder.increaseIndent();
		final WizardGroupAppender qualityAppender = facade.createDefaultGroupAppender(keyChecker,
				CreateList.<Key<?>>of(p.qualityRadius, p.qualityMinCentral, p.qualityMinRegion));
		qualityFiltersBuilder.append(qualityAppender, state, wcf);
		qualityFiltersBuilder.decreaseIndent();

		variantStep.addConstraint(wcf.constraints(state).asSpecified(keyChecker,
				p.secondPageGroup.extractContainedKeys(p.getKeyObjects())));

		// first wizard step
		final WizardStepBuilder settingsStep = builder.prependStep(p.firstPageGroup.getId(),
				p.firstPageGroup.getTitle());
		final WizardState settingState = settingsStep.getState();

		final WizardGroupBuilder variantBuilder = settingsStep.appendGroup(p.conflictResolutionGroup.getTitle());
		variantBuilder.appendWidget(wcf.widget(state, p.conflictResolution));

		final WizardGroupAppender variantAppender = facade.createDefaultGroupAppender(keyChecker,
				CreateList.<Key<?>>of(p.minCoverage, p.minFrequency, p.addConflictAnnotations));
		variantBuilder.append(variantAppender, settingState, wcf);

		// Read Qualtiy filters group
		final WizardGroupBuilder primerBuilder = settingsStep.appendGroup(p.extensionSettingsGroup.getTitle());
		primerBuilder.appendWidget(wcf.widget(settingState, p.extendStartEnd));
		primerBuilder.appendWidget(wcf.widget(settingState, p.minCoverageExtend));
		primerBuilder.appendWidget(wcf.widget(settingState, p.trimPrimers));
		final WizardGroupAppender primerAppender = facade.createDefaultGroupAppender(keyChecker,
				CreateSet.of(p.trimLinkerList));
		primerAppender.append(primerBuilder, settingState, wcf);

		final WizardGroupBuilder indelBuilder = settingsStep.appendGroup(p.inDelResolutionGroup.getTitle());
		final WizardGroupAppender indelAppender = facade.createDefaultGroupAppender(keyChecker,
				CreateList.<Key<?>>of(p.inDelResolution, p.minBreakpoint));
		indelAppender.append(indelBuilder, settingState, wcf);

		settingsStep.addConstraint(wcf.constraints(settingState).asSpecified(keyChecker,
				p.firstPageGroup.extractContainedKeys(p.getKeyObjects())));

		return builder.getFirstStep();
	}

	@Override
	protected AlgoSaveWizardStepModel getAlgoSaveWizardStepModel(final AlgoParameters parameters) {
		final ConsensusInterpreter p = new ConsensusInterpreter(parameters);
		return WizardFacade.getInstance().createDefaultSaveStepModel(p.createKeyChecker(getManager()),
				p.getKeyObjects());
	}

	@Override
	public String getName() {
		return ConsensusAlgo.NAME;
	}

	@Override
	public String getToolTip() {
		return "Extract consensus sequence from read mapping.";
	}

	@Override
	public int getPreferredMenuLocation() {
		return 15;
	}

	@Override
	public ClcIcon createIcon() {
		return new DefaultClcIcon("actions/reassemble");
	}

	@Override
	public MultiSelectRestrictor createRestrictor(final WarningReceptor warningReceptor) {
		return new MultiSelectClassRestrictor(new Class[] { ReadMappingObject.class }, "Select Read Mappings");
	};
}
