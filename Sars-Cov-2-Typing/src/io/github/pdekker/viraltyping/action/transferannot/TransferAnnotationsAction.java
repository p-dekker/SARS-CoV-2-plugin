package io.github.pdekker.viraltyping.action.transferannot;

import com.clcbio.api.base.algorithm.Algo;
import com.clcbio.api.base.algorithm.parameter.AlgoParameters;
import com.clcbio.api.base.algorithm.parameter.keys.KeyChecker;
import com.clcbio.api.clc.gui.wizard.WizardFacade;
import com.clcbio.api.free.actions.framework.ActionGroup;
import com.clcbio.api.free.algorithm.AlgoAction;
import com.clcbio.api.free.algorithm.wizard.AlgoSaveWizardStepModel;
import com.clcbio.api.free.datatypes.ClcObject;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.NucleotideSequence;
import com.clcbio.api.free.gui.components.MultiSelectClassRestrictor;
import com.clcbio.api.free.gui.components.MultiSelectRestrictor;
import com.clcbio.api.free.gui.icon.ClcIcon;
import com.clcbio.api.free.gui.icon.DefaultClcIcon;
import com.clcbio.api.free.wizard.dynamic.ClcWizardStepModel;

import io.github.pdekker.viraltyping.action.actiongroup.ViralTypingActionGroup;
import io.github.pdekker.viraltyping.algo.transferannot.TransferAnnotationInterpreter;
import io.github.pdekker.viraltyping.algo.transferannot.TransferAnnotationsAlgo;

public class TransferAnnotationsAction extends AlgoAction {
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
		return new TransferAnnotationsAlgo(getManager());
	}

	@Override
	public String getClassKey() {
		return TransferAnnotationsAlgo.ID;
	}

	@Override
	protected int getOutputObjectsCount(final AlgoParameters parameters, final ClcObject[] selectedObjects) {
		return 1;
	}

	@Override
	public ClcWizardStepModel getFirstStep(final AlgoParameters parameters, final ClcWizardStepModel nextStep) {
		final TransferAnnotationInterpreter p = new TransferAnnotationInterpreter(parameters);
		final KeyChecker keyChecker = p.createKeyChecker(getManager());
		final WizardFacade facade = WizardFacade.getInstance();
		return facade.createDefaultParameterSteps(keyChecker, p.getKeyObjects(), nextStep);
	}

	@Override
	protected AlgoSaveWizardStepModel getAlgoSaveWizardStepModel(final AlgoParameters parameters) {
		final TransferAnnotationInterpreter p = new TransferAnnotationInterpreter(parameters);
		return WizardFacade.getInstance().createDefaultSaveStepModel(p.createKeyChecker(getManager()),
				p.getKeyObjects());
	}

	@Override
	public String getName() {
		return TransferAnnotationsAlgo.NAME;
	}

	@Override
	public String getToolTip() {
		return "Transfer annotations";
	}

	@Override
	public int getPreferredMenuLocation() {
		return 20;
	}

	@Override
	public ClcIcon createIcon() {
		return new DefaultClcIcon("clcobjects/feature");
	}

	@Override
	public MultiSelectRestrictor createRestrictor(final WarningReceptor warningReceptor) {
		return new MultiSelectClassRestrictor(new Class[] { NucleotideSequence.class }, "Select Sequence");
	};
}
