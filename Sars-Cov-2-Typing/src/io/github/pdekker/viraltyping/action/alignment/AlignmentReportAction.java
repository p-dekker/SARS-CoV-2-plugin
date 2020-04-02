package io.github.pdekker.viraltyping.action.alignment;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.clcbio.api.base.algorithm.Algo;
import com.clcbio.api.base.algorithm.parameter.AlgoParameters;
import com.clcbio.api.base.algorithm.selection.SelectionAddConstraint;
import com.clcbio.api.base.math.misc.Source;
import com.clcbio.api.clc.gui.wizard.WizardBuilder;
import com.clcbio.api.clc.gui.wizard.WizardContentsFactory;
import com.clcbio.api.clc.gui.wizard.WizardFacade;
import com.clcbio.api.clc.gui.wizard.WizardGroupBuilder;
import com.clcbio.api.clc.gui.wizard.WizardState;
import com.clcbio.api.clc.gui.wizard.WizardStepBuilder;
import com.clcbio.api.free.actions.framework.ActionGroup;
import com.clcbio.api.free.algorithm.AlgoAction;
import com.clcbio.api.free.datatypes.ClcObject;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.alignment.Alignment;
import com.clcbio.api.free.gui.components.MultiSelectClassRestrictor;
import com.clcbio.api.free.gui.components.MultiSelectRestrictor;
import com.clcbio.api.free.gui.icon.ClcIcon;
import com.clcbio.api.free.gui.icon.DefaultClcIcon;
import com.clcbio.api.free.wizard.dynamic.ClcWizardStepModel;

import io.github.pdekker.viraltyping.action.actiongroup.ViralTypingActionGroup;
import io.github.pdekker.viraltyping.algo.aligment.AlignmentReportAlgo;
import io.github.pdekker.viraltyping.algo.aligment.AlignmentReportInterpreter;
import io.github.pdekker.viraltyping.algo.aligment.AlignmentReportInterpreter.ReferenceDetermination;

public class AlignmentReportAction extends AlgoAction {
	private static final long serialVersionUID = 1;
	public static final String PLUGIN_GROUP = "free";

	@Override
	protected void addToActionGroup() {
		final ActionGroup ag = manager.getActionManager().findActionGroup(ViralTypingActionGroup.CLASS_KEY);
		if (ag != null) {
			ag.addAction(this);
		}
	}

	@Override
	public Algo createAlgo() {
		return new AlignmentReportAlgo(getManager());
	}

	@Override
	public String getName() {
		return AlignmentReportAlgo.NAME;
	}

	@Override
	public double getVersion() {
		return AlignmentReportAlgo.VERSION;
	}

	@Override
	public String getClassKey() {
		return AlignmentReportAlgo.ID;
	}

	@Override
	public int getPreferredMenuLocation() {
		return 90;
	}

	@Override
	public ClcIcon createIcon() {
		return new DefaultClcIcon("actions/proteinreport");
	}

	@Override
	public String getToolTip() {
		return "Create an alignment report";
	}

	@Override
	public MultiSelectRestrictor createRestrictor(final WarningReceptor warningReceptor) {
		return new MultiSelectClassRestrictor(new Class<?>[] { Alignment.class }, "Select an alignment")
				.allowAddWhen(SelectionAddConstraint.SINGLE_ELEMENT);
	}

	@Override
	public ClcWizardStepModel getFirstStep(AlgoParameters parameters, ClcWizardStepModel nextStep) {
		final AlignmentReportInterpreter p = new AlignmentReportInterpreter(parameters);
		final WizardFacade facade = WizardFacade.getInstance();
		final WizardContentsFactory wcf = facade.getWizardContentsFactory();
		final WizardBuilder builder = facade.createWizardBuilder(parameters, nextStep);

		final WizardStepBuilder step = builder.prependStep(p.firstPageGroup.getId(), p.firstPageGroup.getTitle());
		final WizardState state = step.getState();

		final WizardGroupBuilder referenceGroup = step.appendGroup(p.referenceGroup.getTitle());
		referenceGroup.appendWidget(wcf.widget(state, p.referenceSelection));

		final Source<List<String>> sequenceIds = () -> {
			final Alignment aln = getInputAlignment(state);
			final List<String> result = IntStream.range(0, aln.getSequenceCount())
					.mapToObj(i -> aln.getSequence(i).getId()).collect(Collectors.toList());

			if (!result.contains(p.sequenceName.get())) {
				p.sequenceName.put(result.get(0));
			}
			return result;
		};

		final Source<Boolean> selectSequence = () -> p.referenceSelection.get() == ReferenceDetermination.SELECT;

		referenceGroup.appendWidget(wcf.comboBoxes(state.enabledWhen(selectSequence)).dynamicNamedValues(p.sequenceName,
				sequenceIds, mapper -> {
					final Alignment aln = getInputAlignment(state);
					for (int i = 0; i < aln.getSequenceCount(); i++) {
						if (aln.getSequence(i).getId().equals(mapper)) {
							return aln.getSequence(i).getName();
						}
					}
					throw new AssertionError("Id not present in alignment");
				}));

		final WizardGroupBuilder optionGroup = step.appendGroup(p.optionsGroup.getTitle());
		optionGroup.appendWidget(wcf.widget(state, p.ignoreGapsAtEnd));

		return builder.getFirstStep();
	}

	@Override
	protected boolean canDoBatchExecution() {
		return false;
	}

	private static Alignment getInputAlignment(WizardState state) {
		if (!state.getStepModel().knowsSelectedObjects()) {
			throw new AssertionError("Cannot select annotation types without knowing input");
		}
		final ClcObject[] inputNow = state.getStepModel().getSelectedObjects();
		if (inputNow.length != 1 || !(inputNow[0] instanceof Alignment)) {
			throw new AssertionError("No or invalid input!");
		}
		return (Alignment) inputNow[0];
	}

	@Override
	protected int getOutputObjectsCount(AlgoParameters parameters, ClcObject[] selectedObjects) {
		return 1;

	}

}
