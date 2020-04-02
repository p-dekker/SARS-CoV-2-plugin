package io.github.pdekker.viraltyping.action.reportmerger;

import java.util.Collection;
import java.util.List;
import java.util.Set;

import com.clcbio.api.base.algorithm.Algo;
import com.clcbio.api.base.algorithm.parameter.AlgoParameters;
import com.clcbio.api.base.algorithm.parameter.keys.KeyChecker;
import com.clcbio.api.base.algorithm.selection.MinElementCountDoneConstraint;
import com.clcbio.api.base.math.misc.Source;
import com.clcbio.api.base.util.CreateList;
import com.clcbio.api.base.util.CreateSet;
import com.clcbio.api.base.util.MapperTools;
import com.clcbio.api.clc.datatypes.report.ReportCompositeElement;
import com.clcbio.api.clc.gui.wizard.WizardBuilder;
import com.clcbio.api.clc.gui.wizard.WizardConstraint;
import com.clcbio.api.clc.gui.wizard.WizardContentsFactory;
import com.clcbio.api.clc.gui.wizard.WizardErrorHandler;
import com.clcbio.api.clc.gui.wizard.WizardFacade;
import com.clcbio.api.clc.gui.wizard.WizardGroupBuilder;
import com.clcbio.api.clc.gui.wizard.WizardStepBuilder;
import com.clcbio.api.free.actions.framework.ActionGroup;
import com.clcbio.api.free.algorithm.AlgoAction;
import com.clcbio.api.free.algorithm.wizard.InputError;
import com.clcbio.api.free.datatypes.ClcObject;
import com.clcbio.api.free.datatypes.report.Report;
import com.clcbio.api.free.datatypes.report.ReportElement;
import com.clcbio.api.free.gui.components.MultiSelectClassRestrictor;
import com.clcbio.api.free.gui.components.MultiSelectRestrictor;
import com.clcbio.api.free.gui.icon.ClcIcon;
import com.clcbio.api.free.gui.icon.DefaultClcIcon;
import com.clcbio.api.free.wizard.dynamic.ClcWizardStepModel;

import io.github.pdekker.viraltyping.action.actiongroup.ViralTypingActionGroup;
import io.github.pdekker.viraltyping.algo.reportmerger.ReportMergerAlgo;
import io.github.pdekker.viraltyping.algo.reportmerger.ReportMergerInterpreter;

public class ReportMergerAction extends AlgoAction {
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
		return new ReportMergerAlgo(getManager());
	}

	@Override
	public String getClassKey() {
		return ReportMergerAlgo.ID;
	}

	@Override
	protected int getOutputObjectsCount(final AlgoParameters parameters, final ClcObject[] selectedObjects) {
		return 1;
	}

	@Override
	public String getName() {
		return ReportMergerAlgo.NAME;
	}

	@Override
	public String getToolTip() {
		return "Merge multple reports";
	}

	@Override
	public int getPreferredMenuLocation() {
		return 60;
	}

	@Override
	public ClcIcon createIcon() {
		return new DefaultClcIcon("actions/proteinreport");
	}

	@Override
	public ClcWizardStepModel getFirstStep(final AlgoParameters parameters, final ClcWizardStepModel nextStep) {
		final ReportMergerInterpreter p = new ReportMergerInterpreter(parameters);
		final WizardFacade facade = WizardFacade.getInstance();
		final WizardContentsFactory wcf = facade.getWizardContentsFactory();
		final WizardBuilder builder = facade.createWizardBuilder(parameters, nextStep);

		// first wizards step
		final WizardStepBuilder settingStep = builder.prependStep(p.firstPageGroup.getId(),
				p.firstPageGroup.getTitle());
		settingStep.addConstraint(new WizardConstraint() {
			@Override
			public void checkValidity(WizardErrorHandler errorHandler) {
				if (p.elements.getAll().size() == 0) {
					errorHandler.put(InputError.createInconsistentInputError("Select at least one report element"));
				}
			}

			@Override
			public Set<String> getParameterKeys() {
				return CreateSet.of(p.elements.getKey());
			}
		});
		final WizardGroupBuilder captionGroupBuilder = settingStep.appendGroup(p.inputGroup.getTitle());
		final Source<Collection<String>> sourceElements = new Source<Collection<String>>() {
			@Override
			public Collection<String> get() {
				final List<String> elements = CreateList.of();
				final ClcObject[] inputNow = settingStep.getState().getStepModel().getSelectedObjects();
				for (final ClcObject o : inputNow) {
					if (o instanceof Report) {
						final Report report = (Report) o;
						for (final ReportElement re : report.getReportElements()) {
							addElements(re, elements);
						}
					}
				}
				return elements;
			}

			private void addElements(ReportElement re, List<String> elements) {
				if (re instanceof ReportCompositeElement) {
					elements.add(((ReportCompositeElement) re).getCaption());
//					for (ReportElement inner : ((ReportCompositeElement) re).getReportElements()) {
//						addElements(inner, elements);
//					}
				}
			}
		};
		captionGroupBuilder.appendWidget(wcf.browsers(settingStep.getState()).pickList(p.elements, sourceElements,
				MapperTools.idMapper(String.class)));
		final KeyChecker checker = p.createKeyChecker(getManager());
		final WizardGroupBuilder optionGroupBuilder = settingStep.appendGroup(p.optionGroup.getTitle());
		optionGroupBuilder.append(
				facade.createDefaultGroupAppender(checker, p.optionGroup.extractContainedKeys(p.getKeyObjects())),
				settingStep.getState(), wcf);
		return builder.getFirstStep();
	}

	@Override
	public MultiSelectRestrictor createRestrictor(final WarningReceptor warningReceptor) {
		return new MultiSelectClassRestrictor(new Class[] { Report.class }, "Select Reports")
				.allowDoneWhen(new MinElementCountDoneConstraint(1));
	};
}
