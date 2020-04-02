package io.github.pdekker.viraltyping.algo.aligment;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import com.clcbio.api.base.algo.alignment.AlignmentAlgoFacade;
import com.clcbio.api.base.algo.alignment.AlignmentParameters;
import com.clcbio.api.base.algo.alignment.AlignmentParameters.AlignmentMode;
import com.clcbio.api.base.algo.alignment.AlignmentParameters.EndGapCost;
import com.clcbio.api.base.algorithm.Algo;
import com.clcbio.api.base.algorithm.AlgoException;
import com.clcbio.api.base.algorithm.AlgoHistoryTools;
import com.clcbio.api.base.algorithm.CallableExecutor;
import com.clcbio.api.base.algorithm.ChannelDescription;
import com.clcbio.api.base.algorithm.Multiplicity;
import com.clcbio.api.base.algorithm.OutputHandler;
import com.clcbio.api.base.algorithm.TemporaryObjectManager;
import com.clcbio.api.base.algorithm.parameter.AlgoParameters;
import com.clcbio.api.base.algorithm.resource.NonExclusive;
import com.clcbio.api.base.math.misc.Target;
import com.clcbio.api.base.process.Activity;
import com.clcbio.api.base.process.NullActivity;
import com.clcbio.api.base.session.ApplicationContext;
import com.clcbio.api.base.session.FactoryManager;
import com.clcbio.api.base.translation.Translator;
import com.clcbio.api.base.translation.TranslatorTools;
import com.clcbio.api.base.util.IteratorTools;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.Sequence;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.SequenceSource;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.alignment.Alignment;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.alignment.AlignmentSequenceIndexer;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.alignment.AlignmentTools;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.alignment.AlignmentTools.ComparisonResults;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.alphabet.Alphabet;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.alphabet.AlphabetTools;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.feature.Feature;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.feature.FeatureTypes;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.list.SequenceList;
import com.clcbio.api.free.datatypes.report.Report;

import io.github.pdekker.viraltyping.algo.aligment.AlignmentReportInterpreter.ReferenceDetermination;

@NonExclusive(minThreads = 1, maxThreads = -1)
public class AlignmentReportAlgo extends Algo {

	public static final String ID = "sars_cov2_alignment_report";
	public static final double VERSION = 1.0;
	public static final String NAME = "Alignment report";
	private static final byte GAP = AlphabetTools.getDnaAlphabet().getSymbolIndex('-');

	private final Translator DNA_TRANSLATOR = TranslatorTools.getStandardDnaToProteinTranslator();

	public static final ChannelDescription<Alignment> INPUT_CHANNEL = ChannelDescription.create("Alignment",
			Alignment.class, "alignment", Multiplicity.ONE);
	public static final ChannelDescription<Report> OUTPUT_CHANNEL = ChannelDescription.create("Table", Report.class,
			"table");

	public AlignmentReportAlgo(final ApplicationContext applicationContext) {
		super(applicationContext);
		addInputChannel(INPUT_CHANNEL.createDefaultInputChannel());
		addOutputChannel(OUTPUT_CHANNEL.createDefaultOutputChannel());
	}

	@Override
	public void alignParametersToChannelUse(Set<? extends ChannelDescription<?>> usedChannels,
			Target<String> alignmentProblems) {
		super.alignParametersToChannelUse(usedChannels, alignmentProblems);
	}

	@Override
	public void checkParametersAndInput(Target<String> problems) {
		super.checkParametersAndInput(problems);
	}

	@Override
	public void calculate(final OutputHandler handler, final CallableExecutor objectModificationExecutor)
			throws AlgoException, InterruptedException {

		final TemporaryObjectManager tom = new TemporaryObjectManager();
		try {
			final AlignmentReportInterpreter p = getInterpreter(getParameters());
			final AlignmentReportBuilder reportBuilder = new AlignmentReportBuilder();
			final ReferenceDetermination referenceDetermination = p.referenceSelection.get();
			final boolean ignoreGapAtEnd = p.ignoreGapsAtEnd.get();
			
			final String referenceId = referenceDetermination == ReferenceDetermination.SELECT ? p.sequenceName.get()
					: null;

			final Alignment aln = (Alignment) getInputObjectsIterator().next();
			final int referenceIndex = getReferenceIndex(aln, referenceDetermination, referenceId);

			final int steps = aln.getSequence(referenceIndex).getFeatureCount(FeatureTypes.CDS) + 1;
			
			final Activity child = startActivity(getActivity(), "Processing  " + aln.getName(), 1.0 / steps, handler);


			processAlignment(referenceIndex, aln, ignoreGapAtEnd, reportBuilder);
			endActivity(child);

			Iterator<Feature> it = aln.getSequence(referenceIndex).getFeatureIterator(FeatureTypes.CDS);
			while (it.hasNext() ) {
				Feature cds = it.next();
				final Activity translateActivity = startActivity(getActivity(), "Translating " + cds.getName(), 1.0 / steps,
						handler);
				final Alignment prot = createProteinAlignment(aln, cds, tom, translateActivity);
				if (prot == null || prot.getSequenceCount() != aln.getSequenceCount()) {
					throw new AlgoException("Protein translation failed, alignment did not work");
				}
				processAlignment(referenceIndex, prot, ignoreGapAtEnd, reportBuilder);
			}
			final Report report = reportBuilder.createReport(NullActivity.INSTANCE);
			report.startNoUndoBlock();
			report.setName(aln.getName());
			report.addHistory(AlgoHistoryTools.createEnrichedEntry(report, this));
			report.endNoUndoBlock();

			handler.postOutputObjects(report, this);
			postToChannel(OUTPUT_CHANNEL, report);
		} finally {
			tom.disposeAll();
		}
	}

	private Alignment createProteinAlignment(Alignment aln, Feature cds, TemporaryObjectManager tom, Activity child)
			throws AlgoException, InterruptedException {
		final List<Sequence> list = new ArrayList<>();
	
		for (int i = 0; i < aln.getSequenceCount(); i++) {
			final Sequence cdsSeq = extractCds(aln.getSequence(i), cds);
			tom.registerClcObject(cdsSeq);
			final Sequence translated = DNA_TRANSLATOR.translate(cdsSeq);
			translated.startNoUndoBlock();
			translated.setName("Translation of " + cds.getName() );
			translated.endNoUndoBlock();
			tom.disposeClcObject(cdsSeq);
			tom.registerClcObject(translated);
			list.add(translated);
		}

		final SequenceList input = FactoryManager.getInstance().getSequenceListFactory().createSequenceListWith(cds.getName(),
				IteratorTools.upcast(list.iterator()));
		tom.registerClcObject(input);

		final AlignmentParameters p = new AlignmentParameters(new AlgoParameters());
		p.setToDefault();
		p.alignmentMode.put(AlignmentMode.FAST);
		p.endGapCost.put(EndGapCost.CHEAP);

		final Alignment protAln = AlignmentAlgoFacade.getInstance()
				.align(Collections.<SequenceSource>singletonList(input), p, child);
		tom.disposeClcObject(input);
		tom.registerClcObject(protAln);
		tom.disposeClcObjects(list);
		return protAln;
	}

	private static Sequence extractCds(Sequence seq, Feature cds) throws AlgoException {
		return seq.getSubsequence(cds.getRegion());
	}

	private static String postionToString(int pos, int insertion) {
		if (insertion == 0) {
			return "" + pos;
		}
		return "" + pos + "." + insertion;
	}

	private void processAlignment(int referenceIndex, Alignment aln, boolean ignoreGapAtEnd,
			AlignmentReportBuilder reportBuilder) {
		final Alphabet alphabet = aln.getSequence(0).getAlphabet();
		final AlignmentSequenceIndexer refIndex = new AlignmentSequenceIndexer(aln, referenceIndex);
		final int refStart = refIndex.getSequenceStart();
		final int refEnd = refIndex.getSequenceEnd();
		int humanReadablePostion = 1;
		int humanReadableInsertionPosition = 0;
		reportBuilder.startAlignment(aln, referenceIndex);
		for (int i = refStart; i < refEnd; i++) {
			final char[] symbols = getSequenceSymbols(aln, ignoreGapAtEnd, i, alphabet);
			if (!sameSymbols(symbols, referenceIndex)) {
				final String pos = postionToString(humanReadablePostion, humanReadableInsertionPosition);
				reportBuilder.addMutationData(pos, symbols, referenceIndex);
			}
			if (symbols[referenceIndex] == GAP) {
				humanReadableInsertionPosition++;
			} else {
				humanReadableInsertionPosition = 0;
				humanReadablePostion++;
			}
		}
		reportBuilder.endAlignment();
	}

	private static char getSymbol(Alphabet alphabet, byte index) {
		return alphabet.getSymbol(index).getCharName();
	}

	private static char[] getSequenceSymbols(Alignment aln, boolean ignoreGapAtEnd, int alignmentPos,
			Alphabet alphabet) {
		final int[] pos = aln.getSequencePositions(alignmentPos);
		final char[] symbols = new char[pos.length];
		for (int index = 0; index < pos.length; index++) {
			final int seqPos = pos[index];
			boolean startOrEndGap = false;
			if (seqPos < 0 && ignoreGapAtEnd) {
				final int negPos = ~seqPos - 1;
				if (negPos < 0 || negPos >= aln.getSequence(index).getLength() - 1) {
					startOrEndGap = true;
				}
			}
			if (seqPos < 0) {
				symbols[index] = startOrEndGap ? ' ' : '-';
			} else {
				symbols[index] = getSymbol(alphabet, aln.getSequence(index).getResidueIndexAt(seqPos));
			}
		}
		return symbols;
	}

	private static boolean sameSymbols(char[] symbols, int referenceIndex) {
		if (symbols.length < 2) {
			return true;
		}
		final char ref = symbols[referenceIndex];
		for (int i = 1; i < symbols.length; i++) {
			if (ref != symbols[i] && symbols[i] != ' ') {
				return false;
			}
		}
		return true;
	}

	private int getReferenceIndex(Alignment aln, ReferenceDetermination referenceDetermination, String referenceId) {
		switch (referenceDetermination) {
		case AUTOMATIC:
			return getAutomaticReferenceIndex(aln);
		case SELECT:
			for (int i = 0; i < aln.getSequenceCount(); i++) {
				if (aln.getSequence(i).getId().equals(referenceId)) {
					return i;
				}
			}
			// fall through
		case FIRST:
		default:
			return 0;

		}
	}

	private int getAutomaticReferenceIndex(Alignment aln) {
		if (aln.getSequenceCount() < 3) {
			return 0;
		}
		int maxIdentities = 0;
		int bestReference = 0;
		for (int i = 0; i < aln.getSequenceCount(); i++) {
			int sumIdentities = 0;
			for (int j = 0; j < aln.getSequenceCount(); j++) {
				if (i != j) {
					// TODO: buffer identities;
					final ComparisonResults comparisons = AlignmentTools.getComparisons(aln, i, j);
					sumIdentities += comparisons.getIdentities();
				}
			}
			if (maxIdentities < sumIdentities) {
				maxIdentities = sumIdentities;
				bestReference = i;
			}
		}
		return bestReference;
	}

	@Override
	public String getName() {
		return NAME;
	}

	private Activity startActivity(final Activity act, final String msg, final double d, final OutputHandler handler)
			throws InterruptedException {
		handler.postStatus(msg, this);
		final Activity child = act.getChildActivity(d);
		child.setCurrentActivity(msg);
		child.checkStop();
		return child;
	}

	private void endActivity(final Activity act) throws InterruptedException {
		act.checkStop();
		act.setProgress(1.0);
	}

	@Override
	protected AlignmentReportInterpreter getInterpreter(final AlgoParameters parameters) {
		return new AlignmentReportInterpreter(parameters);
	}

	@Override
	public String getClassKey() {
		return ID;

	}

	@Override
	public double getVersion() {
		return VERSION;
	}
}
