package io.github.pdekker.viraltyping.algo.consensus;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.clcbio.api.base.algorithm.Algo;
import com.clcbio.api.base.algorithm.AlgoException;
import com.clcbio.api.base.algorithm.AlgoHistoryTools;
import com.clcbio.api.base.algorithm.AlgoOutputNamingTools;
import com.clcbio.api.base.algorithm.CallableExecutor;
import com.clcbio.api.base.algorithm.ChannelDescription;
import com.clcbio.api.base.algorithm.OutputHandler;
import com.clcbio.api.base.algorithm.SimpleAlgoExecuter;
import com.clcbio.api.base.algorithm.TemporaryObjectManager;
import com.clcbio.api.base.algorithm.blast.BlastFactory;
import com.clcbio.api.base.algorithm.blast.BlastProgram;
import com.clcbio.api.base.algorithm.blast.parameters.LocalBlastProgramAlgoParameters;
import com.clcbio.api.base.algorithm.blast.parameters.impl.local.LocalBlastnParameters;
import com.clcbio.api.base.algorithm.parameter.AlgoParameters;
import com.clcbio.api.base.algorithm.parameter.interpreters.DeNovoV2Parameters;
import com.clcbio.api.base.algorithm.parameter.interpreters.DeNovoWithMappingParameters;
import com.clcbio.api.base.algorithm.parameter.interpreters.DeNovoWithMappingParameters.MappingMode;
import com.clcbio.api.base.algorithm.parameter.interpreters.TrimInterpreter;
import com.clcbio.api.base.algorithm.resource.NonExclusive;
import com.clcbio.api.base.math.misc.DoubleInt;
import com.clcbio.api.base.math.misc.Target;
import com.clcbio.api.base.process.Activity;
import com.clcbio.api.base.process.NullActivity;
import com.clcbio.api.base.session.ApplicationContext;
import com.clcbio.api.base.session.FactoryManager;
import com.clcbio.api.base.util.CreateList;
import com.clcbio.api.base.util.iterator.BulkByteIterator;
import com.clcbio.api.clc.datatypes.bioinformatics.blast.BlastHit;
import com.clcbio.api.clc.datatypes.bioinformatics.blast.BlastHsp;
import com.clcbio.api.clc.datatypes.bioinformatics.blast.BlastOutput;
import com.clcbio.api.free.datatypes.ClcObject;
import com.clcbio.api.free.datatypes.bioinformatics.gis.track.SymbolTrack;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.BasicSequence;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.NucleotideSequence;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.Sequence;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.SequenceFactory;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.alphabet.AlphabetTools;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.feature.Feature;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.list.SequenceList;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.list.SequenceListBuilder;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.list.SequenceListBuilderByReadGroup;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.readgroup.ReadGroup;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.region.Region;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.symbol.SymbolSource;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.MatchIntersection;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.MatchList;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.MatchListIntersections;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.MatchListSession;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.PairedEndLocalAlignment;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.ReadMapping;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.ReadMappingObject;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.ScatteredLocalAlignment;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.SequenceCluster;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.SequenceClusterList;
import com.clcbio.api.free.datatypes.bioinformatics.sequencecluster.UnalignedEnds;
import com.clcbio.api.free.datatypes.bioinformatics.trim.TrimAdapterList;
import com.clcbio.api.free.datatypes.framework.history.HistoryEntry;
import com.clcbio.api.free.datatypes.report.Report;
import com.clcbio.api.genomics.base.algo.assembly.AssemblyAlgoFacade;
import com.clcbio.api.genomics.base.algo.trim.TrimAlgoFacade;
import com.clcbio.api.genomics.base.algo.trim.Trimmer;

import io.github.pdekker.viraltyping.algo.consensus.ConsensusBuilder.CoverageInformation;
import io.github.pdekker.viraltyping.algo.consensus.ConsensusInterpreter.ConflictResolution;
import io.github.pdekker.viraltyping.algo.consensus.ConsensusInterpreter.PrimerMode;

@NonExclusive(minThreads = 1, maxThreads = 1)
public class ConsensusAlgo extends Algo {

	public static final String ID = "sars_cov2_consensus_creator";
	private static final double VERSION = 1.5;
	public static final String NAME = "Extract Consensus";
	public static long algoVersionUID = 2202964220434324834L;
	private static byte N = AlphabetTools.getDnaAlphabet().getMostAmbiguousSymbolIndex();

	public static final ChannelDescription<ReadMappingObject> INPUT_CHANNEL = new ChannelDescription<ReadMappingObject>(
			"Read Mapping", "Read Mapping", ReadMappingObject.class, "read-mapping");
	public static final ChannelDescription<NucleotideSequence> CONSENSUS_OUTPUT = ChannelDescription
			.create("Consensus", NucleotideSequence.class, "consensus");
	public static final ChannelDescription<SymbolTrack> CONSENSUS_OUTPUT_TRACK = ChannelDescription.create("Consensus",
			SymbolTrack.class, "consensus-track");
	public static final ChannelDescription<Report> CONSENSUS_REPORT = ChannelDescription.create("Report", Report.class,
			"report");

	public ConsensusAlgo(final ApplicationContext applicationContext) {
		super(applicationContext);
		addInputChannel(INPUT_CHANNEL.createDefaultInputChannel());
		addOutputChannel(CONSENSUS_OUTPUT.createDefaultOutputChannel());
		addOutputChannel(CONSENSUS_REPORT.createDefaultOutputChannel());
	}

	@Override
	public void alignParametersToChannelUse(Set<? extends ChannelDescription<?>> usedChannels,
			Target<String> alignmentProblems) {
		super.alignParametersToChannelUse(usedChannels, alignmentProblems);
		final ConsensusInterpreter p = getInterpreter(getParameters());
		p.createReport.put(usedChannels.contains(CONSENSUS_REPORT));
	}

	@Override
	public void checkParametersAndInput(Target<String> problems) {
		super.checkParametersAndInput(problems);
	}

	@Override
	public void calculate(final OutputHandler handler, final CallableExecutor objectModificationExecutor)
			throws AlgoException, InterruptedException {

		final ConsensusInterpreter p = new ConsensusInterpreter(getParameters());

		final TemporaryObjectManager tom = new TemporaryObjectManager();

		final Trimmer trimmer = createTrimmer(p.trimPrimers.get(),
				(TrimAdapterList) p.trimLinkerList.getClcObject(getApplicationContext()));

		final List<String> problematicSegments = CreateList.of();

		try {
			final int minCoverage = p.minCoverage.get();
			final int minCoverageExtend = p.minCoverageExtend.get();
			final double minFrequency = p.minFrequency.get();
			final boolean addConflicts = p.addConflictAnnotations.get();
			final double minProblematicVariantFreq = p.inDelResolution.get() ? p.minBreakpoint.get() : 2.0; // ratio 2.0
																											// is
																											// impossoble
																											// so mean's
																											// we don't
																											// use this
																											// option.
			final ConsensusReportBuilder reportBuilder = p.createReport.get() ? new ConsensusReportBuilder() : null;

			final ConflictResolution conflictResolution = p.conflictResolution.get();

			final SequenceListBuilder listBuilder = FactoryManager.getInstance().getSequenceListFactory()
					.createBuilder();
			tom.registerCleanupable(listBuilder);

			final List<ClcObject> output = CreateList.of();

			final boolean extend = p.extendStartEnd.get();
			final ReadMapping mapping = ((ReadMappingObject) getInputObjectsIterator().next()).asSequenceMapping();
//			for (int index = 0; index < mapping.size(); index++) {
				final SymbolSource mainSequence = mapping.getMainSequence(0);

				final MatchList matches = mapping.getMatchList(0);
				final Activity child = startActivity(getActivity(), "Processing " + mainSequence.getName(),
						1.0 / mapping.size(), handler);
				final ConsensusBuilder builder = ConsensusBuilder.createBuilder(mainSequence, extend, minCoverage,
						minCoverageExtend, minFrequency, addConflicts, conflictResolution, minProblematicVariantFreq);
				tom.registerCleanupable(builder);
				final ConsensusIterator it = new ConsensusIterator(builder, matches, mainSequence, p);

				final DoubleInt extension = extend ? new DoubleInt(0, 0) : null;

				int mainStart = 0;
				int mainEnd = mainSequence.getLength();

				if (extend) {
					mainStart = it.searchStart();
					extension.n1 = it.fixStart(mainStart);
					mainEnd = it.searchEnd();
				}

				it.iterate(mainStart, mainEnd, false, child);

				if (extend) {
					extension.n2 = it.fixEnd(mainEnd);
				}
				builder.finish();

				Sequence cons = builder.getConsensus();
				final Region toBeFixed = ConsensusBuilder.getBreakPointRegion(cons);

				if (toBeFixed != null) {
					handler.postStatus("Running local de novo to improve consensus", this);
					child.setCurrentActivity("Running local de novo to improve consensus");
					final boolean fixed = runLocalDeNovo(cons, toBeFixed, matches, tom, handler);
					if (!fixed) {
						handler.postStatus("Consensus was not updated", this);
						child.setCurrentActivity("Consensus was not updated");
					}

				}

				DoubleInt trimmedBases = null;
				if (trimmer != null) {
					final DoubleInt trimRegion = trimmer.getGoodRegionBounds(cons);
					trimmedBases = new DoubleInt(Math.min(trimRegion.n1, trimRegion.n2),
							cons.getLength() - Math.max(trimRegion.n1, trimRegion.n2));

					if (trimRegion != null) {
						final int start = Math.min(trimRegion.n1, trimRegion.n2);
						final int end = Math.max(trimRegion.n1, trimRegion.n2);

						if (start > 0 || end < cons.getLength()) {
							tom.registerClcObject(cons);
							final Sequence trimmed = cons.getSubsequence(new Region(start, end));
							tom.disposeClcObject(cons);
							cons = trimmed;
						}
					}
				}
				if (reportBuilder != null) {
					final CoverageInformation coverInfo = builder.getCoverageInformation();
					reportBuilder.addCoverageInformation( cons, coverInfo);
					reportBuilder.addFragmentInformation( cons, extension, trimmedBases);
				}
				if (sequenceContainsN(cons)) {
					problematicSegments.add(cons.getName());
				}
				tom.disposeCleanupable(builder);
				endActivity(child);

			final NucleotideSequence result = (NucleotideSequence) cons;
			tom.deregisterCleanupable(listBuilder);

			final HistoryEntry he1 = AlgoHistoryTools.createEnrichedEntry(result, this);
			result.startNoUndoBlock();
			result.addHistory(he1);
			result.setName(AlgoOutputNamingTools.createRetaggedName(mapping.getObject().getName(), "consensus"));
			result.endNoUndoBlock();
			postToChannel(CONSENSUS_OUTPUT, result);
			output.add(result);

			if (reportBuilder != null) {
				final HistoryEntry he2 = AlgoHistoryTools.createEnrichedEntry(result, this);
				final Report report = reportBuilder.createReport(NullActivity.INSTANCE);
				report.startNoUndoBlock();
				report.addHistory(he2);
				report.setName(AlgoOutputNamingTools.createRetaggedName(mapping.getObject().getName(), "report"));
				report.endNoUndoBlock();
				postToChannel(CONSENSUS_REPORT, report);
				output.add(report);
			}

			if (!problematicSegments.isEmpty()) {
				final String msg = problematicSegments.stream()
						.collect(Collectors.joining(", ", "Following segments have N in sequence: ", "."));
				handler.postMessage(msg, this);
				;
			}

			handler.postOutputObjects(output, this);
		} finally {
			tom.disposeAll();
		}
	}

	private boolean sequenceContainsN(Sequence cons) {
		final BulkByteIterator it = cons.getSymbolIterator();
		while (it.hasNext()) {
			if (it.next() == N) {
				return true;
			}
		}
		return false;
	}

	private boolean runLocalDeNovo(Sequence cons, Region toBeFixed, MatchList matches, TemporaryObjectManager tom,
			OutputHandler handler) throws AlgoException, InterruptedException {
		// extract reads;
		final List<SequenceList> reads = extractReads(matches, toBeFixed, tom);
		tom.registerClcObjects(reads);

		if (reads.isEmpty()) {
			return false;
		}
		final ClcObject contigs = runDeNovo(reads, tom);
		tom.disposeClcObjects(reads);
		if (contigs == null) {
			return false;
		}
		//handler.postOutputObjects(contigs, this); // enable for debugging puproses
		final BlastOutput output = runBlast(cons, contigs, tom);
		return updateConsensus(cons, output);
	}

	private boolean updateConsensus(Sequence cons, BlastOutput output) {
		if (output == null) {
			return false;
		}
		final BlastHit[] hits = output.getIteration(0).getHits();
		if (hits == null || hits.length == 0) {
			return false;
		}
		final BlastHsp best = hits[0].getHsp(0);
		if (best.getIdentity() == best.getAlignLen()) {
			// we have no mismatches so we can't update consensus...
			return false;
		}
		final int start = Math.max(0, best.getQueryFrom() - 1); // blast is 1 based
		final int end = best.getQueryTo();
		final String hseq = best.getHseq().replaceAll("-", "");
		final byte[] seqSymbols = AlphabetTools.convertToBytes(hseq, AlphabetTools.getDnaAlphabet());
		cons.startNoUndoBlock();
		cons.replaceSymbols(seqSymbols, start, end);
		cons.addFeature(
				new Feature("Local assembly", new Region(start, start + seqSymbols.length), ConsensusBuilder.FAILURES));
		cons.endNoUndoBlock();
		return true;
	}

	private List<SequenceList> extractReads(MatchList matches, Region toBeFixed, TemporaryObjectManager tom) {
		final SequenceListBuilderByReadGroup builder = FactoryManager.getInstance().getSequenceListFactory()
				.createBuilderByReadGroup("reads");
		tom.registerCleanupable(builder);

		final MatchListSession session = MatchListSession.forMatchList(matches);
		final MatchIntersection intersection = MatchListIntersections.mainRange(session,
				toBeFixed.getFirstPos().getMin(), toBeFixed.getLastPos().getMax());
		final List<ReadGroup> groups = matches.getReadGroups();
		final SequenceFactory seqFac = FactoryManager.getInstance().getSequenceFactory();

		while (intersection.findMatch(UnalignedEnds.EXCLUDE)) {
			final ScatteredLocalAlignment sla = intersection.currentMatch();
			final Region r = toBeFixed.getIntersection(sla.getMainStartPosition(), sla.getMainEndPosition());
			if (r == null) {
				continue;
			}
			final ReadGroup rg = sla.getReadGroupIndex() == 0 ? null : groups.get(sla.getReadGroupIndex() - 1);

			if (sla instanceof PairedEndLocalAlignment) {
				final PairedEndLocalAlignment pla = (PairedEndLocalAlignment) sla;
				final BasicSequence left = seqFac.createBasicSequence(pla.getOriginalLeftSequence());
				final BasicSequence right = seqFac.createBasicSequence(pla.getOriginalRightSequence());
				builder.addPair(rg, left, right);
			} else {
				final BasicSequence read = seqFac.createBasicSequence(sla.getBasicSequence());
				builder.addSingle(rg, read);
			}
		}
		final List<SequenceList> reads = builder.finish();
		tom.deregisterCleanupable(builder);
		return reads;
	}

	private ClcObject runDeNovo(List<SequenceList> reads, TemporaryObjectManager tom) throws AlgoException {
		final AlgoParameters parms = new AlgoParameters();
		final DeNovoV2Parameters p1 = new DeNovoV2Parameters(parms);
		final DeNovoWithMappingParameters p2 = new DeNovoWithMappingParameters(parms);
		p1.setToDefault();
		p2.setToDefault();

		p1.performScaffolding.put(false);
		p2.mapReads.put(MappingMode.MAP_BACK);
		p2.updateContigs.put(true);
		p2.createReport.put(false);

		final Algo algo = AssemblyAlgoFacade.getInstance().createDeNovoAssemblyAlgo(getApplicationContext());
		final SimpleAlgoExecuter executer = new SimpleAlgoExecuter();
		try {
			executer.executeAndThrowExceptions(algo, parms, reads, NullActivity.INSTANCE);
		} catch (final Exception e) {
			throw new AlgoException(e);
		}

		ClcObject result = null;
		for (final ClcObject obj : executer.getOutputObjects()) {
			tom.registerClcObject(obj);
			if (obj instanceof SequenceCluster || obj instanceof SequenceClusterList) {
				result = obj;
			} else {
				tom.disposeClcObject(obj);
			}
		}
		return result;
	}

	private BlastOutput runBlast(Sequence cons, ClcObject contigs, TemporaryObjectManager tom) throws AlgoException {
		final List<Sequence> seqs = extractAndRenameContigs(contigs);
		final AlgoParameters p = new AlgoParameters();

		final Algo algo = BlastFactory.getInstance().getLocalBlastAlgo(getApplicationContext());
		final LocalBlastProgramAlgoParameters parm1 = new LocalBlastProgramAlgoParameters(p, getApplicationContext());
		final LocalBlastnParameters parm2 = new LocalBlastnParameters(p);
		parm2.setToDefault();
		parm1.setToDefault();
		parm1.setBlastProgram(BlastProgram.blastn);
		parm1.setSubject(seqs);
		parm2.setWord_size(13);
		parm2.setEvalue(1.0);
		parm2.setMaxNumHits(1);

		final SimpleAlgoExecuter executer = new SimpleAlgoExecuter();
		try {
			cons.increaseUsage();
			seqs.forEach(ClcObject::increaseUsage);
			contigs.increaseUsage();
			executer.executeAndThrowExceptions(algo, p, java.util.Collections.singletonList(cons), NullActivity.INSTANCE);
			seqs.forEach(ClcObject::decreaseUsage);; // deletes contigs!

		} catch (final Exception e) {
			throw new AlgoException(e);
		}

		BlastOutput output = null;
		for (final ClcObject obj : executer.getOutputObjects()) {
			if (obj instanceof BlastOutput) {
				output = (BlastOutput) obj;
			} else {
				tom.registerClcObject(obj);
				tom.disposeClcObject(obj);
			}
		}
		return output;
	}

	private List<Sequence> extractAndRenameContigs(ClcObject obj) {
		final List<Sequence> results = new ArrayList<>();
		if (obj instanceof SequenceClusterList) {
			final SequenceClusterList sl = (SequenceClusterList) obj;
			int i = 1;
			for (final SequenceCluster c : sl.getClusters()) {
				results.add(extractConsensusSequence(c, i++));
			}
		}
		if (obj instanceof SequenceCluster) {
			results.add(extractConsensusSequence((SequenceCluster) obj, 1));
		}
		return results;
	}

	private Sequence extractConsensusSequence(SequenceCluster sc, int i) {
		Sequence seq = FactoryManager.getInstance().getSequenceFactory().createSequence(sc.getMainSequence());
		seq.setName("contig-" + i);
		return seq;
	}

	private Trimmer createTrimmer(PrimerMode primerMode, TrimAdapterList list) {
		if (primerMode == PrimerMode.IGNORE) {
			return null;
		}
		final TrimInterpreter p = new TrimInterpreter(new AlgoParameters());
		p.setToDefault();
		p.keyAmbigiousTrim.put(false);
		p.keyQualityTrim.put(false);
		p.keyLeftTrim.put(false);
		p.keyRightTrim.put(false);
		p.keyTrimLinkers.putClcObject(list);
		return TrimAlgoFacade.getInstance().createTrimmer(getApplicationContext(), p.getAlgoParameters());
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
	protected ConsensusInterpreter getInterpreter(final AlgoParameters parameters) {
		return new ConsensusInterpreter(parameters);
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
