package io.github.pdekker.viraltyping.algo.transferannot;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
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
import com.clcbio.api.base.algorithm.AlgoOutputNamingTools;
import com.clcbio.api.base.algorithm.CallableExecutor;
import com.clcbio.api.base.algorithm.ChannelDescription;
import com.clcbio.api.base.algorithm.OutputHandler;
import com.clcbio.api.base.algorithm.TemporaryObjectManager;
import com.clcbio.api.base.algorithm.parameter.AlgoParameters;
import com.clcbio.api.base.math.misc.Target;
import com.clcbio.api.base.process.Activity;
import com.clcbio.api.base.process.NullActivity;
import com.clcbio.api.base.session.ApplicationContext;
import com.clcbio.api.base.session.FactoryManager;
import com.clcbio.api.base.translation.Translator;
import com.clcbio.api.base.translation.TranslatorTools;
import com.clcbio.api.base.util.Comprehension;
import com.clcbio.api.base.util.CreateList;
import com.clcbio.api.base.util.IteratorTools;
import com.clcbio.api.base.util.StringTools;
import com.clcbio.api.free.datatypes.ClcObject;
import com.clcbio.api.free.datatypes.bioinformatics.gis.Chromosome;
import com.clcbio.api.free.datatypes.bioinformatics.gis.Genome;
import com.clcbio.api.free.datatypes.bioinformatics.gis.TrackFactory;
import com.clcbio.api.free.datatypes.bioinformatics.gis.track.FeatureTrack;
import com.clcbio.api.free.datatypes.bioinformatics.gis.track.FeatureTrackBuilder;
import com.clcbio.api.free.datatypes.bioinformatics.gis.track.SymbolTrack;
import com.clcbio.api.free.datatypes.bioinformatics.gis.track.SymbolTrackBuilder;
import com.clcbio.api.free.datatypes.bioinformatics.gis.track.Track;
import com.clcbio.api.free.datatypes.bioinformatics.gis.track.VariantTrack;
import com.clcbio.api.free.datatypes.bioinformatics.gis.track.VariantTrackBuilder;
import com.clcbio.api.free.datatypes.bioinformatics.gis.track.VariantTracklet;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.BasicSequence;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.NucleotideSequence;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.Sequence;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.SequenceFactory;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.SequenceSource;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.alignment.Alignment;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.alignment.AlignmentSequenceIndexer;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.feature.Feature;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.feature.FeatureTools;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.feature.FeatureTypes;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.index.BasicIndexer;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.interval.Interval;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.interval.SimpleInterval;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.list.SequenceList;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.position.OpenPosition;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.position.Position;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.position.SimplePosition;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.region.Region;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.region.RegionTools;
import com.clcbio.api.free.datatypes.bioinformatics.variant.Variant;
import com.clcbio.api.free.datatypes.framework.history.HistoryEntry;
import com.clcbio.api.free.datatypes.report.Report;

import io.github.pdekker.viraltyping.algo.consensus.ConsensusBuilder;

public class TransferAnnotationsAlgo extends Algo {

	public static final String ID = "transfer_annotations";
	private static final double VERSION = 2.2;
	public static final String NAME = "Transfer Annotations";

	private final NumberFormat percFormatter = getFormatter();
	
	public final static String matpept = "Mature peptide";
	public final static String gene = FeatureTypes.GENE;

	public static final ChannelDescription<NucleotideSequence> INPUT_CHANNEL = new ChannelDescription<NucleotideSequence>(
			"Sequence", "Sequence", NucleotideSequence.class, "Sequence");
	public static final ChannelDescription<NucleotideSequence> NUCL_OUTPUT = ChannelDescription
			.create("Annotated Sequence", NucleotideSequence.class, "annotated-sequence");

	public static final ChannelDescription<SymbolTrack> SYMBOL_TRACK_OUTPUT = ChannelDescription.create("Genome",
			SymbolTrack.class, "symbol_track");
	public static final ChannelDescription<FeatureTrack> CDS_TRACK_OUTPUT = ChannelDescription.create("CDS Track",
			FeatureTrack.class, "cds_track");

	public static final ChannelDescription<FeatureTrack> FAILURE_TRACK_OUTPUT = ChannelDescription
			.create("Failure Track", FeatureTrack.class, "failure_track");

	public static final ChannelDescription<FeatureTrack> GENE_TRACK_OUTPUT = ChannelDescription.create("Gene Track",
			FeatureTrack.class, "gene_track");

	public static final ChannelDescription<VariantTrack> VARIANT_TRACK_OUTPUT = ChannelDescription
			.create("Variant Track", VariantTrack.class, "variant_track");

	public static final ChannelDescription<Report> REPORT = ChannelDescription.create("Report", Report.class, "report");

	public final SequenceFactory FACTORY = FactoryManager.getInstance().getSequenceFactory();

	private final Translator DNA_TRANSLATOR = TranslatorTools.getStandardDnaToProteinTranslator();

	public TransferAnnotationsAlgo(final ApplicationContext applicationContext) {
		super(applicationContext);
		addInputChannel(INPUT_CHANNEL.createDefaultInputChannel());
		addOutputChannel(NUCL_OUTPUT.createDefaultOutputChannel());
		addOutputChannel(SYMBOL_TRACK_OUTPUT.createDefaultOutputChannel());
		addOutputChannel(CDS_TRACK_OUTPUT.createDefaultOutputChannel());
		addOutputChannel(GENE_TRACK_OUTPUT.createDefaultOutputChannel());
		addOutputChannel(VARIANT_TRACK_OUTPUT.createDefaultOutputChannel());
		addOutputChannel(FAILURE_TRACK_OUTPUT.createDefaultOutputChannel());
		addOutputChannel(REPORT.createDefaultOutputChannel());
	}

	private static NumberFormat getFormatter() {
		final NumberFormat formatter = NumberFormat.getPercentInstance();
		formatter.setMinimumFractionDigits(1);
		return formatter;
	}

	@Override
	public void checkParametersAndInput(Target<String> problems) {
		super.checkParametersAndInput(problems);
		final TransferAnnotationInterpreter p = getInterpreter(getParameters());

		if (getInputObjectsCount() != 1) {
			problems.put("Select exactly one sequence as input!");
		}
		final Sequence input = (Sequence) getInputObjectsIterator().next();
		final Sequence annot = (Sequence) p.reference.getClcObject(getApplicationContext());
		if (!input.getAlphabet().equals(annot.getAlphabet())) {
			problems.put("Don't mix protein and nucleotide sequences");
		}
	}

	@Override
	public void calculate(final OutputHandler handler, final CallableExecutor objectModificationExecutor)
			throws AlgoException, InterruptedException {

		final TransferAnnotationInterpreter p = getInterpreter(getParameters());

		final TemporaryObjectManager tom = new TemporaryObjectManager();

		final TransferAnnotationReportBuilder reportBuilder = p.createReport.get()
				? new TransferAnnotationReportBuilder()
				: null;
		try {

			int failure = 0;

			final CdsFeatureProcessor processor = new CdsFeatureProcessor();

			final NucleotideSequence inputSeq = (NucleotideSequence) getInputObjectsIterator().next();
			final Sequence annotatedInput = (Sequence) p.reference.getClcObject(getApplicationContext());

			final Activity child = startActivity(getActivity(), "Processing " + inputSeq.getName(), 1.0, handler);

			final Sequence annotatedOutput = transferAnnotations(inputSeq, annotatedInput, tom, child, processor);

			final List<ClcObject> output = CreateList.of();

			if (annotatedOutput == null) {
				failure++;
				handler.postMessage(inputSeq.getName() + " did not align, could not transfer annotations", this);
			} else {
				final HistoryEntry he = AlgoHistoryTools.createEnrichedEntry(annotatedInput, this);
				annotatedOutput.startNoUndoBlock();
				annotatedOutput.addHistory(he);
				annotatedOutput.setName(AlgoOutputNamingTools.createRetaggedName(inputSeq.getName(), "annot"));
				annotatedOutput.endNoUndoBlock();
				output.add(annotatedOutput);
				postToChannel(NUCL_OUTPUT, inputSeq);
			}
			endActivity(child);

			final List<Track<?>> tracks = CreateList.of();
			final Genome genome = Genome.fromSequences(Collections.singleton(annotatedOutput));

			final SymbolTrack st = asSymbolTrack(annotatedOutput, genome);
			output.add(st);
			postToChannel(SYMBOL_TRACK_OUTPUT, st);

			final VariantTrack varTrack = asVariantTrack(annotatedOutput, genome);
			if (varTrack != null) {
				postToChannel(VARIANT_TRACK_OUTPUT, varTrack);
				output.add(varTrack);
			}

			//CDS
			final FeatureTrack cdsTrack = asFeatureTrack(annotatedOutput, genome, FeatureTypes.CDS);
			if (reportBuilder != null) {
				reportBuilder.addCds(processor.reportRows);
			}
			if (cdsTrack != null) {
				postToChannel(CDS_TRACK_OUTPUT, cdsTrack);
				output.add(cdsTrack);
			}

			final FeatureTrack failureTrack = asFeatureTrack(annotatedOutput, genome, ConsensusBuilder.FAILURES);
			if (failureTrack != null) {
				postToChannel(FAILURE_TRACK_OUTPUT, failureTrack);
				output.add(failureTrack);
			}

			final FeatureTrack geneTrack = asFeatureTrack(annotatedOutput, genome, gene);
			if (geneTrack != null) {
				postToChannel(GENE_TRACK_OUTPUT, geneTrack);
				output.add(geneTrack);
			}

			output.addAll(tracks);
			if (reportBuilder != null) {
				final Report rep = reportBuilder.createReport(NullActivity.INSTANCE);
				rep.startNoUndoBlock();
				rep.setName(AlgoOutputNamingTools.createRetaggedName(inputSeq.getName(), "Report"));
				rep.addHistory(AlgoHistoryTools.createEnrichedEntry(rep, this));
				rep.endNoUndoBlock();
				output.add(rep);
				postToChannel(REPORT, rep);
			}
			if (processor.warning) {
				handler.postMessage("Please check CDS regions!", this);
			}

			if (failure > 0) {
				handler.postMessage("Some sequence did not align, please check input!", this);
			}
			handler.postOutputObjects(output, this);

		} finally {
			tom.disposeAll();
		}

	}

//	public interface FeatureProcessor extends Mapper<Feature, Feature> {
//		public void init(BasicSequence bs);
//	}

	private final class CdsFeatureProcessor {

		private final List<String[]> reportRows = new ArrayList<>();
		private boolean warning = false;

		public void updateFeature(Feature f, BasicSequence seq) {
			// for cds we have to double check if we still have a valid reading
			// frame.
			final Region r = f.getRegion();
			final String oldName = f.getName();
			final boolean validLength = r.getSize() % 3 == 0;
			final List<String> msg = new ArrayList<>();

			final int expectedLen = getExpectedProtLength(f);
			int translatedLen = expectedLen;
			f.removeAnnotation("translation");

			final String protString = getProteinString(seq, r);

			final String percentageOrf;

			int startPos = r.getFirstPos().getMin();
			int endPos = r.getLastPos().getMax();
			boolean uncertainStart = false;
			boolean uncertainEnd = false;

			boolean warning = true;

			if (validLength && protString.startsWith("M") && protString.indexOf('*') == protString.length() - 1) {
				// everything looks goods, we are ready with this cds...
				f.addAnnotation("translation", protString);
				final double perc = protString.length() * 1.0 / expectedLen;
				percentageOrf = percFormatter.format(perc);
				msg.add("Ok");
				warning = false;

			} else {
				final int stop = protString.indexOf('*');

				final boolean internalStop = stop >= 0 && stop < protString.length() - 1;
				if (!validLength && internalStop) {
					msg.add("Frameshift detected");
					uncertainStart = true;
					uncertainEnd = true;

				} else {
					if (!protString.startsWith("M")) {
						final int end = r.getIntervalAt(0).getLastPos().getMax();
						final int start = end % 3;
						// the M should be between start & end..
						final String aa = getProteinString(seq, new Region(start, end));
						final int lastStop = Math.max(0, aa.lastIndexOf('*')); // converts -1 (not found) to 0
						final int firstMethionine = aa.indexOf('M', lastStop);

						if (firstMethionine == -1) {
							msg.add("Start codon not found");
							uncertainStart = true;
						} else {
							startPos = start + firstMethionine * 3;
							msg.add("Start codon moved");
							translatedLen -= firstMethionine;

						}
					}
					if (!protString.endsWith("*")) {
						int lengthPreviousExons = 0;
						for (int i = 0; i < r.numIntervals() - 1; i++) {
							lengthPreviousExons += r.getIntervalAt(i).getMaxLength();
						}
						final Interval last = r.getIntervalAt(r.numIntervals() - 1);
						final int start = last.getFirstPos().getMin() + lengthPreviousExons % 3; // keep in frame with
																									// previous exons;
						final int length = (seq.getLength() - start) / 3;
						final int end = start + length * 3;
						final String aa = getProteinString(seq, new Region(start, end));
						final int firstStop = aa.indexOf('*');
						if (firstStop == -1) {
							msg.add("Stop codon not found");
							uncertainEnd = true;
						} else {
							endPos = (firstStop + 1) * 3 + start;
							msg.add("Stop codon moved");
							translatedLen += firstStop;
						}
					}
					if (protString.indexOf('*') > 0 && protString.indexOf('*') != protString.length() - 1) {
						msg.add("Internal stop codon found");
						uncertainEnd = true;
						translatedLen -= protString.length() - protString.indexOf('*');
					}
				}

				final List<Interval> copies = new ArrayList<>();
				for (int i = 0; i < r.numIntervals(); i++) {
					final Interval interval = r.getIntervalAt(i);
					final Position start = i == 0 ? getPostion(startPos, uncertainStart) : interval.getFirstPos();
					final Position end = i == r.numIntervals() - 1 ? getPostion(endPos, uncertainEnd)
							: interval.getLastPos();
					copies.add(new SimpleInterval(start, end, interval.getStrand()));
				}
				f.setRegion(new Region(copies));

				final double perc = translatedLen * 1.0 / expectedLen;

				if (uncertainStart || uncertainEnd) {
					f.setName("Invalid");
					f.addAnnotation("Error", "Invalid translated protein");
					percentageOrf = "~" + percFormatter.format(perc);
				} else {
					percentageOrf = percFormatter.format(perc);
					f.addAnnotation("translation", getProteinString(seq, f.getRegion()));
				}
			}
			if (warning) {
				this.warning = true;
			}

			final String[] reportRow = new String[] { seq.getName(), oldName, percentageOrf, String.join(",", msg) };
			reportRows.add(reportRow);
		}

		private Position getPostion(int pos, boolean uncertain) {
			return uncertain ? OpenPosition.valueOf(pos) : SimplePosition.valueOf(pos);
		}

		private int getExpectedProtLength(Feature f) {
			final Object[] annot = f.getAnnotation("translation");
			if (annot != null && annot.length > 0) {
				return annot[0].toString().length();
			}
			return f.getRegion().getSize() / 3;
		}

		private String getProteinString(BasicSequence seq, Region r) {
			final BasicSequence bs = FactoryManager.getInstance().getSequenceFactory().createSubSequence(seq, r);
			final BasicSequence prot = DNA_TRANSLATOR.translateBasicSequence(bs);
			final String result = prot.getString();
			disposeWhenNeeded(bs);
			disposeWhenNeeded(prot);
			return result;
		}

		private void disposeWhenNeeded(BasicSequence bs) {
			if (bs instanceof Sequence) {
				((Sequence) bs).increaseUsage();
				((Sequence) bs).decreaseUsage();
			}
		}
	}

	private FeatureTrack asFeatureTrack(Sequence copy, Genome genome, String... types) throws InterruptedException {

		final Set<String> typeSet = new HashSet<String>(Arrays.asList(types));
		final HashSet<String> foundTypes = new HashSet<String>(FeatureTools.getTypes(copy));
		if (!foundTypes.removeAll(typeSet)) {
			return null;
		}
		final String suffix = StringTools.concatenateWithSeparator(",", types);

		final TrackFactory factory = FactoryManager.getInstance().getTrackFactory();
		final String name = AlgoOutputNamingTools.createRetaggedName(copy.getName(), suffix);
		final FeatureTrackBuilder builder = factory.createFeatureTrackBuilder(name, genome);
		int index = 0;

		final Chromosome chr = genome.chromosome(index++);
		final Iterator<Feature> it = copy.getFeatureIterator(typeSet, 0, copy.getLength());
		final int count = copy.getFeatureCount(typeSet, 0, copy.getLength());

		builder.add(factory.createFeatureTracklet(chr, it, count, NullActivity.INSTANCE));
		return builder.finish(typeSet);
	}

	private VariantTrack asVariantTrack(Sequence copy, Genome genome) throws InterruptedException {
		if (!FeatureTools.getTypes(copy).contains(FeatureTypes.CONFLICT)) {
			return null;
		}
		final TrackFactory factory = FactoryManager.getInstance().getTrackFactory();
		final String name = AlgoOutputNamingTools.createRetaggedName(copy.getName(), "Variant");
		final VariantTrackBuilder builder = factory.createVariantTrackBuilder(name, genome);
		int index = 0;

		final Chromosome chr = genome.chromosome(index++);
		final Iterator<Feature> it = copy.getFeatureIterator(FeatureTypes.CONFLICT);
		final Iterator<Iterator<Variant>> itIt = IteratorTools.getMappedIterator(it,
				source -> ConsensusBuilder.asVariant(source).iterator());

		final VariantTracklet tracklet = factory.createVariantTracklet(chr,
				Comprehension.flattenIteratorOfIterators(itIt), NullActivity.INSTANCE);

		builder.add(tracklet);
		return builder.finish();
	}

	private SymbolTrack asSymbolTrack(Sequence copy, Genome genome) throws InterruptedException {
		final TrackFactory factory = FactoryManager.getInstance().getTrackFactory();
		final String name = AlgoOutputNamingTools.createRetaggedName(copy.getName(), "Genome");
		final SymbolTrackBuilder builder = factory.createSymbolTrackBuilder(name, genome, copy.getAlphabet());

		builder.add(factory.createSymbolTracklet(copy.getName(), copy.getSymbolIterator(), copy.getAlphabet(),
				NullActivity.INSTANCE));
		return builder.finish();
	}

	private Sequence transferAnnotations(Sequence inputSeq, Sequence annotSeq, TemporaryObjectManager tom,
			Activity child, CdsFeatureProcessor processor) throws AlgoException, InterruptedException {

		BasicSequence inputCopy = null;
		BasicSequence annotCopy = null;

		try {
			inputCopy = FactoryManager.getInstance().getSequenceFactory().createBasicSequence(inputSeq);
			if (inputCopy instanceof ClcObject) {
				((ClcObject) inputCopy).increaseUsage();
			}

			annotCopy = FactoryManager.getInstance().getSequenceFactory().createBasicSequence(annotSeq);
			if (annotCopy instanceof ClcObject) {
				((ClcObject) annotCopy).increaseUsage();
			}

			final SequenceList input = FactoryManager.getInstance().getSequenceListFactory()
					.createSequenceListWith("alignment", inputCopy, annotCopy);
			tom.registerClcObject(input);
			final AlignmentParameters p = new AlignmentParameters(new AlgoParameters());
			p.setToDefault();
			p.alignmentMode.put(AlignmentMode.FAST);
			p.endGapCost.put(EndGapCost.CHEAP);

			final Alignment aln = AlignmentAlgoFacade.getInstance()
					.align(Collections.<SequenceSource>singletonList(input), p, child);
			tom.disposeClcObject(input);
			tom.registerClcObject(aln);
			if (aln != null && aln.getSequenceCount() == 2) {
				final int inputIndex;
				final int annotIndex;
				if (aln.getSequence(0).getName().equals(inputCopy.getName())) {
					inputIndex = 0;
					annotIndex = 1;
				} else {
					inputIndex = 1;
					annotIndex = 0;
				}

				if (inputCopy instanceof ClcObject) {
					((ClcObject) inputCopy).startNoUndoBlock();
				}
				final BasicIndexer inputProvider = new AlignmentSequenceIndexer(aln, inputIndex);
				final BasicIndexer annotProvider = new AlignmentSequenceIndexer(aln, annotIndex);

				for (final Feature f : annotCopy.getFeatures()) {
					final Region rAnnot = RegionTools.convertToAlignmentCoordinates(f.getRegion(), annotProvider);
					final Region rInput = RegionTools.convertFromAlignmentCoordinates(rAnnot, inputProvider);
					final Feature feature = f.getFeatureWithRegion(rInput);
					feature.addAnnotation("note", "\"Copied from " + annotCopy.getName() + "\"");
					if (feature.getType().equals(FeatureTypes.CDS)) {
						processor.updateFeature(feature, inputCopy);
					}
					inputCopy.addFeature(feature);
				}
				if (inputCopy instanceof ClcObject) {
					((ClcObject) inputCopy).endNoUndoBlock();
				}
				tom.disposeClcObject(aln);
				return FactoryManager.getInstance().getSequenceFactory().createSequence(inputCopy);
			} else {
				return null;
			}
		} finally {
			if (inputCopy instanceof ClcObject) {
				((ClcObject) inputCopy).decreaseUsage();
			}
			if (annotCopy instanceof ClcObject) {
				((ClcObject) annotCopy).decreaseUsage();
			}
		}
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
	protected TransferAnnotationInterpreter getInterpreter(final AlgoParameters parameters) {
		return new TransferAnnotationInterpreter(parameters);
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
