package io.github.pdekker.viraltyping.algo.aligment;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.clcbio.api.base.algorithm.AlgoOutputNamingTools;
import com.clcbio.api.base.math.misc.SingleInt;
import com.clcbio.api.clc.algorithms.report.AbstractReportCalculator;
import com.clcbio.api.clc.datatypes.report.ClcReportTableModel;
import com.clcbio.api.clc.datatypes.report.ReportCompositeElement;
import com.clcbio.api.clc.datatypes.report.ReportTableElement;
import com.clcbio.api.clc.datatypes.report.ReportTextElement;
import com.clcbio.api.clc.datatypes.report.SimpleReport;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.BasicSequence;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.alignment.Alignment;
import com.clcbio.api.free.datatypes.report.Report;
import com.clcbio.api.free.datatypes.report.ReportElement;

class AlignmentReportBuilder extends AbstractReportCalculator {
	Integer columns = null;
	String alnname = null;
	String[] header = null;
	List<String[]> rows = new ArrayList<>();
	List<ReportElement> elements = new ArrayList<>();

	int lastPos = -1;

	public void startAlignment(final Alignment aln, final int referenceIndex) {
		if (columns != null || alnname != null) {
			throw new IllegalStateException("Previous alignment is not properly finished");
		}
		final BasicSequence ref = aln.getSequence(referenceIndex);

		alnname = AlgoOutputNamingTools.extractBaseString(aln.getName());
		columns = aln.getSequenceCount() + 2;

		int index = 0;
		header = new String[columns];
		header[index++] = "Position";
		header[index++] = AlgoOutputNamingTools.extractBaseString(ref.getName());
		for (int i = 0; i < aln.getSequenceCount(); i++) {
			if (i == referenceIndex) {
				continue;
			}
			header[index++] = AlgoOutputNamingTools.extractBaseString(aln.getSequence(i).getName());
		}
		header[index++] = "total";
	}

	public void addMutationData(String pos, char[] symbols, int referenceIndex) {
		if (alnname == null || columns == null || header == null) {
			throw new IllegalStateException("Data collections is not started");
		}
		final Map<Character, SingleInt> counts = new LinkedHashMap<>();

		final String[] data = new String[columns];
		int index = 0;
		data[index++] = pos;
		final char refSymbol = symbols[referenceIndex];
		data[index++] = "" + refSymbol;
		for (int i = 0; i < symbols.length; i++) {
			final SingleInt count = counts.get(symbols[i]);
			if (count != null) {
				count.n++;
			} else {
				counts.put(symbols[i], new SingleInt(1));
			}
			if (i != referenceIndex) {
				data[index++] = "" + symbols[i];
			}
		}
		data[index++] = counts.entrySet().stream().map(e -> e.getKey() + " (" + e.getValue() + ")")
				.collect(Collectors.joining(", "));
		rows.add(data);
	}

	public void endAlignment() {
		if (alnname == null || columns == null || header == null) {
			throw new IllegalStateException("Data collections is not started");
		}
		final ReportCompositeElement rce = new ReportCompositeElement();
		rce.setCaption("Mutations compared to: " + alnname);
		if (rows.isEmpty()) {
			final ReportTextElement rte = new ReportTextElement("No mutation for this alignment");
			rce.addReportElement(rte);
		} else {
			final ClcReportTableModel model = new ClcReportTableModel(header, rows.toArray(new String[rows.size()][]));
			rce.addReportElement(new ReportTableElement(model));

		}
		elements.add(rce);
		alnname = null;
		columns = null;
		header = null;
		rows.clear();
	}

	@Override
	protected Report createReport(List<ReportElement> elements) throws InterruptedException {
		return new SimpleReport(elements, null);
	}

	@Override
	protected List<ReportElement> createReportElements() throws InterruptedException {
		return elements;
	}

}