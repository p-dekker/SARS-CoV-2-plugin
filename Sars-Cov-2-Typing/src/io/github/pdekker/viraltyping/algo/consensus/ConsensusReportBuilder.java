package io.github.pdekker.viraltyping.algo.consensus;

import java.text.NumberFormat;
import java.util.List;

import com.clcbio.api.base.math.FrequencyDistribution;
import com.clcbio.api.base.math.misc.DoubleInt;
import com.clcbio.api.base.util.CreateList;
import com.clcbio.api.clc.algorithms.report.AbstractReportCalculator;
import com.clcbio.api.clc.datatypes.report.ClcReportTableModel;
import com.clcbio.api.clc.datatypes.report.ReportCompositeElement;
import com.clcbio.api.clc.datatypes.report.ReportTableElement;
import com.clcbio.api.clc.datatypes.report.SimpleReport;
import com.clcbio.api.free.datatypes.bioinformatics.sequence.BasicSequence;
import com.clcbio.api.free.datatypes.report.Report;
import com.clcbio.api.free.datatypes.report.ReportElement;

import io.github.pdekker.viraltyping.algo.consensus.ConsensusBuilder.CoverageInformation;

public class ConsensusReportBuilder extends AbstractReportCalculator {

	public final List<ReportElement> elements = CreateList.of();

	List<String[]> fragmentData = CreateList.of();
	List<String[]> coverageData = CreateList.of();

	public void addCoverageInformation(BasicSequence bs, CoverageInformation coverInfo) {
		final NumberFormat nf = NumberFormat.getNumberInstance();
		nf.setMaximumFractionDigits(1);
		nf.setGroupingUsed(false);

		final FrequencyDistribution fd = coverInfo.fd;
		final String name = bs.getName();
		if (fd.isEmpty()) {
			fd.add(0); // prevents "null" in report
		}
		final String[] row = new String[5];
		row[0] = name;
		row[1] = "" + fd.getMin();
		row[2] = "" + fd.getMax();
		row[3] = nf.format(fd.getAverage()) + " \u00B1 " + nf.format(fd.getStandardDeviation());
		row[4] = coverInfo.lowCoverageRegions();

		coverageData.add(row);

	}

	public void addFragmentInformation( BasicSequence bs, DoubleInt extension, DoubleInt trimRegion) {
		final String[] row = new String[7];
		row[0] = bs.getName();
		row[6] = "" + bs.getLength();
		row[5] = ConsensusBuilder.hasFailures(bs) ? "Yes" : "No";
		if (extension == null) {
			// no extension choosen
			row[1] = "-";
			row[2] = "-";
		} else {
			row[1] = "" + extension.n1;
			row[2] = "" + extension.n2;
		}
		if (trimRegion == null) {
			// no trim choosen.
			row[3] = "-";
			row[4] = "-";
		} else {
			row[4] = "" + trimRegion.n1;
			row[5] = "" + trimRegion.n2;
		}
		fragmentData.add(row);
	}

	@Override
	protected Report createReport(List<ReportElement> elements) throws InterruptedException {
		return new SimpleReport(elements, null);
	}

	@Override
	protected List<ReportElement> createReportElements() throws InterruptedException {
		final ReportCompositeElement rce1 = new ReportCompositeElement();
		rce1.setCaption("Genome information");
		final ReportTableElement tabel1 = asModel(fragmentData, "Name", "Extension left", "Extension right",
				"Primer left", "Primer right", "Problematic regions", "Final length");
		rce1.addReportElement(tabel1);
		elements.add(rce1);

		final ReportCompositeElement rce2 = new ReportCompositeElement();
		rce2.setCaption("Coverage information");
		final ReportTableElement tabel2 = asModel(coverageData, "Name", "Min", "Max", "Mean \u00B1 StdDev",
				"Low Coverage Regions");
		rce2.addReportElement(tabel2);
		elements.add(rce2);

		return elements;
	}

	private static ReportTableElement asModel(List<String[]> data, String... columns) {
		final ClcReportTableModel model = new ClcReportTableModel(columns, data.toArray(new String[data.size()][]));
		return new ReportTableElement(model);
	}

}
