package io.github.pdekker.viraltyping.algo.transferannot;

import java.util.List;

import com.clcbio.api.base.util.CreateList;
import com.clcbio.api.base.util.Null;
import com.clcbio.api.clc.algorithms.report.AbstractReportCalculator;
import com.clcbio.api.clc.datatypes.report.ClcReportTableModel;
import com.clcbio.api.clc.datatypes.report.ReportCompositeElement;
import com.clcbio.api.clc.datatypes.report.ReportFrontPageStringElement;
import com.clcbio.api.clc.datatypes.report.ReportTableElement;
import com.clcbio.api.clc.datatypes.report.ReportTextElement;
import com.clcbio.api.clc.datatypes.report.SimpleReport;
import com.clcbio.api.free.datatypes.report.Report;
import com.clcbio.api.free.datatypes.report.ReportElement;

public class TransferAnnotationReportBuilder extends AbstractReportCalculator {

	public final List<ReportElement> elements = CreateList.of(new ReportFrontPageStringElement());

	@Override
	protected Report createReport(List<ReportElement> elements) throws InterruptedException {
		return new SimpleReport(elements, null);
	}

	public void addCds(List<String[]> reportRows) {
		final ReportCompositeElement rce = new ReportCompositeElement();
		rce.setCaption("CDS Information");
		if (Null.isEmpty(reportRows)) {
			rce.addReportElement(new ReportTextElement("No CDS annotation found"));
		} else {
			final String[] header = new String[] { "Fragment", "Protein", "% expected CDS length", "Result" };
			final ClcReportTableModel model = new ClcReportTableModel(header,
					reportRows.toArray(new String[reportRows.size()][]));
			rce.addReportElement(new ReportTableElement(model));
		}
		elements.add(rce);
	}

	@Override
	protected List<ReportElement> createReportElements() throws InterruptedException {
		return elements;
	}

}
