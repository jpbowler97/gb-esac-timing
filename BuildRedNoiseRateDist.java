

import gb.esac.io.AsciiDataFileWriter;
import gb.esac.montecarlo.TimmerKonig;
import gb.esac.tools.MinMax;
import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.ITree;
import java.text.DecimalFormat;

public class BuildRedNoiseRateDist { 

    private static DecimalFormat twoDigits = new DecimalFormat("#.00");

    public static void main(String[] args)  throws Exception {

	//  Create AIDA factories
	IAnalysisFactory af = IAnalysisFactory.create();
	ITree tree = af.createTreeFactory().create();
	IHistogramFactory hf = af.createHistogramFactory(tree);
	
	double alpha = 2.0;
	double duration = 1e3;
	if ( args.length == 2 ) {
	    alpha = (Double.valueOf(args[0])).doubleValue();
	    duration = (Double.valueOf(args[1])).doubleValue();
	}

	int nBins = (int) Math.pow(2,14);
	double meanRate = 1d;
	double[] rates = TimmerKonig.getTimmerRates(meanRate, duration, alpha, nBins);
	double max = MinMax.getMax(rates);
	int i=0;
	while ( i < 9 ) {
	    rates = TimmerKonig.getTimmerRates(meanRate, duration, alpha, nBins);
	    max = Math.max(max, MinMax.getMax(rates));
	    i++;
	}
	double histoMax = 1.5*max;
	IHistogram1D histoOfRates = hf.createHistogram1D("histoOfRates", nBins, 0, histoMax);
	i=0;
	while ( i < 100 ) {
	    rates = TimmerKonig.getTimmerRates(meanRate, duration, alpha, nBins);
	    for ( int j=0; j < rates.length; j++ ) {
		histoOfRates.fill(rates[j]);
	    }
	    i++;
	}
	AsciiDataFileWriter out = new AsciiDataFileWriter("histoOfRedNoiseRates_index_"+alpha+"_"+duration+".qdp");
	out.writeHisto(histoOfRates, "Poisson Rate Parameter (cps)");

    }
}
