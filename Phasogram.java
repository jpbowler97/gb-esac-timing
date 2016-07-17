package gb.esac.timing;

import gb.esac.eventlist.EventList;
import gb.esac.tools.BasicStats;
import gb.esac.tools.MinMax;
import gb.esac.tools.FitsUtils;
import gb.esac.tools.DataUtils;
import gb.esac.binner.Binner;
import gb.esac.io.AsciiDataFileReader;
import gb.esac.io.AsciiDataFileWriter;
import gb.esac.aida.functions.SineFunction;
import gb.esac.timeseries.*;

import hep.aida.*;
import hep.aida.ref.histogram.*;

import java.util.Arrays;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.text.DecimalFormat;
import nom.tam.fits.*;
import cern.jet.stat.Descriptive;
import cern.colt.list.DoubleArrayList;

import org.apache.log4j.Logger;


/**
 *
 * @author <a href="mailto:guilaume.belanger@esa.int">Guillaume Belanger</a>
 * @version (Feb 2015, ESAC)
 *
 *                Previous: Dec 2007
 */

public class Phasogram {

    private static Logger logger  = Logger.getLogger(Phasogram.class);

    public static void main(String[] args) throws FitsException, IOException, Exception {

	DecimalFormat numberFormat = new DecimalFormat("0.000");

	//  Handle arguments
	String evlistName = null;
	double period = 0;
	double[] periods = new double[1];
	int nPhaseBins = 15;
	boolean fileIsFits = true;
	if ( args.length > 3 ) {
	    periods = new double[args.length-2];
	    evlistName = args[0];
	    for ( int i=1; i < args.length-1; i++ ) {
		periods[i] = (Double.valueOf(args[i])).doubleValue();
	    }
	    nPhaseBins = (Integer.valueOf(args[args.length-1])).intValue();
	}
	else if ( args.length == 3 ) {
	    evlistName = args[0];
	    period = (Double.valueOf(args[1])).doubleValue();
	    nPhaseBins = (Integer.valueOf(args[2])).intValue();
	}
	else if ( args.length == 2 ) {
	    evlistName = args[0];
	    period = (Double.valueOf(args[1])).doubleValue();
	}
	else {
	    logger.error("Usage: Phasogram eventlistName period (nPhaseBins [15])");
	    System.exit(0);
	}


	// Construct Phasogram
	EventList evlist = new EventList(evlistName);
	double[] phases = TimingUtils.getPhases(evlist.getArrivalTimes(), period);
	Arrays.sort(phases);
	double binWidth = 1d/nPhaseBins;
	TimeSeries ts = TimeSeriesMaker.makeTimeSeries(phases, binWidth, 0, 1);
	ts.writeCountsAsQDP("phaso.qdp");




	
	/**  Construct the Sine function to fit the phasogram  **/
// 	IFunction sine = new SineFunction("Sine Function");
// 	double amplitude = maxPhaseHeight;
// 	double xOffset = -0.25;
// 	double yOffset = meanCount;
// 	sine.setParameter("period", 1);
// 	sine.setParameter("amplitude", amplitude);
// 	sine.setParameter("xOffset", xOffset);
// 	sine.setParameter("yOffset", yOffset);


	/**  Set up JAIDA factories  **/
// 	IAnalysisFactory af = IAnalysisFactory.create();
//  	IFitFactory fitF   = af.createFitFactory();
// 	IFitter fitter = fitF.createFitter("Chi2", "jminuit");

	/**  Construct a histogram of the phases from 0 to 1  **/
// 	double[] heights = new double[2*nPhaseBins+2];
// 	double[] errs = new double[2*nPhaseBins+2];
// 	double[] edges = new double[2*nPhaseBins+1];
// 	double width = 1.0/nPhaseBins;
// 	for ( int i=1; i <= 2*nPhaseBins; i++ ) {
// 	    edges[i-1] = width*(i-1);
// 	    heights[i] = orderedPhaseBinHeights[i-1];
// 	    errs[i] = orderedPhaseBinHeightsErr[i-1];
// 	}
// 	edges[2*nPhaseBins] = 2.0;
// 	VariableAxis xaxis = new VariableAxis(edges);
// 	Histogram1D phaso = new Histogram1D("Phaso", "Phaso", xaxis);
// 	phaso.setContents(heights, errs, null, null, null);


	/**  Fit  **/
// 	IFitResult fitResult = fitter.fit(phaso, sine);
// 	double fittedPeriod = fitResult.fittedParameter("period");
// 	double fittedAmp = fitResult.fittedParameter("amplitude");
// 	logger.info("Chi2 = "+fitResult.quality());
// 	logger.info("DOF = "+fitResult.ndf());

	/**  Display the results  **/
// 	IPlotter plotter = af.createPlotterFactory().create("Phasogram");
// 	plotter.createRegion();
// 	plotter.region(0).plot(phaso);
// 	plotter.region(0).plot(fitResult.fittedFunction() );
// 	plotter.region(0).style().statisticsBoxStyle().setVisible(true);	

// 	plotter.show();
	
// 	/**  Display light curve and phasogram  **/
// 	Runtime rt = Runtime.getRuntime();
// 	rt.exec("qdp "+lcFilename);
// 	rt.exec("qdp "+phasoFilename);

    }
}
