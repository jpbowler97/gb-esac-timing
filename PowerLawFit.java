package gb.esac.timing;


import gb.esac.io.AsciiDataFileFormatException;
import gb.esac.io.AsciiDataFileReader;
import gb.esac.tools.LeastSquaresFitter;
import hep.aida.IAnalysisFactory;
import hep.aida.IDataPointSet;
import hep.aida.IDataPointSetFactory;
import hep.aida.IFitData;
import hep.aida.IFitFactory;
import hep.aida.IFitResult;
import hep.aida.IFitter;
import hep.aida.IFunction;
import hep.aida.IFunctionFactory;
import hep.aida.IPlotter;
import hep.aida.ITree;
import java.io.IOException;
import java.text.DecimalFormat;
import org.apache.log4j.Logger;


public class PowerLawFit {

    static Logger logger  = Logger.getLogger(PowerLawFit.class);

    private static double[] x, y, yErr;
    private static double norm, index;
    private static boolean logSpace = false;
    private static boolean isXspecFile = false;
    private static String filename;

    private static DecimalFormat number = new DecimalFormat("0.0000");
    private static DecimalFormat sci = new DecimalFormat("0.0000E00");

    public static void main (String[] args) throws Exception  {

	readArguments(args);
	readDataFile();
	calculateLS();
	displayFitResults();
    }


    private static void readArguments(String[] args) {

	if ( args.length > 3 || args.length < 1 ) {
	    logger.info("Usage: java PowerLawFit dataFilename (xspec) (log)");
	    System.exit(-1);
	}
	else if ( args.length == 1 ) {
	    filename = args[0];
	}
	else if ( args.length == 2 ) {
	    filename = args[0];
	    if ( args[1].equals("xspec") ) {
		isXspecFile = true;
	    }
	    else if ( args[1].equals("log") ) {
		logSpace = true;
	    }
	    else {
		logger.info("Usage: java PowerLawFit dataFilename (xspec) (log)");
		System.exit(-1);
	    }
	}
	else {
	    filename = args[0];
	    if ( args[1].equals("xspec") && args[2].equals("log") ) {
		isXspecFile = true;
		logSpace = true;
	    }
	    else {
		logger.info("Usage: java PowerLawFit dataFilename (xspec) (log)");
		System.exit(-1);
	    }
	}
    }

    private static void readDataFile() throws AsciiDataFileFormatException, IOException  {

	//  Get data from the data file: Assume col 0= frequency, and col 1=power
	AsciiDataFileReader dataFile = new AsciiDataFileReader(filename); 
	if ( isXspecFile ) {
	    double[] energy = dataFile.getDblCol(0);
	    //double[] energyErr = dataFile.getDblCol(1);  do not need this for the fit
	    double[] flux = dataFile.getDblCol(2);
	    double[] fluxErr = dataFile.getDblCol(3);
	    x = energy;
	    y = flux;
	    yErr = fluxErr;
	}
	else {
	    double[] freq = dataFile.getDblCol(0);
	    double[] pow = dataFile.getDblCol(1);
	    int n = freq.length;
	    x = freq;
	    y = pow;
	    yErr = new double[n];
	    for ( int i=0; i < y.length; i++ ) {
		yErr[i] = 1.0E-4;
	    }
	}
    }

    private static void calculateLS() {

	double  err_index, err_norm;
	if ( logSpace ) {

	    //  Transform to log space if requested
	    if ( logSpace ) {
		for ( int i=0; i < x.length; i++ ) {
		    x[i] = Math.log10(x[i]);
		    //yErr[i] = Math.abs(yErr[i]/y[i]*Math.log10(y[i]));
		    y[i] = Math.log10(y[i]);
		}
	    }
	    double[] result = LeastSquaresFitter.leastSquaresFitLine(x, y);
	    index = result[0];
	    err_index = result[1];
	    norm = result[2];
	    err_norm = result[3];
	    logger.info("Fitted function is: f(x) = m*x + b");
	    logger.info("Slope (m) = "+number.format(index)+" +/- "+number.format(err_index));
	    logger.info("Intercet (b) = "+number.format(norm)+" +/- "+number.format(err_norm));

	}
	else {
	    double[] result = LeastSquaresFitter.leastSquaresFitPowerLaw(x, y);
	    index=result[0];
	    norm=result[1];
	    logger.info("Fitted function is f(x) = N*x^b");
	    logger.info("Power (b) = "+number.format(index));
	    logger.info("Norm (N) = "+sci.format(norm));
	}


    }

    private static void displayFitResults() {

	//  Construct factories
	IAnalysisFactory af = IAnalysisFactory.create();
	ITree tree = af.createTreeFactory().create();
	IDataPointSetFactory dpsf = af.createDataPointSetFactory(tree);
 	IFunctionFactory funcF  = af.createFunctionFactory(tree);
  	IFitFactory fitF   = af.createFitFactory();
  	IFitter fitter = fitF.createFitter("BinnedMaximumLikelihood", "jminuit", "noClone=true");
	if ( logSpace )  fitter.setFitMethod("LS");
	else fitter.setFitMethod("BinnedMaximumLikelihood");

	//  Create 2D IDataPointSet
 	IDataPointSet dps = dpsf.create("dps", "Power Spectrum", 2);
	for ( int i = 0; i < x.length; i++ ) {
 	    dps.addPoint();
 	    dps.point(i).coordinate(0).setValue( x[i] );
 	    dps.point(i).coordinate(1).setValue( y[i] );
 	    //yErr[i] = 0.00000001;
 	    dps.point(i).coordinate(1).setErrorPlus( yErr[i] );
 	    dps.point(i).coordinate(1).setErrorMinus( yErr[i] );
	}

	//  Construct the power law function
	IFunction f = null;
	if ( logSpace ) {   //  Straigt line function
	    f = funcF.createFunctionFromScript("func", 1, "norm+index*x[0]", "norm,index","norm+index*x");
	    f.setParameter("norm", norm);
	    f.setParameter("index", index);
	}
	else {   //  Power law function
	    f = funcF.createFunctionFromScript("func", 1, "norm*pow(x[0],index)", "norm,index","norm*pow(x[0], index)"); 
	    f.setParameter("norm", norm);
	    f.setParameter("index", index);
	}

	//  Perform the fit
	IFitData data = fitF.createFitData();
	data.create1DConnection(dps, 0, 1);
 	IFitResult fittedFunction = fitter.fit(data, f);
	
	//  Display the results
	IPlotter plotter = af.createPlotterFactory().create("Plot Power Spectrum");
	plotter.region(0).plot( dps );
	plotter.region(0).plot( fittedFunction.fittedFunction() );
	plotter.show();
	plotter.region(0).style().statisticsBoxStyle().setVisible(true);
    }

}
