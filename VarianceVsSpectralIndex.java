package gb.esac.timing;

import java.util.Arrays;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;

import gb.esac.tools.Stats;
import gb.esac.tools.Converter;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;


public class VarianceVsSpectralIndex {

    static Logger logger  = Logger.getLogger(VarianceVsSpectralIndex.class);


    public static void main (String[] args) throws IOException, TimingException  {

	PropertyConfigurator.configure("/Users/gbelanger/javaProgs/gb/esac/logger.config");

	DecimalFormat num = new DecimalFormat("0.000000");
	DecimalFormat display = new DecimalFormat("0.000");
	DecimalFormat sci = new DecimalFormat("0.0000E0");

	
	/**  Handle arguments  **/
	double duration = 10000;
	double alphaMin = 0;
	double alphaMax = 4;
	int nSims = 100;
	if ( args.length == 4 ) {
	    duration = (Double.valueOf(args[0])).doubleValue();
	    alphaMin = (Double.valueOf(args[1])).doubleValue();
	    alphaMax = (Double.valueOf(args[2])).doubleValue();
	    nSims = (Integer.valueOf(args[3])).intValue();
	}
	else {
	    System.out.println("Usage: java VarianceVsSpectralIndex duration alphaMin alphaMax nSims");
	    System.exit(-1);
	}
	logger.info("Running VarianceVsSpectralIndex");
	logger.info("Duration = "+duration);
	logger.info("Alpha min = "+alphaMin);
	logger.info("Alpha max = "+alphaMax);
	logger.info("Number of simulated event lists = "+nSims);

	/**  Define mean rates (same as in SpectralIndexVsCountRate.java  **/
	double[] countRates = new double[]{0.25, 1, 5, 25, 100, 250};
 	//double[] countRates = new double[]{0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250};
	logger.info("Program will loop over these mean count rates");
	for ( int i=0; i < countRates.length; i++ ) {
	    logger.info("  "+countRates[i]+" cts/s");
	}

	/**  Define step size in alpha  **/
	double stepSize = 0.2;
	double alphaRange = alphaMax - alphaMin;
	int nAlphaBins = (new Double(Math.round(alphaRange/stepSize))).intValue();
	stepSize = alphaRange/nAlphaBins;
	int nAlphas = nAlphaBins + 1;
	double[] alphas = new double[nAlphas];
	//logger.info("Spectral index values are:");
	for ( int i=0; i < nAlphas; i++ ) {
	    alphas[i] = alphaMin + i*stepSize;
	    //logger.info("   "+alphas[i]);
	}

	/**  Define optimal number of bins  **/
	int nbins = 128;

	/**  PrintWriter  **/
	PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("varVsAlpha.qdp")));
	String[] header = new String[] {
	    "DEV /XS",
	    "READ SERR 2 3 4 5 6 7 8 9 10 11", 
	    "TIME OFF", "LINE ON",
	    "LOG Y ON",
	    "LW 4", "CS 1.3", "LAB T", "LAB F",
	    "LAB X Spectral Index",
	    "LAB Y Light Curve Variance (cts s\\u-1\\d)\\u2",	    
	    "VIEW 0.1 0.2 0.9 0.8",
	    "!"
	};
	for ( int k=0; k < header.length; k++ ) {
	    pw.println(header[k]);
	}
	
	
	/**  Start loop on values of alpha  **/
	logger.info("Starting loop on values of alpha");
	for ( int a=0; a < nAlphas; a++ ) {

	    double[] vars = new double[nSims];

	    /**  define alpha for this loop  **/
	    double alpha = alphas[a];
	    pw.print(alpha+"\t");
	    System.out.print("Log  :   Alpha = "+display.format(alpha)+" ... ");

	    /**  Start loop on values of count rate  **/
	    //logger.info("    Starting loop on count rates");
	    for ( int r=0; r < countRates.length; r++ ) {
		
		double meanRate = countRates[r];

		/**  Start loop on nSims  **/
		for ( int i=0; i < nSims; i++ ) {
		    
		    double[] t = ArrivalTimes.generateRedArrivalTimes(meanRate, duration, alpha);
		    double[] lc = (double[]) TimingUtils.makeLC(t, nbins)[1];
		    vars[i] = Stats.getVariance(lc);
		}

		/**  Get the average variance  **/
		double[] aveAndVar = Stats.getRunningAveAndVar(vars);
		double aveVar = aveAndVar[0];
		double varStdDev = Math.sqrt(aveAndVar[1]);
 		//System.out.print("Log  :       Rate "+(r+1)+" = "+meanRate+" cts/s \t");		
 		//System.out.println("Var = "+display.format(aveVar)+" +/- "+display.format(varStdDev));

		/**  Print to file  **/
		pw.print(num.format(aveVar)+"\t"+num.format(varStdDev)+"\t");

	    }
	    System.out.println("done");
	    pw.println();
	    pw.flush();
	}
	pw.println("!");
	pw.println("HARD varVsAlpha.ps/ps");
	pw.close();
	
    }
}
