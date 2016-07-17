package gb.esac.timing;

import java.util.Arrays;
import java.io.*;
import java.text.DecimalFormat;

import hep.aida.ref.histogram.Histogram1D;

import gb.esac.tools.Stats;
import gb.esac.tools.Analysis;
import gb.esac.tools.Binner;
import gb.esac.tools.Complex;
import gb.esac.tools.FFT;
import gb.esac.tools.Converter;
import gb.esac.io.DataFile;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;


public class SpectralIndexVsCountRate {


    static Logger logger  = Logger.getLogger(SpectralIndexVsCountRate.class);

    
    public static void main (String[] args) throws IOException, TimingException {


	DecimalFormat num = new DecimalFormat("0.000000");
	DecimalFormat sci = new DecimalFormat("0.0000E0");

	
	/**  Handle arguments  **/
	double index = 1;
	double duration = 10000;
	int nspecs = 100;
	if ( args.length == 3 ) {
	    index = (Double.valueOf(args[0])).doubleValue();
	    duration = (Double.valueOf(args[1])).doubleValue();
	    nspecs = (Integer.valueOf(args[2])).intValue();
	}
	else {
	    System.out.println("Usage: java SpectralIndexVsCountRate index duration nSpecs");
	    System.exit(-1);
	}
	logger.info("Running CountRateEffect");
	logger.info("Spectral index = "+index);
	logger.info("Duration = "+duration);
	logger.info("Number of spectra = "+nspecs);


	/**  Define count rates and bin times  **/
 	double[] countRates = new double[]{0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250};
 	double[] binTimes = new double[]{9.7, 19.5, 39.0, 78.1, 156.2, 312.4};
// 	double[] countRates = new double[]{0.1, 5, 10};
// 	double[] binTimes = new double[]{19.5, 78.1, 312.4};


	/**  Write the header to the output file  **/
	PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("alphaVsCounts.qdp")));
	String[] header = new String[] {
	    "! QDP file",
	    "READ SERR 2 3 4 5 6 7",
	    "DEV /XS", "VIEW 0.1 0.2 0.9 0.8",
	    "CS 1.3", "LW 4", "LAB t", "LAB f",
	    "TIME OFF",
	    "LAB X \"Count rate (cts/s)\"",
	    "LAB Y \"Estimated Spectral Index\"",
	    "LOG ON", "LINE ON", "MA ON",
	    "!!!!!!!!!!  DEFINE THE COLOURS",
	    "CO 2 ON 2",
	    "CO 3 ON 3",
	    "CO 4 ON 4",
	    "CO 5 ON 5",
	    "CO 6 ON 6",
	    "CO 8 ON 7",
	    "!!!!!!!!!!  DEFINE THE MARKERS",
	    "MA 2 ON 2",
	    "MA 7 ON 3",
	    "MA 4 ON 4",
	    "MA 5 ON 5",
	    "MA 6 ON 6",
	    "MA 12 ON 7",
	    "MA SIZE 2",
	    "!!!!!!!!!!  DEFINE THE LABELS",
	    "LAB 1 VPOS 0.79 0.55 \"bintime\" CS 1.2 JUST RIGHT",
	    "LAB 2 VPOS 0.68 0.30 \"10 s\" CO 2 MA 2 MSIZE 1.5 CS 1.5 JUST LEFT",
	    "LAB 3 VPOS 0.68 0.34 \"20 s\" CO 3 MA 7 MSIZE 1.5 CS 1.5 JUST LEFT",
	    "LAB 4 VPOS 0.68 0.38 \"40 s\" CO 4 MA 4 MSIZE 1.5 CS 1.5 JUST LEFT",
	    "LAB 5 VPOS 0.68 0.42 \"80 s\" CO 5 MA 5 MSIZE 1.5 CS 1.5 JUST LEFT",
	    "LAB 6 VPOS 0.68 0.46 \"160 s\" CO 6 MA 6 MSIZE 1.5 CS 1.5 JUST LEFT",
	    "LAB 7 VPOS 0.68 0.50 \"320 s\" CO 8 MA 12 MSIZE 1.5 CS 1.5 JUST LEFT",
	    "!!!!!!!!!!  DEFINE THE DATA",
	    "! bin times are:  (data)",
	    "! 10 = 9.730 s",
	    "! 20 = 19.529 s",
	    "! 40 = 39.058 s",
	    "! 80 = 78.115 s",
	    "! 160 = 156.230 s",
	    "! 320 = 312.460 s",
	    "!"
	};

	for ( int k=0; k < header.length; k++ ) {
	    pw.println(header[k]);
	}


	/**  Loop on count rates  **/
	double meanRate = 0;
	for ( int c=0; c < countRates.length; c++ ) {

	    meanRate = countRates[c];
	    pw.print(num.format(meanRate)+"\t");
	    logger.info("Count rate = "+meanRate);

	    /**  Loop on bin times  **/
	    double binTime = 0;
	    for ( int t=0; t < binTimes.length; t++ ) {

		/**  Define binTime for this loop  **/
		binTime = binTimes[t];
		logger.info(" Bin time = "+binTime);
		
		/**  Determine power of 2 number of bins  **/
		double n = Math.round(Math.log10(duration/binTime)/Math.log10(2));
		int nTimeBins = (new Double(Math.pow(2, n))).intValue();
		double newBinTime = duration/nTimeBins;
		double nuMin = 1/duration;
		double nuMax = 1/(2*newBinTime);
		double[] testFreqs = TimingUtils.getFFTFrequencies(nTimeBins, duration);
		double[] x2 = Converter.lin2logSpace(testFreqs);
		int nFreqs = testFreqs.length;
		double[] sumOfPows = new double[nFreqs];
		logger.info("  nBins = "+nTimeBins);
		logger.info("  New bin time for event list powspec = "+num.format(newBinTime)+" s");
		logger.info("  Nu min (1/duration) = "+sci.format(nuMin)+" Hz");
		logger.info("  Nu max (1/2*newBinTime) = "+sci.format(nuMax)+" Hz");
		logger.info("  "+nFreqs+" frequencies will be tested");

		/**  Generate powspecs and get LS best fit values for the slopes  **/
		logger.info("  Simulating "+nspecs+" power spectra ...");

		//int nbins = 4096;
		//double[] testFreqs = TimingUtils.getFFTFrequencies(nbins, duration);
		//int nFreqs = testFreqs.length;
		//double[] x1  = Converter.lin2logSpace(testFreqs);


		//  Variables for Timmer-Konig powspecs
		//double powBinTime = duration/nbins;
		//Complex[] fourierComp = null;
		//double[] pow1 = null;
		//double[] y1  = new double[nbins];  // pow in log space
		double[] lsFitResults = new double[4];
		//double[] slopes = new double[nspecs];

		//  Variables for evlist powspecs
		//Complex[] ifft = null;
		//double[] lc = null;
		//Histogram1D lcHisto = null;
		//Histogram1D cdfHisto = null;
		//double[] times = null;
		//double[][] psd = null;
		//double[] freqs = new double[nTimeBins/2];
		//double[] pow2 = new double[nTimeBins/2];
		//double[] x2 = new double[nTimeBins/2];
		//double[] y2 = new double[nTimeBins/2];
		double[] evlistSlopes = new double[nspecs];
		double tzero = 0;

		for ( int i=0; i < nspecs; i++ ) {

// 		    /**  Construct Timmer-Koenig powspec  **/
// 		    fourierComp = TimmerKonig.getTimmerFourierComp(index, duration, nbins);
// 		    pow1 = TimmerKonig.makeTimmerPowspec(fourierComp);
// 		    y1[j] = Converter.lin2logSpace(pow1);

// 		    /**  Fit power spectrum  **/
// 		    lsFitResults = Stats.leastSquaresFitLine(x1, y1);
// 		    /** function is: pow = a + b*testFreqs. Method returns [a][b][err_a][err_b]  **/
// 		    slopes[i] = -lsFitResults[1];

// 		    /**  Inverse transform the fourierComp to get the light curve  **/
// 		    lc = TimmerKonig.makeTimmerLC(fourierComp);

// 		    /**  construct CDF from the light curve  **/
// 		    lcHisto = Converter.array2histo("light curve", tzero, powBinTime, lc);
// 		    cdfHisto = Stats.getCDFHisto(lcHisto);
		    
// 		    /**  Construct event list by drawing from the CDF  **/
// 		    int nevents = (new Double(meanRate*duration)).intValue();
// 		    times = Analysis.getRandom(cdfHisto, nevents);
// 		    Arrays.sort(times);

		    double[] times = ArrivalTimes.generateRedArrivalTimes(meanRate, duration, index);

		    /**  Make PSD from event list  **/
		    double[][] psd = Periodograms.makePSD_FFT(times, nTimeBins, "leahy");
		    double[] pow2 = psd[1];
		    for ( int f=0; f < nFreqs; f++ ) {
			
			sumOfPows[f] += pow2[f];
		    }

		    /**  Fit  **/
		    double[] y2 = Converter.lin2logSpace(pow2);
		    lsFitResults = Stats.leastSquaresFitLine(x2, y2);
		    /** function is: pow = a + b*testFreqs, method returns [a][b][err_a][err_b]  **/
		    evlistSlopes[i] = -lsFitResults[1];
		}

	
		/**  Bin the results **/
		int nHistoBins = 40;
// 		double minSlope1 = Stats.getMin(slopes);
// 		double maxSlope1 = Stats.getMax(slopes);
// 		double slopeBinWidth1 = (maxSlope1 - minSlope1)/nHistoBins;
// 		double[] histoOfSlopes1 = Binner.binData(slopes, minSlope1, maxSlope1, nHistoBins);
		
		//double minSlope2 = Stats.getMin(evlistSlopes);
		//double maxSlope2 = Stats.getMax(evlistSlopes);
		//double slopeBinWidth2 = (maxSlope2 - minSlope2)/nHistoBins;
		//double[] histoOfSlopes2 = Binner.binData(evlistSlopes, minSlope2, maxSlope2, nHistoBins);
		

		/**  Construct the axes for histos and Gaussian function  **/
// 		double[] slopeAveAndVar1 = Stats.getRunningAveAndVar(slopes);
// 		double ave1 = slopeAveAndVar1[0];
// 		double var1 = slopeAveAndVar1[1];
// 		double norm1 = nspecs*slopeBinWidth1;
// 		double[] slopeAxis1 = new double[nHistoBins];
		double[] slopeAveAndVar2 = Stats.getRunningAveAndVar(evlistSlopes);
		//double ave2 = slopeAveAndVar2[0];
		double var2 = slopeAveAndVar2[1];
		//double norm2 = nspecs*slopeBinWidth2;
		//double[] slopeAxis2 = new double[nHistoBins];
// 		double[] gauss1 = new double[nHistoBins];
 		//double[] gauss2 = new double[nHistoBins];
		//for ( int i=0; i < nHistoBins; i++ ) {
// 		    slopeAxis1[i] = minSlope1 + (0.5 + i)*slopeBinWidth1;
// 		    gauss1[i] = norm1*Math.sqrt(1/(2*Math.PI*var1)) * 
// 			Math.exp( -Math.pow((slopeAxis1[i]-ave1), 2)/(2*var1) );
		//slopeAxis2[i] = minSlope2 + (0.5 + i)*slopeBinWidth2;
		//gauss2[i] = norm2*Math.sqrt(1/(2*Math.PI*var2)) * 
		//Math.exp( -Math.pow((slopeAxis2[i]-ave2), 2)/(2*var2) );
		//}


		/**  Fit the average power spectrum to get the mean index  **/
		double[] avePow = new double[nFreqs];
		for ( int f=0; f < nFreqs; f++ ) {
		    avePow[f] = sumOfPows[f]/nspecs;
		}
		double[] y2 = Converter.lin2logSpace(avePow);
		lsFitResults = Stats.leastSquaresFitLine(x2, y2);
		/** function is: pow = a + b*testFreqs, method returns [a][b][err_a][err_b]  **/
		double aveSlope = -lsFitResults[1];
		double aveSlopeErr= lsFitResults[3];

		/**  Write the mean of the  slopes to the output file  **/
		//pw.print(num.format(ave2)+"\t"+num.format(Math.sqrt(var2/nspecs))+"\t");
		pw.print(num.format(aveSlope)+"\t"+num.format(aveSlopeErr)+"\t");
		pw.flush();

		/**  Compute Pearson's Chi2 test  for the Gaussians **/
// 		double chi2_1 = Stats.computeChi2Test(histoOfSlopes1, gauss1);
// 		logger.info("Mean of Timmer-Koenig PSD slopes = "+
// 				   num.format(ave1)+" +/- "+ num.format(Math.sqrt(var1))
// 				   +", chi2 for Gaussian = "+num.format(chi2_1));
		
// 		double chi2_2 = Stats.computeChi2Test(histoOfSlopes2, gauss2);
// 		logger.info("    Mean slope = "+
// 				   num.format(ave2)+" +/- "+ num.format(Math.sqrt(var2))
// 				   +", chi2 for Gaussian = "+num.format(chi2_2));

		
// 		/**  Write the histograms  as QDP files  **/
// 		String mean = num.format(ave1);
// 		String stdDev = num.format(Math.sqrt(var1));
// 		String chi2 = num.format(chi2_1);
// 		String dof = Integer.toString(nHistoBins-1);
// 		String binN = Integer.toString(nspecs);
// 		String[] header = new String[] {
// 		    "DEV /XS",
// 		    "READ 1 2 3", "TIME OFF",
// 		    "LINE STEP ON 2", "LINE ON 3",
// 		    "LW 3", "CS 1.3", "LAB T", "LAB F",
// 		    "LAB X Spectral Index",
// 		    "LAB Y Entries Per Bin",
// 		    "LAB 3 VPOS 0.175 0.7 \"\\gm = "+mean+"\" JUST LEFT",
// 		    "LAB 4 VPOS 0.175 0.67 \"\\gs = "+stdDev+"\" JUST LEFT",
// 		    "LAB 5 VPOS 0.175 0.64 \"N = "+binN+"\" JUST LEFT",
// 		    "LAB 6 VPOS 0.825 0.7 \"\\gx\\u2\\d("+dof+") = "+chi2+"\" JUST RIGHT",
// 		    "VIEW 0.1 0.2 0.9 0.8",
// 		    "!",
// 		    "Mean = "+ave1,
// 		    "Std dev = "+Math.sqrt(var1),
// 		    "Pearson's Chi2 = "+chi2_1,
// 		    "!"
// 		};
// 		DataFile histoOfSlopesFile = new DataFile("histoOfSlopes.qdp");
// 		histoOfSlopesFile.writeData(header, slopeAxis1, histoOfSlopes1, gauss1);
		
		
// 		mean = num.format(ave2);
// 		stdDev = num.format(Math.sqrt(var2));
// 		chi2 = num.format(chi2_2);
// 		binN = Integer.toString(nspecs);
// 		header = new String[] {
// 		    "DEV /XS",
// 		    "READ 1 2 3", "TIME OFF",
// 		    "LINE STEP ON 2", "LINE ON 3",
// 		    "LW 3", "CS 1.3", "LAB T", "LAB F",
// 		    "LAB X Fit Spectral Index",
// 		    "LAB Y Entries Per Bin",	    
// 		    "LAB 3 VPOS 0.175 0.7 \"\\gm = "+mean+"\" JUST LEFT",
// 		    "LAB 4 VPOS 0.175 0.67 \"\\gs = "+stdDev+"\" JUST LEFT",
// 		    "LAB 5 VPOS 0.175 0.64 \"N = "+binN+"\" JUST LEFT",
// 		    "LAB 6 VPOS 0.825 0.7 \"\\gx\\u2\\d("+dof+") = "+chi2+"\" JUST RIGHT",
// 		    "VIEW 0.1 0.2 0.9 0.8",
// 		    "!",
// 		    "Mean = "+ave2,
// 		    "Std dev = "+Math.sqrt(var2),
// 		    "Pearson's Chi2 = "+chi2_2,
// 		    "!"
// 		};
// 		DataFile histoOfErrorsFile = new DataFile("histoOfEvlistSlopes.qdp");
// 		histoOfErrorsFile.writeData(header, slopeAxis2, histoOfSlopes2, gauss2);


	    }
	    pw.println();
	    pw.flush();
	}
	pw.close();

    }
}
