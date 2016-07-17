package gb.esac.timing;

import gb.esac.io.DataFile;
import gb.esac.tools.*;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.BufferedWriter;
import java.io.FileWriter;

import java.text.DecimalFormat;
import java.util.Vector;
import java.util.Arrays;

import nom.tam.fits.*;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;



public class EstimateSpectralIndex {

    static Logger logger  = Logger.getLogger(GenerateEventList.class);

    public static DecimalFormat num = new DecimalFormat("0.00000");
    public static DecimalFormat freq = new DecimalFormat("0.000E00");


    public static void main(String[] args) throws IOException, Exception  {


	PropertyConfigurator.configure("/Users/gbelanger/javaProgs/gb/esac/logger.config");


	String inputFilename = null;
	int nSims = 50;
    
	/**  Read input filename  **/
	if ( args.length == 2 ) {
	    inputFilename = args[0];
	    nSims = (Integer.valueOf(args[1])).intValue();
	}
	else {
	    logger.error("Usage: java EstimateSpectralIndex filename nSims");
	    System.exit(-1);
	}

	/**  Get data  **/
	Object[] lcData = TimingUtils.readLCFile(inputFilename);

	int ncols = ((Integer) lcData[0]).intValue();
	double duration = ((Double) lcData[1]).doubleValue();
	double oldBinTime = ((Double) lcData[2]).doubleValue();
	double[] times = (double[]) lcData[3];
	double[] binTimes = (double[]) lcData[4];
	double[] binEdges = (double[]) lcData[5];
	double[] rate = (double[]) lcData[6];
	double[] error = (double[]) lcData[7];
	double[] signif = (double[]) lcData[8];
	boolean dataTypeIsRate = ((Boolean) lcData[9]).booleanValue();
	boolean useBinEdges = ((Boolean) lcData[10]).booleanValue();


	/**  Determine and print out duration and mean count rate  **/
	System.out.println("Log  : Duration = "+num.format(duration)+" s");
	int nevents = 0;
	double meanRate = 0;
	double wMeanRate = 0;
	double meanIndivError = 0;
	double meanSignalToNoise = 0;
	double equivMeanRate = 0;
	if ( dataTypeIsRate ) {
	    meanRate = Stats.getMean(rate);
	    nevents = (new Double(meanRate*duration)).intValue();
	    //wMeanRate = Stats.getWMean(rate, error);
	    meanIndivError = Stats.getMean(error);
	    meanSignalToNoise = meanRate/meanIndivError;
	    equivMeanRate = Math.pow(meanSignalToNoise, 2)/oldBinTime;
	    System.out.println("Log  : Mean rate = "+num.format(meanRate)+" cts/s");
	    //System.out.println("Log  : Weighted mean rate = "+num.format(wMeanRate)+" cts/s");
	    System.out.println("Log  : Average number of events (from mean) = "+nevents);
	    System.out.println("Log  : Equivalent mean rate = "+num.format(equivMeanRate));
	}
	else {
	    nevents = times.length;
	    meanRate = times.length/duration;
	    equivMeanRate = meanRate;
	    System.out.println("Log  : Mean rate = "+num.format(meanRate)+" cts/s");
	    System.out.println("Log  : Number of events = "+nevents);
	}


	/**  Define max and min bintimes  **/
	double bintimeMax = duration/128;
	double min = duration/256;
	double bintimeMin = Math.max(min, oldBinTime);
	double nuMax = 1/(2*bintimeMin);
	double nuMin = 1/duration;
	int nIFS = TimingUtils.getnIFS(duration, nuMin, nuMax);
	System.out.println("Log  : Min bintime = "+num.format(bintimeMin)+" s");
	System.out.println("Log  : Nu max (1/2*bintimeMin) = "+freq.format(nuMax)+" Hz");
	System.out.println("Log  : nIFS in this range = "+nIFS);

	
	/**  Define the range of bin times between T/128 and 10s  **/
	Vector bintimesVec = new Vector();
	Vector lcBinsVec = new Vector();
	double bintime = bintimeMax;
	while ( bintime >= bintimeMin ) {
	    double nbins = (new Double(Math.ceil(duration/bintime))).intValue();
	    double n = Math.round(Math.log10(nbins)/Math.log10(2));
	    int nTimebins = (new Double(Math.pow(2, n))).intValue();
	    lcBinsVec.add(nTimebins);
	    bintime = duration/nTimebins;
	    bintimesVec.add(bintime);
	    bintime /= 2;
	}
	double[] timeBins = new double[bintimesVec.size()];
	int[] nLCbins = new int[lcBinsVec.size()];
	System.out.println("Log  : Bin times used between "+bintimeMin+" and T/64 s are: ");
	for ( int i=timeBins.length-1; i >= 0; i-- ) {
	    timeBins[i] = ((Double) bintimesVec.elementAt(i)).doubleValue();
	    nLCbins[i] = ((Integer) lcBinsVec.elementAt(i)).intValue();
	    System.out.println("Log  :   "+num.format(timeBins[i])+" s  ("+nLCbins[i]+" bins)");
	}


	/**  Write the header of output QDP file  **/
	double[] alpha = new double[]{1.5, 1.75, 2.0};
	String outName = "alphaVsBintime.qdp";
	PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outName)));
	String[] head = new String[] {
	    "! QDP File",
	    "DEV /XS",
	    "VIEW 0.1 0.2 0.9 0.8",
	    "LW 8 ON 2", "LW 4", "CS 1.3",
	    "TIME OFF", "LAB T", "LAB F",
	    "LINE ON", "LOG X ON",
	    "LAB X Bin Time (s)",
	    "LAB Y Estimated Spectral Index",
	    "!"
	};
	for ( int i=0; i < head.length; i++ ) pw.println(head[i]);
	pw.print("READ SERR ");
	for ( int i=0; i < alpha.length+1; i++ ) {
	    pw.print(+(i+2)+" ");
	}
	pw.println();
	for ( int i=0; i < alpha.length+1; i++ ) {	
	    pw.println("CO "+(i+1)+" ON "+(i+2));
	}
	int[] markers = new int[] {17, 4, 6, 7, 11};
	for ( int i=0; i < alpha.length+1; i++ ) {
	    pw.println("MA "+markers[i]+" ON "+(i+2));
	}
	pw.println("MA SIZE 1.5");
	pw.println("MA SIZE 2 ON 2");
	pw.println("LAB 1 VPOS 0.15 0.75 \"Data\" CO 1 MA "+markers[0]+" MSIZE 1.5 CS 1.3");
	for ( int i=0; i < alpha.length; i++ ) {
	    pw.println("LAB "+(i+2)+" VPOS 0.15 "+(0.75-(i+1)*0.04)
		       +" \"Sim \\ga = "+alpha[i]+"\" CO "+(i+2)
		       +" MA "+markers[i+1]+" MSIZE 1.5 CS 1.3"); 
	}
	pw.println("!");
	pw.flush();


	/**  Loop on bin times to estimate the spectral index   **/
	System.out.println("Log  : Estimating spectral index ... ");
	double[][] psd = null;
	double[][] simPSD = null;
	double[] simFreq = null;
	double[] simPow = null;
	double slope = 0;
	double slopeErr = 0;
	double[] lsFitResults = new double[4];
	double[] simSlopes = new double[nSims];
	double[] simSlopeErrs = new double[nSims];

	for ( int i=nLCbins.length-1; i >= 0; i-- ) {

	    System.out.println("Log  :  Bin time = "+num.format(timeBins[i])+" s ("+nLCbins[i]+" bins)");
	    pw.print(num.format(timeBins[i])+"\t");

	    if ( dataTypeIsRate ) {

		double[] newTimes = new double[nLCbins[i]];
		double halfBinTime = timeBins[i]/2;
		for ( int j=0; j < nLCbins[i]; j++ ) {
		    newTimes[j] = halfBinTime + j*timeBins[i];
		}
 		double[] newRates = null;
		double[] newErrors = null;
		double[][] newRatesAndErrors = null;
		if ( useBinEdges ) {
		    
		    System.out.println("Log  :   Using bin edges");
		    if ( ncols == 4 ) {
			newRatesAndErrors = Binner.rebinRates(rate, error, binEdges, timeBins[i]);
			newRates = newRatesAndErrors[0];
			newErrors = newRatesAndErrors[1];
		    }
		    else {
			newRates = Binner.rebinRates(rate, binEdges, timeBins[i]);
		    }
		}
		else {
		    newRates = Binner.rebinRates(rate, oldBinTime, timeBins[i]);
		}

		psd = Periodograms.makePSD_FFT(newTimes, newRates, nLCbins[i], "leahy");
	    }

	    else {
		psd = Periodograms.makePSD_FFT(times, nLCbins[i], "leahy");
	    }
	    
	    /** Move to log-log space and fit  **/ 
	    double[] freq = psd[0];
	    double[] pow = psd[1];
	    double[] x = new double[freq.length];
	    double[] y = new double[pow.length];
	    for ( int j=0; j < freq.length; j++ ) {
		x[j] = Math.log10(freq[j]);
		y[j] = Math.log10(pow[j]);
	    }
	    /**  y(x) = a + b*testFreqs: method returns [a][b][err_a][err_b]  **/
	    lsFitResults = Stats.leastSquaresFitLine(x, y);
	    slope = -lsFitResults[1];
	    slopeErr = lsFitResults[3];
	    pw.print(slope+"\t"+slopeErr+"\t");

	    System.out.println("Log  :   Comparing with simulations:");
	    for ( int k=0; k < alpha.length; k++ ) {

		System.out.print("Log  :    index = "+alpha[k]+", simulating "+nSims+" event lists ... ");
		simPow = new double[freq.length];
		double[] sumOfPows = new double[freq.length];

		for ( int j=0; j < nSims; j++ ) {

		    double[] t = ArrivalTimes.generateRedArrivalTimes(equivMeanRate, duration, alpha[k]);
		    simPSD = Periodograms.makePSD_FFT(t, nLCbins[i], "leahy");
		    simFreq = simPSD[0];
		    simPow = simPSD[1];
		    x = new double[simFreq.length];
		    y = new double[simPow.length];
		    for ( int m=0; m < simPow.length; m++ ) {
			x[m] = Math.log10(simFreq[m]);
			y[m] = Math.log10(simPow[m]);
			sumOfPows[m] += simPow[m];
		    }
// 		    lsFitResults = Stats.leastSquaresFitLine(x, y);
// 		    simSlopes[j] = -lsFitResults[1];
// 		    simSlopeErrs[j] = lsFitResults[3];
		}
// 		double ave = Stats.getWMean(simSlopes, simSlopeErrs);
// 		double sig = Math.sqrt(Stats.getWVariance(simSlopes, simSlopeErrs)/nSims);

		double[] avePow = new double[freq.length];
		for ( int m=0; m < freq.length; m++ ) {
		    avePow[m] = Math.log10(sumOfPows[m]/nSims);
		}
		lsFitResults = Stats.leastSquaresFitLine(x, avePow);
		double ave = -lsFitResults[1];
		double sig = lsFitResults[3];
		logger.info("index = "+num.format(ave)+" +/- "+ num.format(sig));
		pw.print(num.format(ave)+"\t"+num.format(sig)+"\t");
	    }
	    pw.println();
	    pw.flush();
	}
	pw.close();
	logger.info("Result written to "+outName);

    }

}
