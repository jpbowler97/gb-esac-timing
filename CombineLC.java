package gb.esac.timing;

import java.io.IOException;


import gb.esac.io.DataFile;
import gb.esac.tools.Binner;


public class CombineLC {


    public static void main(String[] args) throws Exception, IOException {


	/**  Handle args  **/
	if ( args.length < 2 ) {
	    System.out.println("Usage java gb.esac.timing.CombineLC lc1.qdp lc2.qdp (lc3.qdp lc4.qdp ...)");
	    System.exit(-1);
	}
	int nfiles = args.length;
	String[] filenames = new String[nfiles];
	for ( int i=0; i < nfiles; i++ )
	    filenames[i] = args[i];

	
	/**  Get the data  **/
	double[] lcTimes = null;
	double tstart = 0;
	double[] startTimes = new double[nfiles];
	double[] oldBintimes = new double[nfiles];
	double newBintime = 0;
	DataFile[] dataFiles = new DataFile[nfiles];
	for ( int i=0; i < nfiles; i++ ) {
	    dataFiles[i] = new DataFile(filenames[i]);
	    lcTimes = dataFiles[i].getDblCol(0);
	    startTimes[i] = lcTimes[0];
	    oldBintimes[i] = lcTimes[1] - lcTimes[0];
	    tstart = Math.max(startTimes[i], lcTimes[0]);
	    newBintime = Math.max(newBintime, (lcTimes[1]-lcTimes[0]));
	}

	boolean startTimesAreEqual = true;
	boolean binTimesAreEqual = true;
	for ( int i=0; i < nfiles; i++ ) {
	    if ( tstart != startTimes[i] ) {
		startTimesAreEqual = false;
	    }
	    if ( newBintime != oldBintimes[i] ) {
		binTimesAreEqual = false;
	    }
	}

	double[][] rates = new double[nfiles][];
 	double[][] errors = new double[nfiles][];
	double[][] rebinnedRates = new double[nfiles][];
	double[][] rebinnedErrors = new double[nfiles][];
	int nCommonBins = Integer.MAX_VALUE;

	if ( ! startTimesAreEqual || ! binTimesAreEqual ) {

	    /**  Rebin the LCs and determine number of common bins  **/
	    for ( int i=0; i < nfiles; i++ ) {
		rates[i] = dataFiles[i].getDblCol(1);
		errors[i] = dataFiles[i].getDblCol(2);
		rebinnedRates[i] = Binner.rebinRatesSimple(rates[i], oldBintimes[i], newBintime, tstart);
		rebinnedErrors[i] = Binner.rebinRateErrors(errors[i], oldBintimes[i], newBintime, tstart);
		nCommonBins = (new Double(Math.min(nCommonBins, rebinnedRates[i].length))).intValue();
		//System.out.println(rates.length+"\t"+errors.length+"\t"+rebinnedRates[i].length+"\t"+rebinnedErrors[i].length);
	    }
	}

	else {

	    System.out.println("Log  : Start times and bin times are equal");

	    /**  Determine the number of common bins  **/
	    for ( int i=0; i < nfiles; i++ ) {
		rates[i] = dataFiles[i].getDblCol(1);
		errors[i] = dataFiles[i].getDblCol(2);
		nCommonBins = (new Double(Math.min(nCommonBins, rates[i].length))).intValue();
		rebinnedRates[i] = rates[i];
		rebinnedErrors[i] = errors[i];
	    }
	}


	/**  Combine the LCs  **/
	double[] combinedRates = new double[nCommonBins];
	double[] combinedErrors = new double[nCommonBins];
	for ( int j=0; j < nCommonBins; j++ ) {
	    double sumOfWeights = 0;
	    double weightedSum = 0;
	    for ( int i=0; i < nfiles; i++ ) {
		double weight = 1/Math.pow(rebinnedErrors[i][j], 2);
		weightedSum += weight*rebinnedRates[i][j];
		sumOfWeights += weight;
	    }
	    combinedRates[j] = weightedSum/sumOfWeights;
	    combinedErrors[j] = 1/Math.sqrt(sumOfWeights);
	}


 	/**  Determine final LC bin centres  **/
	double[] lcBinCentres = new double[nCommonBins];
	double halfBin = 0.5*newBintime;
	for ( int j=0; j < nCommonBins; j++ ) {
	    lcBinCentres[j] = tstart + halfBin + (j*newBintime);
	}


	/**  Write final combined LC  **/
	DataFile lcFile = null;
	String outFilename = "lc-combined.qdp";
	try { lcFile = new DataFile(outFilename); }
	catch (IOException e) {
		System.out.println("Error [CombineLC]: Could not open or create "+outFilename); 
		System.exit(-1);
	}
	String[] header = new String[] {
	    "DEV /XS",
	    "READ SERR 2",
	    "TIME OFF", "LINE ON",
	    "LAB T", "LAB F", "LW 3", "CS 1.3",
	    "LAB Y Rate (cts/s)",
	    "LAB X Time (s)",
	    "VIEW 0.1 0.2 0.9 0.8",
	    "!"
	};
	try { lcFile.writeData(header, lcBinCentres, combinedRates, combinedErrors); }
	catch (IOException e) {
	    System.out.println("Error: Could not write to file "+outFilename);
	    System.exit(-1);
	}


    }

}
