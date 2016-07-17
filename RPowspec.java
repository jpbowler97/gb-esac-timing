package gb.esac.timing;

import gb.esac.io.*;
import gb.esac.tools.*;

import java.io.IOException;
import java.io.FileNotFoundException;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import nom.tam.fits.*;



/**
 *  This class constructs the R2 and Z2 periodograms for an event list
 *  using the arrival times, or for a binned light curve with
 *  Time, Rate and Error columns.
 *
 * @version  21 October 2008 (last modified)
 * @author   Guillaume Belanger (ESAC, Spain)
 *
 */

public class RPowspec {

    public static void main(String[] args) throws Exception {

	DecimalFormat num = new DecimalFormat("0.00000");
	DecimalFormat freq = new DecimalFormat("0.0000E00");

	/***  Handle arguments  **/
	String evlistName = null;
	double pmin = 999;
	double pmax = 999;
	int sampling = 1;
	int nHarm = 1;
	boolean fileIsFits = true;
	if ( args.length == 1 ) {
	    evlistName = args[0];
	}
	else if ( args.length == 4 ) {
	    evlistName = args[0];
	    pmin = (Double.valueOf(args[1])).doubleValue();
	    pmax = (Double.valueOf(args[2])).doubleValue();
	    sampling = (Integer.valueOf(args[3])).intValue();
	}
	else {
	    System.out.println("Usage: java RPowspec evfile (pmin) (pmax) (samplingFactor)");
	    System.exit(-1);
	}

	//  Check if file is ascii or fits
	if ( evlistName.endsWith("fits") || evlistName.endsWith("fits.gz") )
	    fileIsFits = true;
	else if ( evlistName.endsWith("dat") || evlistName.endsWith("txt") || evlistName.endsWith("qdp") ) 
	    fileIsFits = false;
	else {
	    System.out.println("Error: File format not recognized");
	    System.out.println("Log  : filename extension can be 'fits, 'fits.gz', 'txt' or 'dat'");
	    System.exit(-1);
	}
	System.out.println("Log  : Running RPowspec");
	

	/**  Get data  **/
	Fits fitsFile;
	double[] times = null;
	double[] binWidths = null;
	double[] rate = null;
	double[] error = null;
	//double[] signif = null;
	boolean noRateCol = true, noErrorCol = true;
	boolean dataTypeIsRate = false;
	if ( fileIsFits ) {

	    System.out.println("Log  : File is in FITS format");
	    fitsFile = Astro.openFits(evlistName); 
	    BinaryTableHDU hdu = null;
	    try { hdu = Astro.getHDU(evlistName, "EVENTS"); }
	    catch (NullPointerException e) {
		hdu = Astro.getHDU(evlistName, "RATE"); 
	    }
		

	    //  Get TIME
	    try { times = (double[]) hdu.getColumn("TIME"); }
	    catch (FitsException e) {
		try { times = (double[]) hdu.getColumn("time"); }
		catch (FitsException e1) 		    
		    { System.out.println("Error: No TIME column found"); System.exit(-1); }
	    }

	    //  Get RATE
	    try { rate = (double[]) hdu.getColumn("RATE"); noRateCol = false; }
	    catch (FitsException e) { 
		try { rate = (double[]) hdu.getColumn("RATE1"); noRateCol = false; }
		catch ( FitsException e2 ) { noRateCol = true; }
		catch ( ClassCastException e3 ) {
		    float[] rateFlt = (float[]) hdu.getColumn("RATE1"); noRateCol = false;
		    rate = Converter.float2double(rateFlt);
		}
	    }		

	    //  Get ERROR
	    try { error = (double[]) hdu.getColumn("ERROR"); noErrorCol = false; }
	    catch (FitsException e) {
		try { error = (double[]) hdu.getColumn("ERROR1"); noErrorCol = false; }
		catch ( FitsException e2 ) { noErrorCol = true; }
		catch ( ClassCastException e3 ) {
		    float[] errorFlt = (float[]) hdu.getColumn("ERROR1"); noErrorCol = false;
		    error = Converter.float2double(errorFlt);
		}
	    }


	    //  Print out results 
	    if ( noRateCol==true && noErrorCol==true ) {
		System.out.println("Log  : File is event list (no RATE nor ERROR column)");
		dataTypeIsRate = false;
	    }
	    else if ( noRateCol==false && noErrorCol==true ) {
		dataTypeIsRate = true;
		System.out.print("Log  : File is light curve file (RATE column found)");
		System.out.print("Warn : No ERROR column. Constructing error ... ");
		for ( int i=0; i < times.length; i++ )
		    error[i] = Math.sqrt(rate[i]*times[i])/times[i];
		System.out.print("done");
	    }
	    else if ( noRateCol==false && noErrorCol==false ) {
		dataTypeIsRate = true;
		System.out.println("Log  : Using RATE and ERROR columns");
	    }
	    else { System.out.println("Error: Unrecognized format");  System.exit(-1); }

	}

	else {

	    System.out.println("Log  : File is in ASCII format");
	    DataFile dataFile = new DataFile(evlistName);
	    int ncols = dataFile.getNumOfDataCols();

	    if ( ncols == 1 ) {
		System.out.println("Log  : There is only 1 col, reading TIME col");
		//  Get TIME
		try { times = dataFile.getDblCol(0); }
		catch (IOException e)
		    { System.out.println("Error: Cannot read column 0"); System.exit(-1); }
		dataTypeIsRate = false;
	    }
	    if ( ncols == 2 ) {
		dataTypeIsRate = true;
		System.out.println("Log  : There are 2 cols, reading TIME and RATE cols");
		//  Get TIME and RATE
		try { times = dataFile.getDblCol(0); }
		catch (IOException e)
		    { System.out.println("Error: Cannot read col 0"); System.exit(-1); }
		try { rate = dataFile.getDblCol(1); }
		catch (IOException e)
		    { System.out.println("Error: Cannot read column 1"); System.exit(-1); }
		System.out.print("Warn : No ERROR column. Constructing errors ... ");
		error = new double[times.length];
		for ( int i=0; i < times.length; i++ )
		    error[i] = Math.sqrt(rate[i]*times[i])/times[i];
		System.out.print("done");
	    }
	    if ( ncols == 3 ) {
		dataTypeIsRate = true;
		System.out.println("Log  : There are 3 cols, reading TIME, RATE and ERROR cols");
		try { times = dataFile.getDblCol(0); }
		catch (IOException e)
		    { System.out.println("Error: Cannot read column 0"); System.exit(-1); }

		try { rate = dataFile.getDblCol(1); }
		catch (IOException e)
		    { System.out.println("Error: Cannot read column 1"); System.exit(-1); }

		try { error = dataFile.getDblCol(2); noErrorCol = false; }
		catch (IOException e) 
		    { System.out.println("Error: Cannot read column 2"); 
			noErrorCol = true;
			System.exit(-1);
		    }
	    }
	    if ( ncols == 4 ) {
		System.out.println("Log  : There are 4 cols, reading TIME, BIN, RATE and ERROR cols");
		try { times = dataFile.getDblCol(0); }
		catch (IOException e)
		    { System.out.println("Error: Cannot read column 0"); System.exit(-1); }

		try { binWidths = dataFile.getDblCol(1); }
		catch (IOException e)
		    { System.out.println("Error: Cannot read column 1"); System.exit(-1); }

		try { rate = dataFile.getDblCol(2); }
		catch (IOException e)
		    { System.out.println("Error: Cannot read column 2"); System.exit(-1); }

		try { error = dataFile.getDblCol(3); noErrorCol = false; }
		catch (IOException e) 
		    { System.out.println("Error: Cannot read column 3"); 
			noErrorCol = true;
			System.exit(-1);
		    }
		dataTypeIsRate = true;
	    }
	}
	

	/**  Determine and print out duration and mean count rate  **/
	double duration = Stats.getMax(times) - Stats.getMin(times);
	System.out.println("Log  : Duration = "+num.format(duration)+" s");
	double meanRate = 0;
	double wMeanRate = 0;
	if ( !dataTypeIsRate ) {
	    meanRate = times.length/duration;
	    System.out.println("Log  : Mean Rate = "+num.format(meanRate)+" cts/s");
	}
	else {
	    meanRate = Stats.getMean(rate);
	    wMeanRate = Stats.getWMean(rate, error);
	    System.out.println("Log  : Mean Rate = "+num.format(meanRate)+" cts/s");
	    System.out.println("Log  : Weighted Mean Rate = "+num.format(wMeanRate)+" cts/s");
	}

	
	/**  Construct binned light curve if data is event list  **/
	if ( dataTypeIsRate == false ) {
	    double bintime = 120;
	    int nbins = (new Double(Math.ceil(duration/bintime))).intValue();
	    bintime = duration/nbins;
	    double[] binHeights =  Binner.binOrderedData(times, nbins);
	    double[] lc = new double[nbins];
	    double[] lcErr = new double[nbins];
	    double[] t = new double[nbins];
	    double[] tErr = new double[nbins];
	    double maxY = -Double.MAX_VALUE;
	    double minY = Double.MAX_VALUE;
	    for ( int i=0; i < nbins; i++ ) {
		lc[i] = binHeights[i]/bintime;
		lcErr[i] = Math.sqrt(binHeights[i])/bintime;
		t[i] = (0.5+i) * bintime;
		tErr[i] = bintime/2;
		maxY = Math.max(maxY, lc[i]+lcErr[i]);
		minY = Math.min(minY, lc[i]-lcErr[i]);
	    }


	    /**  Write binned lc as QDP file  **/
	    DataFile qdpFile = new DataFile("lc_120s.qdp");
	    String[] header = new String[] {
		"! QDP data file", 
		"DEV /XS", "READ SERR 1 2",
		"LAB TITLE", "LAB FILE", "TIME OFF",
		"LINE ON", "LW 3", "CS 1.3", 
		"R Y "+minY+" "+maxY,
		"LAB X Time (s)",
		"LAB Y Count rate (cts/s)",
		"VIEW 0.1 0.2 0.9 0.8",
		"!"
	    };
	    qdpFile.writeData(header, t, tErr, lc, lcErr);
	}


	/**  Define pmin and pmax if these were not specified in arguments  **/
	if ( pmin == 999 ) pmin = 100;
	if ( pmax == 999 ) pmax = duration/3;
	double nuMin = 1/pmax;
	double nuMax = 1/pmin;
	System.out.println("Log  : Minimum test period = "+ pmin +" s ("+freq.format(nuMax)+" Hz)");
	System.out.println("Log  : Maximum test period = "+ pmax +" s ("+freq.format(nuMin)+" Hz)");


	/**  Define test frequencies  **/
	int nIFS =(new Double(Math.floor(duration*(nuMax - nuMin)))).intValue();
	double ifsOffset = 0;
	double[] trialFreqs = TimingUtils.getFourierTestFrequencies(nuMin, nuMax, duration, sampling, ifsOffset);
	int nTrials = trialFreqs.length;
	System.out.println("Log  : Number of IFS = "+nIFS);
	System.out.println("Log  : Sampling factor = "+ sampling);
	System.out.println("Log  : Total number of trials = "+ nTrials);


	/**  Construct periodogram of data  **/
	double[] z2 = new double[nTrials];
	double[] alpha = new double[nTrials];
	double[] beta = new double[nTrials];
	double[] r2 = new double[nTrials];
	double[] lomb = new double[nTrials];
	double z2max = 0;
	double z2maxFreq = 0;
	int z2maxIndex = 0;
	double r2max = 0;
	double r2maxFreq = 0;
	int r2maxIndex = 0;
	double lombMax = 0;
	double lombMaxFreq = 0;
	int lombMaxIndex = 0;
	double period = 0;

	if ( dataTypeIsRate ) {

	    /**  Modified R-test and Lomb-Scargle  **/

	    for ( int i=0; i < nTrials; i++ ) {		
		period = 1.0/trialFreqs[i];
		r2[i] = Stats.getModRayleighPower(times, rate, error, period);
		lomb[i] = Stats.getLombPower(times, rate, period);
		
	    }	
	    r2max = Stats.getMax(r2);
	    r2maxIndex = Stats.getIndex(r2, r2max);
	    r2maxFreq = trialFreqs[r2maxIndex];

	    lombMax = Stats.getMax(lomb);
	    lombMaxIndex = Stats.getIndex(lomb, lombMax);
	    lombMaxFreq = trialFreqs[lombMaxIndex];

	}
	else {

	    /**  Modified R-test  **/

	    for ( int i=0; i < nTrials; i++ ) {		
		
		period = 1.0/trialFreqs[i];		
		r2[i] = Stats.getModRayleighPower(times, period);

		//alphaData[i] = (Stats.getZ2stats(times, period, nHarm))[0];
		//betaData[i] = (Stats.getZ2stats(times, period, nHarm))[1];
 		//z2[i] = (Stats.getZ2stats(times, period, nHarm))[2];
		z2[i] = Stats.getZ2Power(times, period, nHarm);
		//if ( z2[i] > z2max ) {
		//    z2maxIndex = i;
		//    z2maxFreq = trialFreqs[i];
		//}
		//z2max = Math.max(z2[i], z2max);
	    }
	    r2max = Stats.getMax(r2);
	    r2maxIndex = Stats.getIndex(r2, r2max);
	    r2maxFreq = trialFreqs[r2maxIndex];


	    /**  Print out the results of the Z test  **/
// 	    System.out.println("Log  : Mean of Z2 = " + Stats.getMean(z2));
// 	    System.out.println("Log  : Variance of Z2 = " + Stats.getVariance(z2));
// 	    System.out.println("Log  : Max Z2 = " + num.format(z2max) + " ("+freq.format(z2maxFreq)+" Hz = "
// 			   + num.format(1/z2maxFreq)+" s)");


// 	    /**  Make histogram of Z2 values  **/
// 	    int histoBins = 50;
// 	    double r2_max = Stats.getMax(r2);
// 	    double lomb_max = Stats.getMax(lomb);
// 	    double min = 0;
// 	    double width = (max - min)/histoBins;
// 	    //double[] histoOfZ2 = Binner.binData(z2, min, max, histoBins);
// 	    int nValues = (new Double(Math.ceil(Stats.getSum(histoOfZ2)))).intValue();
// 	    double halfWidth = width/2;
// 	    double[] centres = new double[histoBins]; 
// 	    double[] normPDF = new double[histoBins];
// 	    for ( int i=0; i < histoBins; i++ ) {
// 		normPDF[i] = histoOfZ2[i]/(nValues*width);
// 		centres[i] = min + halfWidth + i*width;
// 	    }
	    
	    
// 	    /**  Calculated the chi squared using f(x) = 1/2 exp(-x/2)  **/
// 	    double expected = 0;
// 	    double pearsonChi2 = 0;
// 	    for ( int i=0; i < histoBins; i++ ) {
// 		expected = 0.5*Math.exp(-centres[i]/2);
// 		pearsonChi2 += 
//	    Math.pow( (histoOfZ2[i] - nValues*width*expected), 2)/(nValues*width*expected);
// 	    }
// 	    int dof = histoBins - 1;
// 	    pearsonChi2 /= dof;	
// 	    System.out.println("Log  : Pearson's Reduced Chi Squared ("+dof+" DOF) = "
// 			       +number.format(pearsonChi2));
	    
	    
	    /**  Write histogram as QDP file  **/
// 	    DataFile qdpFile = new DataFile("z2Histo.qdp");
// 	    String[] header = new String[] {
// 		"! QDP data file", 
// 		"DEV /xs",
// 		"LABEL TITLE", "LABEL FILE", "TIME OFF", 
// 		"LW 3", "CSIZE 1.3", "LINE STEP",
// 		"LAB X Z\\u2\\d Power",
// 		"LAB Y Normalized PDF",
// 		"LOG OFF",
// 		"VIEW 0.1 0.2 0.9 0.8"
// 	    };
// 	    qdpFile.writeData(header, centres, normPDF);
	    
	    
	    /**  Write Z2 periodogram as a QDP file  **/
	    DataFile qdpFile = new DataFile("zPow.qdp");
	    String[] header = new String[] {
		"! QDP data file", 
		"DEV /XS",
		"LABEL TITLE", "LABEL FILE", "TIME OFF", 
		"LINE ON", "LW 3", "CSIZE 1.3",
		"LABEL X Frequency (Hz)", 
		"LABEL Y Z\\u2\\d Power",
		"LOG X ON", "LOG Y OFF",
		"VIEW 0.1 0.2 0.9 0.8"
	    };
	    qdpFile.writeData(header, trialFreqs, z2);
 	}


	/**  Write the results of the Modified Rayleigh Test  **/
	System.out.println("Log  : Mean of R2 = " + Stats.getMean(r2));
	System.out.println("Log  : Variance of R2 = " + Stats.getVariance(r2));
	System.out.println("Log  : Max R2 = " + num.format(r2max) + " ("+freq.format(r2maxFreq)+" Hz = "
			   + num.format(1/r2maxFreq)+" s)");
	
	
	/**  Write modified Rayleigh periodogram as a QDP file  **/
	DataFile qdpFile = new DataFile("rPow.qdp");
	String[] header = new String[] {
	    "! QDP data file", 
	    "DEV /XS",
	    "LABEL TITLE", "LABEL FILE", "TIME OFF", 
	    "LINE ON", "LW 3", "CSIZE 1.3",
	    "LAB X Frequency (Hz)", 
	    "LAB Y Rayleigh Power",
	    "LOG X ON",
	    "VIEW 0.1 0.2 0.9 0.8"
	};
	if ( dataTypeIsRate ) 
	    qdpFile.writeData(header, trialFreqs, r2, lomb);
	else
	    qdpFile.writeData(header, trialFreqs, r2);
	

	/**  Make histogram of R2 and Lomb power values  **/
	int histoBins = 30;
	double min = 0;
	double r2_max = 15; //Math.ceil(Stats.getMax(r2)+0.5);
	double lomb_max = 10; //Math.ceil(Stats.getMax(lomb)+0.5);
	double r2_width = (r2_max - min)/histoBins;
	double lomb_width = (lomb_max - min)/histoBins;
	double[] histoOfR2 = Binner.binData(r2, min, r2_max, histoBins);
	double[] histoOfLomb = Binner.binData(lomb, min, lomb_max, histoBins);
	int nValues = (new Double(Math.ceil(Stats.getSum(histoOfR2)))).intValue();
	double r2_halfWidth = r2_width/2;
	double lomb_halfWidth = lomb_width/2;
	double[] r2_centres = new double[histoBins]; 
	double[] lomb_centres = new double[histoBins]; 
	double[] normPDF_r2 = new double[histoBins];
	double[] normPDF_lomb = new double[histoBins];
	double[] expect_r2 = new double[histoBins];
	double[] expect_lomb = new double[histoBins];
	double pearsonChi2_r2 = 0;
	double pearsonChi2_lomb = 0;

	for ( int i=0; i < histoBins; i++ ) {
	    normPDF_r2[i] = histoOfR2[i]/(nValues*r2_width);
	    r2_centres[i] = r2_halfWidth + i*r2_width;
	    expect_r2[i] = 0.5*Math.exp(-r2_centres[i]/2);
	    pearsonChi2_r2 += Math.pow( (histoOfR2[i] - nValues*r2_width*expect_r2[i]), 2) / 
		(nValues*r2_width*expect_r2[i]);

	    normPDF_lomb[i] = histoOfLomb[i]/(nValues*lomb_width);
	    lomb_centres[i] = lomb_halfWidth + i*lomb_width;
	    expect_lomb[i] = Math.exp(-lomb_centres[i]);
	    pearsonChi2_lomb += Math.pow( (histoOfLomb[i] - nValues*lomb_width*expect_lomb[i]), 2) / 
		(nValues*lomb_width*expect_lomb[i]);
	}
	int dof = histoBins - 1;
	pearsonChi2_r2 /= dof;
	pearsonChi2_lomb /= dof;
	System.out.println("Log  : Pearson's Reduced Chi Squared for ("+dof+" DOF) = "
			   +freq.format(pearsonChi2_r2)+" (R-test) and "+
			   freq.format(pearsonChi2_lomb)+" (Lomb)");
	
	
	/**  Write histogram as QDP file  **/
	if ( dataTypeIsRate ) {
	    PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("r2Histo.qdp")));
	    header = new String[] {
		"! QDP data file", 
		"DEV /XS",
		"LABEL TITLE", "LABEL FILE", "TIME OFF", 
		"LW 3", "CSIZE 1.3", 
		"LINE STEP 1", "LINE ON 2", "LINE STEP 3", "LINE ON 4",
		"CO 2 ON 1", "CO 1 ON 2", "CO 3 ON 3", "CO 1 ON 4",
		"LAB X Power",
		"LAB Y Normalized PDF",
		"LAB  2 COL 2 CS 1.3 JUS Lef",
		"LAB  2 VPOS 0.3 0.35 \"f\\d2\\u(x) = 1/2 e\\u-x/2\\d (Rayleigh: \\gx\\u2\\d\\d2\\u)\"",
		"LAB  3 COL 3 CS 1.3 JUS Lef",
		"LAB  3 VPOS 0.2 0.7 \"f\\d1\\u(x) = e\\u-x\\d (Lomb-Scargle)\"",
		"VIEW 0.1 0.17 0.9 0.8",
		"SKIP SINGLE",
		"R X 0",
		"!"
	    };
	    for ( int i=0; i < header.length; i++ ) pw.println(header[i]);
	    for ( int i=0; i < histoBins; i++ ) {
		pw.println(num.format(r2_centres[i])+"\t"+ normPDF_r2[i]+"\t"+expect_r2[i]);
	    }
	    pw.println("NO NO");
	    for ( int i=0; i < histoBins; i++ ) {
		pw.println(num.format(lomb_centres[i])+"\t"+ normPDF_lomb[i]+"\t"+expect_lomb[i]);
	    }
	    pw.flush();
	    pw.close();
	}
	else {
	    qdpFile = new DataFile("r2Histo.qdp");
	    header = new String[] {
		"! QDP data file", 
		"DEV /XS",
		"LABEL TITLE", "LABEL FILE", "TIME OFF", 
		"LW 3", "CSIZE 1.3", 
		"LINE STEP 2", "LINE ON 3",
		"CO 2 ON 2", "CO 1 ON 3",
		"LAB X Modified Rayleigh Power",
		"LAB Y Normalized PDF",
		"LAB  2 COL 1 CS 1.5 JUS Lef",
		"LAB  2 POS 6 0.35 \"f\\d2\\u(x) = 1/2 e\\u-x/2\\d\"",
		"VIEW 0.1 0.2 0.9 0.8",
		"!"
	    };
	    qdpFile.writeData(header, r2_centres, normPDF_r2, expect_r2);
	}
	
    }  
}

/////  Plotter

// 	/**  Plotter **/
// 	//  Set up plotter (for available IAxisStyle parameters see bottom) 
// 	IAnalysisFactory af = IAnalysisFactory.create();
// 	IPlotterFactory plotterF = af.createPlotterFactory();
// 	IPlotter plotter = plotterF.create();
// 	double x = 0, y1 = 0, y2 = 0.5, w = 1, h = 0.5;
// 	IPlotterRegion region0 = plotter.createRegion(x, y1, w, h);
// 	IPlotterRegion region1 = plotter.createRegion(x, y2, w, h);

// 	//  Define plot style for Light curve
// 	IAxisStyle x0AxisStyle = plotterF.createAxisStyle();
// 	x0AxisStyle.setParameter("label", "Time (s)");
// 	IPlotterStyle plotterStyle0 = plotterF.createPlotterStyle();
// 	plotterStyle0.setAxisStyleX(x0AxisStyle);
// 	region0.applyStyle(plotterStyle0);

// 	//  Define plot style for Periodogram
// 	IAxisStyle x1AxisStyle = plotterF.createAxisStyle();
// 	x1AxisStyle.setParameter("scale", "logarithmic");
// 	x1AxisStyle.setParameter("label", "Period (s)");
// 	IAxisStyle yAxisStyle = plotterF.createAxisStyle();
// 	//yAxisStyle.setParameter("scale", "logarithmic");
// 	IPlotterStyle plotterStyle1 = plotterF.createPlotterStyle();
// 	plotterStyle1.setAxisStyleX(x1AxisStyle);
// 	plotterStyle1.setAxisStyleY(yAxisStyle);
// 	region1.applyStyle(plotterStyle1);
// 	/**  Plotter end  **/

// 	//  Construct and display light curve
// 	double bintime = 100;
// 	int nbins = (new Double(Math.ceil(duration/bintime))).intValue();
// 	double tfirst = 0;
// 	double tzero = Stats.getMin(times);
// 	double tlast = Stats.getMax(times);
// 	IAxis axis = new FixedAxis(nbins, tfirst, tlast);
// 	Histogram1D lc = new Histogram1D("", "Light Curve", axis);
// 	double[] rebinnedRates = null;
// 	double[] rebinnedErrors = null;
// 	if ( dataTypeIsRate ) {
// 	    rebinnedRates = Binner.rebinRates(rate, nbins, bintime);
// 	    rebinnedErrors = Binner.rebinRateErrors(error, nbins, bintime);
// 	    lc = AnalysisTools.array2histo("Light Curve", tfirst, bintime, rebinnedRates, rebinnedErrors);
// 	}
// 	else 
// 	    for ( int i=0; i < times.length; i++ ) {
// 		times[i] -= tzero;
// 		lc.fill(times[i]);
// 	    }
// 	region0.plot(lc);
// 	plotter.show();



/////       IAxisStyle  PARAMETERS

// 	String[] params = axisStyle.availableParameters();
// 	for ( int i=0; i < params.length; i++ ) {
// 	    System.out.print("Parameter: "+params[i]+" ");
// 	    String[] paramOptions = axisStyle.availableParameterOptions(params[i]);
// 	    System.out.print("(");
// 	    for ( int j=0; j < paramOptions.length; j++ )
// 		System.out.print(paramOptions[j]+", ");
// 	    System.out.println(")");
// 	}
// 	// Output of above is:
// 	Parameter: scale (lin, linear, log, logarithmic, )
// 	Parameter: type (double, int, date, string, )
// 	Parameter: label ()
// 	Parameter: allowZeroSuppression (true, false, )
// 	Parameter: yAxis (Y0, Y1, )


// 	//  Fit the light curve with a low degree polynomial
// 	IFitFactory fitf = af.createFitFactory();
// 	IFitter fitter = fitf.createFitter("chi2");
// 	IFitData fitData = fitf.createFitData();
// 	fitData.create1DConnection(lc);
// 	fitData.range(0).excludeAll();
// 	fitData.range(0).include(0, 4000);
// 	IFitResult fitResult = fitter.fit(fitData,"p2", new double[]{-1, 0.025, -6.25e-6});
// 	fitData.range(0).excludeAll();
// 	fitData.range(0).include(4000, 8900);
// 	IFitResult fitResult2 = fitter.fit(fitData,"p2", new double[]{-174.8,0.08,-4.8e-6});
