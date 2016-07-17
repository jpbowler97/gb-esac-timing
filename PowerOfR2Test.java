package gb.esac.timing;

import gb.esac.io.*;
import gb.esac.tools.*;
import gb.esac.aida.functions.ChiSquareFunction;

import java.io.*;
import java.util.Arrays;
import java.text.DecimalFormat;

import cern.jet.random.engine.MersenneTwister64;
import cern.jet.random.Uniform;
import cern.jet.stat.Descriptive;
import cern.jet.stat.Probability;
import cern.colt.list.DoubleArrayList;

import hep.aida.*;

import nom.tam.util.ArrayFuncs;



public class PowerOfR2Test {


    /**  Define Global Variables  **/
    public static int nPulseFreqs = 100;
    public static double duration = 10000;
    public static double meanRate = 0.5;
    public static int nCyclesMin = 3;
    public static int nCyclesMax = 10;
    public static double snrMin = 0;
    public static double snrMax = 10;
    public static double alpha = 0.0;

    
    /**  Number formats  **/
    public static DecimalFormat number = new DecimalFormat("0.0000");
    public static DecimalFormat freq = new DecimalFormat("0.00000E00");


    /**  Start Main  **/
    public static void main (String[] args) throws IOException {

	
	/**  Handle arguments  **/
	handleArguments(args);


	/**  Determine and define pulsed fractions that will be tested  **/
	int nSNRs = 30;
	double step = (snrMax - snrMin)/(nSNRs-1);
	double[] snr = new double[nSNRs];
	double[] pulsedFracs = new double[nSNRs];
	double bigN = duration*meanRate;
	System.out.println("Log  : "+nSNRs+" values of s/n (pulsed fraction) will be tested:");


 	/**  Print out SNRs  **/
 	//System.out.println("Log  :     s/n     pf");
 	for ( int i=0; i < nSNRs; i++ ) {
 	    snr[i] = snrMin +i*step;
 	    pulsedFracs[i] = snr[i]*Math.sqrt(bigN)/bigN;  // snr = p*sqrt(N)
  	    //System.out.println("Log  :    "+number.format(snr[i])+"   "+number.format(pulsedFracs[i]));
 	}


	/**  Define modulation frequencies **/
	double nuMin = nCyclesMin/duration;
	double nuMax = nCyclesMax/duration;
	int nIFS = (new Double(Math.floor(duration*(nuMax - nuMin)))).intValue();
	System.out.println("Log  : Range of modulation frequencies:");
	System.out.println("Log  : nuMin ("+nCyclesMin+" cycles/duration) = "+freq.format(nuMin)+" Hz"); 
	System.out.println("Log  : nuMax ("+nCyclesMax+" cycles/duration) = "+freq.format(nuMax)+" Hz"); 
	System.out.println("Log  : This range contains "+nIFS+" IFS, over the period range ["+
			   number.format(1/nuMax)+", "+number.format(1/nuMin)+"] s");

	/**  Uniformly distributed  **/
// 	System.out.println("Log  : Defining "+nPulseFreqs+" uniformly distributed freqs in ["
// 			   +freq.format(nuMin)+", "+freq.format(nuMax)+"] Hz");
//  	for ( int i=0; i < nPulseFreqs; i ++ ) {
//   	    pulseFreqs[i] = nuMin + i*(nuMax-nuMin)/(nPulseFreqs-1);
// 	}

	/**  Randomly distributed  **/
	System.out.println("Log  : Drawing "+nPulseFreqs+" random, uniformly distributed freqs in ["
			   +freq.format(nuMin)+", "+freq.format(nuMax)+"] Hz");
	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
	Uniform randUni = new Uniform(nuMin, nuMax, engine);
	double[] pulseFreqs = new double[nPulseFreqs];
	for ( int i=0; i < nPulseFreqs; i ++ ) {
  	    pulseFreqs[i] = randUni.nextDouble();
	}
  	Arrays.sort(pulseFreqs);

 	/**  Print out frequencies  **/
// 	System.out.println("Log  : Test frequencies:");
// 	for ( int i=0; i < nPulseFreqs; i ++ )
// 	    System.out.println("Log  : "+(i+1)+") "+freq.format(pulseFreqs[i])+" Hz \t"+
// 			       number.format(1/pulseFreqs[i])+" s");


	/**  Set up JAIDA factories  **/
// 	IAnalysisFactory af = IAnalysisFactory.create();
// 	ITree tree = af.createTreeFactory().create();
// 	IHistogramFactory hf = af.createHistogramFactory(tree);
// 	IFunctionFactory funcF  = af.createFunctionFactory(tree);
//  	IFitFactory fitF   = af.createFitFactory();
// 	IFitter fitter = fitF.createFitter("Chi2", "jminuit");


	/**  Display histograms **/

// 	IPlotter plotter = af.createPlotterFactory().create("Peak Heights and Positions");
// 	plotter.createRegions(2,2);

// 	IPlotterStyle plotterStyle = plotter.style();
// 	IDataStyle dataStyle = plotterStyle.dataStyle();
// 	IFillStyle fillStyle = dataStyle.fillStyle();
// 	ILineStyle lineStyle = dataStyle.lineStyle();
// 	IMarkerStyle markerStyle = dataStyle.markerStyle();
	
// 	lineStyle.setParameter("color", "black");
// 	markerStyle.setParameter("shape", "dot");
// 	markerStyle.setParameter("size", "5");
// 	markerStyle.setParameter("color", "black");
// 	fillStyle.setParameter("color", "white");

// 	dataStyle.setLineStyle(lineStyle);
// 	dataStyle.setFillStyle(fillStyle);
// 	dataStyle.setMarkerStyle(markerStyle);

	/**  Up to here  **/


	/**  Start simulations  **/
	double r = meanRate;
	double t = duration;
// 	int nbins = 256;
	double sampling = 10;


	/**  Define arrays to keep track of peaks  **/
// 	double[] fftPeaks = new double[nPulseFreqs];
//  	double[] fftLocs = new double[nPulseFreqs];
// 	double[] modRayPeaks = new double[nPulseFreqs];
// 	double[] modRayLocs = new double[nPulseFreqs];
// 	double[] lombPeaks = new double[nPulseFreqs];
// 	double[] lombLocs = new double[nPulseFreqs];

	//DoubleArrayList fftPeaks = new DoubleArrayList();
	DoubleArrayList modRayPeaks = new DoubleArrayList();
	//DoubleArrayList lombPeaks = new DoubleArrayList();
	//double[][] fftPSD = null;
	double[][] modRayPSD = null;
	//double[][] lombPSD = null;
	//DoubleArrayList fftNoisePowers = new DoubleArrayList(); 
	DoubleArrayList modRayNoisePowers = new DoubleArrayList();
	//DoubleArrayList lombNoisePowers = new DoubleArrayList();
	

	/**  Create printwriters to write results to text files  **/
	File out = new File("powerOfR2Test.qdp");
	String outFilename = null;
	if ( out.exists() )
	    outFilename = "powerOfR2Test_2.qdp";
	else
	    outFilename = "powerOfR2Test.qdp";

	PrintWriter powerFile = new PrintWriter(new BufferedWriter(new FileWriter(outFilename)));
	double xmin = -0.005;
	double xmax = pulsedFracs[nSNRs-1]+0.001;
	double confLevel = 0.9;
	String[] header = new String[]{
	    "DEV /XS",
	    "READ 2 3",
	    "LINE ON", "LW 4", "CS 1.3",
	    "LAB FILE", "TIME OFF",
	    "VIEW 0.1 0.2 0.9 0.8",
	    "CO 1 ON 2", "CO 2 ON 3", "CO 3 ON 4",
// 	    "MA SIZE 1.2",
// 	    "MA 11 ON 2",
// 	    "MA 6 ON 3",
// 	    "MA 4 ON 4",
// 	    "LAB X Signal-to-Noise Ratio",
	    "LAB X Pulsed Fraction",
	    "LAB Y Detection Probability >= 3\\gs",
	    "! R Y -0.05 1.05",
	    "! R X "+number.format(xmin)+" "+number.format(xmax),
//  	    "LAB 2 VPOS 0.2 0.7 \"Modified Rayleigh\" CO 2 JUST LEFT CS 1.2",
//  	    "LAB 3 VPOS 0.2 0.65 \"Lomb-Scargle\" CO 3 JUST LEFT CS 1.2",
//  	    "LAB 4 VPOS 0.2 0.6 \"FFT\" CO 1 JUST LEFT CS 1.2",
	    "LAB 5 POS "+number.format(xmin)+" "+number.format(confLevel)+" LS 2 LINE 0 0.99 \"",
	    "!"
	};
	for ( int k=0; k < header.length; k++ ) {
	    powerFile.println(header[k]);
	}


	/**   Start loop on SNRs  **/
	System.out.println("Log  : Starting simulations:");
	for ( int j=0; j < nSNRs; j++ ) {

	    double pf = pulsedFracs[j];
 	    System.out.println("Log  : s/n "+(j+1)+" = "+number.format(snr[j])+" (pf="+number.format(pf*100)+" %)");
	    //powerFile.print(snr[j] + "\t");
	    powerFile.print(freq.format(pf) + "\t");
	    System.out.println("Log  :   Time (in) = "+DateUtils.time());

	    /**  Start loop on Test Frequencies  **/
// 	    double nFFTPeaks = 0.0;
// 	    double nModRayPeaks = 0.0;
// 	    double nLombPeaks = 0.0;
	    for ( int i=0; i < nPulseFreqs; i ++ ) {

		double period = 1/pulseFreqs[i];
//  		System.out.println("Log  : Modulation frequency "+(i+1)+" = "+freq.format(pulseFreqs[i])+
//  				   " Hz (period = "+number.format(1/pulseFreqs[i])+" s)");

		/**  Generate arrival times **/
		double nEvents = duration*meanRate;
		double nPulsedEvents = pf*nEvents;
		double pulsedMeanRate = nPulsedEvents/duration;
		double nRedNoiseEvents = nEvents - nPulsedEvents;
		double redNoiseMeanRate = nRedNoiseEvents/duration;
		double[] pulsedTimes = 
		    WhiteNoiseGenerator.generateModulatedArrivalTimes(pulsedMeanRate, duration, period, 0.999, engine);
		double[] redNoiseTimes = 
		    RedNoiseGenerator.generateArrivalTimes(redNoiseMeanRate, duration, alpha, engine);

		
		/**  Make periodograms of red noise in the absence of pulsation  **/
		t = duration;
		//double binTime = t/nbins;
		//double nuMaxFFT = 1/(2*binTime);
		double nuMinFFT = 1/t;
		double ifs = nuMinFFT;
		double searchMin = pulseFreqs[i] - 1.2*ifs;
		double searchMax = pulseFreqs[i] + 1.2*ifs;
		//Object[] lc = TimingUtils.makeLC(redNoiseTimes, nbins);
		//double[] t_lc = (double[]) lc[0];
		//double[] rates = (double[]) lc[1];
		//double[][] fftPSD_noise = Periodograms.makePSD_FFT(redNoiseTimes, nbins, "leahy-like");
		double[] modRayPSD_noise = Periodograms.makePSD_ModRay(redNoiseTimes, searchMin, searchMax, sampling)[1];
		//double[] lombPSD_noise = Periodograms.makePSD_Lomb(t_lc, rates, searchMin, searchMax, sampling)[1];


		/**  Store the powers as reference for the significance of the peaks  **/

		/** FFT  **/
// 		int idx = 0;
// 		while ( fftPSD_noise[0][idx] < searchMin ) idx++;
// 		while ( fftPSD_noise[0][idx] <= searchMax ) {
// 		    fftNoisePowers.add(fftPSD_noise[1][idx]);
// 		    idx++;
// 		}

		/**  Mod Ray  **/
		for ( int p=0; p < modRayPSD_noise.length; p++ ) {
		    //System.out.println(modRayPSD_noise[p]);
		    modRayNoisePowers.add(modRayPSD_noise[p]);
		}

		/**  Lomb  **/
// 		for ( int p=0; p < lombPSD_noise.length; p++ ) {
// 		    lombNoisePowers.add(lombPSD_noise[p]);
// 		}


		/**  Combine and sort the events  **/
		int nTimes = pulsedTimes.length + redNoiseTimes.length - 2;
		double[] times = new double[nTimes];
		if ( pulsedTimes.length <= 2 )
		    times = redNoiseTimes;
		else {
		    int n = 0;
		    for ( int p=1; p < pulsedTimes.length-1; p++ ) {
			times[p-1] = pulsedTimes[p];
			n++;
		    }
		    for ( int p=0; p < redNoiseTimes.length; p++ ) {
			times[p+n] = redNoiseTimes[p];
		    }
		    Arrays.sort(times);
		}
		//double[] times = ArrivalTimes.generateRedModulatedArrivalTimes(r, duration, alpha, period, pf);


		/**  Construct PSDs of the combined event lists  **/
		//lc = TimingUtils.makeLC(times, nbins);
		//t_lc = (double[]) lc[0];
		//rates = (double[]) lc[1];
		//fftPSD = Periodograms.makePSD_FFT(times, nbins, "leahy-like");
		modRayPSD = Periodograms.makePSD_ModRay(times, searchMin, searchMax, sampling);
		//lombPSD = Periodograms.makePSD_Lomb(t_lc, rates, searchMin, searchMax, sampling);


 		/**  Find peak In the FFT PSD **/
// 		double fftPeakHeight = 0;
// 		double fftPeakLoc = 0;
// 		idx=0;
// 		while ( fftPSD[0][idx] < searchMin ) idx++;
// 		while ( fftPSD[0][idx] <= searchMax ) {
// 		    if ( fftPSD[1][idx] >= fftPeakHeight ) {
// 			fftPeakLoc = fftPSD[0][idx];
// 			fftPeakHeight = fftPSD[1][idx];
// 		    }
// 		    idx++;
// 		}
// 		fftPeaks.add(fftPeakHeight);
// 		fftPeaks[i] = fftPeakHeight;
// 		if ( fftPeakHeight >= 13.215 ) {
// 		    System.out.println("yes");
// 		    nFFTPeaks++;
// 		}


		/**  Find peak In the Modified Rayleight PSD (search range is already restricted) **/
		double modRayPeakHeight = Stats.getMax(modRayPSD[1]);
		modRayPeaks.add(modRayPeakHeight);		
// 		idx  = Stats.getIndex(modRayPSD[1], modRayPeakHeight);
// 		double modRayPeakLoc = modRayPSD[0][idx];
// 		modRayPeaks[i] = modRayPeakHeight;
// 		modRayLocs[i] = modRayPeakLoc;
// 		if ( modRayPeakHeight >= 13.215 )
// 		    nModRayPeaks++;


		/**  Find peak In the Lomb PSD (search range is already restricted) **/
// 		double lombPeakHeight = Stats.getMax(lombPSD[1]);
// 		lombPeaks.add(lombPeakHeight);
// 		idx  = Stats.getIndex(lombPSD[1], lombPeakHeight);
// 		double lombPeakLoc = lombPSD[0][idx];
// 		lombPeaks[i] = lombPeakHeight;
// 		lombLocs[i] = lombPeakLoc;
// 		if ( lombPeakHeight >= 6.608 )
// 		    nLombPeaks++;

		
// 		/**  Get median for peak heights  **/
// 		double[] copy = fftPeaks;
// 		Arrays.sort(copy);
// 		double fftPeakMedian = copy[copy.length/2 -1];

// 		copy = modRayPeaks;
// 		Arrays.sort(copy);
// 		double modRayPeakMedian = copy[copy.length/2 -1];


// 		/**  Get mean of Mod Ray locations  **/
// 		double[] aveAndVar = Stats.getRunningAveAndVar(modRayLocs);
// 		double ave = aveAndVar[0];
// 		double sig = Math.sqrt(aveAndVar[1]);


// 		/**  Determine histogram min and max bounds  **/
// 		double fftMinLoc = 999;
// 		double fftMaxLoc = -999;
// 		double modRayMinLoc = ave - 3.5*sig;
// 		double modRayMaxLoc = ave + 3.5*sig;
// 		double fftMinPeak = 0;
// 		double fftMaxPeak = fftPeakMedian + 8*Math.sqrt(fftPeakMedian);
// 		double modRayMinPeak = 0;
// 		double modRayMaxPeak = modRayPeakMedian + 8*Math.sqrt(modRayPeakMedian);
// 		for ( int n=0; n < nPulseFreqs; n++ ) {
// 		    fftMinLoc = Math.min(fftMinLoc, fftLocs[i]);
// 		    fftMaxLoc = Math.max(fftMaxLoc, fftLocs[i]);
// 		}
//   		fftMinLoc = fftMinLoc - 0.2;
//   		fftMaxLoc = fftMaxLoc + 0.2;


// 		/**  Construct JAIDA histograms  **/
// 		int peakBins = 13;
// 		int locBins = 13;
// 		IHistogram1D histoOfFFTPeaks = 
// 		    hf.createHistogram1D("FFT Peak Heights", peakBins, fftMinPeak, fftMaxPeak);
// 		IHistogram1D histoOfFFTLocs = 
// 		    hf.createHistogram1D("FFT Peak Locations", locBins, fftMinLoc, fftMaxLoc);
// 		IHistogram1D histoOfModRayPeaks = 
// 		    hf.createHistogram1D("R2 Peaks Heigts", peakBins, modRayMinPeak, modRayMaxPeak);
// 		IHistogram1D histoOfModRayLocs = 
// 		    hf.createHistogram1D("R2 Peak Locations", locBins, modRayMinLoc, modRayMaxLoc);
		
		
// 		/**  Fill histograms  **/
// 		for ( int n=0; n < nPulseFreqs; n++ ) {

// 		    histoOfFFTPeaks.fill(fftPeaks[i]);
// 		    histoOfFFTLocs.fill(fftLocs[i]);
// 		    histoOfModRayPeaks.fill(modRayPeaks[i]);
// 		    histoOfModRayLocs.fill(modRayLocs[i]);
// 		}


// 		/**  Construct the ChiSquare function to fit the distribution of peaks  **/
// 		IFunction chi2 = new ChiSquareFunction("Chi Square Function");
		

//  		/**  Perform the fit on FFT peak heights  **/
// 		double dof = histoOfFFTPeaks.mean();
// 		double norm = nPulseFreqs*(fftMaxPeak - fftMinPeak)/peakBins;
// 		chi2.setParameter("dof", dof);
// 		chi2.setParameter("norm", norm);
//  		IFitResult fftPeaksFitResult = fitter.fit(histoOfFFTPeaks, chi2);
//  		double fftFittedDOF = fftPeaksFitResult.fittedParameter("dof");
//  		double fftFittedDOF_sig = Math.sqrt(2*fftFittedDOF);

//  		//IFitResult fftLocsFitResult = fitter.fit(histoOfFFTLocs, "g");
//  		//double meanFFTLocs = fftLocsFitResult.fittedParameter("mean");
//  		//double sigmaFFTLocs = fftLocsFitResult.fittedParameter("sigma");


//  		/**  Perform the fit on Mod Ray peak heights  **/
// 		dof = histoOfModRayPeaks.mean();
// 		norm = nPulseFreqs*(modRayMaxPeak - modRayMinPeak)/peakBins;
// 		chi2.setParameter("dof", dof);
// 		chi2.setParameter("norm", norm);
//  		IFitResult modRayPeaksFitResult = fitter.fit(histoOfModRayPeaks, chi2);
//  		double modRayFittedDOF = modRayPeaksFitResult.fittedParameter("dof");
// 		double modRayFittedDOF_sig = Math.sqrt(2*modRayFittedDOF);


//  		/**  Perform the fit on Mod Ray peak location  **/		
//  		IFitResult modRayLocsFitResult = fitter.fit(histoOfModRayLocs, "g");
//  		double modRayFittedLoc = modRayLocsFitResult.fittedParameter("mean");
//  		double modRayFittedLoc_sig = modRayLocsFitResult.fittedParameter("sigma");


	
		/**  Display the results  **/

// 		plotter.clearRegions();

// 		plotter.region(0).plot(histoOfFFTPeaks);
//  		plotter.region(0).plot(fftPeaksFitResult.fittedFunction() );
//  		plotter.region(0).style().statisticsBoxStyle().setVisible(true);
		
// 		plotter.region(1).plot(histoOfFFTLocs);
//  		//plotter.region(1).plot(fftLocsFitResult.fittedFunction() );
//  		plotter.region(1).style().statisticsBoxStyle().setVisible(true);

// 		plotter.region(2).plot(histoOfModRayPeaks);
//  		plotter.region(2).plot(modRayPeaksFitResult.fittedFunction() );
//  		plotter.region(2).style().statisticsBoxStyle().setVisible(true);

// 		plotter.region(3).plot(histoOfModRayLocs);
//  		plotter.region(3).plot(modRayLocsFitResult.fittedFunction() );
//  		plotter.region(3).style().statisticsBoxStyle().setVisible(true);
		
// 		plotter.show();

                                 /**  Up to here  **/
		

// 		/**  Print results to file  **/
// 		fftPeakHeightsFile.print(fftFittedDOF +"\t"+ fftFittedDOF_sig +"\t");
// 		fftPeakLocsFile.print(histoOfFFTLocs.mean() +"\t"+histoOfFFTLocs.rms() +"\t");
// 		modRayPeakHeightsFile.print(modRayFittedDOF +"\t"+modRayFittedDOF_sig +"\t");
// 		modRayPeakLocsFile.print(modRayFittedLoc+"\t"+modRayFittedLoc_sig +"\t");


// 		/**  Strore the relevant results from each set of event lists  **/
// 		/**  Construct data histograms  **/
// 		double[] histoOfFFTPeaks = Binner.binData(fftPeaks, 25);
// 		double[] histoOfFFTLocs = Binner.binData(fftLocs, 25);
// 		double[] histoOfModRayPeaks = Binner.binData(modRayPeaks, 25);
// 		double[] histoOfModRayLocs = Binner.binData(modRayLocs, 25);


	    }
	    /** End of loop on pulsed fraction  **/
	    System.out.println("Log  :   Time (out) = "+DateUtils.time());	    


	    /**  Calculate probability from distribution of noise powers  **/
// 	    fftNoisePowers.sort();
	    modRayNoisePowers.sort();
// 	    lombNoisePowers.sort();
// 	    fftPeaks.sort();
	    modRayPeaks.sort();
// 	    lombPeaks.sort();

	    double threeSigmaProb = 0.00135;
	    double phiPercent = 1.0 - threeSigmaProb;
	    /**  Using the quantile  **/
// 	    double threeSigmaEquiv_fft = Descriptive.quantile(fftNoisePowers, phiPercent);
	    double threeSigmaEquiv_modRay = Descriptive.quantile(modRayNoisePowers, phiPercent);
// 	    double threeSigmaEquiv_lomb = Descriptive.quantile(lombNoisePowers, phiPercent);
	    //System.out.println(threeSigmaEquiv_fft+"\t"+threeSigmaEquiv_modRay+"\t"+threeSigmaEquiv_lomb);

// 	    double fftProb = 1.0 - Descriptive.quantileInverse(fftPeaks, threeSigmaEquiv_fft);
	    double r2Prob = 1.0 - Descriptive.quantileInverse(modRayPeaks, threeSigmaEquiv_modRay);
// 	    double lombProb = 1.0 - Descriptive.quantileInverse(lombPeaks, threeSigmaEquiv_lomb);
	    System.out.println("Log  :   3 sigma equivalent = "+number.format(threeSigmaEquiv_modRay));
	    System.out.println("Log  :   Estimated probability (from "+nPulseFreqs+") = "+number.format(r2Prob));

// 	    double counter = 0.0;
// 	    for ( int p=0; p < nPulseFreqs; p++ ) {
// 		double peak = modRayPeaks.get(p);
// 		if ( peak > threeSigmaEquiv_modRay ) {
// 		    counter++;
// // 		    System.out.println(peak+" > "+threeSigmaEquiv_modRay+" ("+counter+")");
// 		}
// 	    }
// 	    System.out.println("Log  : "+counter+" peaks above 3 sigma");
// 	    double probFromNPeaks = counter/nPulseFreqs;
// 	    System.out.println("Log  :   Estimated probability (from ratio) = "+probFromNPeaks);
			       

	    /**  Print results to file  **/
// 	    powerFile.println(freq.format(fftProb)+"\t"+freq.format(r2Prob)+"\t"+freq.format(lombProb));
	    powerFile.println(freq.format(r2Prob));
	    powerFile.flush();


	    /**  Clear all the DoubleArrayLists  **/
//  	    fftNoisePowers.clear();
	    modRayNoisePowers.clear();
// 	    lombNoisePowers.clear();
// 	    fftPeaks.clear();
	    modRayPeaks.clear();
// 	    lombPeaks.clear();	    

	}
	/**  End of loop on modulation freqs  **/

	powerFile.println("! \n HARD powerOfR2Test.ps/cps");
	powerFile.close();

    }


    public static void handleArguments(String[] args) {

	System.out.println("Usage: java PowerOfR2Test nPulseFreqs duration meanRate nCyclesMin nCyclesMax snrMin snrMax alpha");
                                                                 
	if ( args.length == 0 ) {
	    System.out.println("Log  : Running PowerOfR2Test");
	    System.out.println("Log  : Default nPulseFreqs = "+nPulseFreqs);
	    System.out.println("Log  : Default duration = "+ duration+" s");
	    System.out.println("Log  : Default mean rate = "+ meanRate+" cts/s");
	    System.out.println("Log  : Default nCycles min = "+ nCyclesMin);
	    System.out.println("Log  : Default nCycles max = "+ nCyclesMax);
	    System.out.println("Log  : Default s/n min = "+ snrMin);
	    System.out.println("Log  : Default s/n max = "+ snrMax);

	}
	else if ( args.length == 1 ) {
	    nPulseFreqs = (Integer.valueOf(args[0])).intValue();
	    System.out.println("Log  : Running PowerOfR2Test");
	    System.out.println("Log  : User defined nPulseFreqs ("+nPulseFreqs+")");
	    System.out.println("Log  : Default duration = "+ duration+" s");
	    System.out.println("Log  : Default mean rate = "+ meanRate+" cts/s");
	    System.out.println("Log  : Default nCycles min = "+ nCyclesMin);
	    System.out.println("Log  : Default nCycles max = "+ nCyclesMax);
	    System.out.println("Log  : Default s/n min = "+ snrMin);
	    System.out.println("Log  : Default s/n max = "+ snrMax);
	}
	else if ( args.length == 2 ) {
	    nPulseFreqs = (Integer.valueOf(args[0])).intValue();
	    duration = (Double.valueOf(args[1])).doubleValue();
	    System.out.println("Log  : Running PowerOfR2Test");
	    System.out.println("Log  : User defined nPulseFreqs ("+nPulseFreqs+")");
	    System.out.println("Log  : User defined duration ("+duration+" s)");
	    System.out.println("Log  : Default mean rate = "+ meanRate+" cts/s");
	    System.out.println("Log  : Default nCycles min = "+ nCyclesMin);
	    System.out.println("Log  : Default nCycles max = "+ nCyclesMax);
	    System.out.println("Log  : Default s/n min = "+ snrMin);
	    System.out.println("Log  : Default s/n max = "+ snrMax);
	}
	else if ( args.length == 3 ) {
	    nPulseFreqs = (Integer.valueOf(args[0])).intValue();
	    duration = (Double.valueOf(args[1])).doubleValue();
	    meanRate = (Double.valueOf(args[2])).doubleValue();
	    System.out.println("Log  : Running PowerOfR2Test");
	    System.out.println("Log  : User defined nPulseFreqs ("+nPulseFreqs+")");
	    System.out.println("Log  : User defined duration ("+duration+" s)");
	    System.out.println("Log  : User defined mean rate ("+meanRate+" cts/s)");
	    System.out.println("Log  : Default nCycles min = "+ nCyclesMin);
	    System.out.println("Log  : Default nCycles max = "+ nCyclesMax);
	    System.out.println("Log  : Default s/n min = "+ snrMin);
	    System.out.println("Log  : Default s/n max = "+ snrMax);
	}
	else if ( args.length == 5 ) {
	    nPulseFreqs = (Integer.valueOf(args[0])).intValue();
	    duration = (Double.valueOf(args[1])).doubleValue();
	    meanRate = (Double.valueOf(args[2])).doubleValue();
	    nCyclesMin = (Integer.valueOf(args[3])).intValue();
	    nCyclesMax = (Integer.valueOf(args[4])).intValue();
	    System.out.println("Log  : Running PowerOfR2Test");
	    System.out.println("Log  : User defined nPulseFreqs ("+nPulseFreqs+")");
	    System.out.println("Log  : User defined duration ("+duration+" s)");
	    System.out.println("Log  : User defined mean rate ("+meanRate+" cts/s)");
	    System.out.println("Log  : User defined nCycles min ("+nCyclesMin+")");
	    System.out.println("Log  : User defined nCycles max ("+nCyclesMax+")");
	    System.out.println("Log  : Default s/n min = "+ snrMin);
	    System.out.println("Log  : Default s/n max = "+ snrMax);
	}
	else if ( args.length == 7 ) {
	    nPulseFreqs = (Integer.valueOf(args[0])).intValue();
	    duration = (Double.valueOf(args[1])).doubleValue();
	    meanRate = (Double.valueOf(args[2])).doubleValue();
	    nCyclesMin = (Integer.valueOf(args[3])).intValue();
	    nCyclesMax = (Integer.valueOf(args[4])).intValue();
	    snrMin = (Double.valueOf(args[5])).doubleValue();
	    snrMax = (Double.valueOf(args[6])).doubleValue();
	    System.out.println("Log  : Running PowerOfR2Test");
	    System.out.println("Log  : User defined nPulseFreqs ("+nPulseFreqs+")");
	    System.out.println("Log  : User defined duration ("+duration+" s)");
	    System.out.println("Log  : User defined mean rate ("+meanRate+" cts/s)");
	    System.out.println("Log  : User defined nCycles min ("+nCyclesMin+")");
	    System.out.println("Log  : User defined nCycles max ("+nCyclesMax+")");
	    System.out.println("Log  : User defined s/n min ("+snrMin+")");
	    System.out.println("Log  : User defined s/n max ("+snrMax+")");
	}
	else if ( args.length == 8 ) {
	    nPulseFreqs = (Integer.valueOf(args[0])).intValue();
	    duration = (Double.valueOf(args[1])).doubleValue();
	    meanRate = (Double.valueOf(args[2])).doubleValue();
	    nCyclesMin = (Integer.valueOf(args[3])).intValue();
	    nCyclesMax = (Integer.valueOf(args[4])).intValue();
	    snrMin = (Double.valueOf(args[5])).doubleValue();
	    snrMax = (Double.valueOf(args[6])).doubleValue();
	    alpha = (Double.valueOf(args[7])).doubleValue();
	    if ( alpha < 0.0 || alpha > 4.0 ) {
		System.out.println("Error: Alpha must be between 0 and 3");
		System.exit(-1);
	    }
	    System.out.println("Log  : Running PowerOfR2Test");
	    System.out.println("Log  : User defined nPulseFreqs ("+nPulseFreqs+")");
	    System.out.println("Log  : User defined duration ("+duration+" s)");
	    System.out.println("Log  : User defined mean rate ("+meanRate+" cts/s)");
	    System.out.println("Log  : User defined nCycles min ("+nCyclesMin+")");
	    System.out.println("Log  : User defined nCycles max ("+nCyclesMax+")");
	    System.out.println("Log  : User defined s/n min ("+snrMin+")");
	    System.out.println("Log  : User defined s/n max ("+snrMax+")");
	    System.out.println("Log  : User defined red noise index ("+alpha+")");
	}
	else System.exit(-1);
    }


}

	    
	

