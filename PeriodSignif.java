package gb.esac.timing;

import gb.esac.tools.Analysis;
import gb.esac.tools.Stats;
import gb.esac.tools.Binner;
import gb.esac.tools.Converter;
import gb.esac.tools.Astro;
import gb.esac.tools.Complex;
import gb.esac.io.DataFile;

import java.io.IOException;
import java.io.FileNotFoundException;
import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.Vector;
import java.util.Arrays;

import nom.tam.fits.*;

import hep.aida.*;
import hep.aida.ref.histogram.*;

import cern.jet.stat.Descriptive;
import cern.jet.stat.Probability;
import cern.colt.list.DoubleArrayList;
import cern.jet.random.engine.MersenneTwister64;
import cern.jet.random.Poisson;
import cern.jet.random.Normal;


public class PeriodSignif {

    
    /**  Define global variables  **/
    public static Fits fitsFile;
    public static double oldBinWidth = 0;
    public static int nbins = 0;
    public static double[] times = null;
    public static double[] binWidths = null;
    public static double sumOfBinWidths = 0;
    public static double[] rate = null;
    public static double[] error = null;
    public static double[][] filledRateAndError = null;
    public static double[] signif = null;
    public static boolean noRateCol = true, noErrorCol = true;
    public static boolean dataTypeIsRate;
    public static boolean samplingIsEven;
    public static boolean fileIsFits;
    public static String dataFileName = null;
    public static int nSims = 50;
    public static DecimalFormat num = new DecimalFormat("0.00000");
    public static DecimalFormat freq = new DecimalFormat("0.000E00");


    /**   START  main  **/
    public static void main(String[] args) throws Exception {

	DecimalFormat number = new DecimalFormat("0.0000");
	DecimalFormat prob = new DecimalFormat("0.0000E00");
	DecimalFormat freq = new DecimalFormat("0.000E00");

	Runtime rt = Runtime.getRuntime();

	/**  HANDLE ARGUMENTS **/
	double alpha = 0;
	double pmin = 0;
	double pmax = 0;
	String testType = null;
	int nHarmonics = 1;
	int sampling = 1;
	int nSimu = 0;

	if ( args.length == 8  ) {
	    dataFileName = args[0];
	    alpha = (Double.valueOf(args[1])).doubleValue();
	    pmin = (Double.valueOf(args[2])).doubleValue();
	    pmax = (Double.valueOf(args[3])).doubleValue();
	    testType = args[4];
	    nHarmonics = (Integer.valueOf(args[5])).intValue();
	    sampling = (Integer.valueOf(args[6])).intValue();
	    nSimu = (Integer.valueOf(args[7])).intValue();
	    if ( ! testType.equals("z2") ) {
		System.out.println("Error: There are 8 args. Test type should be 'z2'");
		System.exit(-1);
	    }
	}
	else if ( args.length == 7  ) {
	    dataFileName = args[0];
	    alpha = (Double.valueOf(args[1])).doubleValue();
	    pmin = (Double.valueOf(args[2])).doubleValue();
	    pmax = (Double.valueOf(args[3])).doubleValue();
	    testType = args[4];
	    sampling = (Integer.valueOf(args[5])).intValue();
	    nSimu = (Integer.valueOf(args[6])).intValue();
	    if ( ! testType.equals("r2") && ! testType.equals("lomb") ) {
		System.out.println("Error: There are 7 args. Test type should be 'r2' or 'lomb'");
		System.exit(-1);
	    }
	}
	else if ( args.length == 6 ) {
	    dataFileName = args[0];
	    alpha = (Double.valueOf(args[1])).doubleValue();
	    pmin = (Double.valueOf(args[2])).doubleValue();
	    pmax = (Double.valueOf(args[3])).doubleValue();
	    testType = args[4];
	    nSimu = (Integer.valueOf(args[5])).intValue();
	    if ( ! testType.equals("fft") ) {
		System.out.println("Error: There are 6 args. Test type should be 'fft'");
		System.exit(-1);
	    }
	}
	else {
	    System.out.println
		("Usage: java PeriodSignif file alpha pmin pmax test (nHarmonics) (sampling) nSimu");
	    System.out.println("    z2: 'nHarmonics' and 'sampling' are required (8 args)");
	    System.out.println("    r2 or lomb: 'sampling' is required (7 args)");
	    System.out.println("    fft: no 'sampling' (6 args)");
	    System.exit(-1);
	}


	/**  Get data  **/
	Object[] lcData = TimingUtils.readLCFile(dataFileName);

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
	double equivVar = 0;
	double var = 0;
	double wVar = 0;
	if ( !dataTypeIsRate ) {
	    nevents = times.length;
	    meanRate = times.length/duration;
	    wMeanRate = meanRate;
	    equivMeanRate = meanRate;
	    System.out.println("Log  : Number of events = "+nevents);
	    System.out.println("Log  : Mean rate = "+num.format(meanRate)+" cts/s");
	}
	else {
	    meanRate = Stats.getMean(rate);
	    nevents = (new Double(meanRate*duration)).intValue();
	    wMeanRate = Stats.getWMean(rate, error);
	    meanIndivError = Stats.getMean(error);
	    meanSignalToNoise = meanRate/meanIndivError;
	    equivMeanRate = Math.pow(meanSignalToNoise, 2)/oldBinWidth;
	    var = Stats.getVariance(rate);
	    wVar = Stats.getWVariance(rate, error);
	    equivVar = Math.sqrt(equivMeanRate/wMeanRate) * wVar;
	    System.out.println("Log  : Mean rate = "+num.format(meanRate)+" cts/s");
	    System.out.println("Log  : Weighted mean rate = "+num.format(wMeanRate)+" cts/s");
	    System.out.println("Log  : Average number of events (from mean) = "+nevents);
	    System.out.println("Log  : Variance = "+number.format(var)+" (cts/s)^2");
	    System.out.println("Log  : Weighted variance = "+number.format(wVar)+" (cts/s)^2");
	    System.out.println("Log  : Equivalent mean rate = "+num.format(equivMeanRate)+" cts/s");
	    System.out.println("Log  : Equivalent variance = "+num.format(equivVar)+" (cts/s)^2");
	}


	/**  Define minPeriod and maxPeriod if these were not specified in arguments  **/
	double nuMin = 1/pmax;
	double nuMax = 1/pmin;
	if ( dataTypeIsRate ) 
	    nuMax = Math.min(1/pmin, rate.length/(2*duration));
	System.out.println("Log  : Minimum test period = "+ pmin +" s ("+freq.format(nuMax)+" Hz)");
	System.out.println("Log  : Maximum test period = "+ pmax +" s ("+freq.format(nuMin)+" Hz)");


 	/**  Define FFT parameters for the data  **/
	int fftBins = (new Double(Math.round(duration/(5/meanRate)))).intValue();
	double n = Math.round(Math.log10(fftBins)/Math.log10(2));
	fftBins = (new Double(Math.pow(2, n))).intValue();
	String norm = "leahy";
//  	double[][] powspec = new double[1][1];
//  	double[] fftPow = new double[1];
//    	double[][] powspec = Periodograms.makePSD_FFT(times, fftBins, norm);
// 	double[] fftFreqs = powspec[0];
//   	double[] fftPow = powspec[1];
// 	System.out.println("Log  : The FFT returned "+fftFreqs.length+" tested frequencies.");
	

	/**  Define test frequencies  **/
	int nIFS = (new Double(Math.floor(duration*(nuMax - nuMin)))).intValue();
	double ifsOffset = 0;
	double[] trialFreqs = TimingUtils.getFourierTestFrequencies(nuMin, nuMax, duration, sampling, ifsOffset);
	double[] fftFreqs = TimingUtils.getFFTFrequencies(fftBins, duration);
	int nTrials = trialFreqs.length;
	System.out.println("Log  : Number of IFS = "+nIFS);
	System.out.println("Log  : Sampling factor = "+ sampling);
	System.out.println("Log  : Test frequencies in range = "+ nTrials);


	/**  Reset the times to (0+dt/2)  **/
	if ( dataTypeIsRate ) {
	    double tzero = times[0] - 0.5*binWidths[0];
	    Analysis.resetToZero(times, tzero);
	}


	/**  Construct periodogram of data  **/
	double[] dataPow = new double[nTrials];
	double period = 0;

	if ( testType.equals("fft") ) {
		double[][] fftPSD = Periodograms.makePSD_FFT(times, fftBins, norm);
		trialFreqs = fftPSD[0];
		dataPow = fftPSD[1];
	}
	else {
	    for ( int i=0; i < nTrials; i++ ) {		

		period = 1.0/trialFreqs[i];
		if ( testType.equals("r2") ) {

		    if ( dataTypeIsRate )
			dataPow[i] = Stats.getModRayleighPower(times, rate, error, period);
		    else
			dataPow[i] = Stats.getModRayleighPower(times, period);
		}

		else if ( testType.equals("lomb") ) {
		
		    if ( dataTypeIsRate ) 
			dataPow[i] = Stats.getLombPower(times, rate, period);
		    else {
			System.out.println("Error: Lomb test can only be applied to binned data");
			System.exit(-1);
		    }
		}

		else  { //   if ( testType.equals("z2") )

		    if ( ! dataTypeIsRate ) {
			dataPow[i] = Stats.getZ2Power(times, period, nHarmonics);
		    }
		    else {
			System.out.println("Error: Z-test can only be applied to unbinned data");
			System.exit(-1);
		    }
		}
	    }
	}
	

	/***************************    SIMULATIONS    **************************/

	System.out.println("Log  : Running simulations  ... ");

	DoubleArrayList[] powerAtTestFreq = new DoubleArrayList[nTrials];
 	for ( int i=0; i < nTrials; i++ ) {
	    powerAtTestFreq[i] = new DoubleArrayList();
	}
  	DoubleArrayList peaks = new DoubleArrayList();
  	DoubleArrayList peaksIndices = new DoubleArrayList();
 	double[] t = null;
	double[] simRates = new double[nbins];
	double[] simErrors = new double[nbins];
	double[] binnedTimes = null;
	double[][] psd = null;
	double[] pow = null;

	double marker = 0;
	int intPart = 0;
	double sumOfAve = 0;
	double sumOfVar = 0;
	double maxPeak = 0;
	double maxPeakSigma = 0;
	double maxPeakFreq = 0;
	int maxPeakIndex = 0;

	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
	Poisson poisson = new Poisson(0, engine);
	Normal gaussian = new Normal(0, 1, engine);
	nevents = (new Double(equivMeanRate*duration)).intValue();
	int nTKbins =16384;  // 2^14 bins for the Timmer LC, and thus 2^13 for the Timmer powspec 
	double bintime = duration/nTKbins;
	double[] timmerLC = null;
	Histogram1D lcHisto, cdfHisto = null;
	double tzero = 0;
	

	/**  START FOR: Loop on nSimu  **/
	System.out.println("Log  : Event list "+1+" of "+nSimu);
	for ( int s=0; s < nSimu; s++ ) {

	    /**  Print out iteration  **/
	    marker = (s+1)/500D;
	    intPart = (new Double(Math.floor((s+1)/500D))).intValue();
	    if ( marker == intPart )  {
  		System.out.println("Log  : Event list "+(s+1)+" of "+nSimu);
//  		System.out.println("Log  : Event list "+(s+1)+" of "+nSimu+": Memory used = "+
// 			  number.format(rt.totalMemory()/1e6)+" Mb");
	    }

	    /** Generate red noise arrival times  **/
	    t = ArrivalTimes.generateRedArrivalTimes(equivMeanRate, duration, alpha);


  	    /**  Bootstrap  **/
// 	    t = Analysis.bootstrap(times, events);
	    
	    
	    /**  If input data is a lightcurve, use indentical sampling  **/
	    if ( dataTypeIsRate ) {

		/**  Bin the data according the binning of input light curve  **/
		binnedTimes = (double[]) Binner.binOrderedData(t, times, binWidths)[0];
		for ( int i=0; i < binnedTimes.length; i++ ) {
		    simRates[i] = binnedTimes[i]/binWidths[i];
		    simErrors[i] = Math.sqrt(binnedTimes[i])/binWidths[i];
		}
 		double ave = Stats.getMean(simRates);
 		double variance = Stats.getVariance(simRates);
		sumOfAve += ave;
		sumOfVar += variance;


		/**  Construct periodogram  **/
		if ( testType.equals("fft") ) {
		    System.out.println("Error: cannot use FFT with input light curve");
		    System.exit(-1);
		}
		else if ( testType.equals("r2") ) {
		    psd = Periodograms.makePSD_ModRay(times, simRates, simErrors, nuMin, nuMax, sampling);
		    pow = psd[1];
		}
		else if ( testType.equals("lomb") ) {
		    psd = Periodograms.makePSD_Lomb(times, simRates, nuMin, nuMax, sampling);
		    pow = psd[1];
		}
		else {
		    System.out.println("Error: cannot use Z2 with input light curve");
		    System.exit(-1);
		}

	    }
	    else {
		/**  Construct periodogram  **/
		if ( testType.equals("fft") ) {
		    psd = Periodograms.makePSD_FFT(t, fftBins, norm);
		    pow = psd[1];
		}
		else if ( testType.equals("r2") ) {
		    psd = Periodograms.makePSD_ModRay(t, nuMin, nuMax, sampling);
		    pow = psd[1];
		}
		else if ( testType.equals("lomb") ) {
		    psd = Periodograms.makePSD_Lomb(t, fftBins, nuMin, nuMax, sampling);
		    pow = psd[1];
		}
		else {
		    psd = Periodograms.makePSD_Z2(t, nuMin, nuMax, nHarmonics, sampling);
		    pow = psd[1];
		}
	    }


	    /**  Store power associated with each frequency and max of PSD  **/
	    maxPeak = 0;
	    for ( int i=0; i < pow.length; i++ ) {
		powerAtTestFreq[i].add(pow[i]);
		if ( pow[i] > maxPeak ) {
		    maxPeak = pow[i];
		    maxPeakIndex = i;
		}
	    }
	    peaks.add(maxPeak);
	    peaksIndices.add(maxPeakIndex);


	    /**  START IF:  Write results every 'marker' event lists  **/
	    if ( marker == intPart )  {

		/**  Construct the averge periodogram  **/
		double[] averagePower = new double[nTrials];
		double[] uncertainty = new double[nTrials];
		double[] threeSigmaEnv = new double[nTrials];
		double mean = 0;
		double variance = 0;
		for ( int i=0; i < nTrials; i++ ) {
		    mean = Descriptive.mean(powerAtTestFreq[i]);
		    averagePower[i] = Math.max(mean, 1.0);
		    variance = Descriptive.sampleVariance(powerAtTestFreq[i], mean);
		    uncertainty[i] = Math.sqrt(variance);
		    threeSigmaEnv[i] = mean + 3*uncertainty[i];
		}
		

		/**  Translate the max peak heights into probabilities  **/
		DoubleArrayList peaksInProb = new DoubleArrayList();
		DoubleArrayList peaksInSigma = new DoubleArrayList();
		int idx = 0;
		double peak = 0;
		for ( int i=0; i < s; i++ ) {
		    peak = peaks.get(i);
		    idx = (int) peaksIndices.get(i);
		    peaksInProb.add(Probability.chiSquareComplemented(averagePower[idx], peak));
		    peaksInSigma.add( (peak - averagePower[idx]) / uncertainty[idx] );
		}


		/**  Count the number of peaks that pass the 3 sigma env for each freq  **/
		DoubleArrayList valuesInProb = new DoubleArrayList();
		double[] dataChi2Probs = new double[nTrials];
		double[] dataPowInSigma = new double[nTrials];
		double value = 0;
		double valueInSigma = 0;
		int nPeaksAbove3sigma = 0;
		for ( int i=0; i < nTrials; i++ ) {
		    dataChi2Probs[i] = Probability.chiSquareComplemented(averagePower[i], dataPow[i]);
		    dataPowInSigma[i] = (dataPow[i] - averagePower[i]) / uncertainty[i];
		    for ( int j=0; j < intPart; j++ ) {
			value = powerAtTestFreq[i].get(j);
			valuesInProb.add(Probability.chiSquareComplemented(averagePower[i], value));
		    }
		}

		/**  Calculate probabilities  **/
		double[] probability = new double[nTrials];
		double[] probFromPeaks = new double[nTrials];
		double[] probFromChi2ProbsPeaks = new double[nTrials];
		double[] probFromChi2ProbsAll = new double[nTrials];
		double[] probFromSigmaValue = new double[nTrials];
		double[] chiSquareComplemented = new double[nTrials];
		double meanPeak = Descriptive.mean(peaks);
		peaks.sort();
		peaksInProb.sort();
		peaksInSigma.sort();
		valuesInProb.sort();
		for ( int i=0; i < nTrials; i++ ) {
		    powerAtTestFreq[i].sort();
		    probability[i] = (1.0 - Descriptive.quantileInverse(powerAtTestFreq[i], dataPow[i]));
		    probFromPeaks[i] = (1.0 - Descriptive.quantileInverse(peaks, dataPow[i]));
		    probFromChi2ProbsPeaks[i] = (Descriptive.quantileInverse(peaksInProb, dataChi2Probs[i]));
		    probFromSigmaValue[i] = (1.0- Descriptive.quantileInverse(peaksInSigma, dataPowInSigma[i]));
		    probFromChi2ProbsAll[i] = (Descriptive.quantileInverse(valuesInProb, dataChi2Probs[i]));
		}

	 
		/**  Write Probability vs. Period to file **/
		DataFile probabilityFile = new DataFile("probability.qdp");
		//double xmin = pmin - 50;
		//double xmax = pmax + 50;
		String[] header = new String[] {
		    "! QDP data file",
		    "DEV /XS",
		    "VIEW 0.1 0.2 0.9 0.8",
		    "LAB TITLE", "LAB FILE", "TIME OFF", 
		    "LW 3", "CSIZE 1.3", 
		    "LAB X Frequency (Hz)",
		    "LAB Y Probability",
		    "LOG ON",
		    //"R X "+Double.toString(xmin)+" "+Double.toString(xmax)
		    "VIEW 0.1 0.2 0.9 0.8",
		    "!"
		};
		probabilityFile.writeData(header, trialFreqs, probability, probFromPeaks, probFromChi2ProbsPeaks);
	    }
	    /**  END IF:  Write results every 500 event lists  **/

	}
	/**  END FOR:  Loop on nSimu  **/


	/**  START IF ( dataTypeIsRate )  **/
	if ( dataTypeIsRate ) {

	    /**  Print out average of averages and of variances  **/
	    System.out.println("Log  : Average of averages = "
			       +number.format(sumOfAve/nSimu)+" cts/s"
			       +" (data equivalent: "+number.format(equivMeanRate)+")");
	    
	    System.out.println("Log  : Average of variances = "
			       +number.format(sumOfVar/nSimu)+" (cts/s)^2"
			       +" (data equivalent: "+number.format(equivVar)+")");


	    /**  Write last simulated light curve  **/
	    DataFile sampleLC = new DataFile("sampleLC.qdp");
	    String[] header = new String[] {
		"! QDP data file", 
		"DEV /XS",
		"READ SERR 2", 
		"LAB TITLE", "LAB FILE", "TIME OFF",
		"LINE ON", "MA 17 ON", 
		"LW 3", "CSIZE 1.3", 
		"LAB X Time (s)", 
		"LAB Y Count rate (cts/s)",
		"VIEW 0.1 0.2 0.9 0.8",
		"!"
	    };
	    sampleLC.writeData(header, times, simRates, simErrors);
	}
	/**  END IF ( dataTypeIsRate )  **/


	/**  Construct the averge periodogram and calculate the 3 sigma envelope  **/
	double[] averagePower = new double[nTrials];
	double[] uncertainty = new double[nTrials];
	double[] threeSigmaEnv = new double[nTrials];
	double threeSigmaProb = 0.00135;
	double phiPercent = 1.0 - threeSigmaProb;
	for ( int i=0; i < nTrials; i++ ) {
	    double mean = Descriptive.mean(powerAtTestFreq[i]);
	    double variance = Descriptive.sampleVariance(powerAtTestFreq[i], mean);
	    averagePower[i] = Math.max(mean, 1.0);
	    uncertainty[i] = Math.sqrt(variance/(nSimu-1));
	    threeSigmaEnv[i] = averagePower[i] + 3*uncertainty[i];
	}


	/**  Translate the max peak heights into probabilities  **/
	DoubleArrayList peaksInProb = new DoubleArrayList();
//  	DoubleArrayList peaksInSigma = new DoubleArrayList();
	int idx = 0;
	double peak = 0;
	for ( int i=0; i < nSimu; i++ ) {
	    peak = peaks.get(i);
	    idx = (int) peaksIndices.get(i);
	    peaksInProb.add(Probability.chiSquareComplemented(averagePower[idx], peak));
//  	    peaksInSigma.add( (peak - averagePower[idx]) / uncertainty[idx] );
// 	    System.out.println(idx);
	}
	

	/**  Count the number of peaks that pass the 3 sigma env for each freq  **/
// 	DoubleArrayList valuesInProb = new DoubleArrayList();
	double[] dataChi2Probs = new double[nTrials];
// 	DoubleArrayList valuesInSigma = new DoubleArrayList();
//  	double[] dataPowInSigma = new double[nTrials];
	double value = 0;
	double valueInSigma = 0;
	int nPeaksAbove3sigma = 0;
	for ( int i=0; i < nTrials; i++ ) {
	    dataChi2Probs[i] = Probability.chiSquareComplemented(averagePower[i], dataPow[i]);
// 	    dataPowInSigma[i] = (dataPow[i] - averagePower[i]) / uncertainty[i];
	    for ( int j=0; j < nSimu; j++ ) {
		value = powerAtTestFreq[i].get(j);
// 		valuesInProb.add(Probability.chiSquareComplemented(averagePower[i], value));
// 		valueInSigma = (value - averagePower[i]) / uncertainty[i];
// 		valuesInSigma.add(valueInSigma);
// 		if ( value > threeSigmaEnv[i] ) {
// 		    nPeaksAbove3sigma++;
// 		}
	    }
	}
// 	System.out.println("Log  : There are "+nPeaksAbove3sigma+" peaks above 3 sigma.");
    

	/**  Write data periodogram with 3 sigma envelope **/
	DataFile qdpFile = new DataFile("dataPow.qdp");
	String[] header = new String[] {
	    "! QDP data file", 
	    "DEV /XS",
	    "LABEL TITLE", "LABEL FILE", "TIME OFF", 
	    "LINE ON", "LW 3", "CSIZE 1.3",
	    "LAB X Frequency (Hz)", 
	    "LAB Y Power",
	    "LOG X ON",
	    "VIEW 0.1 0.2 0.9 0.8"
	};
	qdpFile.writeData(header, trialFreqs, dataPow, threeSigmaEnv);


	/**  Write average periodogram  **/
	DataFile averagePSDFile = new DataFile("averagePSD.qdp");
	String yLabel = null;
	if ( testType.equals("r2") ) yLabel = "Modified Rayleigh Power";
	else if ( testType.equals("lomb") ) yLabel = "Lomb Power";
	else if ( testType.equals("z2") ) yLabel = "Z\\u2\\d Power";
	else yLabel = "Leahy Power";
	header = new String[] {
	    "! QDP data file", 
	    "DEV /XS", 
	    "READ SERR 2", 
	    "LAB TITLE", "LAB FILE", "TIME OFF",
	    "LINE ON", 
	    "! MA 17 ON", 
	    "LW 3", "CSIZE 1.3", 
	    "LAB X Frequency (Hz)", 
	    "LAB Y "+yLabel,
	    "LOG X ON",
	    "VIEW 0.1 0.2 0.9 0.8",
	    "!"
 	};
 	averagePSDFile.writeData(header, trialFreqs, averagePower, uncertainty);


	/**  Calculate probabilities  **/
	double[] probability = new double[nTrials];
	double[] probFromPeaks = new double[nTrials];
	double[] probFromChi2ProbsPeaks = new double[nTrials];
// 	double[] probFromChi2ProbsAll = new double[nTrials];
	double[] probFromSigmaValue = new double[nTrials];
	double[] chiSquareComplemented = new double[nTrials];
// 	double[] trialPeriods = new double[nTrials];
	double meanPeak = Descriptive.mean(peaks);
	peaks.sort();
	peaksInProb.sort();
// 	peaksInSigma.sort();
// 	valuesInProb.sort();
	for ( int i=0; i < nTrials; i++ ) {
	    powerAtTestFreq[i].sort();
	    probability[i] = (1.0 - Descriptive.quantileInverse(powerAtTestFreq[i], dataPow[i]));
	    probFromPeaks[i] = (1.0 - Descriptive.quantileInverse(peaks, dataPow[i]));
	    probFromChi2ProbsPeaks[i] = (Descriptive.quantileInverse(peaksInProb, dataChi2Probs[i]));
// 	    probFromSigmaValue[i] = (1.0- Descriptive.quantileInverse(peaksInSigma, dataPowInSigma[i]));
// 	    probFromChi2ProbsAll[i] = (Descriptive.quantileInverse(valuesInProb, dataChi2Probs[i]));
	}

	 
	/**  Write Probability vs. Period to file **/
	DataFile probabilityFile = new DataFile("probability.qdp");
	//double xmin = pmin - 50;
	//double xmax = pmax + 50;
	header = new String[] {
	    "! QDP data file",
	    "DEV /XS",
	    "VIEW 0.1 0.2 0.9 0.8",
	    "LAB TITLE", "LAB FILE", "TIME OFF", 
	    "LW 3", "CSIZE 1.3", 
	    "LAB X Frequency (Hz)",
	    "LAB Y Probability",
	    "LOG ON",
	    //"R X "+Double.toString(xmin)+" "+Double.toString(xmax)
	    "VIEW 0.1 0.2 0.9 0.8",
	    "!"
	};
 	probabilityFile.writeData(header, trialFreqs, probability, probFromPeaks, probFromChi2ProbsPeaks);
    }


//     public static void getData() throws Exception, IOException, FitsException {
	
// 	/**  Get data  **/
// 	if ( fileIsFits ) {        /**  FITS FORMAT  **/
	    
// 	    System.out.println("Log  : File is in FITS format");
// 	    fitsFile = Astro.openFits(dataFileName);
// 	    BinaryTableHDU hdu = null;

// 	    try {   /**  Event file  **/

// 		/**  HDU  **/
// 		hdu = Astro.getHDU(dataFileName, "EVENTS");
// 		System.out.println("Log  : File is event list (HDU name = EVENTS)");
// 		dataTypeIsRate = false;

// 		/**  TIME  **/
// 		try { times = (double[]) hdu.getColumn("TIME"); }
// 		catch ( FitsException e )
// 		    { System.out.println("Error: No TIME column found"); System.exit(-1); }
// 	    }
// 	    catch ( NullPointerException e ) {

// 		try {  /** Light curve  **/

// 		    /**  HDU  **/
// 		    hdu = Astro.getHDU(dataFileName, "RATE"); 
// 		    System.out.println("Log  : File is light curve (HDU name = RATE");
// 		    dataTypeIsRate = true;

// 		    /**  TIME  **/
// 		    try { times = (double[]) hdu.getColumn("TIME"); }
// 		    catch ( FitsException e1 )
// 			{ System.out.println("Error: No TIME column found"); System.exit(-1); }

// 		    /**  RATE  **/
// 		    try { rate = (double[]) hdu.getColumn("RATE"); noRateCol = false; }
// 		    catch (FitsException e1) { 
// 			try { rate = (double[]) hdu.getColumn("RATE1"); noRateCol = false; }
// 			catch ( ClassCastException e2 ) {
// 			    float[] rateFlt = (float[]) hdu.getColumn("RATE1"); noRateCol = false;
// 			    rate = Converter.float2double(rateFlt);
// 			}
// 			catch ( FitsException e3 )
// 			    { System.out.println("Error: No RATE column found"); System.exit(-1); }
// 		    }

// 		    /**  ERROR  **/
// 		    try { error = (double[]) hdu.getColumn("ERROR"); noErrorCol = false; }
// 		    catch (FitsException e1) {
// 			try { error = (double[]) hdu.getColumn("ERROR1"); }
// 			catch ( ClassCastException e2 ) {
// 			    float[] errorFlt = (float[]) hdu.getColumn("ERROR1");
// 			    error = Converter.float2double(errorFlt);

// 			    // Calculate significance in each bin
// 			    signif = new double[error.length];
// 			    for ( int i=0; i < error.length; i++ ) {
// 				signif[i] = rate[i]/error[i];
// 			    }
// 			}
// 			catch ( FitsException e3 ) { 
// 			    System.out.print("Warn : No ERROR column. Constructing error ... ");
// 			    for ( int i=0; i < times.length; i++ )
// 				error[i] = Math.sqrt(rate[i]*times[i])/times[i];
// 			    System.out.print("done");
// 			}
// 		    }
// 		}
// 		catch ( NullPointerException e2 ) {
// 		    System.out.println("Error: Cannot get HDU: no 'EVENTS' nor 'RATE'");
// 		    System.exit(-1);
// 		}
// 	    }
// 	}

// 	else {	    /**  ASCII FORMAT  **/

// 	    System.out.println("Log  : File is in ASCII format");
// 	    DataFile dataFile = new DataFile(dataFileName);

// 	    //  Find out how many cols there are
// 	    int ncols = 0;
// 	    try { 
// 		ncols = dataFile.getNumOfDataCols();
// 		if ( ncols == 1 ) {
// 		    System.out.println("Log  : File is event list (1 col = TIME)");
// 		    dataTypeIsRate = false;
// 		    times = dataFile.getDblCol(0);
// 		}
// 		else if ( ncols == 2 ) {
// 		    System.out.println("Log  : File is light curve (2 cols = TIME and RATE)");
// 		    dataTypeIsRate = true;
// 		    times = dataFile.getDblCol(0);
// 		    rate = dataFile.getDblCol(1);
// 		    System.out.print("Warn : No ERROR column. Constructing error ... ");
// 		    for ( int i=0; i < times.length; i++ )
// 			error[i] = Math.sqrt(rate[i]*times[i])/times[i];
// 		    System.out.print("done");
// 		}
// 		else if ( ncols == 3 ) {
// 		    System.out.println("Log  : File is light curve (3 cols = TIME, RATE, ERROR)");
// 		    dataTypeIsRate = true;
// 		    times = dataFile.getDblCol(0);
// 		    rate = dataFile.getDblCol(1);
// 		    error = dataFile.getDblCol(2);
// 		    nbins = rate.length;
// 		    double fixedBinWidth = times[1] - times[0];
// 		    binWidths = new double[nbins];
// 		    for ( int i=0; i < nbins; i++ ) {
// 			binWidths[i] = fixedBinWidth;
// 		    }
// 		    sumOfBinWidths = nbins*fixedBinWidth;
// 		    oldBinWidth = fixedBinWidth;
// 		}
// 		else if ( ncols == 4 ) {
// 		    System.out.println("Log  : File is light curve (4 cols = TIME, BINWIDTH, RATE, ERROR)");
// 		    dataTypeIsRate = true;
// 		    times = dataFile.getDblCol(0);
// 		    double[] halfBinWidths = dataFile.getDblCol(1);
// 		    nbins = halfBinWidths.length;
// 		    binWidths = new double[nbins];
// 		    for ( int i=0; i < nbins; i++ ) 
// 			binWidths[i] = 2*halfBinWidths[i];

// 		    for ( int i=0; i < nbins-1; i++ ) {
// 			double overlap = (times[i+1]-halfBinWidths[i+1]) - (times[i]+halfBinWidths[i]);
// 			if ( overlap < 0 ) {
// 			    binWidths[i] += overlap;
// 			    binWidths[i+1] += overlap;
// 			}
// 			sumOfBinWidths += binWidths[i];
// 		    }


// 		    /**  Check out if there are data gaps  **/
// 		    double[] deltaT = new double[times.length-1];
// 		    for ( int i=0; i < deltaT.length; i++ )   deltaT[i] = times[i+1] - times[i];
// 		    double dtMin = Stats.getMin(deltaT);
// 		    double dtMax = Stats.getMax(deltaT);
// 		    double dtAve = Stats.getMean(deltaT);
// 		    if ( dtMin != dtMax ) {
// 			System.out.println("Warn : Sampling is NOT even: There are data gaps.");
// 			System.out.println("Log  : Ave delta T = "+num.format(dtAve)+" s");
// 			System.out.println("Log  : Min delta T = "+dtMin+" s");
// 			System.out.println("Log  : Max delta T = "+dtMax+" s");
// 			samplingIsEven = false;
// // 			System.out.println("Log  : Filling data gaps with local running average");
// // 			filledRateAndError = Analysis.fillDataGaps(rate, error);
// 		    }			

// 		    /**  Figure out if bin time is variable  **/
// 		    double aveBinWidth = Stats.getMean(binWidths);
// 		    double minBinWidth = Stats.getMin(binWidths);
// 		    double maxBinWidth = Stats.getMax(binWidths);
// 		    if ( minBinWidth != maxBinWidth ) {
// 			System.out.println("Warn : Bin time is variable.");
// 			System.out.println("Log  : Ave bin time = "+aveBinWidth+" s");
// 			System.out.println("Log  : Min bin time = "+minBinWidth+" s");
// 			System.out.println("Log  : Max bin time = "+maxBinWidth+" s");
// 		    }
// 		    else { System.out.println("Log  : Integration time is uniform: "+maxBinWidth+" s"); }
// 		    System.out.println("Log  : Using ave delta T as old bin time");
// 		    oldBinWidth = dtAve;

// 		    /**  Get the rate and error cols  **/
// 		    rate = dataFile.getDblCol(2);
// 		    error = dataFile.getDblCol(3);

// 		    signif = new double[error.length];
// 		    for ( int i=0; i < error.length; i++ ) {
// 			signif[i] = rate[i]/error[i];
// 		    }

// 		}
// 		else {
// 		    System.out.println("Error: File format not recognised");
// 		    System.exit(-1);
// 		}


// 	    }
// 	    catch ( IOException e ) {
// 		System.out.println("Error: Problem reading ASCII data file");
// 		System.exit(-1);
// 	    }
// 	}

//     }


}
