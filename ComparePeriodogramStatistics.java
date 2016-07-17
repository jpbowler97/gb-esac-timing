package gb.esac.timing;

import java.util.Arrays;
import java.util.Date;

import cern.colt.list.DoubleArrayList;
import cern.jet.random.engine.MersenneTwister64;
import cern.jet.stat.Descriptive;
import gb.esac.eventlist.EventList;
import gb.esac.eventlist.EventListSelector;
import gb.esac.io.AsciiDataFileWriter;
import gb.esac.likelihood.ExponentialLikelihood;
import gb.esac.montecarlo.RedNoiseGenerator;
import gb.esac.montecarlo.WhiteNoiseGenerator;
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.LombScarglePeriodogram;
import gb.esac.periodogram.ModifiedRayleighPeriodogram;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;
import org.apache.log4j.Logger;
import gb.esac.tools.DataUtils;
import gb.esac.binner.Binner;

public class ComparePeriodogramStatistics {

    private static Logger logger  = Logger.getLogger(ComparePeriodogramStatistics.class);
    // Arguments
    public static double duration = 1e4;
    public static double meanRate = 0.5;
    public static double period = duration/3.;
    public static double alpha = 0;
    public static double snrMin = 0;
    public static double snrMax = 10;
    public static double snrStep = 0.5;
    public static int nRuns = 50;
    // calculated from args
    public static double bigN = duration*meanRate;
    // Argument aliases
    private static double t = duration;
    private static double mu = meanRate;
    private static double p = period;
    private static double a = alpha;
    private static double min = snrMin;
    private static double max = snrMax;
    private static double d = snrStep;
    private static int n = nRuns;

    // Define here
    private static int sampling = 21;

    // MAIN
    public static void main(String[] args)  throws Exception {
	handleArguments(args);
	//evaluateDetectProbAsFuncOfSNR();
	evaluateProbOfFalsePositivesAsFuncOfSampling();
    }

    // Global variables used in several methods
    private static DoubleArrayList xAxis = new DoubleArrayList();
    private static DoubleArrayList fft_prob = new DoubleArrayList();
    private static DoubleArrayList fft_lower = new DoubleArrayList();
    private static DoubleArrayList fft_upper = new DoubleArrayList();	
    private static DoubleArrayList ls_prob = new DoubleArrayList();
    private static DoubleArrayList ls_lower = new DoubleArrayList();
    private static DoubleArrayList ls_upper = new DoubleArrayList();	
    private static DoubleArrayList r2_prob = new DoubleArrayList();
    private static DoubleArrayList r2_lower = new DoubleArrayList();
    private static DoubleArrayList r2_upper = new DoubleArrayList();

    private static void evaluateDetectProbAsFuncOfSNR() throws Exception {
	clearLists();
	double snr = snrMin;
	while ( snr <= snrMax ) {
	    double[] probWithBounds = monteCarloProb(snr);
	    addToLists(probWithBounds);
	    xAxis.add(snr);
	    snr += snrStep;
	}
	trimLists();
	String filename = "detectFracVsSNR_t_"+t+"_mu_"+mu+"_P_"+p+"_alpha_"+a+"_snr_"+min+"-"+max+"_step_"+d+"_n_"+n+".qdp";
	String xLabel = "Signal-to-Noise Ratio";
	String yLabel = "Detection Fraction (1-\\gb)";
	writeResultsToFile(filename, xLabel, yLabel);
    }

    private static void evaluateProbOfFalsePositivesAsFuncOfSampling() throws Exception {
	clearLists();
	double snr = 0;
	int sMin = 1;
	int sMax = 21;
	// Looping on sampling factor by modifying the value of this global variable used elsewhere
	sampling = sMin;
	while ( sampling <= sMax ) {
	    double[] probWithBounds = monteCarloProb(snr);
	    addToLists(probWithBounds);
	    xAxis.add(sampling);
	    sampling += 2;
	}	    
	trimLists();
	String filename = "falseDetectFrac_t_"+t+"_mu_"+mu+"_alpha_"+a+"_samp_"+sMin+"-"+sMax+"_n_"+n+".qdp";
	String xLabel = "IFS Sampling Factor";
	String yLabel = "Detection Fraction (1-\\ga)";
	writeResultsToFile(filename, xLabel, yLabel);
    }

    private static double[] monteCarloProb(double snr) throws Exception {
	double pulsedFraction = snr/Math.sqrt(bigN);
	logger.info("SNR value "+snr+" (pulsed fraction = "+pulsedFraction+")");
	double[] fftPeakHeights = new double[nRuns];
	double[] fftPeakLocs = new double[nRuns];
	double[] fftLikelihoodRatios = new double[nRuns];
	double[] lsPeakHeights = new double[nRuns];
	double[] lsPeakLocs = new double[nRuns];
	double[] lsLikelihoodRatios = new double[nRuns];
	double[] r2PeakHeights = new double[nRuns];
	double[] r2PeakLocs = new double[nRuns];
	double[] r2LikelihoodRatios = new double[nRuns];
	for ( int i=0; i < nRuns; i++ ) {
	    logger.info("Run number "+(i+1)+" (of "+nRuns+")");
	    double[] peakHeightsLocationsAndLRs = simulateObservation(duration, meanRate, alpha, period, pulsedFraction);
	    fftPeakHeights[i] = peakHeightsLocationsAndLRs[0];
	    fftPeakLocs[i] = peakHeightsLocationsAndLRs[1];
	    fftLikelihoodRatios[i] = peakHeightsLocationsAndLRs[2];
	    lsPeakHeights[i] = peakHeightsLocationsAndLRs[3];
	    lsPeakLocs[i] = peakHeightsLocationsAndLRs[4];
	    lsLikelihoodRatios[i] = peakHeightsLocationsAndLRs[5];
	    r2PeakHeights[i] = peakHeightsLocationsAndLRs[6];
	    r2PeakLocs[i] = peakHeightsLocationsAndLRs[7];
	    r2LikelihoodRatios[i] = peakHeightsLocationsAndLRs[8];
	}
	// Make histograms of everything for fixed snr of 5
	// if ( snr == 5. ) {
	//     makeHistograms(fftPeakHeights,fftPeakLocs,fftLikelihoodRatios,
	// 		   lsPeakHeights,lsPeakLocs,lsLikelihoodRatios,
	// 		   r2PeakHeights,r2PeakLocs,r2LikelihoodRatios);
	// }
	// Calculate and return detection likelihood intervals
	double[] fft_prob = getLikelihoodIntervals(fftLikelihoodRatios);
	double[] ls_prob = getLikelihoodIntervals(lsLikelihoodRatios);
	double[] r2_prob = getLikelihoodIntervals(r2LikelihoodRatios);
	return new double[] {fft_prob[0], fft_prob[1], fft_prob[2], 
			     ls_prob[0], ls_prob[1], ls_prob[2], 
			     r2_prob[0], r2_prob[1], r2_prob[2]
	};
    }

    private static double[] getLikelihoodIntervals(double[] likelihoods) {
	DoubleArrayList list = new DoubleArrayList(likelihoods);
	list.sort();
	//  These are probability ratios with respect to the mode at zero for the normal distribution
	double twoAndHalfSigma = 4.3937e-2;
	double threeSigma = 1.1109e-2;
	double threeAndHalfSigma = 2.1875e-3;
	double nEqualOrLessLikelyThanThreeSigma = Descriptive.rankInterpolated(list, threeSigma);
	double fractionOfDetections = nEqualOrLessLikelyThanThreeSigma/likelihoods.length;
	double nEqualOrLessLikelyThanThreeAndHalfSigma = Descriptive.rankInterpolated(list, threeAndHalfSigma);
	double fractionThreeAndHalfSigma = nEqualOrLessLikelyThanThreeAndHalfSigma/likelihoods.length;
	double nEqualOrLessLikelyThanTwoAndHalfSigma = Descriptive.rankInterpolated(list, twoAndHalfSigma);
	double fractionTwoAndHalfSigma = nEqualOrLessLikelyThanTwoAndHalfSigma/likelihoods.length;
	return new double[] {fractionOfDetections, fractionThreeAndHalfSigma, fractionTwoAndHalfSigma};
    }

    // Global variables used in simulateObservation()
    private static MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
    private static ExponentialLikelihood expL = new ExponentialLikelihood();
    private static PeriodogramMaker psdMaker = new PeriodogramMaker();
    private static TimeSeriesMaker tsMaker = new TimeSeriesMaker();
    private static WhiteNoiseGenerator whiteNoise = new WhiteNoiseGenerator();
    private static RedNoiseGenerator redNoise = new RedNoiseGenerator();

    private static double[] simulateObservation(double duration,double meanRate,double alpha,double period,double pulsedFraction) throws Exception {
	double t = duration;
	int nbins = 1024;
	double binTime = t/nbins;
	double nuMaxFFT = 1/(2*binTime);
	double nuMinFFT = 1/t;
	double ifs = nuMinFFT;
	double pulseFrequency = 1/period;
	double searchMin = pulseFrequency - 2.5*ifs;
	double searchMax = pulseFrequency + 2.5*ifs;
	double nEvents = duration*meanRate;
	double nPulsedEvents = pulsedFraction*nEvents;
	double pulsedMeanRate = nPulsedEvents/duration;
	int nNoiseEvents = (int)(nEvents - nPulsedEvents);
	// Generate arrival times: background and pulsed events
	double[] noiseTimes = redNoise.generateArrivalTimes(meanRate, duration, alpha, engine);
	double[] pulsedTimes = whiteNoise.generateModulatedArrivalTimes(pulsedMeanRate, duration, period, 0.999, engine);
	// Make periodograms of background
	EventList evlist = new EventList(noiseTimes);
	TimeSeries ts = tsMaker.makeTimeSeries(evlist, nbins);
	FFTPeriodogram fft = psdMaker.makePlainFFTPeriodogram(ts);
	LombScarglePeriodogram ls = psdMaker.makeLombScarglePeriodogram(ts, searchMin, searchMax, sampling);
	// Scale by factor of 2 so that powers are like Rayleigh powers
	//ls = (LombScarglePeriodogram) ls.scale(2.);	
	ModifiedRayleighPeriodogram r2 = psdMaker.makeModifiedRayleighPeriodogram(evlist, searchMin, searchMax, sampling);
	// Get background power values in the test frequency range
	DoubleArrayList lsNoisePowers = new DoubleArrayList();
	DoubleArrayList r2NoisePowers = new DoubleArrayList();
	DoubleArrayList fftNoisePowers = new DoubleArrayList();
	DoubleArrayList fftTestFreqs = new DoubleArrayList();
	double[] f = fft.getFreqs();
	double[] p = fft.getPowers();
	int idx = 0;
	//// skip frequencies below the first test frequency
	while ( f[idx] < searchMin ) idx++;
	//// store each of the frequencies in the test range
	while ( f[idx] <= searchMax ) {
	    fftNoisePowers.add(p[idx]);
	    fftTestFreqs.add(f[idx]);
	    //// find the index of the test frequency in the LS and R2 which have identical frequency sampling
	    double[] freqs = ls.getFreqs();
	    double[] pows = ls.getPowers();
	    //int k = DataUtils.getClosestIndexInSortedData(f[idx], freqs);
	    int k = idx*sampling;
	    //// store the background power values
	    //System.out.println("FFT test frequency and power (idx="+idx+"): "+f[idx]+"\t"+p[idx]);
	    lsNoisePowers.add(pows[k]);
	    //System.out.println("LS test frequency and power (idx="+k+"): "+freqs[k]+"\t"+pows[k]);
	    pows = r2.getPowers();
	    r2NoisePowers.add(pows[k]);
	    //System.out.println("R2 test frequency and power (idx="+k+"): "+freqs[k]+"\t"+pows[k]+"\n");
	    //// go to the next frequency in the FFT test frequencies
	    idx++;
	}
	fftNoisePowers.trimToSize();
	lsNoisePowers.trimToSize();
	r2NoisePowers.trimToSize();
	//// compute the average (the ML estimate) of these background powers at the test frequencies
	double avgFFTNoisePower = Descriptive.mean(fftNoisePowers);
	double avgLombNoisePower = Descriptive.mean(lsNoisePowers);
	double avgR2NoisePower = Descriptive.mean(r2NoisePowers);
	// Merge events lists by replacing background noise events by pulsed events
	double[] someNoiseTimes = EventListSelector.getRandomArrivalTimes(evlist, nNoiseEvents);
	int nTimes = pulsedTimes.length + someNoiseTimes.length;
	nTimes = nTimes - 2; // Because both lists will start at 0 and end at T; we must exclude these from one of the lists
	double[] times = new double[nTimes];
	if ( pulsedTimes.length <= 2 )
	    times = noiseTimes;
	else {
	    int n = 0;
	    for ( int k=1; k < pulsedTimes.length-1; k++ ) {
		times[k-1] = pulsedTimes[k];
		n++;
	    }
	    for ( int k=0; k < someNoiseTimes.length; k++ ) {
		times[k+n] = someNoiseTimes[k];
	    }
	    Arrays.sort(times);
	}
	// Make evlist, ts and periodograms for combined background and pulsed events
	evlist = new EventList(times);
	ts = tsMaker.makeTimeSeries(evlist, nbins);
	fft = psdMaker.makePlainFFTPeriodogram(ts);
	ls = psdMaker.makeLombScarglePeriodogram(ts, searchMin, searchMax, sampling);
	r2 = psdMaker.makeModifiedRayleighPeriodogram(evlist, searchMin, searchMax, sampling);
	// Get the values of the peaks heights, locations and likelihood ratios for each likelihood function
	//// for FFT
	idx = 0;
	f = fft.getFreqs();
	p = fft.getPowers();
	double fftPeakHeight = 0;
	double fftPeakLoc = 0;
	while ( f[idx] < searchMin ) idx++;
	while ( f[idx] <= searchMax ) {
	    if ( p[idx] >= fftPeakHeight ) {
		fftPeakHeight = p[idx];
		fftPeakLoc = f[idx];
	    }
	    idx++;
	}
	double[] backgroundPowers = fftNoisePowers.elements();
	double fftLikelihoodOfBackgroundMLE = expL.getLikelihood(avgFFTNoisePower, backgroundPowers);
	double fftLikelihoodOfPeak = expL.getLikelihood(fftPeakHeight, backgroundPowers);
	double fftLikelihoodRatio = fftLikelihoodOfPeak/fftLikelihoodOfBackgroundMLE;
	double logOfFFTLikelihoodRatio = -2*Math.log(fftLikelihoodRatio);
	//// for LS
	backgroundPowers = lsNoisePowers.elements();
	double lombPeakHeight = ls.getMaxPower();
	double lombPeakLoc = ls.getFreqAtMaxPower();
	double lombLikelihoodOfBackgroundMLE = expL.getLikelihood(avgLombNoisePower, backgroundPowers);
 	double lombLikelihoodOfPeak = expL.getLikelihood(lombPeakHeight, backgroundPowers);;
	double lombLikelihoodRatio = lombLikelihoodOfPeak/lombLikelihoodOfBackgroundMLE;
	double logOfLombLikelihoodRatio = -2*Math.log(lombLikelihoodRatio);
	//// for R2
	backgroundPowers = r2NoisePowers.elements();
	double r2PeakHeight = r2.getMaxPower();
	double r2PeakLoc = ls.getFreqAtMaxPower();
	double r2LikelihoodOfBackgroundMLE = expL.getLikelihood(avgR2NoisePower, backgroundPowers);
 	double r2LikelihoodOfPeak = expL.getLikelihood(r2PeakHeight, backgroundPowers);;
	double r2LikelihoodRatio = r2LikelihoodOfPeak/r2LikelihoodOfBackgroundMLE;
	double logOfR2LikelihoodRatio = -2*Math.log(r2LikelihoodRatio);
	//// Return the results
	return new double[] {fftPeakHeight, fftPeakLoc, fftLikelihoodRatio,
			     lombPeakHeight, lombPeakLoc, lombLikelihoodRatio,
			     r2PeakHeight, r2PeakLoc, r2LikelihoodRatio,
	};
    }

    private static void makeHistograms(double[] fft_heights, double[] fft_locs, double[] fft_lr,
				       double[] ls_heights, double[] ls_locs, double[] ls_lr,
				       double[] r2_heights, double[] r2_locs, double[] r2_lr) throws Exception {
	String[] names = new String[] {"fftPeakHeights","fftPeakLocs", "fftLR", "lsPeakHeight", "lsPeakLocs", "lsLR", "r2PeakHeights", "r2PeakLocs", "r2LR"};
    	int nbins = (int)Math.floor(Math.min(15., nRuns/10.));
	AsciiDataFileWriter[] out = new AsciiDataFileWriter[names.length];
	String[] filenames = new String[names.length];
	String[] xLabels = new String[names.length];
	for ( int i=0; i < names.length; i++ ) {
	    filenames[i] = "histo_"+names[i]+".qdp";
	    xLabels[i] = names[i];
	    out[i] = new AsciiDataFileWriter(filenames[i]);
	}
	out[0].writeHisto(Binner.makeHisto(fft_heights, nbins), xLabels[0]);
	out[1].writeHisto(Binner.makeHisto(fft_locs, nbins), xLabels[1]);
	out[2].writeHisto(Binner.makeHisto(fft_lr, nbins), xLabels[2]);
	out[3].writeHisto(Binner.makeHisto(ls_heights, nbins), xLabels[3]);
	out[4].writeHisto(Binner.makeHisto(ls_locs, nbins), xLabels[4]);
	out[5].writeHisto(Binner.makeHisto(ls_lr, nbins), xLabels[5]);
	out[6].writeHisto(Binner.makeHisto(r2_heights, nbins), xLabels[6]);
	out[7].writeHisto(Binner.makeHisto(r2_locs, nbins), xLabels[7]);
	out[8].writeHisto(Binner.makeHisto(r2_lr, nbins), xLabels[8]);
    }

    private static void writeResultsToFile(String filename, String xLabel, String yLabel) throws Exception {
	AsciiDataFileWriter out = new AsciiDataFileWriter(filename);
	String[] header = new String[] {
	    "DEV /XS",
	    "LAB T", "LAB F",
	    "TIME OFF",
	    "LINE ON",
	    "MA OFF",
	    "MA 17 ON 1", "MA 17 ON 4", "MA 17 ON 7",
	    "LS 4 ON 1", "LS 4 ON 4", "LS 4 ON 7", 
	    "CO 1 ON 1", "CO 1 ON 2", "CO 1 ON 3", 
	    "CO 4 ON 4", "CO 4 ON 5", "CO 4 ON 6", 
	    "CO 2 ON 7", "CO 2 ON 8", "CO 2 ON 9", 
	    "MA SIZE 1",
	    "LW 4", "CS 1.5",
	    "LAB X "+xLabel,
	    "LAB Y "+yLabel,
	    "LAB 2 VPOS 0.25 0.82 \"R\\u2\\d\\d1\\u\" JUST LEFT MA 45 CO 2 CS 1.5",
	    "LAB 3 VPOS 0.25 0.77 \"LS\" JUST LEFT MA 45 CO 4 CS 1.5",
	    "LAB 4 VPOS 0.25 0.72 \"FFT\" JUST LEFT MA 45 CO 1 CS 1.5",
	    "VIEW 0.2 0.1 0.8 0.9",
	    "SKIP SINGLE",
	    "R Y -0.05 1.05",
	    "!"
	};
	out.writeData( header, xAxis.elements(),
		      fft_prob.elements(), fft_lower.elements(), fft_upper.elements(), 
		      ls_prob.elements(), ls_lower.elements(), ls_upper.elements(), 
		      r2_prob.elements(), r2_lower.elements(), r2_upper.elements() );
    }

    private static void addToLists(double[] results) {
    	fft_prob.add(results[0]);
    	fft_lower.add(results[1]);
    	fft_upper.add(results[2]);
    	ls_prob.add(results[3]);
    	ls_lower.add(results[4]);
    	ls_upper.add(results[5]);
    	r2_prob.add(results[6]);
    	r2_lower.add(results[7]);
    	r2_upper.add(results[8]);
    }
 
    private static void trimLists() {
    	fft_prob.trimToSize();
    	fft_lower.trimToSize();
    	fft_upper.trimToSize();
    	ls_prob.trimToSize();
    	ls_lower.trimToSize();
    	ls_upper.trimToSize();
    	r2_prob.trimToSize();
    	r2_lower.trimToSize();
    	r2_upper.trimToSize();
    	xAxis.trimToSize();
    }

    private static void clearLists() {
    	fft_prob.clear();
    	fft_lower.clear();
    	fft_upper.clear();
    	ls_prob.clear();
    	ls_lower.clear();
    	ls_upper.clear();
    	r2_prob.clear();
    	r2_lower.clear();
    	r2_upper.clear();
    	xAxis.clear();
    }

    public static void handleArguments(String[] args) {
	if ( args.length == 1 ) {
	    duration = (Double.valueOf(args[0])).doubleValue();
	    logger.info("Running ComparePeriodogramStatistics");
	    logger.info("User-defined duration "+duration+" s");
	    logger.info("(Default) mean rate = ("+ meanRate+" cts/s)");
	    logger.info("(Default) period = ("+ period+")");
	    logger.info("(Default) alpha = ("+ alpha+")");
	    logger.info("(Default) snrMin = ("+ snrMin+")");
	    logger.info("(Default) snrMax = ("+ snrMax+")");
	    logger.info("(Default) snrStep = ("+ snrStep+")");
	    logger.info("(Default) nRuns = ("+ nRuns+")");
	}
	else if ( args.length == 2 ) {
	    duration = (Double.valueOf(args[0])).doubleValue();
	    meanRate = (Double.valueOf(args[1])).doubleValue();
	    logger.info("Running ComparePeriodogramStatistics");
	    logger.info("User-defined duration "+duration+" s");
	    logger.info("User-defined mean rate "+meanRate+" cts/s");
	    logger.info("(Default) period = ("+ period+")");
	    logger.info("(Default) alpha = ("+ alpha+")");
	    logger.info("(Default) snrMin = ("+ snrMin+")");
	    logger.info("(Default) snrMax = ("+ snrMax+")");
	    logger.info("(Default) snrStep = ("+ snrStep+")");
	    logger.info("(Default) nRuns = ("+ nRuns+")");
	}
	else if ( args.length == 3 ) {
	    duration = (Double.valueOf(args[0])).doubleValue();
	    meanRate = (Double.valueOf(args[1])).doubleValue();
	    period = (Double.valueOf(args[2])).doubleValue();
	    logger.info("Running ComparePeriodogramStatistics");
	    logger.info("User-defined duration "+duration+" s");
	    logger.info("User-defined mean rate "+meanRate+" cts/s");
	    logger.info("User-defined period = "+ period);
	    logger.info("(Default) alpha = ("+ alpha+")");
	    logger.info("(Default) snrMin = ("+ snrMin+")");
	    logger.info("(Default) snrMax = ("+ snrMax+")");
	    logger.info("(Default) snrStep = ("+ snrStep+")");
	    logger.info("(Default) nRuns = ("+ nRuns+")");
	}
	else if ( args.length == 4 ) {
	    duration = (Double.valueOf(args[0])).doubleValue();
	    meanRate = (Double.valueOf(args[1])).doubleValue();
	    period = (Double.valueOf(args[2])).doubleValue();
	    alpha = (Double.valueOf(args[3])).doubleValue();
	    logger.info("Running ComparePeriodogramStatistics");
	    logger.info("User-defined duration "+duration+" s");
	    logger.info("User-defined mean rate "+meanRate+" cts/s");
	    logger.info("User-defined period = "+ period);
	    logger.info("User-defined alpha = "+ alpha);
	    logger.info("(Default) snrMin = ("+ snrMin+")");
	    logger.info("(Default) snrMax = ("+ snrMax+")");
	    logger.info("(Default) snrStep = ("+ snrStep+")");
	    logger.info("(Default) nRuns = ("+ nRuns+")");
	}
	else if ( args.length == 5 ) {
	    duration = (Double.valueOf(args[0])).doubleValue();
	    meanRate = (Double.valueOf(args[1])).doubleValue();
	    period = (Double.valueOf(args[2])).doubleValue();
	    alpha = (Double.valueOf(args[3])).doubleValue();
	    nRuns = (Integer.valueOf(args[4])).intValue();
	    logger.info("Running ComparePeriodogramStatistics");
	    logger.info("User-defined duration "+duration+" s");
	    logger.info("User-defined mean rate "+meanRate+" cts/s");
	    logger.info("User-defined period = "+ period);
	    logger.info("User-defined alpha = "+ alpha);
	    logger.info("(Default) snrMin = ("+ snrMin+")");
	    logger.info("(Default) snrMax = ("+ snrMax+")");
	    logger.info("(Default) snrStep = ("+ snrStep+")");
	    logger.info("User-defined nRuns = "+ nRuns);
	}
	else if ( args.length == 8 ) {
	    duration = (Double.valueOf(args[0])).doubleValue();
	    meanRate = (Double.valueOf(args[1])).doubleValue();
	    period = (Double.valueOf(args[2])).doubleValue();
	    alpha = (Double.valueOf(args[3])).doubleValue();
	    snrMin = (Double.valueOf(args[4])).doubleValue();
	    snrMax = (Double.valueOf(args[5])).doubleValue();
	    snrStep = (Double.valueOf(args[6])).doubleValue();
	    nRuns = (Integer.valueOf(args[7])).intValue();
	    logger.info("Running ComparePeriodogramStatistics");
	    logger.info("User-defined duration "+duration+" s");
	    logger.info("User-defined mean rate "+meanRate+" cts/s");
	    logger.info("User-defined period = "+ period);
	    logger.info("User-defined alpha = "+ alpha);
	    logger.info("User-defined snrMin = "+ snrMin);
	    logger.info("User-defined snrMax = "+ snrMax);
	    logger.info("User-defined snrStep = "+ snrStep);
	    logger.info("User-defined nRuns = "+ nRuns);
	}
	else {
	    logger.error("Usage: java gb.esac.timing.ComparePeriodogramStatistics duration (meanRate [0.5]) (period [duration/3]) (alpha [0]) (snrMin snrMax snrStep [0. 10. 0.5]) (nRuns [50])");
	    System.exit(0);
	}
    }

}


/**


   private static double[] pulsedFracs; = new double[nSNRs];
   private static void definePulsedFractions(int n) {
   pulsedFracs = new double[n];
   double step = (snrMax - snrMin)/(n-1);
   double[] snr = new double[n];
   double bigN = duration*meanRate;
   logger.info(""+n+" values of s/n (pulsed fraction) will be tested:");
   //logger.info("    s/n     pf");
   for ( int i=0; i < n; i++ ) {
   snr[i] = snrMin +i*step;
   pulsedFracs[i] = snr[i]/Math.sqrt(bigN);  // snr = pf*sqrt(N)
   //logger.info("   "+number.format(snr[i])+"   "+number.format(pulsedFracs[i]));
   }
   }

   private static double[] pulseFreqs;
   private static void definePulsedFrequencies(int n) {
   pulseFreqs = new double[n];
   double nuMin = nCyclesMin/duration;
   double nuMax = nCyclesMax/duration;
   int nIFS = (new Double(Math.floor(duration*(nuMax - nuMin)))).intValue();
   logger.info("Range of modulation frequencies:");
   logger.info("nuMin ("+nCyclesMin+" cycles/duration) = "+freq.format(nuMin)+" Hz"); 
   logger.info("nuMax ("+nCyclesMax+" cycles/duration) = "+freq.format(nuMax)+" Hz"); 
   logger.info("contains "+nIFS+" IFS, over the period range ["+number.format(1/nuMax)+", "+number.format(1/nuMin)+"] s");
   ////  Uniformly distributed
   //logger.info("Defining "+nPulseFreqs+" uniformly spaced freqs in ["+freq.format(nuMin)+", "+freq.format(nuMax)+"] Hz");
   // 			   
   //for ( int i=0; i < nPulseFreqs; i ++ ) {
   //    pulseFreqs[i] = nuMin + i*(nuMax-nuMin)/(nPulseFreqs-1);
   //}
   ////  Randomly distributed
   logger.info("Drawing "+nPulseFreqs+" randomly distributed freqs in ["+freq.format(nuMin)+", "+freq.format(nuMax)+"] Hz");
   MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
   Uniform randUni = new Uniform(nuMin, nuMax, engine);
   for ( int i=0; i < nPulseFreqs; i ++ ) {
   pulseFreqs[i] = randUni.nextDouble();
   }
   Arrays.sort(pulseFreqs);
   // 	logger.info("Test frequencies:");
   // 	for ( int i=0; i < nPulseFreqs; i ++ ) {
   // 	    logger.info((i+1)+") "+freq.format(pulseFreqs[i])+" Hz \t"+
   // 			       number.format(1/pulseFreqs[i])+" s");
   //      }
   }

   private static IAnalysisFactory af;
   private static ITree tree;
   private static IHistogramFactory hf;
   private static IFunctionFactory funcF;
   private static IFitFactory fitF;
   private static IFitter fitter;
   private static IPlotter plotter;
   private static void setupJAIDA() {
   af = IAnalysisFactory.create();
   tree = af.createTreeFactory().create();
   hf = af.createHistogramFactory(tree);
   funcF  = af.createFunctionFactory(tree);
   fitF   = af.createFitFactory();
   fitter = fitF.createFitter("Chi2", "jminuit");
   // Plotter
   plotter = af.createPlotterFactory().create("Peak Heights and Positions");
   plotter.createRegions(2,2);
   IPlotterStyle plotterStyle = plotter.style();
   IDataStyle dataStyle = plotterStyle.dataStyle();
   IFillStyle fillStyle = dataStyle.fillStyle();
   ILineStyle lineStyle = dataStyle.lineStyle();
   IMarkerStyle markerStyle = dataStyle.markerStyle();
   lineStyle.setParameter("color", "black");
   markerStyle.setParameter("shape", "dot");
   markerStyle.setParameter("size", "5");
   markerStyle.setParameter("color", "black");
   fillStyle.setParameter("color", "white");
   dataStyle.setLineStyle(lineStyle);
   dataStyle.setFillStyle(fillStyle);
   dataStyle.setMarkerStyle(markerStyle);
   }
	    
	

*/
