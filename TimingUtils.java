package gb.esac.timing;


import gb.esac.binner.BinningException;
import gb.esac.binner.Rebinner;
import gb.esac.eventlist.EventListException;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesException;
import gb.esac.timeseries.TimeSeriesFileException;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.tools.Converter;
import gb.esac.tools.DataUtils;
import nom.tam.fits.*;
import org.apache.log4j.Logger;


/**
 *
 * @version   July 2010 (last modified)
 * @author   Guillaume Belanger (ESAC, Spain)
 *
 */

public final class TimingUtils {

    static Logger logger  = Logger.getLogger(TimingUtils.class);

    public static double[] getFFTFrequencies(int nbins, double duration) throws TimingException {

	//    Check that nLCBins is a power of 2 >= nOldBins 
	double n1 = Math.floor(Math.log10(nbins)/Math.log10(2));
	double n2 = Math.log10(nbins)/Math.log10(2);
	if ( n1 != n2 ) {
	    new TimingException("nbins = "+nbins+" is not a power of 2");
	}

	//    Define the frequencies  
	double nuMin = 1/duration;
 	int nFreqs = nbins/2;
 	double[] frequencies = new double[nFreqs];
	for ( int i=0; i < nFreqs; i++ ) {
	    frequencies[i] = nuMin*(i+1);
	    //System.out.println(frequencies[i]);
	}

	return frequencies;
    }
    

    public static double[] getFourierTestPeriods(double _pmin, double _pmax, double duration, double oversampling, double _ifsOffset) {

	double numin = 1/_pmax;
	double numax = 1/_pmin;
	int nIFS = (new Double(Math.ceil(duration*(numax - numin)))).intValue() +1;
	double[] testFreqs = getFourierTestFrequencies(numin, numax, duration, oversampling, _ifsOffset); 
	int nFreqs = testFreqs.length;
	double[] testPeriods = new double[nFreqs];
	for ( int i=0; i < nFreqs; i++ ) {
	    testPeriods[i] = 1/testFreqs[nFreqs-1-i];
	}

	return testPeriods;
    }


    /**  The offset must be given as a fraction of an IFS in frequency space 
	 i.e offset = 0.5 means that each independent frequency will be 
	 (i/T) + offset*(i/T)
	 
	 Oversampling specifies the factor. For example, oversampling=3, will
	 result in 3 test frequencies per IFS
    **/
    public static double[] getFourierTestFrequencies(double nuMin, double nuMax, double duration, double oversampling, double _ifsOffset) {
	
	//    Determine the number of IFS in interval  
	int nIFS = getnIFS(duration, nuMin, nuMax);
	int nFreqs = (new Double(nIFS*oversampling)).intValue();

	
	//    Define and initialise array of test frequencies  
	double[] testFrequencies = new double[nFreqs];
	double freq = 0;
	int n = 0;
	int i = 1;
	double step = 1/(oversampling*duration);
	while ( n < nFreqs && freq < nuMax ) {
	    freq = (i + _ifsOffset)*step;
	    if ( freq >= nuMin && freq <= nuMax ) {
		testFrequencies[n] = freq;
		//System.out.println("YES \t"+i+"\t"+n+"\t"+testFrequencies[n]);
		n++;
	    }
	    //else 	    System.out.println("NO \t"+i+"\t"+n+"\t"+freq);

	    i++;
	}
	return testFrequencies;
    }


    public static int getnIFS(double duration, double nuMin, double nuMax) {
	
	return (new Double(Math.ceil(duration*(nuMax - nuMin)))).intValue();
    }


    public static double[] getPhases(double[] times, double period) {

	double[] phases = new double[times.length];
	double tOverP = 0;
	for ( int i=0; i < times.length; i++ ) {
	    tOverP = times[i]/period;
	    phases[i] = tOverP - Math.floor(tOverP);
	}
	//Arrays.sort(phases);
	return phases;
    }


    public static double[] getIFSEdgesInFreqSpace(double nuMin, double nuMax, double duration, double oversampling, double _ifsOffset) {

	//    Determine the number of IFS in interval  
	int nIFS = (new Double(Math.round(duration*(nuMax - nuMin)))).intValue();
	int nFreqs = (new Double((nIFS+1)*oversampling)).intValue();

	
	//    Define and initialise array of test frequencies  
	double[] ifsEdgesInFreqSpace = new double[nFreqs+1];
	double freq = 1/(oversampling*duration);
	int n = 0;
	int i = 2;
	while ( n < nFreqs ) {
	    freq = i/(oversampling*duration) + _ifsOffset/duration;
	    if ( freq >= nuMin && freq <= nuMax ) {
		ifsEdgesInFreqSpace[n] = freq + 0.5/(oversampling*duration);
		//System.out.println(n+" "+freq+" "+ifsEdgesInFreqSpace[n]);
		n++;
		i++;
	    }
	    else i++;
	}
	ifsEdgesInFreqSpace[n] = ifsEdgesInFreqSpace[n-1] + 1/(oversampling*duration);
	//System.out.println(n+" "+ifsEdgesInFreqSpace[n]); 

	return ifsEdgesInFreqSpace;
    }


    public static double[] getIFSEdgesInPeriodSpace(double _pmin, double _pmax, double duration, double oversampling, double _ifsOffset) {

	double numin = 1/_pmax;
	double numax = 1/_pmin;
	double[] ifsEdgesInFreqSpace = getIFSEdgesInFreqSpace(numin, numax, duration, oversampling, _ifsOffset);
	int nEdges = ifsEdgesInFreqSpace.length;
	double[] ifsEdgesInPeriodSpace = new double[nEdges];
	for ( int i=0; i < nEdges; i++ ) {
	    ifsEdgesInPeriodSpace[i] = 1/ifsEdgesInFreqSpace[nEdges-1 - i];
	    //System.out.println(i+" "+ifsEdgesInPeriodSpace[i]);
	}
	return ifsEdgesInPeriodSpace;
    }




    public static void applyFracexpCorrection(double[] ratePerFrame, double[] fracexpPerFrame) {

	int nbins = ratePerFrame.length;
	for ( int i=0; i < nbins; i++ ) {
	    ratePerFrame[i] *= fracexpPerFrame[i];
	}
    }


    public static Object[] getFracexpCorrectedLC(BinaryTableHDU eventsHDU, BinaryTableHDU expoHDU, double lcBinTime) throws FitsException, BinningException, EventListException, TimeSeriesFileException, TimeSeriesException  {

	//    Get data from the event file  
	double[] eventTimes = (double[]) eventsHDU.getColumn("TIME");
	double[] frameTime = (double[]) expoHDU.getColumn("TIME");
	float[] frameTimedel_f = (float[]) expoHDU.getColumn("TIMEDEL");
	float[] frameFracexp_f = (float[]) expoHDU.getColumn("FRACEXP");
	double[] frameTimedel = Converter.float2double(frameTimedel_f);
	double[] frameFracexp = Converter.float2double(frameFracexp_f);

	
	//    Get some time stats  
 	double tfirst = frameTime[0] - 0.5*frameTimedel[0];
 	double tlast = frameTime[frameTime.length-1] + 0.5*frameTimedel[frameTimedel.length-1];
 	double totalTime = tlast - tfirst;
	logger.info("tfirst = "+tfirst);
	logger.info("tlast = "+tlast);
	logger.info("total = "+totalTime);


	//    Define the frame edges  
	double[] newFrameTime = DataUtils.resetToZero(frameTime, frameTimedel[0]/2);
	frameTime = newFrameTime;
	int nframes = frameTime.length;
	double[] frameEdges = new double[2*nframes];
	double[] rightEdges = new double[nframes];
	for ( int i=0; i < nframes; i++ ) {
	    frameEdges[2*i] = frameTime[i] - 0.5*frameTimedel[i];
	    frameEdges[2*i+1] = frameTime[i] + 0.5*frameTimedel[i]; 
	    rightEdges[i] = frameEdges[2*i+1];
	}


	//    Calculate the gaps between frames (if any) 
	double[] gaps = new double[nframes];
	gaps[0] = 0;
	double deadTime = 0;
	for ( int i=1; i < nframes; i++ ) {
	    gaps[i] = frameEdges[2*i] - frameEdges[2*i-1];
	    deadTime += gaps[i];
	}
	logger.info("deadtime = "+deadTime+" ("+(deadTime/totalTime)+" %)");


	//  Rebin the fracexp
	boolean directSum = true;
	double[] binnedFracexp = Rebinner.rebinRates(frameFracexp, frameEdges, lcBinTime, directSum);


	//    Bin the arrival times to construct the raw light curve  
	int nbins = binnedFracexp.length;
	TimeSeries rawLC = TimeSeriesMaker.makeTimeSeries(eventTimes, nbins);
	double[] lcBinCentres = rawLC.getBinCentres();
	double[] rates = rawLC.getRates();
	double[] errors = rawLC.getErrorsOnRates();


	//    Apply the fracexp correction  
	double[] fracexpCorrRates = new double[nbins];
	double[] fracexpCorrErrors = new double[nbins];
	for ( int i=0; i < nbins; i ++ ) {
	    fracexpCorrRates[i] = rates[i]/binnedFracexp[i];
	    fracexpCorrErrors[i] = errors[i]/binnedFracexp[i];
	}

	return new Object[] {lcBinCentres, fracexpCorrRates, fracexpCorrErrors};
    }


    public static Object[] getFracexpGTICorrectedLC(BinaryTableHDU eventsHDU, BinaryTableHDU expoHDU, BinaryTableHDU gtiHDU, double lcBinTime) throws FitsException, BinningException, EventListException, TimeSeriesFileException, TimeSeriesException {

	//    Get the fractional exposure corrected  light curve  	
	Object[] lc = getFracexpCorrectedLC(eventsHDU, expoHDU, lcBinTime);
	double[] lcTimes = (double[]) lc[0];
	double[] rates = (double[]) lc[1];
	double[] errors = (double[]) lc[2];

	//    Apply the GTI correction  
	double[] gti_start = (double[]) gtiHDU.getColumn("START");
	double[] gti_stop = (double[]) gtiHDU.getColumn("STOP");
	double[][] gtiCorrected = gtiCorrect(rates, lcTimes, gti_start, gti_stop);

	double[] gtiCorrRates = gtiCorrected[0];
	double[] gtiCorrErrors = gtiCorrected[1];

	return new Object[] {lcTimes, gtiCorrRates, gtiCorrErrors};
    }


    public static double[][] gtiCorrect(double[] _rate, double[] _lcTime, double[] _gtiStart, double[] _gtiStop) {

	double[] corrRate = new double[_rate.length];
	double[] corrError = new double[_rate.length];

	int k = 0;
	int l = 0;
	double bti = 0;

	for ( int j=0; j < _lcTime.length; j++ ) {

	    // initialize eff_lcTimedel as the nominal bin width
	    double lcTimedel = _lcTime[1] - _lcTime[0];
	    double eff_lcTimedel = lcTimedel;
		    

	    // Check if there is a stop in the bin
	    if ( _gtiStop[k] <= _lcTime[j] ) {

		try {

		    // Keep looping while there still are stops or starts in this bin
		    while ( _gtiStart[k+1] <= _lcTime[j] || _gtiStop[k] <= _lcTime[j] && k < _gtiStart.length ) {

			// If the next start after a stop is in the next bin
			// then the bti is from the stop to the end of this bin
			// Otherwise the bti is the diff between the stop and the next start
			if ( _gtiStart[k+1] >= _lcTime[j] ) bti += _lcTime[j] - _gtiStop[k];
			else bti += (_gtiStart[k+1] - _gtiStop[k]);
			k++;
		    }
		}
		catch (ArrayIndexOutOfBoundsException e) {}

		eff_lcTimedel -= bti;
		bti = 0;
	    }

	    if ( j > 0 && _gtiStart[k] >= _lcTime[j] ) 
		bti = _gtiStart[k] - _lcTime[j];


	    // Calculate the ratio effective bin size to nominal bin size
	    double gtiCorrFactor = eff_lcTimedel/lcTimedel;

	    //  Apply correction to rate and calculate error
	    corrRate[j]  = _rate[j]/gtiCorrFactor;
	    corrError[j] = Math.sqrt(_rate[j]/lcTimedel)/gtiCorrFactor;

// 	    //  This is temporary and return rates and errors without corrections
// 	    corrRate[j]  = _rate[j];
// 	    corrError[j] = Math.sqrt(_rate[j]/lcTimedel);


	    // Sum the exposure fraction within the bin
	    // 		    int n = 0;
	    // 		    double expoTimeInBin = 0;
	    // 		    float   sumOfFracexp = 0;
	    // 		    while ( expoTime[l] <= _lcTime[j] && l < expoTime.length ) {
	    // 			if ( l==0 ) expoTimeInBin = 0;
	    // 			else expoTimeInBin += expoTime[l] - expoTime[l-1];
	    // 			sumOfFracexp += fracexp[l];
	    // 			l++;
	    // 			n++; 
	    // 		    }
	    // 		    System.out.print(" sumOfFracexp =  "+sumOfFracexp+" ");	
		

	    // Calculate final correction factor
	    // 		    double avgFracexp = 0;
	    // 		    if ( sumOfFracexp != 0 && sumOfFracexp != Double.NaN ) 
	    // 			avgFracexp = sumOfFracexp/n;

	    // 		    double expoFrac = expoTimeInBin/lcTimedel;
	    // 		    double totalExpoFrac = expoFrac*avgFracexp;
	    // 		    double corrFactor = gtiCorrFactor*totalExpoFrac;

	    //  		    double corrFactor = gtiCorrFactor;

	    // if ( corrFactor < 0.4 ) corrFactor = 0.5;

	}
	double[][] corrRateAndError = new double[][]{corrRate, corrError};
	return corrRateAndError;
    }


    public static double[] gtiCorrRate(double[] _rate, double[] _lcTime, double[] _gtiStart, double[] _gtiStop) {

	double[] corrRate = new double[_rate.length];

	int k = 0;
	int l = 0;
	double bti = 0;


	for ( int j=0; j < _lcTime.length; j++ ) {

	    // initialize eff_lcTimedel as the nominal bin width
	    double lcTimedel = _lcTime[1] - _lcTime[0];
	    double eff_lcTimedel = lcTimedel;
		    

	    // Check if there is a stop in the bin
	    if ( _gtiStop[k] <= _lcTime[j] ) {

		try {

		    // Keep looping while there still are stops or starts in this bin
		    while ( _gtiStart[k+1] <= _lcTime[j] || _gtiStop[k] <= _lcTime[j] && k < _gtiStart.length ) {

			// If the next start after a stop is in the next bin
			// then the bti is from the stop to the end of this bin
			// Otherwise the bti is the diff between the stop and the next start
			if ( _gtiStart[k+1] >= _lcTime[j] ) bti += _lcTime[j] - _gtiStop[k];
			else bti += (_gtiStart[k+1] - _gtiStop[k]);
			k++;
		    }
		}
		catch (ArrayIndexOutOfBoundsException e) {}

		eff_lcTimedel -= bti;
		bti = 0;
	    }

	    if ( j > 0 && _gtiStart[k] >= _lcTime[j] ) bti = _gtiStart[k] - _lcTime[j];


	    // Calculate the ratio effective bin size to nominal bin size
	    double gtiCorrFactor = eff_lcTimedel/lcTimedel;

	    corrRate[j]  = _rate[j]/gtiCorrFactor;

		    
	}
	return corrRate;

    }

    public static double[] gtiCorrError(double[] _rate, double[] _lcTime, double[] _gtiStart, double[] _gtiStop) {

	double[] corrRate = new double[_rate.length];
	double[] corrError = new double[_rate.length];

	int k = 0;
	int l = 0;
	double bti = 0;


	for ( int j=0; j < _lcTime.length; j++ ) {

	    // initialize eff_lcTimedel as the nominal bin width
	    double lcTimedel = _lcTime[1] - _lcTime[0];
	    double eff_lcTimedel = lcTimedel;
		    

	    // Check if there is a stop in the bin
	    if ( _gtiStop[k] <= _lcTime[j] ) {

		try {

		    // Keep looping while there still are stops or starts in this bin
		    while ( _gtiStart[k+1] <= _lcTime[j] || _gtiStop[k] <= _lcTime[j] && k < _gtiStart.length ) {

			// If the next start after a stop is in the next bin
			// then the bti is from the stop to the end of this bin
			// Otherwise the bti is the diff between the stop and the next start
			if ( _gtiStart[k+1] >= _lcTime[j] ) bti += _lcTime[j] - _gtiStop[k];
			else bti += (_gtiStart[k+1] - _gtiStop[k]);
			k++;
		    }
		}
		catch (ArrayIndexOutOfBoundsException e) {}

		eff_lcTimedel -= bti;
		bti = 0;
	    }

	    if ( j > 0 && _gtiStart[k] >= _lcTime[j] ) bti = _gtiStart[k] - _lcTime[j];


	    // Calculate the ratio effective bin size to nominal bin size
	    double gtiCorrFactor = eff_lcTimedel/lcTimedel;
	    
	    corrRate[j]  = _rate[j]/gtiCorrFactor;
	    corrError[j] = Math.sqrt(_rate[j]/lcTimedel)/gtiCorrFactor;


	}
	return corrError;
    }


    public static double getGtiSum(double _tfirst, double _tlast, double[] _gtiStart, double[] _gtiStop) {

	double timebin = _tlast - _tfirst;
	double eff_timebin = timebin;
	double bti = 0;

	int k = 0;
	while ( _gtiStart[k] < _tfirst ) k++;

	if ( _gtiStart[k] < _tfirst ) k--;

	if ( _gtiStop[k] <= _tlast ) {

	    try {
		// Keep looping while there still are stops or starts
		while ( _gtiStart[k+1] <= _tlast || _gtiStop[k] <= _tlast && k < _gtiStart.length ) {    

		    // If the next start after a stop is in the next bin
		    // then the bti is from the stop to the end of this bin
		    // Otherwise the bti is the diff between the stop and the next start
		    if ( _gtiStart[k+1] >= _tlast ) bti += _tlast - _gtiStop[k];
		    else bti += (_gtiStart[k+1] - _gtiStop[k]);
		    k++;
		}
	    }
	    catch (ArrayIndexOutOfBoundsException e) {}
	    eff_timebin -= bti;
	    bti = 0;
	}
	return eff_timebin;
    }


}
