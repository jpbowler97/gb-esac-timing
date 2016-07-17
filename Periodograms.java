package gb.esac.timing;

import gb.esac.tools.Stats;
import gb.esac.binner.Binner;
import gb.esac.binner.BinningException;
import gb.esac.tools.Complex;
import gb.esac.tools.FFT;
import gb.esac.tools.Analysis;
import gb.esac.tools.MyFFT;
import gb.esac.tools.ComplexNumbers;

import org.apache.log4j.Logger;


public final class Periodograms {


    static Logger logger  = Logger.getLogger(Periodograms.class);	

	
    public static double[][] makePSD_FFT(double[] times, int nLCbins, String normType) throws TimingException , BinningException {
		
		
	/**  Check that nLCbins is a power of 2 >= nOldBins **/
	double n1 = Math.floor(Math.log10(nLCbins)/Math.log10(2));
	double n2 = Math.log10(nLCbins)/Math.log10(2);
	if ( n1 != n2 ) {
	    System.out.println("Error [makePSD_FFT(double[], int, String)]: nbins = "+nLCbins+" is not a power of 2");
	    System.exit(-1);
	}
		
		
	/**  Construct light curve **/
	Object[] lc = TimingUtils.makeLC(times, nLCbins);
	double[] t = (double[]) lc[0];
	double[] r = (double[]) lc[1];
		
	return makePSD_FFT(t, r, normType);
    }
	
	
    public static double[][] makePSD_FFT(double[] times, double[] rate, String normType) throws TimingException, BinningException  {
		
		
	int nOldBins = times.length;
	double[] t = new double[nOldBins];
	double[] r = new double[nOldBins];
	for ( int i=0; i < nOldBins; i++ ) {
	    t[i] = times[i];
	    r[i] = rate[i];
	}
		
	/**  Determine closest power of 2 for nNewBins  **/
	double n = Math.floor(Math.log10(nOldBins)/Math.log10(2));
	int nNewBins = (new Double(Math.pow(2, n))).intValue();
	double oldBinTime = times[1] - times[0];
	double duration = nOldBins*oldBinTime;
	double newBinTime = duration/nNewBins;
	// 	System.out.println("Log  : Old bintime (times[1] - times[0]) = "+number.format(oldBinTime)+" s");
	// 	System.out.println("Log  : New bintime (duration/nNewBins) = "+number.format(newBinTime)+" s");
		
	return makePSD_FFT(t, r, nNewBins, normType);
    }
	
	
    public static double[][] makePSD_FFT(double[] times, double[] rates, int nLCbins, String normType) throws TimingException, BinningException  {
		
		
	int nOldBins = times.length;
	double[] t = new double[nOldBins];
	double[] r = new double[nOldBins];
	for ( int i=0; i < nOldBins; i++ ) {
	    t[i] = times[i];
	    r[i] = rates[i];
	}
		
	/**  Check that nLCbins is a power of 2 >= nOldBins **/
	double n1 = Math.floor(Math.log10(nLCbins)/Math.log10(2));
	double n2 = Math.log10(nLCbins)/Math.log10(2);
	if ( n1 != n2 ) {
	    System.out.println("Error [makePSD_FFT]: nLCbins = "+nLCbins+" is not a power of 2");
	    System.exit(-1);
	}
	if ( nLCbins > nOldBins ) {
	    System.out.println("Error: [makePSD_FFT]: nLCbins = "+nLCbins+" > nOldBins ("+nOldBins+")");
	    System.exit(-1);
	}
		
	/**  Check the normalization type  **/
	if ( !normType.equals("leahy") && !normType.equals("miyamoto") 
	     && !normType.equals("variance") && !normType.equals("leahy-like") ) {
	    System.out.println("Error [makePSD_FFT]: Normalisation must be 'leahy', 'miyamoto', 'variance' or 'leahy-like'");
	    System.exit(-1);
	}
		
	/**  Determine nuMin and newBinTime  **/
	// 	System.out.println("Log  : Building Power Spectrum ...");
	double oldBinTime = t[1] - t[0];
	double duration = nOldBins*oldBinTime;
	double nuMin = 1/duration;
	// 	System.out.println("Log  : Min test frequency and step (1/T) = "+ freq.format(nuMin) +" Hz");
	// 	System.out.println("Log  : Original bin time = "+number.format(oldBinTime)+" s");
	double newBinTime = duration/nLCbins;
	// 	double nuMax = 1/(2*newBinTime);
	// 	int nIFS = TimingUtils.getnIFS(nLCbins);
	// 	System.out.println("Log  : Closest power of 2 for power spectral bins = "+nSpecBins);
	// 	System.out.println("Log  : New bin time = "+number.format(newBinTime)+" s");
	// 	System.out.println("Log  : Max test frequency (1/2*newBinTime) = "+freq.format(nuMax)+" Hz");
	// 	System.out.println("Log  : Number of IFS in this range = "+nIFS);
		
		
	/**  Rebin light curve and subtract mean count rate  **/
	LightCurveMaker lcMaker = new LightCurveMaker();
	LightCurve lc = lcMaker.makeLightCurve(t, r);
	double meanCountRate = lc.getMean();
	if ( nLCbins < nOldBins )
	    lc.rebin(nLCbins);
	lc.subtractMean();
	double[] binCentres = lc.getTimes();
	double[] binHeights = lc.getRates();
	double[] binHeightsErr = lc.getErrors();
	if ( lc.thereAreGapsInRates )
	    lc.fillGapsInRates("average");

		
// 	/**  Fill data gaps if any with average from either side of it  **/
// 	for ( int i=0; i < nLCbins; i++ ) {
// 	    if ( Double.isNaN(binHeights[i]) ) {
// 		int j = i+1;
// 		if ( j < nLCbins ) {
// 		    while ( j < nLCbins && Double.isNaN(binHeights[j]) )  j++;
// 		    binHeights[i] = (binHeights[i-1] + binHeights[j])/2;
// 		    binHeightsErr[i] = (binHeightsErr[i-1] + binHeightsErr[j])/2;
// 		}
// 	    }
// 	}
		

	/**  Make FFT with MyFFT.java of binned LC  **/
	int nn = nLCbins;
	double[] fftBinHeights = MyFFT.fft(ComplexNumbers.myComplex(binHeights), nn, +1);
	double[] power = ComplexNumbers.getPower(fftBinHeights);
		

	/**  Drop second half of power spectrum corresponding to negative frequencies **/
	int size = nn/2 +1;
	double[] pow = new double[size];
	for ( int i=0; i < size; i++ ) {
	    pow[i] = power[i];
	}
		

// 	/**  Make FFT with FFT.java of binned LC  **/
//  	Complex[] fftBinHeights2 = FFT.fft(complex);
// 	for ( int i=0; i < size; i++ ) {	
// 	    System.out.println(fftBinHeights[i] - fftBinHeights2[i].real());
// 	}
		
// 	/**  Make Power spectrum **/
// 	int nFreqs = testFreqs.length;
// 	//System.out.println(nFreqs);
// 	double[] pow = new double[nFreqs];
// 	for ( int i=0; i < nFreqs; i++ ) {
// 	    pow[i] = fftBinHeights[i+1].pow();
// 	}
		
		
	/**  Normalise the power  spectrum  
	     Notes:
	     1) In the Miyamoto normalization, the integrated PSD, i.e., 
	     the area under the PSD estimated by the periodogramme in a given
	     frequency range, defines the square of the total rms variability, 
	     i.e. the fractional amount by which the lightcurve is sinusoidally
	     modulated in the given frequency range (RMS = sigma/mean). 
	     The given rms is the square root of the integrated PSD, mutiplied 
	     by the mean flux of the corresponding lightcurve segments.
	**/
		
		
	/**  Calculate the integral of the power spectrum  **/
	//System.out.println("Log  : Normalising PSD ...");
	//double varOfLc = Stats.getVariance(binHeights);
	//double meanIndivErrorSquared = Math.pow(Stats.getMean(binHeightsErr), 2);
	double freqBinWidth = nuMin;
	double totalPower = 0;
	//System.out.println("Log  : Variance of light curve = "+freq.format(varOfLc)+" (cts/s)^2");
	for ( int i=0; i < pow.length; i++ ) {
	    totalPower += pow[i]*freqBinWidth;
	    //System.out.println(nuCentres[i] +"\t"+ power[i]);
	}
	double avePower = totalPower/pow.length;
	//System.out.println("Log  : Total power in spectrum = "+freq.format(totalPower));
		
		
	/**  Define the different  normalisations  **/
	double norm = 0;
	double varNorm = 2.0/freqBinWidth;
	//double varNorm = varOfLc/totalPower;
	//double rmsNorm = (2*duration)/Math.pow(meanCountRate*nLCbins, 2);
	//double miyamotoNorm = ((varOfLc - meanIndivErrorSquared)/Math.pow(meanCountRate, 2)) / totalPower;
	//double miyamotoNorm = (varOfLc/Math.pow(meanCountRate, 2)) / totalPower;
	double miyamotoNorm = 2.0/(freqBinWidth*Math.pow(meanCountRate, 2));
	double leahyNorm = 2.0/(freqBinWidth*meanCountRate);
	double leahyLikeNorm = 2*freqBinWidth/avePower;
	if ( normType.equals("leahy") ) {
	    norm = leahyNorm;
	}
	else if ( normType.equals("miyamoto") ) {
	    norm = miyamotoNorm;
	}
	else if ( normType.equals("variance") ) {
	    norm = varNorm;
	}
	else {
	    norm = leahyLikeNorm;
	}
		
	/**  Apply the normalization  **/
	for ( int i=0; i < pow.length; i++ ) {
	    pow[i] *= norm;
	}
		
	/**  Drop the first term of the spectrum  **/
	int nValues = pow.length -1;
	double[] finalPow = new double[nValues];
	for ( int i=0; i < nValues; i++ )
	    finalPow[i] = pow[i+1];
		
	/**  Determine the test frequencies  **/
	double[] testFreqs = TimingUtils.getFFTFrequencies(nLCbins, duration);
		
	return new double[][] { testFreqs, finalPow};
    }
	
	
    public static double[][] makePSD_ModRay(double[] times, double nuMin, double nuMax, double samplingFactor) {
		

	double[] t = new double[times.length];
	for ( int i=0; i < times.length; i++ )
	    t[i] = times[i];
		

	/**  Check nuMin  **/
	double duration = t[t.length-1] - t[0];
	double ifs = 1/duration;
	if ( nuMin < ifs ) {
	    System.out.println("Error [makePSD_ModRay]: nuMin < 1/duration");
	    System.exit(-1);
	}
		
	/**  Define the test frequencies  **/
	double ifsOffset = 0;
	double[] trialFreqs = TimingUtils.getFourierTestFrequencies(nuMin, nuMax, duration, samplingFactor, ifsOffset);
	int nTrials = trialFreqs.length;
		
		
	/**  Construct periodogram  **/
	double[] r2 = new double[nTrials];
	double period = 0;
	double[][] psd = new double[2][nTrials];
	for ( int i=0; i < nTrials; i++ ) {		
	    period = 1.0/trialFreqs[i];
	    r2[i] = Stats.getModRayleighPower(t, period);
	    psd[0][i] = trialFreqs[i];
	    psd[1][i] = r2[i];
	    //System.out.println((1/psd[0][i])+"\t"+psd[1][i]);
	}
		
	return psd;
    }
	
    public static double[][] makePSD_ModRay(double[] times, double[] rates, double[] errors, double nuMin, double nuMax, double samplingFactor) {
		
	/**  Assumptions
	     
	1) times[0] = firstHalfBin
		 
	**/
		
	double[] t = new double[times.length];
	for ( int i=0; i < times.length; i++ )
	    t[i] = times[i];

		
	/**  Check nuMin  **/
	int nBins = times.length;
	double binTime = t[1] - t[0];
	double duration = nBins*binTime;
	double ifs = 1/duration;
	if ( nuMin < ifs ) {
	    System.out.println("Error [makePSD_ModRay]: nuMin < 1/duration");
	    System.exit(-1);
	}
		
	/**  Define the test frequencies  **/
	double ifsOffset = 0;
	double[] trialFreqs = TimingUtils.getFourierTestFrequencies(nuMin, nuMax, duration, samplingFactor, ifsOffset);
	int nTrials = trialFreqs.length;
		
		
	/**  Construct periodogram  **/
	double[] r2 = new double[nTrials];
	double period = 0;
	double[][] psd = new double[2][nTrials];
	for ( int i=0; i < nTrials; i++ ) {		
	    period = 1.0/trialFreqs[i];
	    r2[i] = Stats.getModRayleighPower(t, rates, errors, period);
	    psd[0][i] = trialFreqs[i];
	    psd[1][i] = r2[i];
	    //System.out.println((1/psd[i][0])+"\t"+psd[i][1]);
	}
		
	return psd;
    }
	
	
    public static double[][] makePSD_Lomb(double[] times, double[] rates, double nuMin, double nuMax, double samplingFactor) {
		
	/**  Assumptions
	     
	1) times[0] = firstHalfBin
		 
	**/
		
	double[] t = times;
		
	/**  Define the test frequencies  **/
	double duration = t[t.length-1] - t[0];
	double ifsOffset = 0;
	double[] trialFreqs = TimingUtils.getFourierTestFrequencies(nuMin, nuMax, duration, samplingFactor, ifsOffset);
	int nTrials = trialFreqs.length;
		
		
	/**  Construct periodogram  **/
	double[] pow = new double[nTrials];
	double period = 0;
	double[][] psd = new double[2][nTrials];
	for ( int i=0; i < nTrials; i++ ) {
	    period = 1.0/trialFreqs[i];
	    pow[i] = Stats.getLombPower(t, rates, period);
	    psd[0][i] = trialFreqs[i];
	    psd[1][i] = pow[i];
	    //System.out.println((1/psd[0][i])+"\t"+psd[1][i]);
	}
		
	return psd;
    }
	
	
    public static double[][] makePSD_Lomb(double[] times, int nbins, double nuMin, double nuMax, double samplingFactor) {
		
	double[] t = Analysis.resetToZero(times);
		
	/**  Construct light curve **/
	double duration = t[t.length-1] - t[0];
	double binTime = duration/(new Double(nbins)).doubleValue();
	double halfBin = binTime/2.0;
	double[] rates = Binner.binData(t, nbins);
	double[] binCentres = new double[nbins];
	for ( int i=0; i < nbins; i++ ) {
	    rates[i] /= binTime;
	}
		
	double[][] psd = makePSD_Lomb(binCentres, rates, nuMin, nuMax, samplingFactor);
	return psd;
    }
	
	
    public static double[][] makePSD_Z2(double[] times, double nuMin, double nuMax, int nHarm, double sampling) {
		
	double[] t = times;
		
	/**  Define the test frequencies  **/
	double duration = t[t.length-1] - t[0];
	double ifsOffset = 0;
	double[] trialFreqs = TimingUtils.getFourierTestFrequencies(nuMin, nuMax, duration, sampling, ifsOffset);
	int nTrials = trialFreqs.length;
		
		
	/**  Construct periodogram  **/
	double[] z2 = new double[nTrials];
	double period = 0;
	double[][] psd = new double[2][nTrials];
	for ( int i=0; i < nTrials; i++ ) {
	    period = 1.0/trialFreqs[i];		
	    z2[i] = Stats.getZ2stats(t, period, nHarm)[2];
	    //z2[i] = Stats.getZ2Power(t, period, nHarm);
	    psd[0][i] = trialFreqs[i];
	    psd[1][i] = z2[i];
	}
		
	return psd;
    }
	
}
