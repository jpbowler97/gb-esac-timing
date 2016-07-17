package gb.esac.timing;

import cern.jet.random.Poisson;
import cern.jet.random.engine.MersenneTwister64;
import gb.esac.binner.Binner;
import gb.esac.binner.BinningException;
import gb.esac.io.AsciiDataFileWriter;
import gb.esac.montecarlo.MonteCarloException;
import gb.esac.montecarlo.TimmerKonig;
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.PeriodogramException;
import gb.esac.tools.BasicStats;
import gb.esac.tools.Complex;
import gb.esac.tools.Converter;
import gb.esac.tools.LeastSquaresFitter;
import gb.esac.tools.MinMax;
import gb.esac.tools.Stats;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Date;
import org.apache.log4j.Logger;


public class TestTimmer {

    static Logger logger  = Logger.getLogger(TestTimmer.class);

    public static void main (String[] args) throws IOException, TimingException, BinningException, MonteCarloException, PeriodogramException {

	DecimalFormat number = new DecimalFormat("0.0000");
	DecimalFormat sci = new DecimalFormat("0.0000E0");
	
	/**  Handle arguments  **/
	double index = 1;
	double duration = 10000;
	double meanRate = 2.5;
	double bintime = 120;
	int nspecs = 500;
	if ( args.length == 5 ) {
	    index = (Double.valueOf(args[0])).doubleValue();
	    duration = (Double.valueOf(args[1])).doubleValue();
	    meanRate = (Double.valueOf(args[2])).doubleValue();
	    bintime = (Double.valueOf(args[3])).doubleValue();
	    nspecs = (Integer.valueOf(args[4])).intValue();
	}
	else {
	    logger.info("Usage: java TestTimmer index duration meanRate bintime nSpecs");
	    System.exit(-1);
	}
	logger.info("Running TestTimmer");
	logger.info("Spectral index = "+index);
	logger.info("Duration = "+duration);
	logger.info("Number of spectra = "+nspecs);

	
	/**  Generate powspecs and get LS best fit values for the slopes  **/
	logger.info("Simulating power spectra ...");
	int nTKbins = 256;
	double tkBinTime = duration/nTKbins;
	double[] testFreqs = TimingUtils.getFFTFrequencies(nTKbins, duration);

	logger.debug("testFreqs.length = "+testFreqs.length);

	int nLCbins = (new Double(Math.pow(2, Math.ceil((Math.log(duration/bintime)/Math.log(2)))))).intValue();
	double newBinTime = duration/nLCbins;
	logger.info("new bin time = "+number.format(newBinTime)+" s");
	double[] newTestFreqs = TimingUtils.getFFTFrequencies(nLCbins, duration);

	logger.debug("newTestFreqs.length = "+newTestFreqs.length);

	double[] lcBinCentres = new double[nLCbins];
	double halfBin = 0.5*newBinTime;
	for ( int i=0; i < nLCbins; i++ )
	    lcBinCentres[i] = halfBin + (new Double(i)).doubleValue()*newBinTime;

	int nTKFreqs = testFreqs.length; // nTKbins/2
	int nLCFreqs = newTestFreqs.length; // nLCbins/2;
	double[] y  = new double[nTKFreqs];
	double[] x1  = new double[nTKFreqs];
	double[] x2  = new double[nLCFreqs];
	x1 = Converter.lin2logSpace(testFreqs);
	x2 = Converter.lin2logSpace(newTestFreqs);
	double[] tkPow = new double[nTKFreqs];
	double[] tkLC = null;
	double[] pow = null;
	double[] lc = new double[nLCbins];
	double[] lcPoisson = null;
	double[] tkSlopes = new double[nspecs];
	double[] tkSlopeErrs = new double[nspecs];
	double[] slopes = new double[nspecs];
	double[] slopeErrs = new double[nspecs];
	double[] lsFitResults = new double[4];
	Complex[] fourierComp = null;
	double[] aveTKPow = new double[nTKFreqs];
	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
	Poisson poisson = new Poisson(0, engine);
	for ( int i=0; i < nspecs; i++ ) {

	    //  Generate the Fourier components

	    logger.debug("nTKbins = "+nTKbins);

	    fourierComp = TimmerKonig.getTimmerFourierComponents(index, duration, nTKbins);

	    logger.debug("fourierComp.length = "+fourierComp.length);

	    tkPow = TimmerKonig.getTimmerPowers(fourierComp);

	    logger.debug("tkPow.length = "+tkPow.length);

	    for ( int p=0; p < nTKFreqs; p++ ) {
		aveTKPow[p] += tkPow[p];
		logger.debug(testFreqs[p] +"\t"+ tkPow[p]);
	    }
	    //  Fit in log space and store the results. NB:  y = a + b*x, method returns [a][b][err_a][err_b]

	    

	    y = Converter.lin2logSpace(tkPow);

	    logger.debug(x1.length+"\t"+y.length);

	    lsFitResults = LeastSquaresFitter.leastSquaresFitLine(x1, y);
	    tkSlopes[i] = -lsFitResults[1];
	    tkSlopeErrs[i] = lsFitResults[3];

// 	    //  Make LC and rebin
// 	    tkLC = TimmerKonig.getTimmerRates(fourierComp);
// 	    lc = Binner.rebinRatesSimple(tkLC, tkBinTime, newBinTime);


// 	    //  Make corresponding powspec  **/
// 	    try { 
// 		pow = Periodograms.makePSD_FFT(lcBinCentres, lc, "leahy")[1];
// 	    }
// 	    catch ( TimingException e ) {
// 		System.out.println("Error [TestTimmer]: Cannot make periodogram");
// 		System.exit(-1);
// 	    }

	    FFTPeriodogram psd = TimmerKonig.getTimmerFFTPeriodogram(fourierComp, duration);
	    pow = psd.getPowers();

	    /**  Fit in log space and store the results  **/
 	    y = Converter.lin2logSpace(pow);
	    lsFitResults = LeastSquaresFitter.leastSquaresFitLine(x2, y);
	    slopes[i] = -lsFitResults[1];
	    slopeErrs[i] = lsFitResults[3];
	}

	for ( int p=0; p < nTKFreqs; p++ ) 
	    aveTKPow[p] /= nspecs;
	
	/**  Write power spectrum  **/
	AsciiDataFileWriter outputFile1 = null;
	try { outputFile1 = new AsciiDataFileWriter("aveTKPow.qdp"); }
	catch (IOException e) {}
	String[] header_pow = new String[] {
	    "! QDP File",
	    "DEV /XS", "READ 1 2", "TIME OFF",
	    "LAB T", "LAB F", "LW 3", "CS 1.3",
	    "LOG ON",
	    "LAB X Frequency (Hz)", "LAB Y Power",
	    "VIEW 0.1 0.2 0.9 0.8",
	    "!"
	};
	try { outputFile1.writeData(header_pow, testFreqs, aveTKPow); }
	catch (IOException e) {}	 	



	/**  Bin the results and write the histos  **/
	double[][] slopesAndErrors = new double[][] { tkSlopes, tkSlopeErrs, slopes, slopeErrs };
	String[] names = new String[] {"tkSlopes", "tkSlopeErrs", "slopes", "slopeErrs"};
	int nHistoBins = 40;

	for ( int i=0; i < 2; i++ ) {

	    double[] s = slopesAndErrors[2*i];
	    double[] e = slopesAndErrors[2*i+1];

	    double minSlope = MinMax.getMin(s);
	    double maxSlope = MinMax.getMax(s);
	    double slopeBinWidth = (maxSlope - minSlope)/nHistoBins;
	    double[] histoOfSlopes = Binner.binData(s, nHistoBins)[0];
	    
	    double minErr = MinMax.getMin(e);
	    double maxErr = MinMax.getMax(e);
	    double errorBinWidth = (maxErr - minErr)/nHistoBins;
	    double[] histoOfErrors = Binner.binData(e, minErr, maxErr, nHistoBins)[0];


	    /**  Construct the axes for histos and Gaussian function  **/
	    double varS = BasicStats.getVariance(s);
	    double varE = BasicStats.getVariance(e);
	    double aveS = BasicStats.getMean(s);
	    double aveE = BasicStats.getMean(e);
	    double normS = nspecs*slopeBinWidth;
	    double normE = nspecs*errorBinWidth;
	    double[] slopeGauss = new double[nHistoBins];
	    double[] errorGauss = new double[nHistoBins];
	    double[] slopeAxis = new double[nHistoBins];
	    double[] errorAxis = new double[nHistoBins];
	    for ( int j=0; j < nHistoBins; j++ ) {
		slopeAxis[j] = minSlope + (0.5 + j)*slopeBinWidth;
		errorAxis[j] = minErr + (0.5 + j)*errorBinWidth;
		slopeGauss[j] = normS*Math.sqrt(1/(2*Math.PI*varS)) * 
		    Math.exp( -Math.pow((slopeAxis[j]-aveS), 2)/(2*varS) );
		errorGauss[j] = normE*Math.sqrt(1/(2*Math.PI*varE)) * 
		    Math.exp( -Math.pow((errorAxis[j]-aveE), 2)/(2*varE) );
	    }


	    /**  Compute Pearson's Chi2 test  for the Gaussians **/
	    logger.info("Mean of slopes = "+number.format(aveS));
	    logger.info("Standard deviation of slopes = "+number.format(Math.sqrt(varS)));
	    double chi2S = Stats.computeChi2Test(histoOfSlopes, slopeGauss);
	    
	    logger.info("Mean of uncertainties on slopes = "+number.format(aveE));
	    logger.info("Standard deviation of uncertainties = "+sci.format(varE));
	    double chi2E = Stats.computeChi2Test(histoOfErrors, errorGauss);


	    /**  Write the histograms  as QDP files  **/
	    String mean = number.format(aveS);
	    String stdDev = number.format(Math.sqrt(varS));
	    String chi2 = number.format(chi2S);
	    String dof = Integer.toString(nHistoBins-1);
	    String binN = Integer.toString(nspecs);
	    String[] header = new String[] {
		"DEV /XS",
		"READ 1 2 3",
		"TIME OFF",
		"LINE STEP ON 2",
		"LINE ON 3",
		"LW 3", "CS 1.3",
		"LAB T", "LAB F",
		"LAB X LS Fit Spectral Index",
		"LAB Y Entries Per Bin",
		"LAB 3 VPOS 0.175 0.7 \"\\gm = "+mean+"\" JUST LEFT",
		"LAB 4 VPOS 0.175 0.67 \"\\gs = "+stdDev+"\" JUST LEFT",
		"LAB 5 VPOS 0.175 0.64 \"N = "+binN+"\" JUST LEFT",
		"LAB 6 VPOS 0.825 0.7 \"\\gx\\u2\\d("+dof+") = "+chi2+"\" JUST RIGHT",
		"VIEW 0.1 0.2 0.9 0.8",
		"!",
		"! Mean = "+aveS,
		"! Std dev = "+Math.sqrt(varS),
		"! Pearson's Chi2 = "+chi2,
		"!"
	    };
	    AsciiDataFileWriter histoOfSlopesFile = new AsciiDataFileWriter("histoOf"+names[2*i]+".qdp");
	    histoOfSlopesFile.writeData(header, slopeAxis, histoOfSlopes, slopeGauss);
	
	    mean = number.format(aveE);
	    stdDev = number.format(Math.sqrt(varE));
	    chi2 = number.format(chi2E);
	    binN = Integer.toString(nspecs);
	    header = new String[] {
		"DEV /XS",
		"READ 1 2",
		"TIME OFF",
		"LINE STEP ON 2",
		"LINE ON 3",
		"LW 3",
		"CS 1.3",
		"LAB T",
		"LAB F",
		"LAB X LS Fit Spectral Index Standard Error",
		"LAB Y Entries Per Bin",	    
		"LAB 3 VPOS 0.175 0.7 \"\\gm = "+mean+"\" JUST LEFT",
		"LAB 4 VPOS 0.175 0.67 \"\\gs = "+stdDev+"\" JUST LEFT",
		"LAB 5 VPOS 0.175 0.64 \"N = "+binN+"\" JUST LEFT",
		"LAB 6 VPOS 0.825 0.7 \"\\gx\\u2\\d("+dof+") = "+chi2+"\" JUST RIGHT",
		"VIEW 0.1 0.2 0.9 0.8",
		"!",
		"! Mean = "+aveE,
		"! Std dev = "+Math.sqrt(varE),
		"! Pearson's Chi2 = "+chi2,
		"!"
	    };
	    AsciiDataFileWriter histoOfErrorsFile = new AsciiDataFileWriter("histoOf"+names[2*i+1]+".qdp");
	    histoOfErrorsFile.writeData(header, errorAxis, histoOfErrors, errorGauss);


	}
    }
}
