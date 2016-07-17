package gb.esac.timing;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.Arrays;

import gb.esac.tools.Stats;
import gb.esac.tools.Binner;
import gb.esac.tools.Analysis;
import gb.esac.tools.Complex;
import gb.esac.tools.Converter;
import gb.esac.io.DataFile;

import cern.jet.random.engine.MersenneTwister64;
import cern.jet.random.Poisson;

import hep.aida.ref.histogram.*;


public class CompareTKImplementation {

    public static DecimalFormat number = new DecimalFormat("0.0000");
    public static DecimalFormat sci = new DecimalFormat("0.0000E0");

    public static void main(String[] args) throws IOException, TimingException  {

	/** Handle args  **/
	double spectralIndex=0;
	double duration=0;
	double meanRate=1;
	double bintime=0;
	int nLC=0;
	if ( args.length != 4 ) {
	    System.out.println("Usage: java CompareTKImplementation spectralIndex duration bintime nLC");
	    System.exit(-1);
	}
	else {
	    spectralIndex = (Double.valueOf(args[0])).doubleValue();
	    duration = (Double.valueOf(args[1])).doubleValue();
	    bintime = (Double.valueOf(args[2])).doubleValue();
	    nLC = (Integer.valueOf(args[3])).intValue();
	}

	/**  Determine closest power of 2 larger than nbins  **/
	int nLCbins = (new Double(Math.pow(2, Math.ceil((Math.log(duration/bintime)/Math.log(2)))))).intValue();
	double newBinTime = duration/nLCbins;
	bintime = newBinTime;
	System.out.println("Log  : bintime = "+bintime);
	

	/**  Initialise some variables **/
	int nTKbins = 16384;
	double lc_TK_bintime = duration/nTKbins;
	//System.out.println("Log  : nTKbins = "+nTKbins);
	//System.out.println("Log  : lc_TK_bintime = "+lc_TK_bintime);
	double[] slopes_scaled = new double[nLC];
	double[] slopeErrs_scaled = new double[nLC];
	double[] slopes_poisson = new double[nLC];
	double[] slopeErrs_poisson = new double[nLC];
	double[] slopes_evlist = new double[nLC];
	double[] slopeErrs_evlist = new double[nLC];
	double[] lsFitResults = new double[4];
	int nevents = (new Double(meanRate*duration)).intValue();
	//System.out.println("Log  : Average number of events = "+nevents);
	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
	Poisson poisson = new Poisson(nevents, engine);


	// Local to the for-loop
	double[] lc_TK, arrivalTimes, t = null;
	double tzero = 0;
	Histogram1D lcHisto, cdfHisto = null;
	Object[] lc_evlist = null;
	double[] lc_binCentres, lc_binHeights = null;
	double[] meanAndVar = null;
	double lc_evlist_mean, lc_evlist_var, lc_TK_mean, lc_TK_var = 0;
	double[] lc_TK_rebinned = null;
	double[] lc_TK_scaled = new double[nLCbins];
	double[] lc_TK_poisson = new double[nLCbins];
	double lc_TK_rebinned_mean, lc_TK_rebinned_var = 0; 
	double lc_TK_scaled_mean, lc_TK_scaled_var = 0;
	double lc_TK_poisson_mean, lc_TK_poisson_var = 0;
	double[][] lcScaled_powspec, lcPoisson_powspec, evlist_powspec = null;
	double[] freqs, lcScaled_pow, lcPoisson_pow, evlist_pow, x, y1, y2, y3 = null;



	//double[] countRates = new double[] {0.2};
	double[] countRates = new double[] {0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100};
	PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("compareTKImplementations.qdp")));
	String[] header = new String[] {
	    "DEV /XS",
	    "VIEW  0.1 0.2 0.9 0.8",
	    "READ SERR 2 3 4",
	    "TIME OFF", "LAB T", "LAB F",
	    "LINE STEP ON 2 4",
	    "LINE ON", "LW 3", "CS 1.1",
	    "MARK 4 ON 2",
	    "MARK 6 ON 3",
	    "MARK 7 ON 4",
	    "MARK SIZE 2.0",
	    "LAB 2 VPOS 0.6 0.4 MA 4 \"scaled\" CO 2 CS 1.3 JUST LEFT",
	    "LAB 3 VPOS 0.6 0.35 MA 6 \"scaled + poisson\" CO 3 CS 1.3 JUST LEFT",
	    "LAB 4 VPOS 0.6 0.3 MA 7 \"event list\" CO 4 CS 1.3 JUST LEFT",
	    "LAB  X  Count Rate (cts/s)",
	    "LAB  Y  Estimated Spectral Index",
	    "R X 0.08 130",
	    "R Y   -0.3",
	    "LOG  X ON 1",
	    "!"
	};
	for ( int i=0; i < header.length; i++ ) 
	    pw.println(header[i]);


	/**  START: Loop on count rates  **/

	for ( int r=0; r < countRates.length; r++ ) {

	    meanRate = countRates[r];
	    nevents = (new Double(meanRate*duration)).intValue();
	    nevents = poisson.nextInt(nevents);
	    System.out.println("\n"+"Log  : Rate = "+meanRate+" cps");

	    /**  Print to file compareTKImplementation.qdp  **/
	    pw.print(meanRate +"\t");

	    /**   START: Loop on nLC  **/
	    double[] chi2_compareScaledAndEvlist = new double[nLC];
	    double[] chi2_compareScaledAndPoisson = new double[nLC];

	    for ( int i=0; i < nLC; i++ ) {

	    
		/**  Generate TK light curve  **/
		lc_TK = TimmerKonig.makeTimmerLC(spectralIndex, duration, nTKbins); 
		//System.out.println("Log  : lc_TK.length = "+ lc_TK.length);
	    

		/** Convert lc to CDF **/
		lcHisto = Converter.array2histo("light curve", tzero, lc_TK_bintime, lc_TK);
		cdfHisto = Stats.getCDFHisto(lcHisto);


		/**  Draw events from CDF  **/
		arrivalTimes = Analysis.getRandom(cdfHisto, nevents);
		Arrays.sort(arrivalTimes);
		t = ArrivalTimes.adjustTimes(arrivalTimes, duration);
		//System.out.println("Log  : arrivalTimes.length = "+arrivalTimes.length);


		/**  Bin events into a LC and get mean and var **/
		lc_evlist = TimingUtils.makeLC(t, nLCbins);
		lc_binCentres = (double[]) lc_evlist[0];
		lc_binHeights = (double[]) lc_evlist[1];
// 		lc_evlist_mean = Stats.getMean(lc_binHeights);
// 		lc_evlist_var = Stats.getVariance(lc_binHeights);
// 		System.out.println("Log  : Arrival times binned: lc_binHeights.length = "+lc_binHeights.length);
// 		System.out.println("Log  : lc_evlist_mean = "+lc_evlist_mean);
// 		System.out.println("Log  : lc_evlist_var = "+lc_evlist_var);


		/**  Scale TK lc with same variance as LC from events and rebin **/
		lc_TK_rebinned = Binner.rebinRatesSimple(lc_TK, lc_TK_bintime, bintime);
// 		System.out.println("Log  : TK lc rebinned: lc_TK_rebinned.length = "+lc_TK_rebinned.length+" "+nLCbins);
 		lc_TK_rebinned_mean = Stats.getMean(lc_TK_rebinned);
 		lc_TK_rebinned_var = Stats.getVariance(lc_TK_rebinned);
		for ( int k=0; k < nLCbins; k++ ) {
// 		    lc_TK_scaled[k] = (lc_TK_rebinned[k] - lc_TK_rebinned_mean)*Math.sqrt(lc_evlist_var/lc_TK_rebinned_var);
		    lc_TK_scaled[k] = lc_TK_rebinned[k]*meanRate;
		}
//   		lc_TK_scaled_mean = Stats.getMean(lc_TK_scaled);
//   		lc_TK_scaled_var = Stats.getVariance(lc_TK_scaled);


		/**  Compare scaled and evlist LC  **/
		chi2_compareScaledAndEvlist[i] = Stats.computeChi2Test(lc_binHeights, lc_TK_scaled, nLCbins);


		/**  Add Poisson noise **/
		lc_TK_poisson = new double[lc_TK_scaled.length];
		for ( int p=0; p < lc_TK_poisson.length; p++ ) {
		    //lc_TK_scaled[p] = lc_TK_scaled[p] - lc_TK_scaled_mean + lc_evlist_mean;
		    double eventsInBin = lc_TK_scaled[p]*bintime;
		    int events = poisson.nextInt(eventsInBin);
		    lc_TK_poisson[p] = events/newBinTime;
		}
// 		System.out.println("Log  : Scaled: lc_TK_scaled.length = "+lc_TK_scaled.length);
// 		System.out.println("Log  : lc_TK_scaled_mean = "+lc_TK_scaled_mean);
// 		System.out.println("Log  : lc_TK_scaled_var = "+lc_TK_scaled_var);


		/**  Compare scaled and Poisson LC  **/
		chi2_compareScaledAndPoisson[i] = Stats.computeChi2Test(lc_TK_poisson, lc_TK_scaled, nLCbins);


		/**  Rescale the variance and mean to those of the event list light curve  **/
// 		lc_TK_poisson_var = Stats.getVariance(lc_TK_poisson);
// 		for ( int k=0; k < lc_TK_poisson.length; k++ ) {
// 		    lc_TK_poisson[k] *= Math.sqrt(lc_evlist_var/lc_TK_poisson_var);
// 		}
// 		lc_TK_poisson_mean = Stats.getMean(lc_TK_poisson);
// 		for ( int k=0; k < lc_TK_poisson.length; k++ ) {
// 		    lc_TK_poisson[k] = lc_TK_poisson[k] - lc_TK_poisson_mean + lc_evlist_mean;
// 		}

// 		lc_TK_poisson_mean = Stats.getMean(lc_TK_poisson);
// 		lc_TK_poisson_var = Stats.getVariance(lc_TK_poisson);

// 		System.out.println("Log  : Poisson: lc_TK_poisson.length = "+lc_TK_poisson.length);
// 		System.out.println("Log  : lc_TK_poisson_mean = "+lc_TK_poisson_mean);
// 		System.out.println("Log  : lc_TK_poisson_var = "+lc_TK_poisson_var);

// 		System.out.println("E: "+lc_evlist_mean+"\t"+lc_TK_scaled_mean+"\t"+lc_TK_poisson_mean);
// 		System.out.println("V: "+lc_evlist_var+"\t"+lc_TK_scaled_var+"\t"+lc_TK_poisson_var);


		/**  Make power spectra  **/
		lcScaled_powspec = Periodograms.makePSD_FFT(lc_binCentres, lc_TK_scaled, "leahy");
		lcPoisson_powspec = Periodograms.makePSD_FFT(lc_binCentres, lc_TK_poisson, "leahy");
		evlist_powspec = Periodograms.makePSD_FFT(lc_binCentres, lc_binHeights, "leahy");
		
		freqs = evlist_powspec[0];
		lcScaled_pow = lcScaled_powspec[1];
		lcPoisson_pow = lcPoisson_powspec[1];
		evlist_pow = evlist_powspec[1];
		
		//System.out.println("Log  : evlist_pow.length = "+evlist_pow.length);
		//System.out.println("Log  : evlist_pow.length = "+evlist_pow.length);


		/**  Transform powspec to Log-space  **/
		x = new double[freqs.length];
		y1 = new double[evlist_pow.length];
		y2 = new double[evlist_pow.length];
		y3 = new double[evlist_pow.length];
		for ( int k=0; k < evlist_pow.length; k ++ ) {
		    x[k] = Math.log10(freqs[k]);
		    y1[k] = Math.log10(lcScaled_pow[k]);
		    y2[k] = Math.log10(lcPoisson_pow[k]);
		    y3[k] = Math.log10(evlist_pow[k]);
		}
		lsFitResults = Stats.leastSquaresFitLine(x, y1);
		/** function is: pow = a + b*testFreqs, method returns [a][b][err_a][err_b]  **/
		slopes_scaled[i] = -lsFitResults[1];
		slopeErrs_scaled[i] = lsFitResults[3];

		lsFitResults = Stats.leastSquaresFitLine(x, y2);
		slopes_poisson[i] = -lsFitResults[1];
		slopeErrs_poisson[i] = lsFitResults[3];

		lsFitResults = Stats.leastSquaresFitLine(x, y3);
		slopes_evlist[i] = -lsFitResults[1];
		slopeErrs_evlist[i] = lsFitResults[3];


	    }

	    /** END: loop on nLC  **/
	    System.out.println("Log  : Reduced chi2 between expected signal and:");
	    System.out.println("Log  :    event list LC = "+number.format(Stats.getMean(chi2_compareScaledAndEvlist)));
	    System.out.println("Log  :    poisson LC = "+number.format(Stats.getMean(chi2_compareScaledAndPoisson)));



	    /**  Bin the results **/
	    int nHistoBins = 25;
	    meanAndVar = Stats.getRunningAveAndVar(slopes_scaled);
	    double ave_scaled = meanAndVar[0];
	    double var_scaled = meanAndVar[1];
	    double sig_scaled = Math.sqrt(var_scaled);
	    double minSlope_scaled = ave_scaled - 3.2*sig_scaled;
	    double maxSlope_scaled = ave_scaled + 3.2*sig_scaled;
	    double slopeBinWidth_scaled = (maxSlope_scaled - minSlope_scaled)/nHistoBins;
	    double[] histoOfSlopes_scaled = 
		Binner.binData(slopes_scaled, minSlope_scaled, maxSlope_scaled, nHistoBins);
	    double sum = Stats.getSum(histoOfSlopes_scaled);
	    int n_scaled = (new Double(sum)).intValue();
	    Analysis.normalise(histoOfSlopes_scaled, 1.0/(sum));

	    meanAndVar = Stats.getRunningAveAndVar(slopes_poisson);
	    double ave_poisson = meanAndVar[0];
	    double var_poisson = meanAndVar[1];
	    double sig_poisson = Math.sqrt(var_poisson);
	    double minSlope_poisson = ave_poisson - 3.2*sig_poisson;
	    double maxSlope_poisson = ave_poisson + 3.2*sig_poisson;
	    double slopeBinWidth_poisson = (maxSlope_poisson - minSlope_poisson)/nHistoBins;
	    double[] histoOfSlopes_poisson = 
		Binner.binData(slopes_poisson, minSlope_poisson, maxSlope_poisson, nHistoBins);
	    sum = Stats.getSum(histoOfSlopes_poisson);
	    int n_poisson = (new Double(sum)).intValue();
	    Analysis.normalise(histoOfSlopes_poisson, 1.0/(sum));

	    meanAndVar = Stats.getRunningAveAndVar(slopes_evlist);
	    double ave_evlist = meanAndVar[0];
	    double var_evlist = meanAndVar[1];
	    double sig_evlist = Math.sqrt(var_evlist);
	    double minSlope_evlist = ave_evlist - 3.2*sig_evlist;
	    double maxSlope_evlist = ave_evlist + 3.2*sig_evlist;
	    double slopeBinWidth_evlist = (maxSlope_evlist - minSlope_evlist)/nHistoBins;
	    double[] histoOfSlopes_evlist = 
		Binner.binData(slopes_evlist, minSlope_evlist, maxSlope_evlist, nHistoBins);
	    sum = Stats.getSum(histoOfSlopes_evlist);
	    int n_evlist = (new Double(sum)).intValue();
	    Analysis.normalise(histoOfSlopes_evlist, 1.0/(sum));


	    /**  Construct the axes for histos and Gaussian function  **/
	    double norm_scaled = slopeBinWidth_scaled; // n_scaled*slopeBinWidth_scaled;
	    double norm_poisson = slopeBinWidth_poisson;  // n_poisson*slopeBinWidth_poisson;
	    double norm_evlist = slopeBinWidth_evlist;  // n_evlist*slopeBinWidth_evlist;
	    double[] slopeGauss_scaled = new double[nHistoBins];
	    double[] slopeGauss_poisson = new double[nHistoBins];
	    double[] slopeGauss_evlist = new double[nHistoBins];
	    double[] slopeAxis_scaled = new double[nHistoBins];
	    double[] slopeAxis_poisson = new double[nHistoBins];
	    double[] slopeAxis_evlist = new double[nHistoBins];
	    double[] diffWithModel_scaled = new double[nHistoBins];
	    double[] diffWithModel_poisson = new double[nHistoBins];
	    double[] diffWithModel_evlist = new double[nHistoBins];
	    for ( int i=0; i < nHistoBins; i++ ) {

		slopeAxis_scaled[i] = minSlope_scaled + (0.5 + i)*slopeBinWidth_scaled;
		slopeAxis_poisson[i] = minSlope_poisson + (0.5 + i)*slopeBinWidth_poisson;
		slopeAxis_evlist[i] = minSlope_evlist + (0.5 + i)*slopeBinWidth_evlist;

		slopeGauss_scaled[i] = norm_scaled*Math.sqrt(1/(2*Math.PI*var_scaled)) * 
		    Math.exp( -Math.pow((slopeAxis_scaled[i]-ave_scaled), 2)/(2*var_scaled) );
		slopeGauss_poisson[i] = norm_poisson*Math.sqrt(1/(2*Math.PI*var_poisson)) * 
		    Math.exp( -Math.pow((slopeAxis_poisson[i]-ave_poisson), 2)/(2*var_poisson) );
		slopeGauss_evlist[i] = norm_evlist*Math.sqrt(1/(2*Math.PI*var_evlist)) * 
		    Math.exp( -Math.pow((slopeAxis_evlist[i]-ave_evlist), 2)/(2*var_evlist) );

		diffWithModel_scaled[i] = slopeGauss_scaled[i] - histoOfSlopes_scaled[i];
		diffWithModel_poisson[i] = slopeGauss_poisson[i] - histoOfSlopes_poisson[i];
		diffWithModel_evlist[i] = slopeGauss_evlist[i] - histoOfSlopes_evlist[i];
	    }


	    /**  Compute Pearson's Chi2 test  for the Gaussians **/
	    System.out.println("Log  : Mean of scaled LC slopes = "+number.format(ave_scaled));
	    System.out.println("Log  : Standard deviation of scaled LC slopes = "
			       +number.format(Math.sqrt(var_scaled)));

	    System.out.println("Log  : Mean of poisson LC slopes = "+number.format(ave_poisson));
	    System.out.println("Log  : Standard deviation of poisson LC slopes = "
			       +number.format(Math.sqrt(var_poisson)));

	    System.out.println("Log  : Mean of evlist slopes = "+number.format(ave_evlist));
	    System.out.println("Log  : Standard deviation of evlist slopes = "+sci.format(var_evlist));

// 	    double chi2_scaled = Stats.computeChi2Test(histoOfSlopes_scaled, slopeGauss_scaled, n_scaled);
// 	    double chi2_poisson = Stats.computeChi2Test(histoOfSlopes_poisson, slopeGauss_poisson, n_poisson);
// 	    double chi2_evlist = Stats.computeChi2Test(histoOfSlopes_evlist, slopeGauss_evlist, n_evlist);


	    /**  Write the histogram of scaled LC slopes  **/
 	    String mean = number.format(ave_scaled);
 	    String stdDev = number.format(Math.sqrt(var_scaled));
//  	    String chi2 = number.format(chi2_scaled);
// 	    String dof = Integer.toString(nHistoBins-1);
// 	    String bigN = Integer.toString(n_scaled);
// 	    String lineStart = number.format(slopeAxis_scaled[0]);
// 	    String lineLength = number.format(slopeAxis_scaled[nHistoBins-1] - slopeAxis_scaled[0]);
// 	    header = new String[] {
// 		"DEV /XS",
// 		"TIME OFF", "LAB T", "LAB F",
// 		"LINE STEP ON 2 4",
// 		"LINE ON 3", "LW 3", "CS 1.1",
// 		"LAB 3 VPOS 0.175 0.7 \"\\gm = "+mean+"\" JUST LEFT",
// 		"LAB 4 VPOS 0.175 0.67 \"\\gs = "+stdDev+"\" JUST LEFT",
// 		"LAB 5 VPOS 0.175 0.64 \"N = "+bigN+"\" JUST LEFT",
// 		"LAB 6 VPOS 0.825 0.7 \"\\gx\\u2\\d("+dof+") = "+chi2+"\" JUST RIGHT",
// 		"WIN 1",
// 		"YPLOT 2 3",
// 		"LOC  0.05 0.302 1 0.8",
// 		"LAB NX OFF",
// 		"LAB Y Normalized PDF",
// 		"WIN 2",
// 		"YPLOT 4",
// 		"LOC  0.05 0.2 1 0.368",
// 		"LAB 10 POS "+lineStart+" 0 LINE 0 1 LS 2 \"",
// 		"LAB X Estimated Spectral Index",
// 		"LAB Y Residues",
// 		"!",
// 		"! Mean = "+ave_scaled,
// 		"! Std dev = "+sig_scaled,
// 		"! Pearson's Chi2 = "+chi2_scaled,
// 		"!"
// 	    };
//  	    DataFile histoOfSlopesFile = new DataFile("histoOfSlopes_scaled_"+meanRate+"cps.qdp");
// 	    histoOfSlopesFile.writeData(header, slopeAxis_scaled, histoOfSlopes_scaled, slopeGauss_scaled, diffWithModel_scaled);


	    /**  Print to file compareTKImplementation.qdp  **/
	    pw.print(mean +"\t"+ stdDev+"\t");


// 	    /**  Write the histogram of evlist slopes **/
 	    mean = number.format(ave_poisson);
 	    stdDev = number.format(Math.sqrt(var_poisson));
//  	    chi2 = number.format(chi2_poisson);
// 	    dof = Integer.toString(nHistoBins-1);
// 	    bigN = Integer.toString(n_poisson);
// 	    lineStart = number.format(slopeAxis_poisson[0]);
// 	    lineLength = number.format(slopeAxis_poisson[nHistoBins-1] - slopeAxis_poisson[0]);
// 	    header = new String[] {
// 		"DEV /XS",
// 		"READ 1 2 3",
// 		"TIME OFF", "LAB T", "LAB F",
// 		"LINE STEP ON 2 4",
// 		"LINE ON 3", "LW 3", "CS 1.1",
// 		"LAB 3 VPOS 0.175 0.7 \"\\gm = "+mean+"\" JUST LEFT",
// 		"LAB 4 VPOS 0.175 0.67 \"\\gs = "+stdDev+"\" JUST LEFT",
// 		"LAB 5 VPOS 0.175 0.64 \"N = "+bigN+"\" JUST LEFT",
// 		"LAB 6 VPOS 0.825 0.7 \"\\gx\\u2\\d("+dof+") = "+chi2+"\" JUST RIGHT",
// 		"WIN 1",
// 		"YPLOT 2 3",
// 		"LOC  0.05 0.302 1 0.8",
// 		"LAB NX OFF",
// 		"LAB Y Normalized PDF",
// 		"WIN 2",
// 		"YPLOT 4",
// 		"LOC  0.05 0.2 1 0.368",
// 		"LAB 10 POS "+lineStart+" 0 LINE 0 1 LS 2 \"",
// 		"LAB X Estimated Spectral Index",
// 		"LAB Y Residues",
// 		"!",
// 		"! Mean = "+ave_poisson,
// 		"! Std dev = "+sig_poisson,
// 		"! Pearson's Chi2 = "+chi2_poisson,
// 		"!"
// 	    };
// 	    histoOfSlopesFile = new DataFile("histoOfSlopes_poisson_"+meanRate+"cps.qdp");
// 	    histoOfSlopesFile.writeData(header, slopeAxis_poisson, histoOfSlopes_poisson, slopeGauss_poisson, diffWithModel_poisson);


	    /**  Print to file compareTKImplementation.qdp  **/
	    pw.print(mean +"\t"+ stdDev+"\t");


// 	    /**  Write the histogram of evlist slopes **/
 	    mean = number.format(ave_evlist);
 	    stdDev = number.format(Math.sqrt(var_evlist));
//  	    chi2 = number.format(chi2_evlist);
// 	    dof = Integer.toString(nHistoBins-1);
// 	    bigN = Integer.toString(n_evlist);
// 	    lineStart = number.format(slopeAxis_evlist[0]);
// 	    lineLength = number.format(slopeAxis_poisson[nHistoBins-1] - slopeAxis_evlist[0]);
// 	    header = new String[] {
// 		"DEV /XS",
// 		"READ 1 2 3",
// 		"TIME OFF", "LAB T", "LAB F",
// 		"LINE STEP ON 2 4",
// 		"LINE ON 3", "LW 3", "CS 1.1",
// 		"LAB 3 VPOS 0.175 0.7 \"\\gm = "+mean+"\" JUST LEFT",
// 		"LAB 4 VPOS 0.175 0.67 \"\\gs = "+stdDev+"\" JUST LEFT",
// 		"LAB 5 VPOS 0.175 0.64 \"N = "+bigN+"\" JUST LEFT",
// 		"LAB 6 VPOS 0.825 0.7 \"\\gx\\u2\\d("+dof+") = "+chi2+"\" JUST RIGHT",
// 		"WIN 1",
// 		"YPLOT 2 3",
// 		"LOC  0.05 0.302 1 0.8",
// 		"LAB NX OFF",
// 		"LAB Y Normalized PDF",
// 		"WIN 2",
// 		"YPLOT 4",
// 		"LOC  0.05 0.2 1 0.368",
// 		"LAB 10 POS "+lineStart+" 0 LINE 0 1 LS 2 \"",
// 		"LAB X Estimated Spectral Index",
// 		"LAB Y Residues",
// 		"!",
// 		"! Mean = "+ave_evlist,
// 		"! Std dev = "+sig_evlist,
// 		"! Pearson's Chi2 = "+chi2_evlist,
// 		"!"
// 	    };
// 	    histoOfSlopesFile = new DataFile("histoOfSlopes_evlist_"+meanRate+"cps.qdp");
// 	    histoOfSlopesFile.writeData(header, slopeAxis_evlist, histoOfSlopes_evlist, slopeGauss_evlist, diffWithModel_evlist);

	    /**  Print to file compareTKImplementation.qdp  **/
	    pw.println(mean +"\t"+ stdDev);
	    pw.flush();

	}

	/**  END: Loop on count rates  **/

	/**  Close file compareTKImplementation.qdp  **/
	pw.flush();
	pw.close();

    }

}
