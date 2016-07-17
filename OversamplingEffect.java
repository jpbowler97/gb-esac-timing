package gb.esac.timing;

import gb.esac.io.DataFile;
import gb.esac.tools.Analysis;
import gb.esac.tools.Stats;
import gb.esac.tools.Binner;
import gb.esac.aida.functions.ChiSquareFunction;

import java.io.PrintWriter;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import java.text.DecimalFormat;

import hep.aida.*;

import cern.jet.random.engine.MersenneTwister64;
import cern.jet.random.Uniform;



public class OversamplingEffect {

    public static void main(String[] args) throws IOException {

	DecimalFormat num = new DecimalFormat("0.00000");
	DecimalFormat sci = new DecimalFormat("0.0000E0");

	/**  Handle arguments  **/
	double duration = 10000;
	double meanRate = 1;
	double binTime = 20;
	int nEvlists = 1000;
	if ( args.length == 4 ) {
	    duration = (Double.valueOf(args[0])).doubleValue();
	    meanRate = (Double.valueOf(args[1])).doubleValue();
	    binTime = (Double.valueOf(args[2])).doubleValue();
	    nEvlists = (Integer.valueOf(args[3])).intValue();
	}
	else {
	    System.out.println("Usage: java DOFVsOversampling duration meanRate binTime nEvlists");
	    System.exit(-1);
	}
	System.out.println("Log  : Running VarianceVsSpectralIndex");
	System.out.println("Log  : Duration = "+duration);
	System.out.println("Log  : Mean rate = "+meanRate);
	System.out.println("Log  : Specified bin time = "+binTime);
	System.out.println("Log  : Number of event lists = "+nEvlists);


	/**  Determine power of 2 number of bins  **/
	double n = Math.ceil(Math.log10(duration/binTime)/Math.log10(2));
	int nTimeBins = (new Double(Math.pow(2, n))).intValue();
	double newBinTime = duration/nTimeBins;
	double nuMin = 1/duration;
	double nuMax = 1/(2*newBinTime);
	int nIFS = (new Double(Math.floor(duration*(nuMax - nuMin)))).intValue();
	System.out.println("Log  : New bin time for event list powspec = "+num.format(newBinTime)+" s");
	System.out.println("Log  : Min test frequency (1/duration) = "+sci.format(nuMin)+" Hz");
	System.out.println("Log  : Max test frequency (1/2*newBinTime) = "+sci.format(nuMax)+" Hz");
	System.out.println("Log  : Number of IFS in this range = "+nIFS);


	/**  Set up JAIDA factories  **/
	IAnalysisFactory af = IAnalysisFactory.create();
	ITree tree = af.createTreeFactory().create();
	IHistogramFactory hf = af.createHistogramFactory(tree);
	IFunctionFactory funcF  = af.createFunctionFactory(tree);
 	IFitFactory fitF   = af.createFitFactory();
	IFitter fitter = fitF.createFitter("Chi2", "jminuit");


	/**  Create the IFunctions to fit the different distributions **/	
	IFunction chi2 = new ChiSquareFunction("Chi Square Function");
	double dof = 0;
	double norm = 0;
	double expMean = 0;
	IFunction exp = funcF.createFunctionByName("Exponential", "e");
	IFunction gaus = funcF.createFunctionByName("Gaussian", "g");


	/**  Set up variables for simulations  **/
	double[] arrivalTimes = null;
	double[] lc = new double[nTimeBins];
	double[] lcTimes = new double[nTimeBins];
	double[] rates = new double[nTimeBins];
	double[] rateErrors = new double[nTimeBins];
	double[] r2Pow = null;
	double[] r2BinnedPow = null;
	double[] lombPow = null;
	int histoBins = 30;
	double histoMin = 0;
	double r2Max = 0;
	double lombMax = 0;
	

	/**  Set up printerWriter for output file **/
	PrintWriter pwR2 = new PrintWriter(new BufferedWriter(new FileWriter("meanVsSampling-R2.qdp")));
	PrintWriter pwR2Binned = new PrintWriter(new BufferedWriter(new FileWriter("meanVsSampling-R2Binned.qdp")));	
	PrintWriter pwLomb = new PrintWriter(new BufferedWriter(new FileWriter("meanVsSampling-Lomb.qdp")));

	String[] header = new String[] {
	    "! QDP data file", 
	    "DEV /XS",
	    "READ SERR 2 3 4",
	    "LABEL TITLE", "LABEL FILE", "TIME OFF", 
	    "LW 3", "CSIZE 1.3", "LINE ON",
	    "LAB X Oversampling Factor",
	    "LAB Y Mean of Distribution",
	    "VIEW 0.1 0.12 0.9 0.8",
	    "SKIP SINGLE",
	    "!"
	};
	for ( int i=0; i < header.length; i++ ) {
	    pwR2.println(header[i]);
	    pwR2Binned.println(header[i]);
	    pwLomb.println(header[i]);
	}
	pwR2.flush();
	pwR2Binned.flush();
	pwLomb.flush();

	
	/**  Get the time at the centre of each LC bin  **/
	for ( int i=0; i < nTimeBins; i++ ) {
	    lcTimes[i] = newBinTime/2.0 + i*newBinTime;
	}


	/**  Start simulations  **/
	double[] spectralIndex = new double[] {0.0, 1.0, 2.0, 3.0};
	double[] samplingFactor = new double[] {1.0, 5.0, 8.0, 10.0};
	double sampling = 0;
	int nSamplings = samplingFactor.length;	

	for ( int idx=0; idx < spectralIndex.length; idx++ ) {

	    double alpha = spectralIndex[idx];
	    System.out.println("Log  :");
	    System.out.println("Log  : Spectral index = "+alpha);

	    /**  Define variables to calculate running ave and var  **/
	    double r2_dof = 0;
	    double r2Binned_dof = 0;
	    double lomb_mean = 0;

	    double[] runAve_r2 = new double[nSamplings];
	    double[] runVar_r2 = new double[nSamplings];
	    double[] sum_r2 = new double[nSamplings];
	    double[] ave_r2 = new double[nSamplings];
	    double[] var_r2 = new double[nSamplings];

	    double[] runAve_r2Binned = new double[nSamplings];
	    double[] runVar_r2Binned = new double[nSamplings];
	    double[] sum_r2Binned = new double[nSamplings];
	    double[] ave_r2Binned = new double[nSamplings];
	    double[] var_r2Binned = new double[nSamplings];

	    double[] runAve_lomb = new double[nSamplings];
	    double[] runVar_lomb = new double[nSamplings];
	    double[] sum_lomb = new double[nSamplings];
	    double[] ave_lomb = new double[nSamplings];
	    double[] var_lomb = new double[nSamplings];

	    for ( int j=0; j < nEvlists; j++ ) {
		
		System.out.println("Log  : Event list "+(j+1)+":");

		/**  Generate arrival times **/
		arrivalTimes = ArrivalTimes.generateRedArrivalTimes(meanRate, duration, alpha);
		

		/**  Construct LC  **/
		lc = Binner.binOrderedData(arrivalTimes, nTimeBins);
		for ( int i=0; i < lc.length; i++ ) {
		    rates[i] = lc[i]/newBinTime;
		    rateErrors[i] = Math.sqrt(lc[i])/newBinTime;
		}


		/**  Loop on sampling factors  **/
		for ( int fac=0; fac < samplingFactor.length; fac++ ) {

		    sampling = samplingFactor[fac];
		    System.out.println("Log  : Oversampling factor = "+sampling);


		    /**  Make PSDs and keep the power values only ([0] = freq, [1] = pow)  **/
		    r2Pow = Periodograms.makePSD_ModRay(arrivalTimes, nuMin, nuMax, sampling)[1];
		    r2BinnedPow = Periodograms.makePSD_ModRay(lcTimes, rates, rateErrors, nuMin, nuMax, sampling)[1];
		    lombPow = Periodograms.makePSD_Lomb(lcTimes, rates, nuMin, nuMax, sampling)[1];


		    /**  Set the max for the histos  **/
		    r2Max = Math.ceil(2*Stats.getMean(r2BinnedPow) +1);
		    lombMax = Math.ceil(2*Stats.getMean(lombPow) +1);
		    //System.out.println(r2Max+"\t"+lombMax);


		    /**  Create and fill histos  **/
		    IHistogram1D r2 = hf.createHistogram1D("R2 histo", histoBins, histoMin, r2Max);
		    IHistogram1D r2Binned = hf.createHistogram1D("R2 histo", histoBins, histoMin, r2Max);
		    IHistogram1D lomb = hf.createHistogram1D("R2 histo", histoBins, histoMin, lombMax);
		    int nPowValues = r2BinnedPow.length;
		    for ( int m=0; m < nPowValues; m++ ) {
			r2.fill(r2Pow[m]);
			r2Binned.fill(r2BinnedPow[m]);
			lomb.fill(lombPow[m]);
		    }


		    /**  Fit each histo and store the DOF  **/
		    //  R2
// 		    dof = r2.mean();
// 		    norm = nPowValues*(r2Max - histoMin)/histoBins;
// 		    chi2.setParameter("dof", dof);
// 		    chi2.setParameter("norm", norm);
// 		    IFitResult fitResult = fitter.fit(r2, chi2);
// 		    r2_dof = fitResult.fittedParameter("dof");
		    r2_dof = Stats.getMean(r2Pow);
		    System.out.println("Log  :  r2 histo DOF = "+num.format(r2_dof));
// 				       +" (quality = "+num.format(fitResult.quality())+")");

		    //  R2 binned
// 		    fitResult = fitter.fit(r2Binned, chi2);
// 		    r2Binned_dof = fitResult.fittedParameter("dof");
		    r2Binned_dof = Stats.getMean(r2BinnedPow);
		    System.out.println("Log  :  r2Binned histo DOF = "+num.format(r2Binned_dof));
// 				       +" (quality = "+num.format(fitResult.quality())+")");

		    //  Lomb
		    expMean = lomb.mean();
		    norm = 0.5*nPowValues*(r2Max - histoMin)/histoBins;
		    exp.setParameter("exponent", -1/expMean);
		    exp.setParameter("amplitude", norm);
		    norm = nPowValues*(lombMax - histoMin)/histoBins;
		    IFitResult fitResult = fitter.fit(lomb, exp);
		    lomb_mean = -fitResult.fittedParameter("exponent");
// 		    lomb_mean = Stats.getMean(lombPow);
		    System.out.println("Log  :  Lomb histo mean = "+num.format(lomb_mean)
 				       +" (quality = "+num.format(fitResult.quality())+")");


		    /**  Calculate the running average and variance  **/
		    ave_r2[fac] = sum_r2[fac]/(j+1);
		    var_r2[fac] += Math.pow((r2_dof - ave_r2[fac]), 2);
		    
		    ave_r2Binned[fac]= sum_r2Binned[fac]/(j+1);
		    var_r2Binned[fac]+= Math.pow((r2Binned_dof - ave_r2Binned[fac]), 2);
		    
		    ave_lomb[fac] = sum_lomb[fac]/(j+1);
		    var_lomb[fac] += Math.pow((lomb_mean - ave_lomb[fac]), 2);


		    runAve_r2[fac] = ave_r2[fac] + (r2_dof - ave_r2[fac])/(j+1);
		    runVar_r2[fac] = var_r2[fac] + (r2_dof - runAve_r2[fac])*(r2_dof - ave_r2[fac]);
		    sum_r2[fac] += r2_dof;

		    runAve_r2Binned[fac]= ave_r2Binned[fac]+ (r2Binned_dof - ave_r2Binned[fac])/(j+1);
		    runVar_r2Binned[fac]= var_r2Binned[fac] + 
			(r2Binned_dof - runAve_r2Binned[fac])*(r2Binned_dof - ave_r2Binned[fac]);
		    sum_r2Binned[fac]+= r2Binned_dof;

		    runAve_lomb[fac] = ave_lomb[fac] + (lomb_mean - ave_lomb[fac])/(j+1);
		    runVar_lomb[fac] = var_lomb[fac] + (lomb_mean - runAve_lomb[fac])*(lomb_mean - ave_lomb[fac]);
		    sum_lomb[fac] += lomb_mean;

		}  /**  End loop on sampling  **/

		
	    } /**  End loop on event lists  **/


	    /**  Finalise the calculation of the variance and write results **/
	    System.out.println("Log  : Done loop on "+nEvlists+" event lists");
	    System.out.println("Log  : OVERALL averages for alpha = "+spectralIndex[idx]+":");
	    for ( int fac=0; fac < samplingFactor.length; fac++ ) {	    

		sampling = samplingFactor[fac];
		System.out.println("Log  :  Oversampling factor = "+sampling);

		/**  Variance  **/
		runVar_r2[fac] /= nEvlists;
		runVar_r2Binned[fac] /= nEvlists;
		runVar_lomb[fac] /= nEvlists;

		/**  Print to file  **/
		pwR2.println(sampling+"\t"+runAve_r2[fac]+"\t"+Math.sqrt(runVar_r2[fac]/nEvlists)+"\t");
		pwR2Binned.println(sampling+"\t"+runAve_r2Binned[fac]+"\t"+Math.sqrt(runVar_r2Binned[fac]/nEvlists)+"\t");
		pwLomb.println(sampling+"\t"+runAve_lomb[fac]+"\t"+Math.sqrt(runVar_lomb[fac]/nEvlists)+"\t");

		/**  Print to screen  **/
		System.out.println("Log  :   R2 mean = "
				   +num.format(runAve_r2[fac])+" +/- "
				   +num.format(Math.sqrt(runVar_r2[fac]/nEvlists)));
		System.out.println("Log  :   R2-binned mean = "
				   +num.format(runAve_r2Binned[fac])+" +/- "
				   +num.format(Math.sqrt(runVar_r2Binned[fac]/nEvlists)));
		System.out.println("Log  :   Lomb mean = "
				   +num.format(runAve_lomb[fac])+" +/- "
				   +num.format(Math.sqrt(runVar_lomb[fac]/nEvlists)));

	    }

	    pwR2.println("NO NO");
	    pwR2Binned.println("NO NO");
	    pwLomb.println("NO NO");

	    pwR2.flush();
	    pwR2Binned.flush();
	    pwLomb.flush();	    

	}  /**  End loop on spectral index  **/

	pwR2.close();
	pwR2Binned.close();
	pwLomb.close();
    }
}
