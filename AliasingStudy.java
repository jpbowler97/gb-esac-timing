package gb.esac.timing;

import cern.jet.random.engine.MersenneTwister64;
import gb.esac.montecarlo.RedNoiseGenerator;
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.periodogram.PeriodogramUtils;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.timeseries.TimeSeriesUtils;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Date;
import org.apache.log4j.Logger;


public class AliasingStudy {

    private static Logger logger  = Logger.getLogger(AliasingStudy.class);

    public static void main(String[] args) throws Exception {

	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());

	double duration = 1e5;
	double mean = 5;
	double period = 100;
	int nSpecs = 50;

	System.out.println("Generating red noise event lists with");
	System.out.println("  Duration = "+duration+" s");
	System.out.println("  Mean = "+mean+" cps");
	//System.out.println("  Period = "+period+" s");

	//  Define spectral indices to use
	double[] alphas = new double[] {0.5, 1, 2};
	System.out.println("Red noise spectral indices used are:");
	for ( int i=0; i < alphas.length; i++ ) {
	    System.out.println("  "+alphas[i]);
	}

	//  Define Nyquist bin time
	double nyquistBinTime = 1/(2*mean);
	int nBins = (int) Math.ceil(duration/nyquistBinTime);
	double n = Math.ceil(Math.log10(nBins)/Math.log10(2));
	int nNewBins = (int) Math.pow(2, n);
	double minBinTime = duration/nNewBins;

	//  Define the different numbers of bins to use
	int[] nTimeBins = new int[] {nNewBins, nNewBins/128};
	double[] binTimes = new double[nTimeBins.length];
 	System.out.println("Bintimes used are:");
	for ( int j=0; j < binTimes.length; j++ ) {
	    binTimes[j] = duration/nTimeBins[j];
 	    System.out.println("  "+binTimes[j]);
	}

	//  Loop on spectral indices
	for ( int i=0; i < alphas.length; i++ ) {

	    System.out.println("Running simulations for alpha = "+alphas[i]);
	    double alpha = alphas[i];

	    //  Initialise the average PSD for each bin time
	    int k=0;
	    System.out.println("  Generating event list "+(k+1)+" (of "+nSpecs+")");
	    FFTPeriodogram[] sum = new FFTPeriodogram[binTimes.length];
	    for ( int j=0; j < binTimes.length; j++ ) {

		double[] arrivalTimes = RedNoiseGenerator.generateArrivalTimes(mean, duration, alpha, nTimeBins[0], engine);
		//double[] arrivalTimes = RedNoiseGenerator.generateModulatedArrivalTimes(mean, duration, alpha, period, 0.1);

		sum[j] = PeriodogramMaker.makeFFTPeriodogram(TimeSeriesMaker.makeTimeSeries(arrivalTimes, binTimes[j]));
		//avg[j] = PeriodogramMaker.makeFFTPeriodogram(TimeSeriesUtils.kalmanFilter(TimeSeriesMaker.makeTimeSeries(arrivalTimes, 10d)));
	    }
	    k++;

	    //  Loop on nSpecs to build the average PSDs
	    FFTPeriodogram[] psd = new FFTPeriodogram[binTimes.length];
	    while ( k < nSpecs ) {

		System.out.println("  Generating event list "+(k+1)+" (of "+nSpecs+")");
		for ( int j=0; j < binTimes.length; j++ ) {

		    double[] arrivalTimes = RedNoiseGenerator.generateArrivalTimes(mean, duration, alpha, nTimeBins[0], engine);
		    //arrivalTimes = RedNoiseGenerator.generateModulatedArrivalTimes(mean, duration, alpha, period, 0.1);

		    psd[j] = PeriodogramMaker.makeFFTPeriodogram(TimeSeriesMaker.makeTimeSeries(arrivalTimes, binTimes[j]));
		    //psd[j] = PeriodogramMaker.makeFFTPeriodogram(TimeSeriesUtils.kalmanFilter(TimeSeriesMaker.makeTimeSeries(arrivalTimes, 10d)));

		    sum[j] = (FFTPeriodogram) sum[j].add(psd[j]);
		}
		k++;
	    }
	    

	    //  Make a plot of the results for each alpha
	    String type = "binning";
	    PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("aliasingForAlpha-"+alpha+"-"+type+".qdp")));

	    String[] header = new String[] {
		"DEV /XS",
		"READ 1 2",
		"LAB T", "LAB F",
		"TIME OFF",
		"LINE STEP",
		"LOG ON",
		"LW 4", "CS 1.3",
		"LAB X Frequency (Hz)",
		"LAB Y Power",
		"LAB 1 VPOS 0.83 0.73 \"\\ga = "+alpha+"\"",
		"LAB 1 JUST RIGHT CSIZE 1.4",
		"VIEW 0.1 0.15 0.9 0.85",
		"SKIP SINGLE",
		"R Y 1.001",
		"!"
	    };
	    for ( int m=0; m < header.length; m++ ) {
		pw.println(header[m]);
	    }
	    
	    for ( int m=0; m < sum.length; m++ ) {
		
		FFTPeriodogram avg = (FFTPeriodogram) sum[m].scale(1d/nSpecs);
		//FFTPeriodogram reb = (FFTPeriodogram) avg.rebin(10, "papadakis");
		
		double[] freqs = avg.getFreqs();
		double[] pow = avg.getPowers();
		k=0;
		while ( k < freqs.length ) {		    
		    pw.println(freqs[k]+"\t"+pow[k]);
		    k++;
		}
		pw.println("NO NO");
		pw.flush();
		
	    }
	    pw.println("HARD aliasingForAlpha-"+alpha+"-"+type+".ps/cps");
	    pw.close();

	}
	
    }
}
