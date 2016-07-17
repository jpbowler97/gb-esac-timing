package gb.esac.timing;


import cern.colt.list.DoubleArrayList;
import cern.jet.stat.Descriptive;
import cern.jet.stat.Probability;
import gb.esac.eventlist.EventList;
import gb.esac.io.AsciiDataFileWriter;
import gb.esac.montecarlo.RedNoiseGenerator;
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.ModifiedRayleighPeriodogram;
import gb.esac.periodogram.Periodogram;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.periodogram.PeriodogramUtils;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.tools.MinMax;
import java.text.DecimalFormat;
import org.apache.log4j.Logger;


public class PeriodogramProbabilities {

    private static Logger logger  = Logger.getLogger(PeriodogramProbabilities.class);

    private static DecimalFormat num = new DecimalFormat("0.0000");
    private static DecimalFormat sci = new DecimalFormat("0.0000E0");

    public static void main(String[] args) throws Exception {


	String filename = "/Users/gbelanger/javaProgs/flare_m1m2pn.fits";
	filename = "/Users/gbelanger/javaProgs/flare.fits";
	String psdType = "r2";
	double alpha = 2;
	int n = 100;
	
	if ( args.length == 4 ) {
	    filename = args[0];
	    alpha = (Double.valueOf(args[1])).doubleValue();
	    psdType = args[2];
	    n = (Integer.valueOf(args[3])).intValue();
	}
	else if ( args.length == 3 ) {
	    filename = args[0];
	    alpha = (Double.valueOf(args[1])).doubleValue();
	    psdType = args[2];
	}
	else if ( args.length == 2 ) {
	    filename = args[0];
	    alpha = (Double.valueOf(args[1])).doubleValue();
	}
	else {
	    logger.error("Usage: java gb.esac.timing.PeriodogramProbabilities filename alpha (psdType [fft | r2]) (nSimulations [100])");
	    System.exit(-1);
	}

	//  Read event list
	EventList evlist = new EventList(filename);

// 	//  Make LC for display with 100 s bins
// 	TimeSeries lc = TimeSeriesMaker.makeTimeSeries(evlist, 100d);
// 	lc.writeCountsAsQDP("lc.qdp");


	//  Make the FFT periodogram and fit the power-law index
	FFTPeriodogram fft = PeriodogramMaker.makePlainFFTPeriodogram(TimeSeriesMaker.makeTimeSeries(evlist, 128));
	double[] fitRes = PeriodogramUtils.fitPowerLawInLogSpace(fft);
	double estimatedAlpha = fitRes[0];
	double error = fitRes[1];

	logger.info("Estimated spectral index is = "+num.format(estimatedAlpha)+" +/- "+num.format(error));
	double countRate = evlist.meanCountRate();
	logger.info("Count rate is ="+ num.format(countRate));


	//  Look up in spectral index surface and change 
	estimatedAlpha = 2;


	// Make Periodogram from evlist
	double duration = evlist.duration();
	int startFreq = 2;
	double nuMin = (startFreq+1)/duration;
	double nuMax = 0.0128;
	nuMax = 5e-3;
	int sampling = 10;
	Periodogram dataPSD = null;
	if ( psdType.equals("fft") ) {
	    dataPSD = (FFTPeriodogram) PeriodogramMaker.makePlainFFTPeriodogram(evlist);
	}
	else {
	    dataPSD =  (ModifiedRayleighPeriodogram) PeriodogramMaker.makeModifiedRayleighPeriodogram(evlist, nuMin, nuMax, sampling);
	}
	dataPSD.writeAsQDP("dataPSD.qdp");


	//  Initialise the Array lists that will store the power values
	double[] dataPow = dataPSD.getPowers();
	double[] dataFreqs = dataPSD.getFreqs();
	int nFreqs = 0;
	if ( psdType.equals("fft") ) {
	    nFreqs = (int) Math.round((nuMax-nuMin)*duration) +1;  // take the number of IFS betwen nuMin and nuMax
	}
	else {
	    nFreqs = dataPow.length;
	}
	DoubleArrayList[] powerAtThisFreq = new DoubleArrayList[nFreqs];
 	for ( int i=0; i < nFreqs; i++ ) {
	    powerAtThisFreq[i] = new DoubleArrayList();
	}
	DoubleArrayList peaks = new DoubleArrayList();


	//  Run the simulation by looping over the segments of a long light curve

// 	int k=0;
// 	double[] times = RedNoiseGenerator.generateArrivalTimes(countRate, n*duration, estimatedAlpha);
// 	EventList simEvlist = new EventList(times);
// 	while ( k < n ) {

// 	    double from = k*duration;
// 	    double to = (k+1)*duration;
// 	    double[] segment = simEvlist.getArrivalTimesFromTo(from, to);
// 	    double[] powers = null;
// 	    if ( psdType.equals("fft") ) {
// 		FFTPeriodogram psd = PeriodogramMaker.makePlainFFTPeriodogram(new EventList(segment));
// 		double[] fftPowers = psd.getPowers();
// 		powers = new double[nFreqs];
// 		for ( int i=0; i < nFreqs; i++ ) {
// 		    powers[i] = fftPowers[i+startFreq]; //  skip to the first test freq
// 		}
// 	    }
// 	    else {
// 		ModifiedRayleighPeriodogram psd = PeriodogramMaker.makeModifiedRayleighPeriodogram(new EventList(segment), nuMin, nuMax, sampling);
// 		powers = psd.getPowers();
// 	    }

// 	    int nValues = Math.min(nFreqs, powers.length);
// 	    for ( int i=0; i < nValues; i++ ) {
// 		powerAtThisFreq[i].add(powers[i]);
// 	    }
// 	    peaks.add(MinMax.getMax(powers));
// 	    k++;
// 	}
// 	peaks.trimToSize();
// 	peaks.sort();


	//  Run the simulation by looping over n separate short light curves

	int k=0;
	while ( k < n ) {

	    double[] segment = RedNoiseGenerator.generateArrivalTimes(countRate, duration, estimatedAlpha);
	    double[] powers = null;
	    if ( psdType.equals("fft") ) {
		FFTPeriodogram psd = PeriodogramMaker.makePlainFFTPeriodogram(new EventList(segment));
		double[] fftPowers = psd.getPowers();
		powers = new double[nFreqs];
		for ( int i=0; i < nFreqs; i++ ) {
		    powers[i] = fftPowers[i+startFreq]; //  skip to the first test freq
		}
	    }
	    else {
		ModifiedRayleighPeriodogram psd = PeriodogramMaker.makeModifiedRayleighPeriodogram(new EventList(segment), nuMin, nuMax, sampling);
		powers = psd.getPowers();
	    }

	    int nValues = Math.min(nFreqs, powers.length);
	    for ( int i=0; i < nValues; i++ ) {
		powerAtThisFreq[i].add(powers[i]);
	    }
	    peaks.add(MinMax.getMax(powers));
	    k++;
	}
	peaks.trimToSize();
	peaks.sort();

	//  Derive the probabilities
	double[] local = new double[nFreqs];
	double[] global = new double[nFreqs];
	double[] testFreqs = new double[nFreqs];
	double[] probFromAvg = new double[nFreqs];
	double[] avgPower = new double[nFreqs];
	double[] avgPowerEnv = new double[nFreqs];
	if ( psdType.equals("r2") ) startFreq=0;
	for ( int i=0; i < nFreqs; i++ ) {
	    testFreqs[i] = dataFreqs[i+startFreq];
	    powerAtThisFreq[i].trimToSize();
	    powerAtThisFreq[i].sort();

	    local[i] = (1.0 - Descriptive.quantileInverse(powerAtThisFreq[i], dataPow[i+startFreq]));
	    global[i] = (1.0 - Descriptive.quantileInverse(peaks, dataPow[i+startFreq]));
	    avgPower[i] = Descriptive.mean(powerAtThisFreq[i]);
	    double variance = Descriptive.sampleVariance(powerAtThisFreq[i], avgPower[i]);
	    //avgPowerEnv[i] = avgPower[i] + Math.sqrt(variance);
	    //double dof = avgPowerEnv[i];
	    //probFromAvg[i] = Probability.chiSquareComplemented(dof, dataPow[i+startFreq]);
	}

	//  Write Probability 
	AsciiDataFileWriter out = new AsciiDataFileWriter("prob.qdp");
	String[] header = new String[] {
	    "! QDP data file",
	    "DEV /XS",
	    "VIEW 0.1 0.2 0.9 0.8",
	    "LAB TITLE", "LAB FILE", "TIME OFF", 
	    "LINE STEP",
	    "LW 4", "CSIZE 1.5", 
	    "LAB X Frequency (Hz)",
	    "LAB Y Probability",
	    "LOG ON",
	    "VIEW 0.2 0.1 0.8 0.9",
	    "!"
	};
	out.writeData(header, testFreqs, global);
	//out.writeData(header, testFreqs, global, probFromAvg);
	//out.writeData(header, testFreqs, dataPow, avgPower);

    }


}
