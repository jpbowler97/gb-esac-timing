package gb.esac.timing;

import gb.esac.binner.BinningUtils;
import gb.esac.binner.Resampler;
import gb.esac.io.AsciiDataFileWriter;
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.periodogram.WindowFunction;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.tools.BasicStats;
import gb.esac.tools.Complex;
import gb.esac.tools.FFT;
import gb.esac.tools.Stats;


public class LeakageFunction {

    public static void main(String[] args) throws Exception {

	//  Generate pure sinusoidal signal
	double duration = 500;
	double period = 0.01d;
	double freq = 1d/period;
	freq = 500;
	double omega = 2*Math.PI*freq;
	int nPointsPerCycle = 50;
	double dt = period/nPointsPerCycle;
	int nCycles = (int) Math.round(duration*period);
	int nPoints = nPointsPerCycle*nCycles;
	double[] t = new double[nPoints];
	double[] sineFunction = new double[nPoints];
	double sumOfBinHeights= 0;
	for ( int i=0; i < nPoints; i++ ) {
	    t[i] = i*dt;
	    sineFunction[i] = Math.sin(omega*t[i]) + 2;
	}
	

	int startPoint = (int) Math.round(nPoints/2d*Math.random());
	startPoint = 0;
	int n = 200;
	double[] segment_sine = new double[n];
	double[] segment_t = new double[n];
	for ( int i=0; i < n; i++ ) {
	    segment_sine[i] = sineFunction[i+startPoint];
	    segment_t[i] = t[i+startPoint];
	}
	double[] oldBinEdges = BinningUtils.getBinEdges(segment_t[0], segment_t[segment_t.length-1], n);
	double nn = Math.floor(Math.log(n)/Math.log(2));
	int nNewBins = (int) Math.pow(2, nn);
	double[] newBinEdges = BinningUtils.getBinEdges(segment_t[0], segment_t[segment_t.length-1], nNewBins);
	double[] resampledSegment = Resampler.resample(segment_sine, oldBinEdges, newBinEdges);

	TimeSeries ts = TimeSeriesMaker.makeTimeSeries(newBinEdges, resampledSegment);
	String normName = "leahy";
	String[] windowNames = (new WindowFunction("rectangular")).getAvailableFunctions();
	int sampling = 32;
	String[] header =  new String[] {
	    "DEV /XS",
	    "READ 1 2",
	    "LAB T", "LAB F",
	    "TIME OFF",
	    "LINE STEP",
	    "LOG Y ON",
	    "LOG X OFF",
	    "LW 4", "CS 1.5",
	    "LAB X Frequency (Hz)",
	    "LAB Y Power",
	    "VIEW 0.2 0.1 0.8 0.9",
	    "SKIP SINGLE",
	    "R X 0 1000",
	    "R Y 1.1e-9",
	    "!"
	};
	for ( int i=0; i < windowNames.length; i++ ) {
	    String windowName = windowNames[i];
	    FFTPeriodogram per = PeriodogramMaker.makeOversampledWindowedFFTPeriodogram(ts, windowName, normName, sampling);
	    per.writeAsQDP("powspec-"+windowName+".qdp", header);
	}

    }

}
