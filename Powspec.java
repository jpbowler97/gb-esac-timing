package gb.esac.timing;


import gb.esac.binner.BinningException;
import gb.esac.eventlist.EventList;
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.ModifiedRayleighPeriodogram;
import gb.esac.periodogram.PeriodogramException;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.periodogram.PeriodogramUtils;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesFileException;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.timeseries.TimeSeriesResampler;
import gb.esac.tools.BasicStats;
import gb.esac.tools.MinMax;
import java.io.IOException;
import java.text.DecimalFormat;
import org.apache.log4j.Logger;



/**
 *  <code>Powspec</code> makes the FFT periodogram of an event list or light curve.
 *
 * @author <a href="mailto:">Guillaume Belanger</a>
 * @date Dec 15, 2010 (last modified)
 * @version 1.0
 */


public class Powspec {

    public static Logger logger = Logger.getLogger(Powspec.class);

    private static DecimalFormat number = new DecimalFormat("0.0000");
    private static DecimalFormat freq = new DecimalFormat("0.0000E00");
    
    private static String filename = null;
    private static double binTime = 0;
    private static String normName = "leahy";
    private static String windowName = "rectangular";
    private static int samplingFactor = 1;
    private static int rebinFactor = 1;
    private static TimeSeries ts;


    public static void main(String[] args) throws Exception {

	readArguments(args);
	makeTimeSeriesFromInputFile();
	makePowspec();
    }



    private static void readArguments(String[] args) {

	if ( args.length == 1 ) {
	    filename = args[0];
	}
	else if ( args.length == 2 ) {
	    filename = args[0];
	    binTime = (Double.valueOf(args[1])).doubleValue();
	}
	else if ( args.length == 3 ) {
	    filename = args[0];
	    binTime = (Double.valueOf(args[1])).doubleValue();
	    windowName = args[2];
	}
	else if ( args.length == 4 ) {
	    filename = args[0];
	    binTime = (Double.valueOf(args[1])).doubleValue();
	    windowName = args[2];
	    normName = args[3];
	}
	else if ( args.length == 5 ) {
	    filename = args[0];
	    binTime = (Double.valueOf(args[1])).doubleValue();
	    windowName = args[2];
	    normName = args[3];
	    samplingFactor = (Integer.valueOf(args[4])).intValue();
	}
	else if ( args.length == 6 ) {
	    filename = args[0];
	    binTime = (Double.valueOf(args[1])).doubleValue();
	    windowName = args[2];
	    normName = args[3];
	    samplingFactor = (Integer.valueOf(args[4])).intValue();
	    rebinFactor = (Integer.valueOf(args[5])).intValue();
	}
	else {
	    logger.info("Usage: java Powspec evfile (binTime) (window [rectangular]) (norm [leahy]) (sampling [1]) (rebinFactor [1])");
	    System.exit(-1);
	}
	logger.info("Running Powspec: Assuming time scale in seconds");
    }

    private static void makeTimeSeriesFromInputFile() throws IOException, TimeSeriesFileException, BinningException {

	ts = TimeSeriesMaker.makeTimeSeries(filename);
	double origBinTime = MinMax.getMax(ts.getBinWidths());
	if ( binTime != 0 ) {   // binTime was specifies as an argument
	    if ( binTime < origBinTime ) {
		logger.error("Specified bintime ("+binTime+" s) must be larger than original bintime ("+origBinTime+" s)");
		System.exit(-1);
	    }
	    else if ( binTime != origBinTime ) {
		ts = TimeSeriesResampler.resample(ts, binTime);
	    }
	}
    }

    private static void makePowspec() throws PeriodogramException, BinningException {

 	FFTPeriodogram fft = PeriodogramMaker.makeOversampledWindowedFFTPeriodogram(ts, windowName, normName, samplingFactor);
	if ( rebinFactor > 1 ) {
	    fft = (FFTPeriodogram) fft.rebin(rebinFactor, "papadakis");
	}
  	fft.writeAsQDP("powspec.qdp");
    }


}
