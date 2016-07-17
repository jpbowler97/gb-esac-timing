package gb.esac.timing;

import gb.esac.eventlist.EventList;
import gb.esac.montecarlo.RedNoiseGenerator;
import gb.esac.montecarlo.WhiteNoiseGenerator;
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.periodogram.PeriodogramUtils;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.tools.FitsUtils;
import java.io.FileOutputStream;
import java.text.DecimalFormat;
import nom.tam.fits.BinaryTable;
import nom.tam.fits.BinaryTableHDU;
import nom.tam.fits.Fits;
import nom.tam.fits.FitsFactory;
import nom.tam.fits.Header;
import nom.tam.util.BufferedDataOutputStream;
import org.apache.log4j.Logger;



/**
 * Describe class <code>GenerateEventList</code> here.
 *
 * @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>
 * @version 1.0 (May 2010, ESAC)
 */
public class GenerateEventList {


    static Logger logger  = Logger.getLogger(GenerateEventList.class);

    public static DecimalFormat number = new DecimalFormat("0.000");


    public static void main(String[] args) throws Exception {


	/**  Handle arguments  **/
	double duration = 0;
	double mean = 0;
	double index = 0;
	double period = 0;
	double pulsedFrac = 0;
	boolean whiteNoise = false;
	boolean redNoise = false;
	boolean whitePulse = false;
	boolean redPulse = false;

	if ( args.length == 2 ) {
	    duration = (Double.valueOf(args[0])).doubleValue();
	    mean = (Double.valueOf(args[1])).doubleValue();
	    whiteNoise = true;
	}
	else if ( args.length == 3 ) {
	    duration = (Double.valueOf(args[0])).doubleValue();
	    mean = (Double.valueOf(args[1])).doubleValue();
	    index = (Double.valueOf(args[2])).doubleValue();
	    if ( index < 0 || index > 4 ) {
		logger.error("Error: index must be between 0 and 4");
		System.exit(-1);
	    }
	    redNoise = true;
	}
	else if ( args.length == 4 ) {
	    duration = (Double.valueOf(args[0])).doubleValue();
	    mean = (Double.valueOf(args[1])).doubleValue();
	    period = (Double.valueOf(args[2])).doubleValue();
	    pulsedFrac = (Double.valueOf(args[3])).doubleValue();
	    if ( pulsedFrac > 1.0 ) {
		logger.error("Error: Pulsed fraction must be smaller than 1.0");
		System.exit(-1);
	    }
	    whitePulse = true;
	}
	else if ( args.length == 5 ) {
	    duration = (Double.valueOf(args[0])).doubleValue();
	    mean = (Double.valueOf(args[1])).doubleValue();
	    index = (Double.valueOf(args[2])).doubleValue();
	    if ( index < 0 || index > 4 ) {
		logger.error("Error: index must be between 0 and 4");
		System.exit(-1);
	    }
	    period = (Double.valueOf(args[3])).doubleValue();
	    pulsedFrac = (Double.valueOf(args[4])).doubleValue();
	    if ( pulsedFrac > 1.0 ) {
		logger.error("Error: Pulsed fraction must be smaller than 1.0");
		System.exit(-1);
	    }
	    redPulse = true;
	}
	else {
	    logger.info("Usage: java GenerateEventList duration mean [index] [period] [pulsedFrac]");
	    logger.info("  2 args = white noise");
	    logger.info("  3 args = red noise");
	    logger.info("  4 args = white noise + sinusoid");
	    logger.info("  5 args = red noise + sinusoid"); 
	    System.exit(-1);
	}


	//  Generate arrival times
	logger.info("Running GenerateEventList");
	double[] times = null;
	if ( whiteNoise ) {
	    times = WhiteNoiseGenerator.generateArrivalTimes(mean, duration);
	}
	else if ( redNoise ) {
 	    times = RedNoiseGenerator.generateArrivalTimes(mean, duration, index);
	}
	else if ( whitePulse ) {
	    times = WhiteNoiseGenerator.generateModulatedArrivalTimes(mean, duration, period, pulsedFrac);
	}
	else {
	    times = RedNoiseGenerator.generateModulatedArrivalTimes(mean, duration, index, period, pulsedFrac);
	}
	EventList list = new EventList(times);


	//  Make and write TimeSeries
	int nBins = (int) Math.round(duration*mean/50);
	TimeSeries lc = TimeSeriesMaker.makeTimeSeries(list, nBins);
	String lcName = "simEvlist-lc.qdp";
	lc.writeCountsAsQDP(lcName);


	//  Fit the periodogram made from the TimeSeries
// 	FFTPeriodogram psd = PeriodogramMaker.makePlainFFTPeriodogram(lc);
// 	double index1 = PeriodogramUtils.fitPowerLawInLinearSpace(psd)[0];
// 	double index2 = PeriodogramUtils.fitPowerLawInLogSpace(psd)[0];
// 	logger.info("Index in linear space = "+index1);
// 	logger.info("Index in log space = "+index2);
	

	//  Write event list as FITS file
	FitsFactory.setUseAsciiTables(false);
	String filename = "/Users/gbelanger/javaProgs/.template.fits";
	BinaryTable binTable = new BinaryTable(new Object[]{list.getArrivalTimes()});
	Fits f = FitsUtils.openFits(filename);
	Header tableHead = f.getHDU(1).getHeader();
  	double telapse = list.duration();
	tableHead.addValue("NAXIS2", list.nEvents(), "number of rows in table");
	tableHead.addValue("TSTART", list.tStart(), "Start time of first frame");
	tableHead.addValue("TSTOP", list.tStop(), "End time of last frame");
	tableHead.addValue("TELAPSE", telapse, "Full time interval for the exposure");
	tableHead.addValue("DURATION", telapse, "[s] Duration of observation");
	tableHead.addValue("ONTIME", telapse, "Sum of GTIs for central CCD");
	tableHead.addValue("LIVETIME", telapse, "Live time for the central CCD");
	BinaryTableHDU eventsHDU = new BinaryTableHDU(tableHead, binTable);
	BufferedDataOutputStream dos = new BufferedDataOutputStream(new FileOutputStream("simEvlist.fits"));
	Fits evlist = new Fits();
	evlist.addHDU(f.getHDU(0));
	evlist.addHDU(eventsHDU);
	evlist.write(dos);

    }

}
