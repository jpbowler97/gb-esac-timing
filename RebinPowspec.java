package gb.esac.timing;


import gb.esac.tools.Analysis;
import gb.esac.tools.Binner;
import gb.esac.io.DataFileReader;
import gb.esac.io.DataFileWriter;

import java.io.IOException;

import nom.tam.util.ArrayFuncs;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;


public class RebinPowspec {


    static Logger logger  = Logger.getLogger(RebinPowspec.class);


    public static void main (String[] args) throws Exception {

	PropertyConfigurator.configure("/Users/gbelanger/javaProgs/gb/esac/logger.config");


	/**  Handle arguments  **/
	DataFileReader dataFile = new DataFileReader("pow.qdp");
	int rebinNpoints = 10;
	String binType = "papadakis";
	if ( args.length > 3 || args.length < 1 ) {
	    logger.error("Usage: java RebinPowspec filename (rebinNpoints [10]) (binType: lin, log, [papadakis])");
	    System.exit(-1);
	}
	else if ( args.length == 3 ) {
	    dataFile = new DataFileReader(args[0]);
	    rebinNpoints = (Integer.valueOf(args[1])).intValue();
	    binType = args[2];
	    if ( ! binType.equals("papadakis") && ! binType.equals("log") && ! binType.equals("lin") ) {

		throw new TimingException("binType must be 'papadakis', 'log' or 'lin'");
	    }
	}
	else if ( args.length == 2 ) {
	    dataFile = new DataFileReader(args[0]);
	    rebinNpoints = (Integer.valueOf(args[1])).intValue();
	}
	else if ( args.length == 1 )  dataFile = new DataFileReader(args[0]);
	else  {
	    logger.info("Usage: java RebinPowspec filename (rebinNpoints [10]) (binType [papadakis])");
	    logger.warn("Using defaults");
	}
	logger.info("Data filename = "+args[0]);
	logger.info("Number of points to rebin = "+rebinNpoints);
	logger.info("Bin type = "+binType);


	/**  Get the data  **/
	logger.info("Getting the data from "+args[0]);
	double[] freq = dataFile.getDblCol(0);
	double[] pow = dataFile.getDblCol(1);
	int nDataBins = freq.length;


	/**  Construct the old bin edges assuming adjacent bins  **/
	double oldBinWidth = freq[1] - freq[0];
	double halfBinWidth = 0.5*oldBinWidth;
	double[] oldBinEdges = new double[2*nDataBins];
	for ( int i=0; i < nDataBins; i++ ) {
	    oldBinEdges[2*i] = freq[i] - halfBinWidth;
	    oldBinEdges[2*i+1] = freq[i] + halfBinWidth;
	}


	/**  Rebin the spectrum  **/
	logger.info("Constructing bins and rebinning");
	double delta = freq[1] - freq[0];
	double xmin = freq[0] - delta/2.0;
	double xmax = freq[nDataBins-1] + delta/2.0;
	int nNewBins = 0;
	double[] newFreq = null;
	double[] newPow = null;
	double rebinFactor = rebinNpoints;
	double[] binEdges = null;
	double[][] papadakis = null;
	double binSize = 0;
	boolean directSum = true;
	if ( binType.equals("log") ) {
	    binSize = Math.log(freq[1]/freq[0]) + Math.log(rebinFactor);
	    nNewBins = (new Double(Math.floor( (Math.log(xmax) - Math.log(xmin))/binSize ))).intValue();
	    newFreq = Binner.getMidPoints(xmin, xmax, nNewBins, binType);
	    binEdges = Binner.makeBins(xmin, xmax, nNewBins, binType);
   	    newPow = Binner.rebinData(freq, pow, binEdges);
       	}
	else if ( binType.equals("lin") ) {
	    binSize = (freq[1] - freq[0])*rebinFactor;
	    nNewBins = (new Double(Math.floor( nDataBins/rebinFactor ))).intValue();
	    newFreq = Binner.getMidPoints(xmin, xmax, nNewBins, binType);
	    binEdges = Binner.makeBins(xmin, xmax, nNewBins, binType);
//  	    newPow = Binner.rebinData(freq, pow, binEdges);
 	    newPow = Binner.rebinRatesSimple(pow, oldBinEdges, binSize, directSum);
	}
	else {
	    papadakis = Binner.rebinAsPapadakis(freq, pow, rebinNpoints);
	    newFreq = papadakis[0];
	    newPow = papadakis[1];
	}

	
	/**  Write the result  **/
	String outFilename = "rebPow.qdp";
	DataFileWriter output = new DataFileWriter(outFilename);
	String[] header = new String[] {
		"! QDP File",
		"DEV /XS",
		"READ 1 2",
		"TIME OFF",
		"LAB T",
		"LAB F",
		"LW 3",
		"CS 1.3",
		"LOG ON",
		"LAB X Frequency (Hz)",
		"LAB Y Power (Leahy)",
		"VIEW 0.1 0.2 0.9 0.8",
		"!"
	    };
	try { output.writeData(header, newFreq, newPow); }
	catch (IOException e) { System.out.println("Error: Cannot write output file"); }
	logger.info("Result written to "+outFilename);

    }


}
