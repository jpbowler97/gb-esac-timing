package gb.esac.timing;

import java.io.File;
import java.io.IOException;

import gb.esac.io.AsciiDataFileWriter;
import gb.esac.tools.FitsUtils;
import gb.esac.tools.DataUtils;

import nom.tam.fits.BinaryTableHDU;
import nom.tam.fits.FitsException;


public class MakeLC {

    public static void main(String[] args) throws Exception, FitsException {


	/**  Handle args  **/
	String srcFilename = null;
	String bkFilename = null;
	File srcFile, bkFile = null;
	double bintime = 0;
	boolean resetToZero = false;
	double startTime = 0;

	if ( args.length == 2 ) {

	    srcFilename = args[0];
	    srcFile = new File(srcFilename);
	    if ( ! srcFile.exists() ) {
		System.out.println("Error: File "+srcFilename+" does not exist");
		System.exit(-1);
	    }
	    try {
		bintime = (Double.valueOf(args[1])).doubleValue();
	    }
	    catch ( NumberFormatException e ) {
		System.out.println("Usage: java gb.esac.timing.MakeLC srcEvlist.fits (bkEvlist.fits) bintime (startTime [0]) (resetToZero, y or [n])");
		System.exit(-1);
	    }
	}

	else if ( args.length == 3 ) {

	    srcFilename = args[0];
	    srcFile = new File(srcFilename);
	    if ( ! srcFile.exists() ) {
		System.out.println("Error: File "+srcFilename+" does not exist");
		System.exit(-1);
	    }
	    bkFilename = args[1];
	    bkFile = new File(args[1]);
	    if ( ! bkFile.exists() ) {
		try {
		    bintime = (Double.valueOf(args[1])).doubleValue();
		}
		catch ( NumberFormatException e ) {
		    System.out.println("Usage: java gb.esac.timing.MakeLC srcEvlist.fits (bkEvlist.fits) bintime (startTime) (resetToZero, y or [n])");
		    System.exit(-1);
		}
		try {
		    startTime = (Double.valueOf(args[2])).doubleValue();
		}
		catch ( NumberFormatException e ) {
		    if ( args[2].equals("y") ) {
			resetToZero = true;
		    }
		    else if ( args[2].equals("n") ) {
			resetToZero = false;
		    }
		    else {
			System.out.println("Error: resetToZero must be 'y' or 'n'");
			System.exit(-1);
		    }

		}
	    }
	    else {
		try {
		    bintime = (Double.valueOf(args[2])).doubleValue();
		}
		catch ( NumberFormatException e ) {
		    System.out.println("Usage: java gb.esac.timing.MakeLC srcEvlist.fits (bkEvlist.fits) bintime (startTime) (resetToZero, y or [n])");
		    System.exit(-1);
		}
	    }
	    
	}

	else if ( args.length == 4 ) {

	    srcFilename = args[0];
	    srcFile = new File(srcFilename);
	    if ( ! srcFile.exists() ) {
		System.out.println("Error: File "+srcFilename+" does not exist");
		System.exit(-1);
	    }
	    bkFilename = args[1];
	    bkFile = new File(args[1]);
	    if ( ! bkFile.exists() ) {
		try {
		    bintime = (Double.valueOf(args[1])).doubleValue();
		}
		catch ( NumberFormatException e ) {
		    System.out.println("Usage: java gb.esac.timing.MakeLC srcEvlist.fits (bkEvlist.fits) bintime (startTime) (resetToZero, y or [n])");
		    System.exit(-1);
		}
		try {
		    startTime = (Double.valueOf(args[2])).doubleValue();
		}
		catch ( NumberFormatException e ) {
		    System.out.println("Usage: java gb.esac.timing.MakeLC srcEvlist.fits (bkEvlist.fits) bintime (startTime) (resetToZero, y or [n])");
		    System.exit(-1);
		}
		if ( args[3].equals("y") )
		    resetToZero = true;
		else if ( args[3].equals("n") )
		    resetToZero = false;
		else {
		    System.out.println("Error: resetToZero must be 'y' or 'n'");
		    System.exit(-1);
		}
	    }
	    else {
		try {
		    bintime = (Double.valueOf(args[2])).doubleValue();
		}
		catch ( NumberFormatException e ) {
		    System.out.println("Usage: java gb.esac.timing.MakeLC srcEvlist.fits (bkEvlist.fits) bintime (startTime) (resetToZero, y or [n])");
		    System.exit(-1);
		}
		try {
		    startTime = (Double.valueOf(args[3])).doubleValue();
		}
		catch ( NumberFormatException e ) {
		    if ( args[3].equals("y") )
			resetToZero = true;
		    else if ( args[3].equals("n") )
			resetToZero = false;
		    else {
			System.out.println("Error: resetToZero must be 'y' or 'n'");
			System.exit(-1);
		    }

		}
	    }
	}	

	else if ( args.length == 5 ) {

	    srcFilename = args[0];
	    bkFilename = args[1];
	    bintime = (Double.valueOf(args[2])).doubleValue();
	    startTime = (Double.valueOf(args[3])).doubleValue();
	    resetToZero = (Boolean.valueOf(args[4])).booleanValue();
	}
	else {
	    System.out.println("Usage: java gb.esac.timing.MakeLC srcEvlist.fits (bkEvlist.fits) bintime (startTime) (resetToZero [false])");
	    System.exit(-1);
	}


	/**  Get the data  **/
	BinaryTableHDU srcEventsHDU = null;
	double[] srcEventTimes = null;
	BinaryTableHDU bkEventsHDU = null;
	double[] bkEventTimes = null;
	srcEventsHDU = FitsUtils.getEventsHDU(srcFilename);
	srcEventTimes = (double[]) srcEventsHDU.getColumn("TIME");

	if ( bkFile.exists() ) {
	    srcEventsHDU = FitsUtils.getEventsHDU(srcFilename);
	    srcEventTimes = (double[]) srcEventsHDU.getColumn("TIME");
	    double tmpMax = Math.max(bkEventTimes[0], srcEventTimes[0]);
	    startTime = Math.max(startTime, tmpMax);
	}
	if ( startTime == 0 ) {
	    startTime = Math.max(startTime, srcEventTimes[0]);
	}


	/**  Make the LC  **/
	Object[] srcLC = TimingUtils.makeLC(srcEventTimes, bintime, startTime);
	double[] lcBinCentres = (double[]) srcLC[0];
	double[] rates = (double[]) srcLC[1];
	double[] errors = (double[]) srcLC[2];
	double[] bkSubRates = rates;
	double[] bkSubErrors = errors;
	int nCommonBins = rates.length;
	if ( bkFile.exists() ) {
	    Object[] bkLC = TimingUtils.makeLC(bkEventTimes, bintime, startTime);
	    double[] bkBinCentres = (double[]) bkLC[0];
	    double[] bkRates = (double[]) bkLC[1];
	    double[] bkErrors = (double[]) bkLC[2];
	    nCommonBins = (new Double(Math.min(rates.length, bkRates.length))).intValue();
	    for ( int i=0; i < nCommonBins; i++ ) {
		bkSubRates[i] = rates[i] - bkRates[i];
		bkSubErrors[i] = Math.sqrt( Math.pow(errors[i], 2) + Math.pow(bkErrors[i], 2));
	    }
	}


	/**  Reset LC bin centre times to zero  **/
	double[] binCentres = lcBinCentres;
	if ( resetToZero ) {
	    binCentres = DataUtils.resetToZero(lcBinCentres, bintime/2);
	}


	/**  Write LC  **/
	AsciiDataFileWriter lcFile = null;
	String outFilename = "lc.qdp";
	try { lcFile = new AsciiDataFileWriter(outFilename); }
	catch (IOException e) {
		System.out.println("Error [MakeLC]: Could not open or create "+outFilename); 
		System.exit(-1);
	}
	String[] header = new String[] {
	    "DEV /XS",
	    "READ SERR 2",
	    "TIME OFF", "LINE ON",
	    "LAB T", "LAB F", "LW 3", "CS 1.3",
	    "LAB Y Rate (cts/s)",
	    "LAB X Time (s)",
	    "VIEW 0.1 0.2 0.9 0.8",
	    "!"
	};
	try { lcFile.writeData(header, binCentres, bkSubRates, bkSubErrors); }
	catch (IOException e) {
	    System.out.println("Error: Could not write to file "+outFilename);
	    System.exit(-1);
	}

    }
}
