package gb.esac.timing;

import gb.esac.io.DataFile;

import java.io.IOException;
import java.text.DecimalFormat;


public class AveragePowspec {

    public static DecimalFormat sci = new DecimalFormat("0.000E0");
    public static DecimalFormat flt = new DecimalFormat("0.000");

    public static void main (String[] args) throws Exception {


	/**  Handle arguments  **/
	DataFile[] dataFiles = null;
	if ( args.length < 2 ) {
	    System.out.println("Usage: java AveragePowspec file1 file2 ... fileN");
	    System.exit(-1);
	}
	else {
	    System.out.println("Log  : Running AveragePowspec");
	    dataFiles = new DataFile[args.length];
	    for ( int i=0; i < args.length; i++ )  {
		DataFile file = new DataFile(args[i]); 
		if ( file.exists() ) {
		    dataFiles[i] = file;
		}
		else {
		    System.out.println("Error: File not found: "+args[i]);
		    System.exit(-1);
		}
	    }
	}
	int nfiles = dataFiles.length;
	System.out.println("Log  : Found "+(nfiles)+" files");

	
	/**  Combine the power spectra  **/
	System.out.println("Log  : Combining power spectra ... ");
	System.out.println("Log  : File 1: "+dataFiles[0].getName());
	double[] refFreq = dataFiles[0].getDblCol(0);
	int nfreqs = refFreq.length;
	double[] pow = dataFiles[0].getDblCol(1);
	if ( pow.length != nfreqs ) {
	    System.out.println("\n"+"Error: All files must have the same number of points");
	    System.exit(-1);
	}
	double[] sum = new double[nfreqs];
	double[] ave = new double[nfreqs];
	double[] var = new double[nfreqs];
	double[] runAve = new double[nfreqs];
	double[] runVar = new double[nfreqs];
	for ( int k=0; k < nfreqs; k++ ) {
	    sum[k] = pow[k];
	    ave[k] = pow[k];
	    var[k] = Math.pow((pow[k] - ave[k]), 2);
	    runAve[k] = 0;
	    runVar[k] = 0;
	}
	double[] freq = null;

	for ( int i=1; i < nfiles; i++ ) {

	    System.out.println("Log  : File "+(i+1)+": "+dataFiles[i].getName());

	    /**  Get frequencies and check if they are the same as first file  **/
	    freq = dataFiles[i].getDblCol(0);
	    if ( freq.length != nfreqs ) {
		System.out.println("\n"+"Error: All files must have the same number of points");
		System.exit(-1);
	    }
// 	    for ( int k=0; k < nfreqs; k++ ) {
// 		if ( freq[k] != refFreq[k] ) {
// 		    System.out.println("\n"+"Error: Frequencies must be indentical in all files");
// 		    System.exit(-1);
// 		}
// 	    }

	    /**  Get the next file's powers  **/
	    pow = dataFiles[i].getDblCol(1);
	    if ( pow.length != nfreqs ) {
		System.out.println("\n"+"Error: All files must have the same number of points");
		System.exit(-1);
	    }

	    /**  Loop through the frequencies  **/
	    for ( int k=0; k < nfreqs; k++ ) {

		/**  Calculate the runAve and runVar from previous values  **/
		runAve[k] = ave[k] + (pow[k] - ave[k])/(i+1);
		runVar[k] = var[k] + (pow[k] - runAve[k])*(pow[k] - ave[k]);
		
		/**  Calculate new values of ave and var  **/
		sum[k] += pow[k];
		ave[k] = sum[k]/(i+1);
		var[k] += Math.pow((pow[k] - ave[k]), 2);

	    }

	}


	/**  Calculate the standard deviation and relative error at each frequency  **/
	double[] sigma = new double[nfreqs];
	double[] relSigma = new double[nfreqs];
	for ( int k=0; k < nfreqs; k++ ) {
	    runVar[k] /= (nfiles-1);
	    sigma[k] = Math.sqrt(runVar[k]/(nfiles));
	    //relSigma[k] = (runAve[k] - sigma[k])/runAve[k];
	}
	
	
	/**  Write the result  **/
	DataFile output = new DataFile("avePow.qdp");
	String[] header = new String[] {
	    "! QDP File",
	    "! col1=freq, col2=avePow, col3=sigma, col4=relSigma",
	    "DEV /XS", 
	    "READ SERR 2",
	    //"PLOT VERT",
	    "TIME OFF","LAB T","LAB F",
	    "LINE ON", "LW 3","CS 1.1","LOG ON",
	    "LOG OFF 4", "LOG X4", "R Y4 -1 1",
	    "LAB X Frequency (Hz)","LAB Y Power (Leahy)",
	    //"VIEW 0.2 0.1 0.8 0.9","!"
	    "VIEW 0.1 0.2 0.9 0.8","!"
	};
	//try { output.writeData(header, refFreq, runAve, sigma, relSigma); }
	try { output.writeData(header, refFreq, runAve, sigma); }
	catch (IOException e) { System.out.println("Error: Cannot write output file"); }
	System.out.println("Log  : Average power spectrum written to "+output.getName());
	
    }
}
