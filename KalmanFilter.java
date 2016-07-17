package gb.esac.timing;

import gb.esac.io.AsciiDataFileReader;
import gb.esac.io.AsciiDataFileWriter;



public class KalmanFilter {

    public static void main(String[] args) throws Exception {

	//  Handle argument
	AsciiDataFileReader inputFile = null;
	if ( args.length != 1 ) {
	    System.out.println("Usage: java KalmanFilter filename");
	    System.exit(-1);
	}
	else {
	    inputFile = new AsciiDataFileReader(args[0]);
	}


	//  Get data from input file 
	int nrows = inputFile.getNDataRows();
	int ncols = inputFile.getNDataCols();
	System.out.println("Log  : Data file has nrows="+nrows+", ncols="+ncols);

	double[] t = inputFile.getDblCol(0);
	double[] sigt = inputFile.getDblCol(1);
	double[] z = inputFile.getDblCol(2);
	double[] sigz = inputFile.getDblCol(3);

	//  Define Kalman variables
	double[] xhat = new double[nrows];
	double[] sigx = new double[nrows];
	double k = 0;

	//  Apply Kalman filter
	xhat[0] = z[0];
	sigx[0] = sigz[0];
	for ( int i=0; i < nrows-1; i++ ) {
	    k = Math.pow(sigz[i], 2) / ( Math.pow(sigz[i],2) + Math.pow(sigz[i+1],2) );
	    xhat[i+1] = xhat[i] + k*(z[i+1] - xhat[i]);
	    sigx[i+1] = Math.sqrt( Math.pow(sigx[i],2) - k*Math.pow(sigx[i],2) );
	}

	//  Write Kalman filtered LC as QDP file 
	AsciiDataFileWriter outputFile = new AsciiDataFileWriter("kalmanLC.qdp");
	String[] header = new String[] {
	  "! QDP File", 
	  "DEV /XS", 
	  "READ SERR 1 2",
	  "LINE ON",
	  "LAB T ", 
	  "LAB F ", 
	  "TIME OFF", 
	  "LAB X Time (s)", 
	  "LAB Y Count Rate (cts/s)",
	  "VIEW 0.1 0.2 0.9 0.8",
	  "LW 3",
	  "CS 1.3",
	  "!" 
	};
	outputFile.writeData(header, t, sigt, xhat, sigx);


// 	//  Write result to QDP file showing both the original and filtered lc
// 	DataFile outputFile = new DataFile("filtData.qdp");
// 	String[] header = new String[] {
// 	  "! QDP File", 
// 	  "DEV /XS", 
// 	  "READ SERR 3",
// 	  "LAB T ", 
// 	  "LAB F ", 
// 	  "TIME OFF", 
// 	  "LAB X Time (s)", 
// 	  "LAB Y Count Rate (cts/s)",
// 	  "VIEW 0.1 0.15 0.9 0.85",
// 	  "LW 3",
// 	  "CS 1.4",
// 	  "PLOT VERT",
// 	  "!" 
// 	};
// 	//outputFile.writeData(header, t, sigt, xhat, sigx);
// 	outputFile.writeData(header, t, xhat, z, sigz);

    }
}
