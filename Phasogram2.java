package gb.esac.timing;

import gb.esac.tools.Stats;
import gb.esac.tools.Astro;
import gb.esac.tools.Analysis;
import gb.esac.tools.Binner;
import gb.esac.io.DataFile;

import java.io.IOException;
import java.io.FileNotFoundException;
import java.text.DecimalFormat;
import java.util.Vector;
import nom.tam.fits.*;
import hep.aida.*;
import hep.aida.ref.histogram.*;
import cern.jet.stat.Descriptive;
import cern.colt.list.DoubleArrayList;


public class Phasogram2 {

    public static void main(String[] args) throws Exception {

	DecimalFormat numberFormat = new DecimalFormat("0.000");

	/**   Handle arguments  **/
	String evlistName = null;
	double period = 0;
	int nPhaseBins = 15;
	boolean fileIsFits = true;
	if ( args.length == 3 ) {
	    evlistName = args[0];
	    period = (Double.valueOf(args[1])).doubleValue();
	    nPhaseBins = (Integer.valueOf(args[2])).intValue();
	}
	else if ( args.length == 2 ) {
	    evlistName = args[0];
	    period = (Double.valueOf(args[1])).doubleValue();
	}
	else {
	    System.out.println("Error: 2 args required, (1 optional)");
	    System.out.println("Log  : arg 1 = event list filename ");
	    System.out.println("Log  : arg 2 = period (s)");
	    System.out.println("Log  : (arg 3) = number of phase bins [15]");
	    System.exit(-1);
	}


	/**   Check if file is ascii or fits  **/
	if ( evlistName.endsWith("fits") || evlistName.endsWith("fits.gz") )
	    fileIsFits = true;
	else if ( evlistName.endsWith("dat") || evlistName.endsWith("txt") )
	    fileIsFits = false;
	else {
	    System.out.println("Error: File format not recognized");
	    System.out.println("Log  : filename extension can be 'fits, 'fits.gz', 'txt' or 'dat'");
	    System.exit(-1);
	}
	

	/**  Get time data from event list  **/
	Fits fitsFile;
	double[] times = null;
	short[] energies = null;
	if ( fileIsFits ) {
	    fitsFile = Astro.openFits(evlistName); 
	    BinaryTableHDU hdu = Astro.getEventsHDU(evlistName);
	    times = (double[]) hdu.getColumn("TIME");
	    energies = (short[]) hdu.getColumn("PI");
	}
	else {
	    DataFile dataFile = new DataFile(evlistName);
	    times = dataFile.getDblCol(0);
	}

	
	/**  Select low energy (2--4keV) and high energy (4--10keV) times  **/
	Vector timesLowE = new Vector();
	Vector timesHighE = new Vector();
	int nRejected = 0;
	int emin = 2000;
	int emid = 4000;
	int emax = 10000;
	for ( int i=0; i < times.length; i++ ) {
	    if ( energies[i] >= emin && energies[i] < emid )  timesLowE.addElement(times[i]);
	    else if ( energies[i] >= emid && energies[i] <= emax )  timesHighE.addElement(times[i]);
	    else  nRejected++;
	}
	System.out.println("Log  : "+nRejected+" event are outside the range 2--10 keV");
	

	/**  Convert the Vectors to Arrays  **/
	double[] timesLE = new double[timesLowE.size()];
	for ( int i=0; i < timesLE.length; i++ ) 
	    timesLE[i] = ((Double) timesLowE.elementAt(i)).doubleValue();
	double[] timesHE = new double[timesHighE.size()];	
	for ( int i=0; i < timesHE.length; i++ ) 
	    timesHE[i] = ((Double) timesHighE.elementAt(i)).doubleValue();


	/**  Determine mean count rate  **/
	double duration = Stats.getRange(times);
	double countRate = times.length/duration;
	double countRateLE = timesLE.length/duration;
	double countRateHE = timesHE.length/duration;
	System.out.println("Log  : Duration = "+numberFormat.format(duration)+" s");
	System.out.println("Log  : Mean count rates are: ");
	System.out.println("Log  : Global (2 - 10 keV): "+numberFormat.format(countRate)+" cts/s");
	System.out.println("Log  : Low energy (2 - 4 keV): "+numberFormat.format(countRateLE)+" cts/s");
	System.out.println("Log  : High energy (4 - 10 keV): "+numberFormat.format(countRateHE)+" cts/s");


	/**  Construct light curve  **/
	double bintime = 120D;
	int nbins = (new Double(Math.ceil(duration/bintime))).intValue();
	double lengthOfLcBin = duration/(new Double(nbins)).doubleValue();
	double tfirst = 0;
	double tzero = times[0];
	double tlast = times[times.length-1] - tzero;
	double[] binnedTimes = Binner.binOrderedData(times, nbins);
	double[] lcHeights = new double[nbins];
	double[] lcHeightsErr = new double[nbins];
	double[] lcValues = new double[nbins];
	double[] lcValuesErr = new double[nbins];
	for ( int i=0; i < nbins; i++ ) {
	    lcHeightsErr[i] = Math.sqrt(binnedTimes[i])/lengthOfLcBin;
	    lcHeights[i] = binnedTimes[i]/lengthOfLcBin;
	    lcValues[i] = i*(duration/(new Double(nbins)).doubleValue());
	    lcValuesErr[i] = duration/(new Double(2*nbins)).doubleValue();
	}
	

	/**  Write light curve as QDP file  **/
	String lcFilename = "lightcurve.qdp";
	DataFile lightcurveFile = new DataFile(lcFilename);
	double minLcHeight = Stats.getMin(lcHeights);
	double minLcHeightErr = lcHeightsErr[Stats.getIndex(lcHeights, minLcHeight)];
	double maxLcHeight = Stats.getMax(lcHeights);
	double maxLcHeightErr = lcHeightsErr[Stats.getIndex(lcHeights, maxLcHeight)];
	String ylowLC = (new Double(minLcHeight - 2*minLcHeightErr)).toString();
	String yhighLC = (new Double(maxLcHeight + 2*maxLcHeightErr)).toString();
	String[] headerLC = new String[] {
	    "! QDP data file : "+lcFilename, 
	    "DEV /XS",
	    "READ SERR 1 2",
	    "LABEL TITLE", 
	    "LABEL FILE",
	    "TIME OFF", 
	    "LWIDTH 3", 
	    "CSIZE 1.2",
	    "LINE ON",
	    "LABEL X Time (s)",
	    "LABEL Y Count Rate (cts/s)",
	    "RANGE Y "+ylowLC+" "+yhighLC,
	    "VIEW 0.1 0.2 0.9 0.8"
	};
	lightcurveFile.writeData(headerLC, lcValues, lcValuesErr, lcHeights, lcHeightsErr);


	/**  Construct  Phasograms  **/
	double[] phases = TimingUtils.getPhases(times, period);
	double[] phasesLE = TimingUtils.getPhases(timesLE, period);
	double[] phasesHE = TimingUtils.getPhases(timesHE, period);
	double[] binnedPhases = Binner.binData(phases, nPhaseBins);
	double[] binnedPhasesLE = Binner.binData(phasesLE, nPhaseBins);
	double[] binnedPhasesHE = Binner.binData(phasesHE, nPhaseBins);

	double xmin = 0;
	double xmax = 1;
	FixedAxis phaseAxis = new FixedAxis(nPhaseBins, xmin, xmax);
	Histogram1D phaseHisto = new Histogram1D("", "Global Phasogram", phaseAxis);
	Histogram1D phaseHistoLE = new Histogram1D("", "Low Energy Phasogram", phaseAxis);
	Histogram1D phaseHistoHE = new Histogram1D("", "High Energy Phasogram", phaseAxis);
	Vector[] energiesInBin = new Vector[nPhaseBins];
	for ( int i=0; i < nPhaseBins; i++ )  energiesInBin[i] = new Vector();
	double binWidth = 1D/nPhaseBins;
 	int binIndex = 0;
 	for ( int i=0; i < phases.length; i++ ) {
 	    phaseHisto.fill(phases[i]);
 	    binIndex = (new Double(Math.floor(phases[i]/binWidth))).intValue();
 	    energiesInBin[binIndex].addElement(energies[i]/1000D);  // convert to keV
  	}
 	for ( int i=0; i < timesLE.length; i++ ) phaseHistoLE.fill(phasesLE[i]);
 	for ( int i=0; i < timesHE.length; i++ ) phaseHistoHE.fill(phasesHE[i]);



	/**  Get phasogram info from phase Histos
	 *    Also get mean energies and errors for each bin  **/
	double[] phaseHeights = new double[nPhaseBins*2];
	double[] phaseHeightsErr = new double[nPhaseBins*2];
	double[] phaseHeightsLE = new double[nPhaseBins*2];
	double[] phaseHeightsErrLE = new double[nPhaseBins*2];
	double[] phaseHeightsHE = new double[nPhaseBins*2];
	double[] phaseHeightsErrHE = new double[nPhaseBins*2];
	double[] phaseValues = new double[nPhaseBins*2];
	double[] phaseValuesErr = new double[nPhaseBins*2];

	double[] meanEnergies = new double[nPhaseBins*2];
	double[] errorsOnMeanEnergies = new double[nPhaseBins*2];
	double lengthOfPhaseBin = duration/(new Double(nPhaseBins)).doubleValue();
// 	double meanPhaseHeight = Stats.getMean(phases);
// 	int idxOfMeanPhaseHeight = Stats.getIndex(phases, meanPhaseHeight);
// 	double phaseOfMeanHeight = phases[idxOfMeanPhaseHeight];


	double[] phaso = new double[2*nPhaseBins];
	double[] phasoErr = new double[2*nPhaseBins];
	double[] phasoLE = new double[2*nPhaseBins];
	double[] phasoErrLE = new double[2*nPhaseBins];
	double[] phasoHE = new double[2*nPhaseBins];
	double[] phasoErrHE = new double[2*nPhaseBins];

	for ( int i=0; i < nPhaseBins; i++ ) {
	    phaso[i] = binnedPhases[i]/lengthOfPhaseBin;
	    phasoErr[i] = Math.sqrt(binnedPhases[i])/lengthOfPhaseBin;
	    phasoLE[i] = binnedPhasesLE[i]/lengthOfPhaseBin;
	    phasoErrLE[i] = Math.sqrt(binnedPhasesLE[i])/lengthOfPhaseBin;
	    phasoHE[i] = binnedPhasesHE[i]/lengthOfPhaseBin;
	    phasoErrHE[i] = Math.sqrt(binnedPhasesHE[i])/lengthOfPhaseBin;

	    phaso[i+nPhaseBins] = binnedPhases[i]/lengthOfPhaseBin;
	    phasoErr[i+nPhaseBins] = Math.sqrt(binnedPhases[i])/lengthOfPhaseBin;
	    phasoLE[i+nPhaseBins] = binnedPhasesLE[i]/lengthOfPhaseBin;
	    phasoErrLE[i+nPhaseBins] = Math.sqrt(binnedPhasesLE[i])/lengthOfPhaseBin;
	    phasoHE[i+nPhaseBins] = binnedPhasesHE[i]/lengthOfPhaseBin;
	    phasoErrHE[i+nPhaseBins] = Math.sqrt(binnedPhasesHE[i])/lengthOfPhaseBin;
	}


	for ( int i=0; i < nPhaseBins; i++ ) {
	    
	    /**  Global 2 - 10 keV phase histo  **/
	    phaseHeights[i] = phaseHisto.binHeight(i)/lengthOfPhaseBin;
	    phaseHeights[i+nPhaseBins] = phaseHeights[i];
	    phaseHeightsErr[i] = phaseHisto.binError(i)/lengthOfPhaseBin;
	    phaseHeightsErr[i+nPhaseBins] = phaseHeightsErr[i];
// 	    phaseValues[i] = i/(new Double(nPhaseBins)).doubleValue() - phaseOfMeanHeight;
// 	    if ( phaseValues[i] < 0 ) phaseValues[i] += 1;
	    phaseValues[i] = i/(new Double(nPhaseBins)).doubleValue();
	    phaseValues[i+nPhaseBins] = 1 + phaseValues[i];
	    phaseValuesErr[i] = 1/(new Double(2*nPhaseBins)).doubleValue();
	    phaseValuesErr[i+nPhaseBins] = phaseValuesErr[i];

	    /**  Low Energy  **/
	    phaseHeightsLE[i] = phaseHistoLE.binHeight(i)/lengthOfPhaseBin;
	    phaseHeightsLE[i+nPhaseBins] = phaseHeightsLE[i];
	    phaseHeightsErrLE[i] = phaseHistoLE.binError(i)/lengthOfPhaseBin;
	    phaseHeightsErrLE[i+nPhaseBins] = phaseHeightsErrLE[i];

	    /**  High Energy  **/
	    phaseHeightsHE[i] = phaseHistoHE.binHeight(i)/lengthOfPhaseBin;
	    phaseHeightsHE[i+nPhaseBins] = phaseHeightsHE[i];
	    phaseHeightsErrHE[i] = phaseHistoHE.binError(i)/lengthOfPhaseBin;
	    phaseHeightsErrHE[i+nPhaseBins] = phaseHeightsErrHE[i];

	    //System.out.println((phaseHeights[i]-binnedPhases[i]/lengthOfPhaseBin) +"\t"+
	    //		       (phaseHeightsLE[i]-binnedPhasesLE[i]/lengthOfPhaseBin));

	    /**  Mean Energy as a function of phase  **/
	    meanEnergies[i] = Stats.getMean(energiesInBin[i]);
	    meanEnergies[i+nPhaseBins] = meanEnergies[i];
	    errorsOnMeanEnergies[i] = Stats.getErrOnMean(energiesInBin[i]);
	    errorsOnMeanEnergies[i+nPhaseBins] = errorsOnMeanEnergies[i];
	}

	
	/**  Determine Y-axis bounds for plot  **/
	String phasoFilename = "phasogram.qdp";
	DataFile phasogramFile = new DataFile(phasoFilename);
	double minPhaseHeight = Stats.getMin(phaseHeights);
	double minPhaseHeightErr = phaseHeightsErr[Stats.getIndex(phaseHeights, minPhaseHeight)];
	double maxPhaseHeight = Stats.getMax(phaseHeights);
	double maxPhaseHeightErr = phaseHeightsErr[Stats.getIndex(phaseHeights, maxPhaseHeight)];
	String ylow = (new Double(minPhaseHeight - 2*minPhaseHeightErr)).toString();
	String yhigh = (new Double(maxPhaseHeight + 2*maxPhaseHeightErr)).toString();


// 	/**  Re-order so that phase=0 corresponds to the minimum height  **/
//  	int indexOfMinHeight = Stats.getIndex(phaseHeights, minPhaseHeight);
// 	double[] orderedPhaseHeights = new double[phaseHeights.length];
// 	double[] orderedPhaseHeightsErr = new double[phaseHeights.length];
// 	double[] orderedPhaseHeightsHE = new double[phaseHeights.length];
// 	double[] orderedPhaseHeightsErrHE = new double[phaseHeights.length];
// 	double[] orderedPhaseHeightsLE = new double[phaseHeights.length];
// 	double[] orderedPhaseHeightsErrLE = new double[phaseHeights.length];
// 	for ( int i=0; i < phaseHeights.length; i++ ) {
// 	    int newIndex = indexOfMinHeight + i;
// 	    if ( newIndex > (phaseHeights.length-1) )  newIndex -= (phaseHeights.length-1);
// 	    orderedPhaseHeights[i] = phaseHeights[newIndex];
// 	    orderedPhaseHeightsErr[i] = phaseHeightsErr[newIndex];
// 	    orderedPhaseHeightsHE[i] = phaseHeightsHE[newIndex];
// 	    orderedPhaseHeightsErrHE[i] = phaseHeightsErrHE[newIndex];
// 	    orderedPhaseHeightsLE[i] = phaseHeightsLE[newIndex];
// 	    orderedPhaseHeightsErrLE[i] = phaseHeightsErrLE[newIndex];
// 	}


	/**  Write 'Count Rate vs. Phase' as QDP file  **/
	String[] header = new String[] {
	    "! QDP data file : "+phasoFilename, 
	    "DEV /XS",
	    "READ SERR 1 2 3 4",
	    "LABEL TITLE", 
	    "LABEL FILE",
	    "TIME OFF", 
	    "LWIDTH 3",
	    "CSIZE 1.2",
	    "LINE STEP",
	    "LABEL X Phase",
	    "RANGE X -0.05 2.05",
	    "RANGE Y "+ylow+" "+yhigh,
	    "LABEL Y3 Count Rate (cts/s)",
	    //"VIEW 0.1 0.1 0.9 0.9",
	    "PLOT VERT"
	};
	//phasogramFile.writeData(header, phaseValues, phaseValuesErr, phaso, phasoErr, phasoHE, phasoErrHE, phasoLE, phasoErrLE);
 	phasogramFile.writeData(header, phaseValues, phaseValuesErr, phaseHeights, phaseHeightsErr, phaseHeightsHE, phaseHeightsErrHE, phaseHeightsLE, phaseHeightsErrLE);
// 	phasogramFile.writeData(header, phaseValues, phaseValuesErr, orderedPhaseHeights, orderedPhaseHeightsErr, orderedPhaseHeightsHE, orderedPhaseHeightsErrHE, orderedPhaseHeightsLE, orderedPhaseHeightsErrLE);


	/**  Write mean energy Vs. phase as QDP file  **/
	String energyPhasoFilename = "energyPhaso.qdp";
	DataFile energyPhasoFile = new DataFile(energyPhasoFilename);
	minPhaseHeight = Stats.getMin(meanEnergies);
	minPhaseHeightErr = errorsOnMeanEnergies[Stats.getIndex(meanEnergies, minPhaseHeight)];
	maxPhaseHeight = Stats.getMax(meanEnergies);
	maxPhaseHeightErr = errorsOnMeanEnergies[Stats.getIndex(meanEnergies, maxPhaseHeight)];
	ylow = (new Double(minPhaseHeight - 2*minPhaseHeightErr)).toString();
	yhigh = (new Double(maxPhaseHeight + 2*maxPhaseHeightErr)).toString();
	header = new String[] {
	    "! QDP data file: "+energyPhasoFilename,
	    "DEV /XS",
	    "READ SERR 1 2",
	    "LABEL TITLE", 
	    "LABEL FILE",
	    "TIME OFF", 
	    "LWIDTH 3",
	    "CSIZE 1.2",
	    "LINE STEP",
	    "LABEL X Phase",
	    "LABEL Y Mean Photon Energy (keV)",
	    "RANGE X -0.05 2.05",
	    "RANGE Y "+ylow+" "+yhigh,
	    "VIEW 0.1 0.2 0.9 0.8"
	};
	energyPhasoFile.writeData(header, phaseValues, phaseValuesErr, meanEnergies, errorsOnMeanEnergies);
    }
}
