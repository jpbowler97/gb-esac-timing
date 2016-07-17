package gb.esac.timing;

import cern.jet.random.engine.MersenneTwister64;
import gb.esac.eventlist.EventList;
import gb.esac.montecarlo.RedNoiseGenerator;
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.periodogram.PeriodogramUtils;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.tools.BasicStats;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Date;


public class SpectralIndexSurface {


    public static void main(String[] args) throws Exception {

	DecimalFormat number = new DecimalFormat("0.0000");
	DecimalFormat fit = new DecimalFormat("0.00");
	DecimalFormat sci = new DecimalFormat("0.0#E0");

	double[] durations = new double[] {1e4};
	for ( int m=0; m < durations.length; m++ ) {

	    double duration = durations[m];
	    double lengthOfSegment = duration;
	    System.out.println("Duration = "+sci.format(duration));

	    int nEventLists = 30;
	    double[] indexes = new double[nEventLists];

	    //  Define the values of the true index that will be used in the simulations
	    double[] trueIndexes = new double[] {0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5};
	    //trueIndexes = new double[] {1.0, 1.5, 2.0, 2.5, 3.0};

	    //  Define log-spaced count rates
	    double[] countRates = new double[] {0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50, 100};
	    //countRates = new double[] {0.2, 1, 5, 20, 100};
	    
	    String outName = "xyz_T_"+sci.format(duration)+"_alpha_"+trueIndexes[0]+"-"+trueIndexes[trueIndexes.length-1]+"_cps_"+countRates[0]+"-"+countRates[countRates.length-1]+".dat";
	    PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outName)));
	
	    //  Loop on true spectral indices
	    for ( int i=0; i < trueIndexes.length; i++ ) {
		    
		double trueIndex = trueIndexes[i];
		System.out.println("Index = "+trueIndex);
	    
		//  Loop on the count rates
		MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
		for ( int j=0; j < countRates.length; j++ ) {
		
		    double countRate = countRates[j];
		    System.out.print("   Count rate = "+countRate+" ... generating event list");

		    //  Generate a single long list of arrival times
// 		    double longDuration = duration*nEventLists;
// 		    double[] arrivalTimes = RedNoiseGenerator.generateArrivalTimes(countRate, longDuration, trueIndex, engine);
// 		    EventList evlist = new EventList(arrivalTimes);
// 		    arrivalTimes = null;

		    //  Loop on the segments
		    int k=0;
		    while ( k < nEventLists ) {

			double[] times = RedNoiseGenerator.generateArrivalTimes(countRate, duration, trueIndex, engine);

// 			double from = k*lengthOfSegment;
// 			double to = (k+1)*lengthOfSegment;
// 			double[] times = evlist.getArrivalTimesFromTo(from, to);

 			TimeSeries ts = TimeSeriesMaker.makeTimeSeries(times, 128);
 			FFTPeriodogram psd = PeriodogramMaker.makeFFTPeriodogram(ts);
 			double[] fitRes = PeriodogramUtils.fitPowerLawInLinearSpace(psd);
 			indexes[k] = fitRes[0];
 			k++;

		    }
		    double[] avgAndVar = BasicStats.getRunningAveAndVar(indexes);
		    double avgIndex = -avgAndVar[0];
		    double error = Math.sqrt(avgAndVar[1]);

		    System.out.println(":  Estimated index = "+number.format(avgIndex)+" +/- "+number.format(error));
		    pw.println(trueIndex+"\t"+countRate+"\t"+number.format(avgIndex)+"\t"+number.format(error));
		    pw.flush();

		}
		pw.println();
		pw.flush();
	    
	    }
	    pw.close();
	}

    }

}
