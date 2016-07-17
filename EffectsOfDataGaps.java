package gb.esac.timing;

import gb.esac.periodogram.ModifiedRayleighPeriodogram;
import gb.esac.periodogram.LombScarglePeriodogram;
import gb.esac.timeseries.TimeSeries;
import gb.esac.eventlist.EventListSelector;
import cern.colt.list.DoubleArrayList;
import gb.esac.eventlist.EventList;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.timeseries.TimeSeriesResampler;
import gb.esac.periodogram.PeriodogramUtils;

public class EffectsOfDataGaps {

    public static void main(String[] args) throws Exception {
	String filename = "oversamplingArtefacts_evlist.fits";
	filename = "fftAndR2_evlist.qdp";
	EventList evlist = new EventList(filename);
	double duration = evlist.duration();
	double gapLength = 0.1*duration;
	int nCases = 3;
	for ( int i=0; i <= nCases; i++ ) {
	    int nGaps = i;
	    makePeriodogram(evlist, nGaps, gapLength);
	}
    }

    private static void makePeriodogram(EventList evlist, int nGaps, double gapLength) throws Exception {
	int nSegments = nGaps + 1;
	double duration = evlist.duration();
	double segLength = (duration - nGaps*gapLength)/nSegments;
	double start = 0;
	double end = segLength;
	DoubleArrayList timeList = new DoubleArrayList();
	for ( int i=0; i < nSegments; i++ ) {
	    double[] times = EventListSelector.getArrivalTimesFromTo(evlist, start, end);
	    for ( int j=0; j < times.length; j++ ) {
		timeList.add(times[j]);
	    }
	    start = end + gapLength;
	    end = start + segLength;
	}
	timeList.trimToSize();
	EventList events = new EventList(timeList.elements());
	int nbins = 512;
	TimeSeries ts = TimeSeriesMaker.makeTimeSeries(events, nbins);
	int sampling = 21;
	double nuMin = 1e-4;
	double nuMax = 0.0256;
	LombScarglePeriodogram ls = PeriodogramMaker.makeLombScarglePeriodogram(ts, nuMin, nuMax, sampling);
	// scale LS so that white noise is normalised equally as for R2
	ls = (LombScarglePeriodogram) ls.scale(2.);
	ModifiedRayleighPeriodogram r2 = PeriodogramMaker.makeModifiedRayleighPeriodogram(events, nuMin, nuMax, sampling);
	ls.writeAsQDP(r2.getPowers(), "r2AndLS_"+nGaps+"gapsOf"+((int)gapLength)+"s.qdp");
	// scale each periodogram so that integrated power equals 1
	ls = (LombScarglePeriodogram) ls.scale(1/PeriodogramUtils.getIntegratedPower(ls));
	r2 = (ModifiedRayleighPeriodogram) r2.scale(1/PeriodogramUtils.getIntegratedPower(r2));
	r2.writeAsQDP(ls.getPowers(), "r2AndLS_scaled_"+nGaps+"gapsOf"+((int)gapLength)+"s.qdp");
	ts = TimeSeriesResampler.resample(ts, nbins/2);
	ts.writeCountsAsQDP("ts_"+(nbins/2)+"bins_"+nGaps+"gapsOf"+((int)gapLength)+"s.qdp");
    }
}