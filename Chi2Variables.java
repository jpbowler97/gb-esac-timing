package gb.esac.timing;


import cern.jet.random.Normal;
import cern.jet.random.Uniform;
import cern.jet.random.engine.MersenneTwister64;
import gb.esac.aida.functions.ChiSquareFunction;
import gb.esac.io.AsciiDataFileWriter;
import hep.aida.*;
import hep.aida.IAnalysisFactory;
import hep.aida.IFitFactory;
import hep.aida.IFitResult;
import hep.aida.IFitter;
import hep.aida.IFunction;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.ITree;
import hep.aida.ref.histogram.FixedAxis;
import hep.aida.ref.histogram.Histogram1D;
import java.io.*;
import java.io.IOException;

public class Chi2Variables {

    public static void main (String[] args) {

	
	/**  Set up JAIDA factories  **/
	IAnalysisFactory af = IAnalysisFactory.create();
	ITree tree = af.createTreeFactory().create();
	IHistogramFactory hf = af.createHistogramFactory(tree);
 	IFitFactory fitF   = af.createFitFactory();
	IFitter fitter = fitF.createFitter("Chi2", "jminuit");
	
	MersenneTwister64 engine = new MersenneTwister64();
	Normal normal = new Normal(-Math.PI, Math.PI, engine);
	Uniform uniform = new Uniform(-Math.PI, Math.PI, engine);


	int gbins = 50;
	double gmin = -5;
	double gmax = 5;
	int pbins = 50;
	double pmin = 0;
	double pmax_2 = 15;
	double pmax_10 = 30;
	IHistogram1D histoGauss = hf.createHistogram1D("Histo Gaussian", gbins, gmin, gmax);
	IHistogram1D histoChi2_1 = hf.createHistogram1D("Histo Chi2", pbins, pmin, 3);
	IHistogram1D histoChi2_2 = hf.createHistogram1D("Histo Chi2", pbins, pmin, pmax_2);
	IHistogram1D histoChi2_10 = hf.createHistogram1D("Histo Chi2", pbins, pmin, pmax_10);
	IHistogram1D histoChi2_scaled = hf.createHistogram1D("Histo Chi2", pbins, pmin, pmax_10);
	int nevents = 100000;
	for ( int i=0; i < nevents; i++ ) {
	    double phase = uniform.nextDouble();
	    double cos = Math.cos(phase);
	    double sin = Math.sin(phase);
	    double cos2 = cos*cos;
	    double sin2 = sin*sin;

	    histoGauss.fill(sin);
	    histoChi2_1.fill(sin2);
	    histoChi2_2.fill(cos2+sin2);
	    histoChi2_10.fill(5*(cos2+sin2));

// 	    histoGauss.fill(normal.nextDouble());
// 	    histoChi2_1.fill(Math.pow(normal.nextDouble(), 2));
// 	    histoChi2_2.fill(Math.pow(normal.nextDouble(), 2) + Math.pow(normal.nextDouble(), 2));
// 	    histoChi2_10.fill(Math.pow(normal.nextDouble(), 2) + Math.pow(normal.nextDouble(), 2) + Math.pow(normal.nextDouble(), 2) + Math.pow(normal.nextDouble(), 2) + Math.pow(normal.nextDouble(), 2) + Math.pow(normal.nextDouble(), 2) + Math.pow(normal.nextDouble(), 2) + Math.pow(normal.nextDouble(), 2) + Math.pow(normal.nextDouble(), 2) + Math.pow(normal.nextDouble(), 2));
// 	    histoChi2_scaled.fill(5*(Math.pow(normal.nextDouble(), 2) + Math.pow(normal.nextDouble(), 2)));
	}


	/**  Scale the histograms to make them PDFs  **/
	double binWidth_2 = (pmax_2 - pmin)/pbins;
	double binWidth_10 = (pmax_10 - pmin)/pbins;
	histoGauss.scale(1.0/(histoGauss.sumAllBinHeights()*(gmax-gmin)/gbins));
	histoChi2_1.scale(1.0/(histoChi2_1.sumAllBinHeights()*(3-pmin)/pbins));
	histoChi2_2.scale(1.0/(histoChi2_2.sumAllBinHeights()*binWidth_2));
	histoChi2_10.scale(1.0/(histoChi2_10.sumAllBinHeights()*binWidth_10));
       	histoChi2_scaled.scale(1.0/(histoChi2_scaled.sumAllBinHeights()*binWidth_10));

	/**  Perform the fits  **/
	IFitResult gaussFitResult = fitter.fit(histoGauss, "g");
	IFunction chi2 = new ChiSquareFunction("Chi Square Function");

	double dof = 1;
	double norm = 0.5;
	chi2.setParameter("dof", dof);
	chi2.setParameter("norm", norm);
	IFitResult chi2_1FitResult = fitter.fit(histoChi2_1, chi2);
	dof = histoChi2_2.mean();
	norm = histoChi2_2.sumAllBinHeights()*binWidth_2;
	chi2.setParameter("dof", dof);
	chi2.setParameter("norm", norm);
	IFitResult chi2_2FitResult = fitter.fit(histoChi2_2, chi2);
	dof = histoChi2_10.mean();
	norm = histoChi2_10.sumAllBinHeights()*binWidth_10;
	chi2.setParameter("dof", dof);
	chi2.setParameter("norm", norm);
	IFitResult chi2_10FitResult = fitter.fit(histoChi2_10, chi2);


	/**  Write the results as QDP files  **/
	double[] g = new double[gbins];
	double[] c1 = new double[gbins];
	double[] c2 = new double[gbins];
	double[] c10 = new double[gbins];
	double[] cScaled = new double[gbins];
	double[] gfunc = new double[gbins];
	double[] c1func = new double[gbins];
	double[] c2func = new double[gbins];
	double[] c10func = new double[gbins];
	double[] cScaledfunc = new double[gbins];
	double[] x1 = new double[gbins];
	double[] x2 = new double[gbins];
	double[] x3 = new double[gbins];
	double[] x4 = new double[gbins];
	FixedAxis axis1 = (FixedAxis) histoGauss.axis();
	FixedAxis axis2 = (FixedAxis) histoChi2_1.axis();
	FixedAxis axis3 = (FixedAxis) histoChi2_2.axis();
	FixedAxis axis4 = (FixedAxis) histoChi2_10.axis();
	for ( int i=0; i < gbins; i++ ) {
	    g[i] = histoGauss.binHeight(i);
	    x1[i] = axis1.binCenter(i);
	    gfunc[i] = gaussFitResult.fittedFunction().value(new double[] {x1[i]});

	    c1[i] = histoChi2_1.binHeight(i);
	    x2[i] = axis2.binCenter(i);
	    c1func[i] = chi2_1FitResult.fittedFunction().value(new double[] {x2[i]});

	    c2[i] = histoChi2_2.binHeight(i);
	    x3[i] = axis3.binCenter(i);
	    c2func[i] = chi2_2FitResult.fittedFunction().value(new double[] {x3[i]});

	    c10[i] = histoChi2_10.binHeight(i);
	    x4[i] = axis4.binCenter(i);
	    c10func[i] = chi2_10FitResult.fittedFunction().value(new double[] {x4[i]});

	    cScaled[i] = histoChi2_scaled.binHeight(i);
	    cScaledfunc[i] = 0.1*Math.exp(-0.1*x4[i]);
	}
	String[] header = new String[] {
	    "! QDP File",
	    "DEV /XS",
	    "TIME OFF", "LAB T", "LAB F", "LW 3", "CS 1.3",
	    "LAB X Value", "LAB Y PDF", 
	    "LINE STEP ON 2", "LINE ON 3",
	    "VIEW 0.1 0.1 0.9 0.9",
	    "!"
	};
	try { 
	    AsciiDataFileWriter file = new AsciiDataFileWriter("normal.qdp");
	    file.writeData(header, x1, g, gfunc);
	    file = new AsciiDataFileWriter("chi2_1Histo.qdp");
	    file.writeData(header, x2, c1, c1func);
	    file = new AsciiDataFileWriter("chi2_2Histo.qdp");
	    file.writeData(header, x3, c2, c2func);
	    file = new AsciiDataFileWriter("chi2_10Histo.qdp");
	    file.writeData(header, x4, c10, c10func);
	    file = new AsciiDataFileWriter("chi2_2scaledBy5Histo.qdp");
	    file.writeData(header, x4, cScaled, cScaledfunc);	    
	}
	catch ( IOException e ) {System.out.println("Warn : Could not write QDP file");}


	/**  Display the results 
	IPlotter plotter = af.createPlotterFactory().create("Plot Power Spectrum");
	plotter.createRegions(2,2);
	
	plotter.region(0).plot(histoGauss);
 	plotter.region(0).plot(gaussFitResult.fittedFunction());
	plotter.region(0).style().statisticsBoxStyle().setVisible(true);
	
	plotter.region(1).plot(histoChi2_2);
 	plotter.region(1).plot(chi2_2FitResult.fittedFunction());
	plotter.region(1).style().statisticsBoxStyle().setVisible(true);

	plotter.region(2).plot(histoChi2_1);
 	plotter.region(2).plot(chi2_1FitResult.fittedFunction());
	plotter.region(2).style().statisticsBoxStyle().setVisible(true);

	plotter.region(3).plot(histoChi2_10);
 	plotter.region(3).plot(chi2_10FitResult.fittedFunction());
	plotter.region(3).style().statisticsBoxStyle().setVisible(true);

	plotter.show();  **/



    }
}
