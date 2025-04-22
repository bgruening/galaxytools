/**
 * This script runs in Fiji
 *
 * It needs the following Fiji update sites:
 * - IJPB-Plugins (MorpholibJ)
 */


import fiji.threshold.Auto_Threshold
import ij.IJ
import ij.ImagePlus
import ij.gui.Roi
import ij.measure.Measurements
import ij.measure.ResultsTable
import ij.plugin.FolderOpener
import ij.plugin.ImageCalculator
import ij.plugin.filter.Analyzer
import ij.plugin.frame.RoiManager
import inra.ijpb.binary.BinaryImages
import inra.ijpb.morphology.Morphology
import inra.ijpb.morphology.Reconstruction
import inra.ijpb.morphology.strel.SquareStrel
import inra.ijpb.segment.Threshold
import net.imagej.ImageJ

// INPUT UI AND CLI
//
#@ File (label="Input directory", style="directory") inputDir
#@ String (label="Dataset ID") datasetId
#@ Double (label="CoV threshold (-1: auto)", default="-1") threshold
#@ Boolean (label="Run headless", default="false") headless
#@ Boolean (label="Save results", default="true") saveResults
#@ String (label="Output directory name (will be created next to input directory)", default="analysis") outDirName

// INPUT FOR TESTING WITHIN IDE
//
//def inputDir = new File("/Users/tischer/Documents/wound-healing-htm-screen/data/input")
//def datasetId = "A3ROI2_Slow"; // C4ROI1_Fast A3ROI2_Slow
//def outDirName = "analysis"
//def threshold = (double) -1.0 // auto
//def headless = false;
//def saveResults = false;
//new ij.ImageJ().setVisible(true)

// FIXED PARAMETERS
//
def cellDiameter = 20
def scratchDiameter = 500
def binningFactor = 2

// DERIVED PARAMETERS
//
def cellFilterRadius = cellDiameter/binningFactor
def scratchFilterRadius = scratchDiameter/binningFactor

println("Cell filter radius: " + cellFilterRadius)
println("Scratch filter radius: " + scratchFilterRadius)

// CODE
//

// open the images
//
IJ.run("Close All");
println("Opening: " + datasetId)
def imp = FolderOpener.open(inputDir.toString(), " filter=(.*"+datasetId+".*)")
println("Number of slices: " + imp.getNSlices())

if ( imp == null  || imp.getNSlices() == 0 ) {
    println("Could not find any files matching the pattern!")
    System.exit(1)
}

// process the images to enhance regions with cells
// using the fact the the cell regions have a higher local variance
//
println("Process images to enhance regions containing cells...")
// remove scaling to work in pixel units
IJ.run(imp,"Properties...", "pixel_width=1 pixel_height=1 voxel_depth=1");
// bin to save compute time
IJ.run(imp, "Bin...", "x=" + binningFactor + " y=" + binningFactor + " z=1 bin=Average");
def binnedImp = imp.duplicate() // keep for saving
// enhance cells
IJ.run(imp, "32-bit", "");
def sdevImp = imp.duplicate()
sdevImp.setTitle(datasetId + " sdev" )
IJ.run(sdevImp, "Find Edges", "stack"); // removes larger structures, such as dirt in the background
IJ.run(sdevImp, "Variance...", "radius=" + cellFilterRadius + " stack");
IJ.run(sdevImp, "Square Root", "stack");
// mean
def meanImp = imp.duplicate()
meanImp.setTitle(datasetId + " mean")
IJ.run(meanImp, "Mean...", "radius=" + cellFilterRadius + " stack");
// cov
def covImp = ImageCalculator.run(sdevImp, meanImp, "Divide create 32-bit stack");
IJ.run(covImp, "Enhance Contrast", "saturated=0.35");
IJ.run(covImp, "8-bit", ""); // otherwise the thresholding does not seem to work
covImp.setTitle(datasetId + " cov" )
if (!headless) covImp.duplicate().show();

// create binary image (cell-free regions are foreground)
//
println("Creating binary image of cell-free regions...")
// configure black background
IJ.run("Options...", "iterations=1 count=1 black");
// determine threshold in first frame, because there we are
// closest to a 50/50 occupancy of the image with signal,
// which is best for most auto-thresholding algorithms
covImp.setPosition(1)
def histogram = covImp.getProcessor().getHistogram()
// Auto threshold with Huang method and
// multiply the threshold with a fixed factor (as is done in CellProfiler),
// based on the observation that the threshold is consistently
// a bit too high, which may be due to
// the fact that the majority of the image is foreground
if ( threshold == -1 )
  threshold = Auto_Threshold.Huang(histogram) * 0.8
println("Threshold: " + threshold)
// create binary image of whole movie,
// using the threshold of the first image
// defining the cell free regions as foreground
def binaryImp = (ImagePlus) Threshold.threshold(covImp, 0.0, threshold)
// dilate the cell free regions, because due to the cell filter radius
// the cell sizes are over estimated (blurred into cell free regions)
binaryImp = Morphology.dilation(binaryImp, SquareStrel.fromRadius((int) cellFilterRadius))
binaryImp.setTitle(datasetId + " binary")
if(!headless) binaryImp.duplicate().show()

// create scratch ROIs
//
println("Creating scratch ROI...")
binaryImp.setPosition(1) // scratch is most visible in first frame
def scratchIp = binaryImp.crop("whole-slice").getProcessor().duplicate();
// identify largest cell free region as scratch region
scratchIp = BinaryImages.keepLargestRegion(scratchIp)
// remove cells inside scratch region
scratchIp = Reconstruction.fillHoles(scratchIp)
if(!headless) new ImagePlus("Scratch", scratchIp.duplicate()).show()
// disconnect from cell free regions outside scratch
scratchIp = Morphology.opening(scratchIp, SquareStrel.fromRadius((int)(scratchFilterRadius/20)))
// in case the morphological opening cut off some cell free
// areas outside the scratch we again only keep the largest region
scratchIp = BinaryImages.keepLargestRegion(scratchIp)
// smoothen scratch edges
scratchIp = Morphology.closing(scratchIp, SquareStrel.fromRadius((int)(scratchFilterRadius/2)))

// convert binary image to ROI, which is handy for measurements
def scratchImp = new ImagePlus("Finale scratch", scratchIp)
if(!headless) scratchImp.show()
IJ.run(scratchImp, "Create Selection", "");
def scratchROI = scratchImp.getRoi()

// measure occupancy of scratch ROI
// `area_fraction` returns the fraction of foreground pixels
// (cell free area) within the measurement ROI
println("Performing measurements...")
IJ.run("Set Measurements...", "area bounding area_fraction redirect=None decimal=2");
def rt = RoiManager.multiMeasure(binaryImp, new Roi[]{scratchROI}, false)

// show results
//
if ( !headless ) {
    rt.show("Results")
    binnedImp.show()
    binnedImp.setTitle(datasetId + " binned")
    binnedImp.setRoi(scratchROI, true)
    binaryImp.setPosition(1)
    binaryImp.show()
    binaryImp.setTitle(datasetId + " binary")
    binaryImp.setRoi(scratchROI, true)
}

// save results
//
if ( saveResults ) {
    // create output directory next to input directory
    def outputDir = new File(inputDir.getParent(), outDirName);
    println("Ensuring existence of output directory: " + outputDir)
    outputDir.mkdir()
    // save table
    rt.save(new File(outputDir, datasetId + ".csv").toString());
    // save binned image with ROI
    binnedImp.setRoi(scratchROI, false)
    IJ.save(binnedImp, new File(outputDir, datasetId + ".tiff").toString());
}

println("Analysis of "+datasetId+" is done!")
if ( headless ) System.exit(0)

// FUNCTIONS
//

// copied from ImageJ RoiManager because
// https://forum.image.sc/t/make-multimeasure-public-in-roimanager/69273
private static ResultsTable multiMeasure(ImagePlus imp, Roi[] rois) {
    int nSlices = imp.getStackSize();
    Analyzer aSys = new Analyzer(imp); // System Analyzer
    ResultsTable rtSys = Analyzer.getResultsTable();
    ResultsTable rtMulti = new ResultsTable();
    rtMulti.showRowNumbers(true);
    rtSys.reset();
    int currentSlice = imp.getCurrentSlice();
    for (int slice=1; slice<=nSlices; slice++) {
        int sliceUse = slice;
        if (nSlices==1) sliceUse = currentSlice;
        imp.setSliceWithoutUpdate(sliceUse);
        rtMulti.incrementCounter();
        if ((Analyzer.getMeasurements()& Measurements.LABELS)!=0)
            rtMulti.addLabel("Label", imp.getTitle());
        int roiIndex = 0;
        for (int i=0; i<rois.length; i++) {
            imp.setRoi(rois[i]);
            roiIndex++;
            aSys.measure();
            for (int j=0; j<=rtSys.getLastColumn(); j++){
                float[] col = rtSys.getColumn(j);
                String head = rtSys.getColumnHeading(j);
                String suffix = ""+roiIndex;
                Roi roi = imp.getRoi();
                if (roi!=null) {
                    String name = roi.getName();
                    if (name!=null && name.length()>0 && (name.length()<9||!Character.isDigit(name.charAt(0))))
                        suffix = "("+name+")";
                }
                if (head!=null && col!=null && !head.equals("Slice"))
                    rtMulti.addValue(head+suffix, rtSys.getValue(j,rtSys.getCounter()-1));
            }
        }
    }
    return rtMulti;
}

