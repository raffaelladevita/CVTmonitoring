package analysis;

import java.io.File;
import java.util.ArrayList;
import objects.Event;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;

/**
 *
 * @author devita
 */
public class Module {    
    
    private final String          moduleName;
    private Map<String,DataGroup> moduleGroup  = new LinkedHashMap<>();
    private EmbeddedCanvasTabbed  moduleCanvas = null;
    private List<String> canvasNames = new ArrayList<>();
            
    private int nevents;
    private boolean cosmics;

    public Module(String name, boolean cosmics){                               
        this.moduleName = name;
        this.cosmics    = cosmics;
        this.init();
    }

    public void analyzeHistos() {
        // analyze the histograms at the end of the file processing
    }

    
    public void createHistos() {
        // create histograms
    }
    
    public void testHistos() {
        // run tests on the filled histograms
    }
    
    public void fillHistos(Event event) {
        // fill the histograms
    }

    public final String getName() {
        return moduleName;
    }

    public final int getNevents() {
        return nevents;
    }

    public boolean isCosmics() {
        return cosmics;
    }

    public EmbeddedCanvasTabbed getCanvas() {
        return moduleCanvas;
    }
    
    public EmbeddedCanvas getCanvas(String name) {
        return moduleCanvas.getCanvas(name);
    }
    
    public List<String> getCanvasNames() {
        return canvasNames;
    }
    
    public Map<String,DataGroup> getHistos() {
        return moduleGroup;
    }
    
    public H1F histo1D(String name, String xTitle, String yTitle, int nbins, double min, double max, int color) {
        H1F histo = new H1F(name, "", nbins, min, max);
        histo.setTitleX(xTitle);
        histo.setTitleY(yTitle);
        histo.setFillColor(color);
        return histo;
    }

    public H2F histo2D(String name, String xTitle, String yTitle, int xBins, double xMin, double xMax, int yBins, double yMin, double yMax) {
        H2F histo = new H2F(name, "", xBins, xMin, xMax, yBins, yMin, yMax);
        histo.setTitleX(xTitle);
        histo.setTitleY(yTitle);
        return histo;
    }

    public final void init() {
        this.nevents = 0;
        createHistos();
    }
    
    public final void processEvent(Event event) {
        // process event
        this.nevents++;
        this.fillHistos(event);
    }

    
    public final EmbeddedCanvasTabbed plotHistos() {
        this.analyzeHistos();
        this.drawHistos();
        return this.moduleCanvas;
    }

    public void drawHistos() {
        for(String key : moduleGroup.keySet()) {            
            this.addCanvas(key);
            this.moduleCanvas.getCanvas(key).draw(moduleGroup.get(key));
            this.setPlottingOptions(key);
            this.moduleCanvas.getCanvas(key).setGridX(false);
            this.moduleCanvas.getCanvas(key).setGridY(false);
        }
    }
    
    public final void addCanvas(String name) {
        if(this.moduleCanvas==null) this.moduleCanvas = new EmbeddedCanvasTabbed(name);
        else                        this.moduleCanvas.addCanvas(name);
        this.canvasNames.add(name);
    }
    
    public final void addCanvas(String... names) {
        for(String name : names) {
            this.addCanvas(name);
        }
    }
    
    public final void setHistos(Map<String,DataGroup> group) {
        this.moduleGroup = group;
    }
    
    public void setPlottingOptions(String name) {
        this.getCanvas().getCanvas(name).setGridX(false);
        this.getCanvas().getCanvas(name).setGridY(false);        
    }

    public void printHistos(String figures) {
        File theDir = new File(figures);
        // if the directory does not exist, create it
        if (!theDir.exists()) {
            boolean result = false;
            try{
                theDir.mkdir();
                result = true;
            } 
            catch(SecurityException se){
                //handle it
            }        
            if(result) {    
            System.out.println(">>>>> Created directory " + figures);
            }
        }
        for(String cname : canvasNames) {
            this.moduleCanvas.getCanvas(cname).save(figures + "/" + this.getName() + "_" + cname + ".png");
        }
    }
        
    public final void readDataGroup(TDirectory dir) {
        for(String key : moduleGroup.keySet()) {
            String folder = this.getName() + "/" + key + "/";
            System.out.println("Reading from: " + folder);
            DataGroup group = this.moduleGroup.get(key);
            int nrows = group.getRows();
            int ncols = group.getColumns();
            int nds   = nrows*ncols;
            DataGroup newGroup = new DataGroup(ncols,nrows);
            for(int i = 0; i < nds; i++){
                List<IDataSet> dsList = group.getData(i);
                for(IDataSet ds : dsList){
                    System.out.println("\t --> " + ds.getName());
                    newGroup.addDataSet(dir.getObject(folder, ds.getName()),i);
                }
            }            
            this.moduleGroup.replace(key, newGroup);
        }
    }
    
    public final void writeDataGroup(TDirectory dir) {
        String folder = "/" + this.getName();
        System.out.println(this.getName());
        dir.mkdir(folder);
        dir.cd(folder);
        for(String key : moduleGroup.keySet()) {
            String subfolder = key + "/";
            dir.mkdir(subfolder);
            dir.cd(subfolder);
            DataGroup group = this.moduleGroup.get(key);
            int nrows = group.getRows();
            int ncols = group.getColumns();
            int nds   = nrows*ncols;
            for(int i = 0; i < nds; i++){
                List<IDataSet> dsList = group.getData(i);
                for(IDataSet ds : dsList){
                    System.out.println("\t --> " + ds.getName());
                    dir.addDataSet(ds);
                }
            }
            dir.cd(folder);
        }
    }

    public void fitDataGroup(DataGroup dg) {
        int nx = dg.getColumns();
        int ny = dg.getRows();
        for(int i=0; i<nx*ny; i++) {
            List<IDataSet> ds = dg.getData(i);
            for(IDataSet d : ds) {
                if(d instanceof H1F)
                    this.fitGauss((H1F) d);
            }
        }
    }
    
    public void fitGauss(H1F histo) {
        double tmp_Mean = histo.getMean();
        int Max_Bin = histo.getMaximumBin();
        double tmp_Amp = histo.getBinContent(Max_Bin);
        double tmp_sigma = histo.getRMS();
        //System.out.println(tmp_Amp);
        F1D f1 = new F1D("f1", "[amp]*gaus(x,[mean],[sigma])", 0, 50.0);
        f1.setParameter(0, tmp_Amp);
        f1.setParameter(1, tmp_Mean);
        f1.setParameter(2, tmp_sigma / 2);
        f1.setLineColor(5);
        f1.setLineWidth(7);
        f1.setOptStat(111110);
        DataFitter.fit(f1, histo, "Q");
    }    
} 
