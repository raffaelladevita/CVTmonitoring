package modules;

import analysis.Constants;
import java.util.ArrayList;
import java.util.List;
import objects.Cluster;
import objects.Event;
import analysis.Module;
import objects.Hit;
import org.jlab.groot.data.H1F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.math.F1D;

/**
 *
 * @author devita
 */
public class ClusterModule extends Module {
    
    public ClusterModule() {
        super("Clusters", false);
    }
    
    public DataGroup clusterGroup(int col) {
        String[] names = new String[]{"SVT", "BMTC", "BMTZ"};
        double[] EMAX = {2000, 10000, 10000};
        DataGroup dg = new DataGroup(3,3);
        for(int i=0; i<names.length; i++) {
            String name = names[i];
            H1F hi_size    = histo1D("hi_size_" + name,   name + " Cluster Size",   "Counts",  20, 0, 20, col);
            H1F hi_energy  = histo1D("hi_energy_" + name, name + " Cluster Energy", "Counts", 100, 0, EMAX[i], col);
            H1F hi_time    = histo1D("hi_time_" + name,   name + " Cluster Time",   "Counts", 440, 0, 440, col);

            dg.addDataSet(hi_size,   0 + i*3);
            dg.addDataSet(hi_energy, 1 + i*3);
            dg.addDataSet(hi_time,   2 + i*3);
        }
        return dg;
    }

    public DataGroup occSVTGroup(int col) {
        DataGroup dg = new DataGroup(3,2);
        for(int il=0; il<Constants.SVTLAYERS; il++) {
            int layer = il+1;
            String name = "layer"+layer;
            int nbins = Constants.SVTSTRIPS*Constants.SVTSECTORS[il];
            H1F hi_occ  = histo1D("hi_occ_" + name, "L" + layer + "Strip",   "Counts", nbins,   0, nbins, col);           
            dg.addDataSet(hi_occ, 0 + il);
        }
        return dg;
    }

    public DataGroup occBMTGroup(String name, int col) {
        DataGroup dg = new DataGroup(3,3);
        for(int ir=0; ir<Constants.BMTREGIONS; ir++) {
            for(int is=0; is<Constants.BMTSECTORS; is++) {
                int sector  = is+1;
                int layer   = Constants.BMTCLAYERS[ir];
                int nstrips = Constants.BMTCSTRIPS[ir];
                if(name.equals("BMTZ")) {
                    layer = Constants.BMTZLAYERS[ir];
                    nstrips = Constants.BMTZSTRIPS[ir];
                }
                H1F hi_occ  = histo1D("hi_occ_" +layer + sector, "L" + layer + "S" + sector + " Strip",   "Counts", nstrips+1,   0, nstrips+1, col);
                dg.addDataSet(hi_occ, is + ir*Constants.BMTSECTORS);
            }
        }
        return dg;
    }   
    
    public DataGroup sizeBMTGroup(String name, int col) {
        DataGroup dg = new DataGroup(3,3);
        for(int ir=0; ir<Constants.BMTREGIONS; ir++) {
            for(int is=0; is<Constants.BMTSECTORS; is++) {
                int sector  = is+1;
                int layer   = Constants.BMTCLAYERS[ir];
                int nstrips = Constants.BMTCSTRIPS[ir];
                if(name.equals("BMTZ")) {
                    layer = Constants.BMTZLAYERS[ir];
                    nstrips = Constants.BMTZSTRIPS[ir];
                }
                H1F hi_clus  = histo1D("hi_clus_" +layer + sector, "L" + layer + "S" + sector + " Strip",   "Size", nstrips+1,   0, nstrips+1, col);
                H1F hi_size  = histo1D("hi_size_" +layer + sector, "L" + layer + "S" + sector + " Strip",   "Size", nstrips+1,   0, nstrips+1, col);
                dg.addDataSet(hi_size, is + ir*Constants.BMTSECTORS);
            }
        }
        return dg;
    }

    @Override
    public void createHistos() {
        this.getHistos().put("Clusters",        this.clusterGroup(44));
        this.getHistos().put("ClustersOnTrack", this.clusterGroup(3));
        this.getHistos().put("ClustersNotOnTrack", this.clusterGroup(44));
        this.getHistos().put("SVT",            this.occSVTGroup(44));
        this.getHistos().put("SVTOnTrack",     this.occSVTGroup(3));
        this.getHistos().put("BMTC",           this.occBMTGroup("BMTC", 44));
        this.getHistos().put("BMTCOnTrack",    this.occBMTGroup("BMTC", 3));
        this.getHistos().put("BMTCsize",       this.sizeBMTGroup("BMTC", 3));
        this.getHistos().put("BMTZ",           this.occBMTGroup("BMTZ", 44));
        this.getHistos().put("BMTZOnTrack",    this.occBMTGroup("BMTZ", 3));
        this.getHistos().put("BMTZsize",       this.sizeBMTGroup("BMTZ", 3));
    }
    
    @Override
    public void fillHistos(Event event) {
        List<Cluster> clustersOnTrack    = new ArrayList<>();
        List<Cluster> clustersNotOnTrack = new ArrayList<>();
        for(Cluster cluster : event.getClusters()) {
            if(cluster.getTrackId()>0)
                clustersOnTrack.add(cluster);
            else
                clustersNotOnTrack.add(cluster);
        }
        this.fillGroup(this.getHistos().get("Clusters"),event.getClusters());
        this.fillGroup(this.getHistos().get("ClustersOnTrack"),clustersOnTrack);
        this.fillGroup(this.getHistos().get("ClustersNotOnTrack"),clustersNotOnTrack);
        this.fillOccSVTGroup(this.getHistos().get("SVT"),event.getClusters());
        this.fillOccSVTGroup(this.getHistos().get("SVTOnTrack"),clustersOnTrack);
        this.fillOccBMTGroup(this.getHistos().get("BMTC"), "BMTC", event.getClusters());
        this.fillOccBMTGroup(this.getHistos().get("BMTCOnTrack"), "BMTC", clustersOnTrack);
        this.fillOccBMTGroup(this.getHistos().get("BMTZ"), "BMTZ", event.getClusters());
        this.fillOccBMTGroup(this.getHistos().get("BMTZOnTrack"), "BMTZ", clustersOnTrack);
        this.fillSizeBMTGroup(this.getHistos().get("BMTCsize"), "BMTC", clustersOnTrack);
        this.fillSizeBMTGroup(this.getHistos().get("BMTZsize"), "BMTZ", clustersOnTrack);
    }
    
    public void fillGroup(DataGroup group, List<Cluster> clusters) {
        for(Cluster cluster : clusters) {
            group.getH1F("hi_size_" + cluster.getName()).fill(cluster.getSize());
            group.getH1F("hi_energy_" + cluster.getName()).fill(cluster.getEnergy());
            group.getH1F("hi_time_" + cluster.getName()).fill(cluster.getTime());            
        }
    }
    
    public void fillOccSVTGroup(DataGroup group, List<Cluster> clusters) {
        for(Cluster cluster : clusters) {
            if(cluster.getName().equals("SVT")) 
                group.getH1F("hi_occ_layer"+cluster.getLayer()).fill(cluster.getSeedStrip()+Constants.SVTSTRIPS*(cluster.getSector()-1));  
        }
    }
    
    public void fillOccBMTGroup(DataGroup group, String name, List<Cluster> clusters) {
        for(Cluster cluster : clusters) {
            if(cluster.getName().equals(name)) {
                group.getH1F("hi_occ_"+cluster.getLayer()+cluster.getSector()).fill(cluster.getSeedStrip()); 
            }
        }
    }
    
    public void fillSizeBMTGroup(DataGroup group, String name, List<Cluster> clusters) {
        for(Cluster cluster : clusters) {
            if(cluster.getName().equals(name)) {
                 group.getH1F("hi_size_"+cluster.getLayer()+cluster.getSector()).fill(cluster.getSeedStrip(),cluster.getSize()); 
            }
        }
    }
    
    @Override
    public void drawHistos() {
        this.setCanvas(new EmbeddedCanvasTabbed("Clusters", "SVT", "BMTC", "BMTZ", "BMTCsize", "BMTZsize"));
        this.getCanvas("Clusters").draw(this.getHistos().get("Clusters"));
        this.getCanvas("Clusters").draw(this.getHistos().get("ClustersNotOnTrack"));
        this.getCanvas("Clusters").draw(this.getHistos().get("ClustersOnTrack"));
        this.getCanvas("SVT").draw(this.getHistos().get("SVT"));
        this.getCanvas("SVT").draw(this.getHistos().get("SVTOnTrack"));
        this.getCanvas("BMTC").draw(this.getHistos().get("BMTC"));
        this.getCanvas("BMTC").draw(this.getHistos().get("BMTCOnTrack"));
        this.getCanvas("BMTZ").draw(this.getHistos().get("BMTZ"));
        this.getCanvas("BMTZ").draw(this.getHistos().get("BMTZOnTrack"));
        this.getCanvas("BMTCsize").draw(this.getHistos().get("BMTCsize"));
        this.getCanvas("BMTZsize").draw(this.getHistos().get("BMTZsize"));
        this.setPlottingOptions("Clusters");
        this.setPlottingOptions("SVT");
        this.setPlottingOptions("BMTC");
        this.setPlottingOptions("BMTZ");
        super.setPlottingOptions("BMTCsize");
        super.setPlottingOptions("BMTZsize");
    }
       
    @Override
    public void setPlottingOptions(String name) {
        for(EmbeddedPad pad : this.getCanvas(name).getCanvasPads())
            pad.getAxisY().setLog(true);
    }

    @Override
    public void analyzeHistos() {
        this.analyzeSizeGroup("BMTC", false);        
        this.analyzeSizeGroup("BMTZ", true);        
    }
       
    public void analyzeSizeGroup(String name, boolean doFit){
        for(int ir=0; ir<Constants.BMTREGIONS; ir++) {
            for(int is=0; is<Constants.BMTSECTORS; is++) {
                int sector  = is+1;
                int layer   = Constants.BMTCLAYERS[ir];
                if(name.equals("BMTZ")) {
                    layer = Constants.BMTZLAYERS[ir];
                }
                H1F h1 = this.getHistos().get(name+"OnTrack").getH1F("hi_occ_"+layer+sector);
                H1F h2 = this.getHistos().get(name+"size").getH1F("hi_size_"+layer+sector);
                if(h1.integral()<h2.integral())
                    h2.divide(h1);
                if(doFit) this.fitSize(h2);
            }
        }
    }
    
    private void fitSize(H1F hi) {
        double range = hi.getXaxis().max();
        F1D fsize = new F1D("fsize", "[p0]+[p1]*cos(x*[f1])+[p2]*cos(x*[f2])+[p3]*cos(x*[f3])+[p4]*cos(x*[f4])",range*0.05, range*0.95);
        double ave   = (double) hi.getIntegral()/hi.getDataSize(0);
        fsize.setParameter(0, ave);
        fsize.setParameter(1, 0);
        fsize.setParameter(2, 3*Math.toRadians(120)/range);
        fsize.setParameter(3, 0);
        fsize.setParameter(4, 6*Math.toRadians(120)/range);
        fsize.setParameter(5, 0);
        fsize.setParameter(6, 9*Math.toRadians(120)/range);
        fsize.setParameter(7, 0);
        fsize.setParameter(8, 12*Math.toRadians(120)/range);
        fsize.setOptStat("1111111111");
        fsize.setLineColor(2);
        fsize.setLineWidth(2);
        DataFitter.fit(fsize, hi, "Q");
    }
}
