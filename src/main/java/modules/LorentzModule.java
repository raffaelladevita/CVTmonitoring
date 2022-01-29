package modules;

import analysis.Constants;
import java.util.ArrayList;
import java.util.List;
import objects.Cluster;
import objects.Event;
import analysis.Module;
import objects.Track;
import objects.Trajectory;
import org.jlab.detector.base.DetectorType;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;

/**
 *
 * @author devita
 */
public class LorentzModule extends Module {
    
    public LorentzModule() {
        super("Lorentz", false);
    }
    
    public DataGroup sizeGroup(int col) {
        DataGroup dg = new DataGroup(3,3);
        for(int ir=0; ir<Constants.BMTREGIONS; ir++) {
            for(int is=0; is<Constants.BMTSECTORS; is++) {
                int sector  = is+1;
                int layer   = Constants.BMTZLAYERS[ir];
                H1F hi_clus  = histo1D("hi_clust_" +layer + sector, "L" + layer + "S" + sector + " #phi (deg)", "Size", 200, -60, 60, col);
                H1F hi_size  = histo1D("hi_value_" +layer + sector, "L" + layer + "S" + sector + " #phi (deg)", "Size", 200, -60, 60, col);
                dg.addDataSet(hi_size, is + ir*Constants.BMTSECTORS);
                dg.addDataSet(hi_clus, is + ir*Constants.BMTSECTORS);
            }
        }
        return dg;
    }

    public DataGroup size2DGroup() {
        DataGroup dg = new DataGroup(3,3);
        for(int ir=0; ir<Constants.BMTREGIONS; ir++) {
            for(int is=0; is<Constants.BMTSECTORS; is++) {
                int sector  = is+1;
                int layer   = Constants.BMTZLAYERS[ir];
                H2F hi_clus  = histo2D("hi_clust_" +layer + sector, "L" + layer + "S" + sector + " #phi (deg)", "z (cm)", 200, -60, 60, 50, -25, 25);
                H2F hi_size  = histo2D("hi_value_" +layer + sector, "L" + layer + "S" + sector + " #phi (deg)", "z (cm)", 200, -60, 60, 50, -25, 25);
                dg.addDataSet(hi_clus, is + ir*Constants.BMTSECTORS);
                dg.addDataSet(hi_size, is + ir*Constants.BMTSECTORS);
            }
        }
        return dg;
    }

    public DataGroup langleGroup(int col) {
        DataGroup dg = new DataGroup(2,2);
        H1F hi_langle     = histo1D("hi_langle", "Local Angle (deg)", "Counts",           200, -60, 60, col);
        H2F hi_langlesize = histo2D("hi_langlesize", "Local Angle (deg)", "Cluster Size", 200, -60, 60, 15, 1, 16);
        H2F hi_langlenev  = histo2D("hi_langlenev",  "#phi (deg)", "Local Angle (deg)",   100, -50, 50, 100, -50, 50);
        H2F hi_langlephi  = histo2D("hi_langlephi",  "#phi (deg)", "Local Angle (deg)",   100, -50, 50, 100, -50, 50);
        dg.addDataSet(hi_langle,     0);
        dg.addDataSet(hi_langlenev,  1);
        dg.addDataSet(hi_langlephi,  1);
        dg.addDataSet(hi_langlesize, 2);
        return dg;
    }

    public DataGroup langlePhiGroup(int col) {
        DataGroup dg = new DataGroup(3,3);
        for(int ir=0; ir<Constants.BMTREGIONS; ir++) {
            for(int is=0; is<Constants.BMTSECTORS; is++) {
                int sector  = is+1;
                int layer   = Constants.BMTZLAYERS[ir];
                H1F hi_clus  = histo1D("hi_clust_" +layer + sector, "L" + layer + "S" + sector + " #phi (deg)", "Local Angle (deg)", 200, -60, 60, col);
                H1F hi_angl  = histo1D("hi_value_" +layer + sector, "L" + layer + "S" + sector + " #phi (deg)", "Local Angle (deg)", 200, -60, 60, col);
                dg.addDataSet(hi_angl, is + ir*Constants.BMTSECTORS);
                dg.addDataSet(hi_clus, is + ir*Constants.BMTSECTORS);
            }
        }
        return dg;
    }

    public DataGroup langle2DGroup() {
        DataGroup dg = new DataGroup(3,3);
        for(int ir=0; ir<Constants.BMTREGIONS; ir++) {
            for(int is=0; is<Constants.BMTSECTORS; is++) {
                int sector  = is+1;
                int layer   = Constants.BMTZLAYERS[ir];
                H2F hi_clus  = histo2D("hi_clust_" +layer + sector, "L" + layer + "S" + sector + " #phi (deg)", "Local Angle (deg)", 100, -60, 60, 50, -50, 50);
                H2F hi_angl  = histo2D("hi_value_" +layer + sector, "L" + layer + "S" + sector + " #phi (deg)", "Local Angle (deg)", 100, -60, 60, 50, -50, 50);
                dg.addDataSet(hi_clus, is + ir*Constants.BMTSECTORS);
                dg.addDataSet(hi_angl, is + ir*Constants.BMTSECTORS);
            }
        }
        return dg;
    }

    @Override
    public void createHistos() {
       this.getHistos().put("SizeNeg",   this.sizeGroup(44));
       this.getHistos().put("SizePos",   this.sizeGroup(42));
       this.getHistos().put("Size",      this.sizeGroup(43));
       this.getHistos().put("Size2D",    this.size2DGroup());
       this.getHistos().put("Langle",    this.langleGroup(3));
       this.getHistos().put("LanglePhi", this.langlePhiGroup(3));
       this.getHistos().put("Langle2D",  this.langle2DGroup());       
    }
    
    @Override
    public void fillHistos(Event event) {
        for(Cluster cluster : event.getClusters()) {
            if(cluster.getTrackId()>0) {
                String detector = cluster.getName();                               
                Track track = event.getTracks().get(event.getTrackMap().get(cluster.getTrackId()));
                if(detector.equals("BMTZ")) {
                    Trajectory traj = null;
                    for(Trajectory t : event.getTrajectories(cluster.getTrackId())) {
                        if(t.getDetector()==DetectorType.CVT && 
                           t.getLayer()==cluster.getLayer()+Constants.SVTLAYERS && 
                           t.getSector()==cluster.getSector()) {
                            traj = t;
                        }
                    }
                    if(traj==null) {
//                        System.out.println("Event " + event.getEvent());
//                        System.out.println(" cluster" + cluster.getName() + " " + cluster.getLayer() + " " + cluster.getSector());
//                        for(Trajectory t : event.getTrajectories(cluster.getTrackId())) 
//                            System.out.println(" traj " + t.getDetector().getName() + " " + t.getLayer() + " " + t.getSector());
                        return;
                    }
                    this.fillSizeGroup(this.getHistos().get("Size"), cluster);
                    this.fillSize2DGroup(this.getHistos().get("Size2D"), cluster, traj);
                    this.fillLangleGroup(this.getHistos().get("Langle"), cluster, traj);
                    this.fillLangle2DGroup(this.getHistos().get("Langle2D"), cluster, traj);
                    if(track.charge()>0) {
                        this.fillSizeGroup(this.getHistos().get("SizePos"), cluster);
                    }
                    else {
                        this.fillSizeGroup(this.getHistos().get("SizeNeg"), cluster);
                        this.fillLanglePhiGroup(this.getHistos().get("LanglePhi"), cluster, traj);
                    }
                }
            }
        }
    }
    
    
    public void fillSizeGroup(DataGroup group, Cluster cluster) {
        int il = cluster.getLayer()-1;
        int ir = il/2;
        double phi = Math.toDegrees((cluster.getCentroid()-Constants.BMTZSTRIPS[ir]/2)*Constants.BMTZPITCH[ir]*1E-4/Constants.BMTRADIUS[il]);
        group.getH1F("hi_clust_"+cluster.getLayer()+cluster.getSector()).fill(phi); 
        group.getH1F("hi_value_"+cluster.getLayer()+cluster.getSector()).fill(phi,cluster.getSize()); 
    }
        
    public void fillSize2DGroup(DataGroup group, Cluster cluster, Trajectory traj) {
        int il = cluster.getLayer()-1;
        int ir = il/2;
        double phi = Math.toDegrees((cluster.getCentroid()-Constants.BMTZSTRIPS[ir]/2)*Constants.BMTZPITCH[ir]*1E-4/Constants.BMTRADIUS[il]);
        group.getH2F("hi_clust_"+cluster.getLayer()+cluster.getSector()).fill(phi, traj.z()); 
        group.getH2F("hi_value_"+cluster.getLayer()+cluster.getSector()).fill(phi, traj.z(), cluster.getSize()); 
    }
        
    public void fillLangleGroup(DataGroup group, Cluster cluster, Trajectory traj) {
        int il = cluster.getLayer()-1;
        int ir = il/2;
        double phi = Math.toDegrees((cluster.getCentroid()-Constants.BMTZSTRIPS[ir]/2)*Constants.BMTZPITCH[ir]*1E-4/Constants.BMTRADIUS[il]);
        group.getH1F("hi_langle").fill(Math.toDegrees(traj.phi())); 
        group.getH2F("hi_langlesize").fill(Math.toDegrees(traj.phi()), cluster.getSize());
        group.getH2F("hi_langlenev").fill(phi, Math.toDegrees(traj.phi()));
        group.getH2F("hi_langlephi").fill(phi, Math.toDegrees(traj.phi()), cluster.getSize());
    }
            
    public void fillLanglePhiGroup(DataGroup group, Cluster cluster, Trajectory traj) {
        int layer  = cluster.getLayer();
        int il = layer-1;
        int ir = il/2;
        double phi = Math.toDegrees((cluster.getCentroid()-Constants.BMTZSTRIPS[ir]/2)*Constants.BMTZPITCH[ir]*1E-4/Constants.BMTRADIUS[il]);
        group.getH1F("hi_clust_"+cluster.getLayer()+cluster.getSector()).fill(phi); 
        group.getH1F("hi_value_"+cluster.getLayer()+cluster.getSector()).fill(phi, -Math.toDegrees(traj.phi()));
    }

    public void fillLangle2DGroup(DataGroup group, Cluster cluster, Trajectory traj) {
            int layer  = cluster.getLayer();
            int il = layer-1;
            int ir = il/2;
            double phi = Math.toDegrees((cluster.getCentroid()-Constants.BMTZSTRIPS[ir]/2)*Constants.BMTZPITCH[ir]*1E-4/Constants.BMTRADIUS[il]);
            group.getH2F("hi_clust_"+cluster.getLayer()+cluster.getSector()).fill(phi, Math.toDegrees(traj.phi())); 
            group.getH2F("hi_value_"+cluster.getLayer()+cluster.getSector()).fill(phi, Math.toDegrees(traj.phi()), cluster.getSize());
    }

    @Override
    public void setPlottingOptions(String name) {
        if(name.equals("Langle")) {
            this.getCanvas(name).getPad(1).getAxisZ().setLog(true);
        }
        else if(name.equals("Langle2D")) {
            for(EmbeddedPad pad : this.getCanvas(name).getCanvasPads())
                pad.getAxisZ().setLog(true);
        }
    }
    
    @Override
    public void analyzeHistos() {
        this.analyzeGroup("Size", true);        
        this.analyzeGroup("SizePos", true);        
        this.analyzeGroup("SizeNeg", true);        
        this.analyzeGroup("Size2D", false);  
        this.analyzeGroup("LanglePhi", false);
        this.analyzeGroup("Langle2D", false);
        this.analyzeLangle("Langle");
    }
     
    public void analyzeLangle(String name) {
        DataGroup dg = this.getHistos().get(name);
        H2F h1 = dg.getH2F("hi_langlenev");
        H2F h2 = dg.getH2F("hi_langlephi");  
        if(h1.getEntries()>0) {
            h2.divide(h1); 
            h1.reset();
        }
        GraphErrors gr_langle = new GraphErrors();
        gr_langle.setTitleX("Local Angle (deg)");
        gr_langle.setTitleY("Average cluster size");
        gr_langle.setMarkerColor(1);
        gr_langle.setMarkerSize(5);
        dg.addDataSet(gr_langle, 3);
        H2F hsize = dg.getH2F("hi_langlesize"); 
        for(int i=0; i<hsize.getSlicesX().size(); i++) {
            double x = hsize.getDataX(i);
            double y = hsize.getSlicesX().get(i).getMean();
            if(hsize.getSlicesX().get(i).integral()>50)
                gr_langle.addPoint(x, y, 0, 0);          
        }
        H1F[] hbins = new H1F[10];
        H2F   hphi  = dg.getH2F("hi_langlephi");
        for(int i=0; i<hphi.getSlicesX().size(); i++) {
            int j  = i/hbins.length;
            H1F hi = hphi.getSlicesX().get(i);
            if(hbins[j]==null) {
                hbins[j] = new H1F("hbins"+j, "","",hi.getDataSize(0), hi.getXaxis().min(), hi.getXaxis().max());
            }
            hbins[j].add(hphi.getSlicesX().get(i));
        }
        for(int i=0; i<hbins.length; i++) {
            switch (i) {
                case 0:
                case 9:
                    hbins[i].setLineColor(2);
                    break;
                case 4:
                case 5:
                    hbins[i].setLineColor(3);
                    break;
                default:
                    hbins[i].setLineColor(4);
                    break;
            }
            hbins[i].divide(hphi.getSlicesX().size()/hbins.length);
            dg.addDataSet(hbins[i], 3);
        }
    } 
    
    public void analyzeGroup(String name, boolean doFit){
        DataGroup dg = this.getHistos().get(name);
        for(int ir=0; ir<Constants.BMTREGIONS; ir++) {
            for(int is=0; is<Constants.BMTSECTORS; is++) {
                int sector  = is+1;
                int layer   = Constants.BMTZLAYERS[ir];
                if(dg.getData("hi_value_"+layer+sector) instanceof H1F) {
                    H1F h1 = dg.getH1F("hi_clust_"+layer+sector);
                    H1F h2 = dg.getH1F("hi_value_"+layer+sector);  
                    if(h1.getEntries()!=0) {
                        h2.divide(h1);
                        h1.reset();
                    }
                    if(doFit) this.fitSize(h2);
                }
                else {
                    H2F h1 = dg.getH2F("hi_clust_"+layer+sector);
                    H2F h2 = dg.getH2F("hi_value_"+layer+sector);  
                    if(h1.getEntries()!=0) {
                        h2.divide(h1);
                        h1.reset();
                    }
                }
            }
        }
    }
    
    private void fitSize(H1F hi) {
        double min = hi.getXaxis().min();
        double max = hi.getXaxis().max();
        double range = max - min;
        F1D fsize = new F1D("fsize", "[p0]+[p1]*cos(x*[f1])+[p2]*cos(x*[f2])+[p3]*cos(x*[f3])+[p4]*cos(x*[f4])",min+range*0.1, max-range*0.05);
        double ave   = (double) hi.getIntegral()/hi.getDataSize(0);
        fsize.setParameter(0, ave);
        fsize.setParameter(1, 0);
        fsize.setParameter(2, 3*Math.PI/180);
        fsize.setParameter(3, 0);
        fsize.setParameter(4, 6*Math.PI/180);
        fsize.setParameter(5, 0);
        fsize.setParameter(6, 9*Math.PI/180);
        fsize.setParameter(7, 0);
        fsize.setParameter(8, 12*Math.PI/180);
        fsize.setOptStat("1111111111");
        fsize.setLineColor(2);
        fsize.setLineWidth(2);
        DataFitter.fit(fsize, hi, "Q");
    }
}
