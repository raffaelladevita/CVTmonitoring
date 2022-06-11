package modules;

import analysis.Constants;
import objects.Cluster;
import objects.Event;
import analysis.Module;
import objects.CVTType;
import objects.Track;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;

/**
 *
 * @author devita
 */
public class ResidualModule extends Module {
    
    private static final double BMAX = 1500;
    private static final double SMAX =  500;
    private static final CVTType[] types = new CVTType[]{CVTType.SVT, CVTType.BMTC, CVTType.BMTZ};
    private int scale = 1;
    
    public ResidualModule(int residualScale) {
        super("Residuals",false);
        this.scale = residualScale;
    }
    
    public DataGroup svtLayerGroup(int layer) {
        int ny = 2;
        if(layer>4) ny = 3;
        DataGroup dg = new DataGroup(Constants.SVTSECTORS[layer-1]/ny, ny);
        for (int is = 0; is < Constants.SVTSECTORS[layer-1]; is++) {
            String name = "L" + layer + "S" + (is + 1);
            H1F hi_res = histo1D("hi_res_" + name, name + " Centroid Residual (um)", "Counts", 100, -SMAX, SMAX, 3);
            dg.addDataSet(hi_res, is);
            
        }
        return dg;
    }

    public DataGroup svtGraphGroup() {
        int ny = 2;
        DataGroup dg = new DataGroup(Constants.SVTLAYERS/ny, ny);
        for (int il = 0; il < Constants.SVTLAYERS; il++) {
            int layer = il+1;
            GraphErrors graph = new GraphErrors("gr_L" + layer);
            graph.setTitleX("Sector (Layer " + layer + ")");
            graph.setTitleY("Mean Centroid Residual (um)");
            graph.setMarkerColor(3);
            graph.setMarkerStyle(0);
            graph.setMarkerSize(6);
            dg.addDataSet(graph, (int) (il/2) + (il%2)*3);  
        }
        return dg;
    }
    public DataGroup bmtLayerGroup(CVTType detector) {
        DataGroup dg = new DataGroup(Constants.BMTREGIONS, Constants.BMTSECTORS);
        for (int ir = 0; ir < Constants.BMTREGIONS; ir++) {
            for (int is = 0; is < Constants.BMTSECTORS; is++) {
                String name = "L" + Constants.BMTCLAYERS[ir] + "S" + (is + 1);
                if(detector==CVTType.BMTZ)
                    name = "L" + Constants.BMTZLAYERS[ir] + "S" + (is + 1);
                H1F hi_res = histo1D("hi_res_" + name, name + " Centroid Residual (um)", "Counts", 100, -BMAX, BMAX, 3);
                dg.addDataSet(hi_res, ir * Constants.BMTSECTORS + is);
            }
        }
        return dg;
    }

    public DataGroup bmt2DLayerGroup(CVTType detector, String parameter, int nbin, double min, double max) {
        DataGroup dg = new DataGroup(Constants.BMTREGIONS, Constants.BMTSECTORS);
        for (int ir = 0; ir < Constants.BMTREGIONS; ir++) {
            for (int is = 0; is < Constants.BMTSECTORS; is++) {
                String name = "L" + Constants.BMTCLAYERS[ir] + "S" + (is + 1);
                if(detector==CVTType.BMTZ)
                    name = "L" + Constants.BMTZLAYERS[ir] + "S" + (is + 1);
                if(parameter.contains("phi")) {
                    min = Constants.BMTMEANPHI[is]-80;
                    max = Constants.BMTMEANPHI[is]+80;                    
                }
                else if(parameter.equals("strip")) {
                    if(detector==CVTType.BMTZ) max = Constants.BMTZSTRIPS[ir];
                    else if(detector==CVTType.BMTC) max = Constants.BMTCSTRIPS[ir];
                }
                H2F hi_res = histo2D("hi_res_" + name, parameter, name + " Residual (um)", nbin, min, max, 100, -BMAX, BMAX);
                dg.addDataSet(hi_res, ir * Constants.BMTSECTORS + is);
            }
        }
        return dg;
    }

    public DataGroup sumGroup(double max) {
        DataGroup dg = new DataGroup(2,2);
        H1F hi_res       = histo1D("hi_res", "Centroid Residual (um)", "Counts", 100, -max, max, 3);
        H2F hi_res_p     = histo2D("hi_res_p",     "p (GeV)",      "Centroid Residual (um)", 100,  0.0, 2.0, 100, -max, max);
        H2F hi_res_theta = histo2D("hi_res_theta", "#theta (deg)", "Centroid Residual (um)", 100,   20, 140, 100, -max, max);
        H2F hi_res_phi   = histo2D("hi_res_phi",   "#phi (deg)",   "Centroid Residual (um)", 100, -180, 180, 100, -max, max);
        dg.addDataSet(hi_res,       0);
        dg.addDataSet(hi_res_p,     1);
        dg.addDataSet(hi_res_theta, 2);
        dg.addDataSet(hi_res_phi,   3);
        return dg;
    }

    @Override
    public void createHistos() {
        for(int i=0; i<types.length; i++) {
            if(i==0) {
                for(int il=0; il<Constants.SVTLAYERS; il++) {
                    int layer = il+1;
                    this.getHistos().put(types[i] + "L" + layer, this.svtLayerGroup(layer));
                }
                this.getHistos().put(types[i].getName(), this.svtGraphGroup());                                
                this.getHistos().put(types[i].getName()+"sum", this.sumGroup(SMAX));                                
            }
            else {
                if(types[i]==CVTType.BMTZ) {
                    this.getHistos().put(types[i].getName()+"strip" + "Pos", this.bmt2DLayerGroup(types[i], "strip", 100, 0, 800));                
                    this.getHistos().put(types[i].getName()+"strip" + "Neg", this.bmt2DLayerGroup(types[i], "strip", 100, 0, 800));                                    
                }
                this.getHistos().put(types[i].getName()+"strip", this.bmt2DLayerGroup(types[i], "strip", 100, 0, 800));                
                this.getHistos().put(types[i].getName()+"phi",   this.bmt2DLayerGroup(types[i], "#phi (deg)",100, -180, 180));                
                this.getHistos().put(types[i].getName()+"theta", this.bmt2DLayerGroup(types[i], "#theta (deg)", 100, 20, 140));                
                this.getHistos().put(types[i].getName()+"p",     this.bmt2DLayerGroup(types[i], "p (GeV)", 100, 0.0, 2.0));                
                this.getHistos().put(types[i].getName()+"size",  this.bmt2DLayerGroup(types[i], "size", 20, 0, 20));                
                this.getHistos().put(types[i].getName(), this.bmtLayerGroup(types[i]));                
                this.getHistos().put(types[i].getName()+"sum", this.sumGroup(BMAX));                
            }
        }
    }
    
    @Override
    public void fillHistos(Event event) {
        for(Cluster cluster : event.getClusters()) {
            if(cluster.getTrackId()>0) {
                CVTType detector = cluster.getType();
                
                int layer  = cluster.getLayer();
                int sector = cluster.getSector();
                String name = "L"+layer+"S"+sector;
                if(!event.getTrackMap().containsKey(cluster.getTrackId())) {
                    System.out.println(cluster.getType() + " " + cluster.getTrackId());
                    for(Track t : event.getTracks()) System.out.println(t.getId());
                }
                Track track = event.getTracks().get(event.getTrackMap().get(cluster.getTrackId()));
                double residual = cluster.getCentroidResidual()*1E4/this.scale;
                if(detector==CVTType.SVT)
                    this.getHistos().get(detector.getName() + "L" + layer).getH1F("hi_res_"+name).fill(residual);
                else {
                    double phi = Math.toDegrees(track.phi());
                    if(Math.abs(phi-Constants.BMTMEANPHI[cluster.getSector()-1])>180) 
                        phi -= Math.signum(Math.abs(phi-Constants.BMTMEANPHI[cluster.getSector()-1]))*360;
                    this.getHistos().get(detector.getName()).getH1F("hi_res_"+name).fill(residual);
                    if(detector==CVTType.BMTZ) {
                        if(track.charge()>0)
                            this.getHistos().get(detector.getName()+"stripPos").getH2F("hi_res_"+name).fill(cluster.getCentroid(),residual);
                        else
                            this.getHistos().get(detector.getName()+"stripNeg").getH2F("hi_res_"+name).fill(cluster.getCentroid(),residual);
                    }
                    this.getHistos().get(detector.getName()+"strip").getH2F("hi_res_"+name).fill(cluster.getCentroid(),residual);
                    this.getHistos().get(detector.getName()+"p").getH2F("hi_res_"+name).fill(track.p(),residual);
                    this.getHistos().get(detector.getName()+"phi").getH2F("hi_res_"+name).fill(phi,residual);
                    this.getHistos().get(detector.getName()+"theta").getH2F("hi_res_"+name).fill(Math.toDegrees(track.theta()),residual);
                    this.getHistos().get(detector.getName()+"size").getH2F("hi_res_"+name).fill(cluster.getSize(),residual);
                }
                this.getHistos().get(detector.getName() + "sum").getH1F("hi_res").fill(residual);
                this.getHistos().get(detector.getName() + "sum").getH2F("hi_res_p").fill(track.p(), residual);
                this.getHistos().get(detector.getName() + "sum").getH2F("hi_res_theta").fill(Math.toDegrees(track.theta()), residual);
                this.getHistos().get(detector.getName() + "sum").getH2F("hi_res_phi").fill(Math.toDegrees(track.phi()), residual);
            }
        }
    }
    

    
    @Override
    public void analyzeHistos() {
        for(int i=0; i<types.length; i++) {
            if(i==0) {
                for(int il=0; il<Constants.SVTLAYERS; il++) {
                    int layer = il+1;
                    GraphErrors graph = this.getHistos().get(types[i].getName()).getGraph("gr_L" + layer);
                    for(int is=0; is<Constants.SVTSECTORS[il]; is++) {
                        int sector = is+1;
                        String name = "L"+layer+"S"+sector;
                        H1F hi = this.getHistos().get(types[i].getName() + "L" + layer).getH1F("hi_res_"+name);
                        this.fitResiduals(hi);
                        double mean  = hi.getFunction().getParameter(1);
                        double sigma = hi.getFunction().getParameter(2);
                        graph.addPoint(sector, mean, 0, sigma);
                    }
                }
            }
            else {
                for(int ir=0; ir<Constants.BMTREGIONS; ir++) {
                    int region = ir+1;
                    for(int is=0; is<Constants.BMTSECTORS; is++) {
                        int sector = is+1;
                        String name = "L" + Constants.BMTCLAYERS[ir] + "S" + (is + 1);
                        if(types[i]==CVTType.BMTZ)
                            name = "L" + Constants.BMTZLAYERS[ir] + "S" + (is + 1);
                        this.fitResiduals(this.getHistos().get(types[i].getName()).getH1F("hi_res_"+name));
                    }
                }
            }
            this.fitResiduals(this.getHistos().get(types[i].getName() + "sum").getH1F("hi_res"));
        }
    }

    @Override
    public void setPlottingOptions(String name) {
        super.setPlottingOptions(name);
        if(name.equals("BMTZsize")) {
            this.setLogZ(name);
        }
    }

    private void fitResiduals(H1F hi) {
        double mean = hi.getDataX(hi.getMaximumBin());
        double amp  = hi.getBinContent(hi.getMaximumBin());
        double rms  = hi.getRMS();
        double sigma = rms/3;
        String name = hi.getName().split("hi")[1];
        F1D f1 = new F1D("f1" + name,"[amp]*gaus(x,[mean],[sigma])",-0.3, 0.3);
        f1.setLineColor(2);
        f1.setLineWidth(2);
        f1.setOptStat("1111");
        f1.setParameter(0, amp);
        f1.setParameter(1, mean);
        f1.setParameter(2, sigma);
        double rmax = mean + 3.0 * Math.abs(sigma);
        double rmin = mean - 3.0 * Math.abs(sigma);
        f1.setRange(rmin, rmax);
        hi.setFunction(f1);
        if(amp>30) 
            DataFitter.fit(f1, hi, "Q"); //No options uses error for sigma 
        else
            f1.setParameter(0, 0);
    }
}
