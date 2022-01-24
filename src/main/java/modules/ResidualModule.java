package modules;

import analysis.Constants;
import objects.Cluster;
import objects.Event;
import analysis.Module;
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
    private static final String[] names = new String[]{"SVT", "BMTC", "BMTZ"};

    
    public ResidualModule() {
        super("Residuals",false);
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
    public DataGroup bmtLayerGroup(String detector) {
        DataGroup dg = new DataGroup(Constants.BMTREGIONS, Constants.BMTSECTORS);
        for (int ir = 0; ir < Constants.BMTREGIONS; ir++) {
            for (int is = 0; is < Constants.BMTSECTORS; is++) {
                String name = "L" + Constants.BMTCLAYERS[ir] + "S" + (is + 1);
                if(detector.equals("BMTZ"))
                    name = "L" + Constants.BMTZLAYERS[ir] + "S" + (is + 1);
                H1F hi_res = histo1D("hi_res_" + name, name + " Centroid Residual (um)", "Counts", 100, -BMAX, BMAX, 3);
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
        for(int i=0; i<names.length; i++) {
            if(i==0) {
                for(int il=0; il<Constants.SVTLAYERS; il++) {
                    int layer = il+1;
                    this.getHistos().put(names[i] + "L" + layer, this.svtLayerGroup(layer));
                }
                this.getHistos().put(names[i], this.svtGraphGroup());                                
                this.getHistos().put(names[i]+"sum", this.sumGroup(SMAX));                                
            }
            else {
                this.getHistos().put(names[i], this.bmtLayerGroup(names[i]));                
                this.getHistos().put(names[i]+"sum", this.sumGroup(BMAX));                
            }
        }
    }
    
    @Override
    public void fillHistos(Event event) {
        for(Cluster cluster : event.getClusters()) {
            if(cluster.getTrackId()>0) {
                String detector = cluster.getName();
                
                int layer  = cluster.getLayer();
                int sector = cluster.getSector();
                String name = "L"+layer+"S"+sector;
                if(!event.getTrackMap().containsKey(cluster.getTrackId())) {
                    System.out.println(cluster.getName() + " " + cluster.getTrackId());
                    for(Track t : event.getTracks()) System.out.println(t.getId());
                }
                Track track = event.getTracks().get(event.getTrackMap().get(cluster.getTrackId()));
                if(detector.equals("SVT"))
                    this.getHistos().get(detector + "L" + layer).getH1F("hi_res_"+name).fill(cluster.getCentroidResidual()*1E4);
                else
                    this.getHistos().get(detector).getH1F("hi_res_"+name).fill(cluster.getCentroidResidual()*1E4);
                this.getHistos().get(detector + "sum").getH1F("hi_res").fill(cluster.getCentroidResidual()*1E4);
                this.getHistos().get(detector + "sum").getH2F("hi_res_p").fill(track.p(), cluster.getCentroidResidual()*1E4);
                this.getHistos().get(detector + "sum").getH2F("hi_res_theta").fill(Math.toDegrees(track.theta()), cluster.getCentroidResidual()*1E4);
                this.getHistos().get(detector + "sum").getH2F("hi_res_phi").fill(Math.toDegrees(track.phi()), cluster.getCentroidResidual()*1E4);
            }
        }
    }
    

    
    @Override
    public void analyzeHistos() {
        for(int i=0; i<names.length; i++) {
            if(i==0) {
                for(int il=0; il<Constants.SVTLAYERS; il++) {
                    int layer = il+1;
                    GraphErrors graph = this.getHistos().get(names[i]).getGraph("gr_L" + layer);
                    for(int is=0; is<Constants.SVTSECTORS[il]; is++) {
                        int sector = is+1;
                        String name = "L"+layer+"S"+sector;
                        H1F hi = this.getHistos().get(names[i] + "L" + layer).getH1F("hi_res_"+name);
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
                        if(names[i].equals("BMTZ"))
                            name = "L" + Constants.BMTZLAYERS[ir] + "S" + (is + 1);
                        this.fitResiduals(this.getHistos().get(names[i]).getH1F("hi_res_"+name));
                    }
                }
            }
            this.fitResiduals(this.getHistos().get(names[i] + "sum").getH1F("hi_res"));
        }
    }

    private void fitResiduals(H1F hi) {
        double mean = hi.getDataX(hi.getMaximumBin());
        double amp  = hi.getBinContent(hi.getMaximumBin());
        double rms  = hi.getRMS();
        double sigma = rms/2;
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
        if(amp>50) 
            DataFitter.fit(f1, hi, "Q"); //No options uses error for sigma 
        else
            f1.setParameter(0, 0);
    }
}
