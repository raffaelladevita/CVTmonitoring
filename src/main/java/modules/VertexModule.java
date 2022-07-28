package modules;

import java.util.ArrayList;
import java.util.List;
import objects.Track;
import objects.Event;
import analysis.Module;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;

/**
 *
 * @author devita
 */
public class VertexModule extends Module {
    
    private final double PMIN = 0.0;
    private final double PMAX = 2.0;
    private final double PHIMIN = -180.0;
    private final double PHIMAX = 180.0;
    private final double THETAMIN = 40.0;
    private final double THETAMAX = 100.0;
    private final double VXYMIN = -0.5;//-10;
    private final double VXYMAX = 0.5;//10;
    private final double VZMIN = -8; //-26;
    private final double VZMAX =  8;//26;
    
    private final double CHI2PIDCUT = 10;
    
    
    public VertexModule() {
        super("Vertex");
    }
    
    public DataGroup createVertexGroup(int col) {
        H1F hi_d0       = histo1D("hi_d0", "d0 (cm)", "Counts", 100, VXYMIN, VXYMAX, col);
        H2F hi_d0phi    = histo2D("hi_d0phi", "#phi (deg)", "d0 (cm)", 30, PHIMIN, PHIMAX, 100, VXYMIN, VXYMAX);
        GraphErrors gr  = new GraphErrors("gr_d0phi");
        gr.setTitleX("#phi (deg)");
        gr.setTitleY("d0 (cm)");
        gr.setMarkerColor(2);
        H1F hi_vz       = histo1D("hi_vz", "vz (cm)", "Counts", 100, VZMIN, VZMAX, col);
        H2F hi_vxy      = histo2D("hi_vxy", "vx (cm)", "vy (cm)", 100, VXYMIN, VXYMAX, 100, VXYMIN, VXYMAX);
        H1F hi_xb       = histo1D("hi_xb", "xb (cm)", "Counts", 10000, VXYMIN, VXYMAX, col);
        H1F hi_yb       = histo1D("hi_yb", "yb (cm)", "Counts", 10000, VXYMIN, VXYMAX, col);
        H2F hi_vxphi    = histo2D("hi_vxphi", "#phi (deg)", "vx (cm)", 100, PHIMIN, PHIMAX, 100, VXYMIN, VXYMAX);
        H2F hi_vyphi    = histo2D("hi_vyphi", "#phi (deg)", "vy (cm)", 100, PHIMIN, PHIMAX, 100, VXYMIN, VXYMAX);
        H2F hi_vzphi    = histo2D("hi_vzphi", "#phi (deg)", "vz (cm)", 100, PHIMIN, PHIMAX, 100, VZMIN, VZMAX);

        DataGroup dgVertex = new DataGroup(4,2);
        dgVertex.addDataSet(hi_d0,       0);
        dgVertex.addDataSet(hi_d0phi,    1);
        dgVertex.addDataSet(gr,          2);
        dgVertex.addDataSet(hi_vz,       3);
        dgVertex.addDataSet(hi_xb,       4);
        dgVertex.addDataSet(hi_yb,       4);
        dgVertex.addDataSet(hi_vxy,      4);
        dgVertex.addDataSet(hi_vxphi,    5);
        dgVertex.addDataSet(hi_vyphi,    6);
        dgVertex.addDataSet(hi_vzphi,    7);
        return dgVertex;
    }

    public DataGroup createMomentsGroup(int icol) {
        H2F hi_d0sinphi = histo2D("hi_d0sinphi", "#phi (deg)", "d0 (cm)", 30, PHIMIN, PHIMAX, 100, 2*VXYMIN, 2*VXYMAX);
        H2F hi_d0cosphi = histo2D("hi_d0cosphi", "#phi (deg)", "d0 (cm)", 30, PHIMIN, PHIMAX, 100, 2*VXYMIN, 2*VXYMAX);
        H1F hi_msinphi  = histo1D("hi_msinphi", "vx (cm)", "Counts", 100, 2*VXYMIN, 2*VXYMAX, icol);
        H1F hi_mcosphi  = histo1D("hi_mcosphi", "vy (cm)", "Counts", 100, 2*VXYMIN, 2*VXYMAX, icol);
        H2F hi_tsinphi  = histo2D("hi_tsinphi", "#theta (deg)", "Counts", 30, THETAMIN, THETAMAX, 100, 2*VXYMIN, 2*VXYMAX);
        H2F hi_tcosphi  = histo2D("hi_tcosphi", "#theta (deg)", "Counts", 30, THETAMIN, THETAMAX, 100, 2*VXYMIN, 2*VXYMAX);
        GraphErrors gr_tsinphi  = new GraphErrors("gr_tsinphi");
        gr_tsinphi.setTitleX("#theta (deg)");
        gr_tsinphi.setTitleY("x0 (cm)");
        gr_tsinphi.setMarkerColor(2);        
        GraphErrors gr_tcosphi  = new GraphErrors("gr_tcosphi");
        gr_tcosphi.setTitleX("#theta (deg)");
        gr_tcosphi.setTitleY("y0 (cm)");
        gr_tcosphi.setMarkerColor(2);        
        
        DataGroup dgMoments = new DataGroup(4,2);
        dgMoments.addDataSet(hi_d0sinphi, 0);
        dgMoments.addDataSet(hi_msinphi,  1);
        dgMoments.addDataSet(hi_tsinphi,  2);
        dgMoments.addDataSet(gr_tsinphi,  3);
        dgMoments.addDataSet(hi_d0cosphi, 4);
        dgMoments.addDataSet(hi_mcosphi,  5);
        dgMoments.addDataSet(hi_tcosphi,  6);
        dgMoments.addDataSet(gr_tcosphi,  7);
        return dgMoments;
    }

    @Override
    public void createHistos() {
        this.getHistos().put("UNegatives", this.createVertexGroup(43));
        this.getHistos().put("UPositives", this.createVertexGroup(47));
        this.getHistos().put("Negatives",  this.createVertexGroup(44));
        this.getHistos().put("Positives",  this.createVertexGroup(42));
        this.getHistos().put("NMoments",   this.createMomentsGroup(44));
        this.getHistos().put("PMoments",   this.createMomentsGroup(42));
    }
    
    @Override
    public void fillHistos(Event event) {
        List<Track> trackPos = new ArrayList<>();
        List<Track> trackNeg = new ArrayList<>();
        List<Track> utrackPos = new ArrayList<>();
        List<Track> utrackNeg = new ArrayList<>();
        for(Track track : event.getTracks()) {
            if(Math.abs(track.getChi2pid())<CHI2PIDCUT) {
                if(track.charge()>0) trackPos.add(track);
                else                 trackNeg.add(track);
            }
        }
        for(Track track : event.getUTracks()) {
            if(Math.abs(track.getChi2pid())<CHI2PIDCUT) {
                if(track.charge()>0) utrackPos.add(track);
                else                 utrackNeg.add(track);
            }
        }
        this.fillGroup(this.getHistos().get("Positives"), trackPos);
        this.fillGroup(this.getHistos().get("Negatives"), trackNeg);
        this.fillGroup(this.getHistos().get("UPositives"), utrackPos);
        this.fillGroup(this.getHistos().get("UNegatives"), utrackNeg);
        this.fillMomentsGroup(this.getHistos().get("NMoments"), utrackNeg);
        this.fillMomentsGroup(this.getHistos().get("PMoments"), utrackPos);
    }
    
    public void fillGroup(DataGroup group, List<Track> tracks) {
        for(Track track : tracks) {
            if(track.getNDF()<2 || track.getChi2()/track.getNDF()>30 || track.pt()<0.2) continue;
            group.getH1F("hi_d0").fill(track.d0());
            group.getH2F("hi_d0phi").fill(Math.toDegrees(track.phi()),track.d0());
            group.getH1F("hi_vz").fill(track.vz());
            group.getH2F("hi_vxy").fill(track.vx(),track.vy());
            group.getH1F("hi_xb").fill(track.xb());
            group.getH1F("hi_yb").fill(track.yb());
            group.getH2F("hi_vxphi").fill(Math.toDegrees(track.phi()),track.vx());
            group.getH2F("hi_vyphi").fill(Math.toDegrees(track.phi()),track.vy());
            group.getH2F("hi_vzphi").fill(Math.toDegrees(track.phi()),track.vz());
        }
    }
    
    public void fillMomentsGroup(DataGroup group, List<Track> tracks) {
        for(Track track : tracks) {
            if(track.getNDF()<2 || track.getChi2()/track.getNDF()>30 || track.pt()<0.2) continue;
            group.getH2F("hi_d0sinphi").fill(Math.toDegrees(track.phi()),-2*track.d00()*Math.sin(track.phi()));
            group.getH2F("hi_d0cosphi").fill(Math.toDegrees(track.phi()), 2*track.d00()*Math.cos(track.phi()));
            group.getH1F("hi_msinphi").fill(-2*track.d00()*Math.sin(track.phi()));
            group.getH1F("hi_mcosphi").fill( 2*track.d00()*Math.cos(track.phi()));
            group.getH2F("hi_tsinphi").fill(Math.toDegrees(track.theta()),-2*track.d00()*Math.sin(track.phi()));
            group.getH2F("hi_tcosphi").fill(Math.toDegrees(track.theta()), 2*track.d00()*Math.cos(track.phi()));
        }
    }
    
    @Override
    public void setPlottingOptions(String name) {
        this.getCanvas(name).setGridX(false);
        this.getCanvas(name).setGridY(false);
        this.setLogZ(name);
        if(name.endsWith("tives")) {
            EmbeddedPad pad = this.getCanvas(name).getCanvasPads().get(4);
            pad.getAxisX().setRange(VXYMIN, VXYMAX);
            pad.getAxisY().setRange(VXYMIN, VXYMAX);
        }
        else if(name.endsWith("ments")) {
            this.getCanvas(name).getCanvasPads().get(3).getAxisY().setRange(VXYMIN/2, VXYMAX/2);
            this.getCanvas(name).getCanvasPads().get(7).getAxisY().setRange(VXYMIN/2, VXYMAX/2);
        }
    }

    @Override
    public void analyzeHistos() {
        this.analyzeGroup("Positives");
        this.analyzeGroup("Negatives");
        this.analyzeGroup("UPositives");
        this.analyzeGroup("UNegatives");
        this.analyzeMomentsGroup("NMoments");
        this.analyzeMomentsGroup("PMoments");
    }
    
    private void analyzeGroup(String name) {
        H2F h2         = this.getHistos().get(name).getH2F("hi_d0phi");
        GraphErrors gr = this.getHistos().get(name).getGraph("gr_d0phi");
        F1D f1 = this.fitD0Phi(h2, gr);
        
        System.out.printf("\nAnalyzing Vertex group: " + name +"\n");
        System.out.printf("d0(phi) = p0 sin(p1 x + p2):\n");
        for(int i=0; i<f1.getNPars(); i++)
            System.out.printf("\t p%d = (%.4f +/- %.4f)\n", i, f1.getParameter(i), f1.parameter(i).error());
        double xb =  10*this.getHistos().get(name).getH1F("hi_xb").getMean();
        double yb =  10*this.getHistos().get(name).getH1F("hi_yb").getMean();
        double dx = -10*f1.getParameter(0)*Math.cos(f1.getParameter(2));
        double dy =  10*f1.getParameter(0)*Math.sin(f1.getParameter(2));
        double edx = 10*Math.sqrt(Math.pow(f1.parameter(0).error()*Math.cos(f1.getParameter(2)),2)+
                                  Math.pow(f1.getParameter(0)*Math.sin(f1.getParameter(2))*f1.parameter(2).error(),2));
        double edy = 10*Math.sqrt(Math.pow(f1.parameter(0).error()*Math.sin(f1.getParameter(2)),2)+
                                  Math.pow(f1.getParameter(0)*Math.cos(f1.getParameter(2))*f1.parameter(2).error(),2));
        System.out.printf("x_offset: (%2.3f +/- %2.3f) mm, y_offset: (%2.3f +/- %2.3f) mm\n", dx, edx, dy, edy); // convert to mm        
        System.out.printf("  with respect to beam spot read from banks: (%2.3f, %2.3f) mm\n", xb, yb);      
        System.out.printf("Update the beam (x,y) position to: (%2.3f, %2.3f) mm\n", xb+dx, yb+dy);       
        System.out.printf("or shift the detector position by: (%2.3f, %2.3f) mm\n", -dx, -dy); 
    }
    
    private void analyzeMomentsGroup(String name) {
        this.getHistos().get(name).getH1F("hi_msinphi").setOptStat("1100");
        this.getHistos().get(name).getH1F("hi_mcosphi").setOptStat("1100");
        this.fitMTheta(this.getHistos().get(name).getH2F("hi_tsinphi"), this.getHistos().get(name).getGraph("gr_tsinphi"));
        this.fitMTheta(this.getHistos().get(name).getH2F("hi_tcosphi"), this.getHistos().get(name).getGraph("gr_tcosphi"));
    }
    
    private F1D fitD0Phi(H2F h2, GraphErrors gr) {
        this.fitSlices(h2, gr);
        F1D f1 = new F1D("f1","[p0]*sin([p1]*x+[p2])", PHIMIN, PHIMAX);
        f1.setParameter(0, (gr.getMax()-gr.getMin())/2.0);
        f1.setParameter(1, Math.PI/180);
        f1.setParLimits(1, Math.PI/180*0.99, Math.PI/180*1.01);
        DataFitter.fit(f1, gr, "Q");
        return f1;
    }
    
    private void fitMTheta(H2F h2, GraphErrors gr) {
        ArrayList<H1F> hslice = h2.getSlicesX();
        for(int i=0; i<hslice.size(); i++) {
            double  x = h2.getXAxis().getBinCenter(i);
            double ex = 0;
            double  y = hslice.get(i).getMean();
            double ey = 0;
            if(hslice.get(i).getIntegral()>200) gr.addPoint(x, y, ex, ey);
        }
    }
    
    public void fitSlices(H2F h2, GraphErrors gr) {
        gr.reset();
        ParallelSliceFitter psf = new ParallelSliceFitter(h2);
        psf.setBackgroundOrder(0);
        this.toDevNull();
        psf.fitSlicesX();
        this.restoreStdOutErr();
        
        GraphErrors mean  = psf.getMeanSlices();
        GraphErrors sigma = psf.getSigmaSlices();
        for(int i=0; i<mean.getDataSize(0); i++) {
            if(Math.abs(sigma.getDataY(i))<h2.getSlicesX().get(i).getRMS())
                gr.addPoint(mean.getDataX(i), mean.getDataY(i), 0, sigma.getDataY(i));
        }
        if(gr.getDataSize(0)<2) {
            System.out.println(gr.getName());
            gr.addPoint(PHIMIN, VXYMIN, 0, 0);
            gr.addPoint(PHIMAX, VXYMAX, 0, 0);
        }
//        gr = mean;
//        ArrayList<H1F> hslice = h2.getSlicesX();
//        for(int i=0; i<hslice.size(); i++) {
//            double  x = h2.getXAxis().getBinCenter(i);
//            double ex = 0;
//            double  y = hslice.get(i).getRMS();
//            double ey = 0;
//            double mean  = hslice.get(i).getDataX(hslice.get(i).getMaximumBin());
//            double amp   = hslice.get(i).getBinContent(hslice.get(i).getMaximumBin());
//            double sigma = hslice.get(i).getRMS()/5;
//            F1D f1 = new F1D("f1_dpdhi_phi","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
//            f1.setParameter(0, amp);
//            f1.setParameter(1, mean);
//            f1.setParameter(2, sigma);
//            DataFitter.fit(f1, hslice.get(i), "Q"); //No options uses error for sigma 
//            if(amp>10) gr.addPoint(x, f1.getParameter(1), ex, f1.getParameter(2));
//        }
    }
}
