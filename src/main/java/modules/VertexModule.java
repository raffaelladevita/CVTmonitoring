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
    private final double THETAMIN = 20.0;
    private final double THETAMAX = 140.0;
    private final double VXYMIN = -0.5;//-10;
    private final double VXYMAX = 0.5;//10;
    private final double VZMIN = -8; //-26;
    private final double VZMAX =  8;//26;
    
    private final double CHI2PIDCUT = 10;
    
    
    
    public VertexModule() {
        super("Vertex");
    }
    
    public DataGroup createGroup(int col) {
        H1F hi_d0       = histo1D("hi_d0", "d0 (cm)", "Counts", 100, VXYMIN, VXYMAX, col);
        H2F hi_d0phi    = histo2D("hi_d0phi", "#phi (deg)", "d0 (cm)", 30, PHIMIN, PHIMAX, 100, VXYMIN, VXYMAX);
        GraphErrors gr  = new GraphErrors("gr_d0phi");
        gr.setTitleX("#phi (deg)");
        gr.setTitleY("d0 (cm)");
        gr.setMarkerColor(2);
        H2F hi_vxy      = histo2D("hi_vxy", "vx (cm)", "vy (cm)", 100, VXYMIN, VXYMAX, 100, VXYMIN, VXYMAX);
        H2F hi_vxphi    = histo2D("hi_vxphi", "#phi (deg)", "vx (cm)", 100, PHIMIN, PHIMAX, 100, VXYMIN, VXYMAX);
        H2F hi_vyphi    = histo2D("hi_vyphi", "#phi (deg)", "vy (cm)", 100, PHIMIN, PHIMAX, 100, VXYMIN, VXYMAX);

        DataGroup dgTrack = new DataGroup(3,2);
        dgTrack.addDataSet(hi_d0,       0);
        dgTrack.addDataSet(hi_d0phi,    1);
        dgTrack.addDataSet(gr,          2);
        dgTrack.addDataSet(hi_vxy,      3);
        dgTrack.addDataSet(hi_vxphi,    4);
        dgTrack.addDataSet(hi_vyphi,    5);
        return dgTrack;
    }

    @Override
    public void createHistos() {
        this.getHistos().put("BeamSpot", this.createGroup(46));
    }
    
    @Override
    public void fillHistos(Event event) {
        List<Track> trackC2pid = new ArrayList<>();
        for(Track track : event.getTracks()) {
            if(Math.abs(track.getChi2pid())<CHI2PIDCUT) trackC2pid.add(track);
        }
        this.fillGroup(this.getHistos().get("BeamSpot"),event.getTracks());
    }
    
    public void fillGroup(DataGroup group, List<Track> tracks) {
        for(Track track : tracks) {
            if(track.getChi2()>3 || track.getNDF()<2 || track.pt()<0.2) continue;
            group.getH1F("hi_d0").fill(track.d0());
            group.getH2F("hi_d0phi").fill(Math.toDegrees(track.phi()),track.d0());
            group.getH2F("hi_vxy").fill(track.vx(),track.vy());
            group.getH2F("hi_vxphi").fill(Math.toDegrees(track.phi()),track.vx());
            group.getH2F("hi_vyphi").fill(Math.toDegrees(track.phi()),track.vy());
        }
    }
    
    @Override
    public void setPlottingOptions(String name) {
        this.getCanvas(name).setGridX(false);
        this.getCanvas(name).setGridY(false);
        this.setLogZ(name);
    }

    @Override
    public void analyzeHistos() {
        H2F h2         = this.getHistos().get("BeamSpot").getH2F("hi_d0phi");
        GraphErrors gr = this.getHistos().get("BeamSpot").getGraph("gr_d0phi");
        this.fitSlices(h2, gr);
        F1D f1 = new F1D("f1","[a]*sin([b]*x+[c])", PHIMIN, PHIMAX);
        f1.setParameter(0, (gr.getMax()-gr.getMin())/2.0);
        f1.setParameter(1, Math.PI/180);
        f1.setParLimits(1, Math.PI/180*0.99, Math.PI/180*1.01);
        DataFitter.fit(f1, gr, "Q");
        double a = f1.getParameter(0);
        double b = f1.getParameter(1);
        double c = f1.getParameter(2);
        double aerr = 10*f1.parameter(0).error(); // convert to mm
        double cerr = 10*f1.parameter(2).error();
        System.out.println(a + " * sin(" + b + " * x + " + c + ")");
        
        //        System.out.println("x_beam: " + (-10*a*Math.cos(c)) + " mm, y_beam: " + (10*a*Math.sin(c)) + " mm");
        System.out.printf("x_beam: %2.3f mm,  y_beam: %2.3f mm\n", (-10*a*Math.cos(c)), (10*a*Math.sin(c))); // convert to mm
        System.out.printf("Rb_err: %2.3f mm, phi0_err: %2.3f deg\n", aerr, cerr);
        
    }
    
    public void fitSlices(H2F h2, GraphErrors gr) {
        gr.reset();
        ArrayList<H1F> hslice = h2.getSlicesX();
        for(int i=0; i<hslice.size(); i++) {
            double  x = h2.getXAxis().getBinCenter(i);
            double ex = 0;
            double  y = hslice.get(i).getRMS();
            double ey = 0;
            double mean  = hslice.get(i).getDataX(hslice.get(i).getMaximumBin());
            double amp   = hslice.get(i).getBinContent(hslice.get(i).getMaximumBin());
            double sigma = hslice.get(i).getRMS()/5;
            F1D f1 = new F1D("f1_dpdhi_phi","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
            f1.setParameter(0, amp);
            f1.setParameter(1, mean);
            f1.setParameter(2, sigma);
            DataFitter.fit(f1, hslice.get(i), "Q"); //No options uses error for sigma 
            if(amp>10) gr.addPoint(x, f1.getParameter(1), ex, f1.parameter(1).error());
        }
    }
}