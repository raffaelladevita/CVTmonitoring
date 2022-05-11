package modules;

import objects.Track;
import objects.Event;
import analysis.Module;
import org.jlab.clas.pdg.PhysicsConstants;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.graphics.IDataSetPlotter;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;

/**
 *
 * @author devita
 */
public class DoubletsModule extends Module {
    
    private final double DP = 1.0;
    private final double PCUT = 1.0;
    private final double PMIN = 0.4;
    private final double PMAX = 2.0;
    private final double PHIMIN = -180.0;
    private final double PHIMAX = 180.0;
    private final double DPHI = 20.0;
    private final double THETAMIN = 30.0;
    private final double THETAMAX = 44.0;
    private final double DTHETA   = 15.0;
    private final double THETACUT = 34.0;
    private final double VZMIN = -10; 
    private final double VZMAX =  10;
    
    private String[] charges;

    public DoubletsModule(double ebeam, boolean lund) {
        super("Doublets", false, ebeam, lund);
    }
    
    private DataGroup selectionGroup() {
        DataGroup dg = new DataGroup(4,2);
        for(int i=0; i<2; i++) {
            H2F hi_dtheta_dphi  = histo2D("hi_dtheta_dphi_"  + charges[i], "#theta_f-#theta_c (deg)", "#phi_f-#phi_c (deg)", 100, -DTHETA, DTHETA, 100, -DPHI, DPHI); 
            H2F hi_dtheta_dp    = histo2D("hi_dtheta_dp_"    + charges[i], "#theta_f-#theta_c (deg)", "p_f-p_c (GeV)",       100, -DTHETA, DTHETA, 100, -DP, DP); 
            H2F hi_theta_dtheta = histo2D("hi_theta_dtheta_" + charges[i], "#theta_f (deg)", "#theta_f-#theta_c (deg)",      100, THETAMIN, THETAMAX, 100, -DTHETA, DTHETA); 
            H2F hi_theta_dphi   = histo2D("hi_theta_dphi_"   + charges[i], "#theta_f (deg)", "#phi_f-#phi_c (deg)",          100, THETAMIN, THETAMAX, 100, -DPHI, DPHI); 
            dg.addDataSet(hi_dtheta_dphi,  0+i*4);
            dg.addDataSet(hi_dtheta_dp,    1+i*4);
            dg.addDataSet(hi_theta_dtheta, 2+i*4);
            dg.addDataSet(hi_theta_dphi,   3+i*4);
        }
        return dg;
    }

    private DataGroup resolutionGroup() {
        DataGroup dg = new DataGroup(4,2);
        for(int i=0; i<2; i++) {
            H1F hi_dphi   = histo1D("hi_dphi_"  + charges[i], "#phi_f-#phi_c (deg)",    "Counts", 100, -DPHI, DPHI, 0); 
            H1F hi_dtheta = histo1D("hi_dtheta_"+ charges[i], "#theta_f-#theta_c (deg)","Counts", 100, -DTHETA, DTHETA, 0); 
            H1F hi_dz     = histo1D("hi_dz_"    + charges[i], "z_f-z_c (cm)",           "Counts", 100,  VZMIN, VZMAX, 0); 
            H1F hi_dp     = histo1D("hi_dp_"    + charges[i], "(p_f-p_c)/p_c",          "Counts", 100, -DP, DP, 0); 
            dg.addDataSet(hi_dphi,   0+i*4);
            dg.addDataSet(hi_dtheta, 1+i*4);
            dg.addDataSet(hi_dz,     2+i*4);
            dg.addDataSet(hi_dp,     3+i*4);
        }
        return dg;
    }
    
    private DataGroup pGroup() {
        DataGroup dg = new DataGroup(4,2);
        for(int i=0; i<2; i++) {
            H2F hi_p_dphi   = histo2D("hi_p_dphi_"  + charges[i], "p_c (deg)", "#phi_f-#phi_c (deg)",     100, PMIN, PMAX, 100, -DPHI, DPHI); 
            H2F hi_p_dtheta = histo2D("hi_p_dtheta_"+ charges[i], "p_c (deg)", "#theta_f-#theta_c (deg)", 100, PMIN, PMAX, 100, -DTHETA, DTHETA); 
            H2F hi_p_dz     = histo2D("hi_p_dz_"    + charges[i], "p_c (deg)", "z_f-z_c (cm)",            100, PMIN, PMAX, 100,  VZMIN, VZMAX); 
            H2F hi_p_dp     = histo2D("hi_p_dp_"    + charges[i], "p_c (deg)", "(p_f-p_c)/p_c",           100, PMIN, PMAX, 100, -DP, DP); 
            dg.addDataSet(hi_p_dphi,   0+i*4);
            dg.addDataSet(hi_p_dtheta, 1+i*4);
            dg.addDataSet(hi_p_dz,     2+i*4);
            dg.addDataSet(hi_p_dp,     3+i*4);
        }
        return dg;
    }
    
    private DataGroup phiGroup() {
        DataGroup dg = new DataGroup(4,2);
        for(int i=0; i<2; i++) {
            H2F hi_phi_dphi   = histo2D("hi_phi_dphi_"  + charges[i], "#phi_c (deg)", "#phi_f-#phi_c (deg)",     100, PHIMIN, PHIMAX, 100, -DPHI, DPHI); 
            H2F hi_phi_dtheta = histo2D("hi_phi_dtheta_"+ charges[i], "#phi_c (deg)", "#theta_f-#theta_c (deg)", 100, PHIMIN, PHIMAX, 100, -DTHETA, DTHETA); 
            H2F hi_phi_dz     = histo2D("hi_phi_dz_"    + charges[i], "#phi_c (deg)", "z_f-z_c (cm)",            100, PHIMIN, PHIMAX, 100,  VZMIN, VZMAX); 
            H2F hi_phi_dp     = histo2D("hi_phi_dp_"    + charges[i], "#phi_c (deg)", "(p_f-p_c)/p_c",           100, PHIMIN, PHIMAX, 100, -DP, DP); 
            dg.addDataSet(hi_phi_dphi,   0+i*4);
            dg.addDataSet(hi_phi_dtheta, 1+i*4);
            dg.addDataSet(hi_phi_dz,     2+i*4);
            dg.addDataSet(hi_phi_dp,     3+i*4);
        }
        return dg;
    }
    
    private DataGroup thetaGroup() {
        DataGroup dg = new DataGroup(4,2);
        for(int i=0; i<2; i++) {
            H2F hi_theta_dphi   = histo2D("hi_theta_dphi_"  + charges[i], "#theta_c (deg)", "#phi_f-#phi_c (deg)",     100, THETAMIN, THETAMAX, 100, -DPHI, DPHI); 
            H2F hi_theta_dtheta = histo2D("hi_theta_dtheta_"+ charges[i], "#theta_c (deg)", "#theta_f-#theta_c (deg)", 100, THETAMIN, THETAMAX, 100, -DTHETA, DTHETA); 
            H2F hi_theta_dz     = histo2D("hi_theta_dz_"    + charges[i], "#theta_c (deg)", "z_f-z_c (cm)",            100, THETAMIN, THETAMAX, 100,  VZMIN, VZMAX); 
            H2F hi_theta_dp     = histo2D("hi_theta_dp_"    + charges[i], "#theta_c (deg)", "(p_f-p_c)/p_c",           100, THETAMIN, THETAMAX, 100, -DP, DP); 
            dg.addDataSet(hi_theta_dphi,   0+i*4);
            dg.addDataSet(hi_theta_dtheta, 1+i*4);
            dg.addDataSet(hi_theta_dz,     2+i*4);
            dg.addDataSet(hi_theta_dp,     3+i*4);
        }
        return dg;
    }

    
    @Override
    public void createHistos() {
        charges = new String[]{"pos", "neg"};
        this.getHistos().put("Selection",  this.selectionGroup());
        this.getHistos().put("Resolution", this.resolutionGroup());
        this.getHistos().put("P",          this.pGroup());
        this.getHistos().put("Phi",        this.phiGroup());
        this.getHistos().put("Theta",      this.thetaGroup());
    }
    
    private Track getElectron(Event event) {
        if(!event.getParticles().isEmpty()) {
            Track t = event.getParticles().get(0);
            if(t.pid()==11 && t.getDetector()==2 && t.p()>2.5)
                return t;
        }
        return null;
    }
    
    private Track getFDTrack(Event event) {
        if(event.getParticles().size()>1) {
            for(int i=1; i< event.getParticles().size(); i++) {
                Track t = event.getParticles().get(i);
                if(t.getDetector()==2 && 
                   t.charge()!=0 &&
                   t.p()>PMIN &&
                   Math.toDegrees(t.theta())>THETAMIN && Math.toDegrees(t.theta())<THETAMAX &&
                   Math.abs(t.getChi2pid())<3 &&
                   t.getSector()!=event.getParticles().get(0).getSector()     )
                    return t;
            }
        }
        return null;
    }
    
    private Track getCDTrack(Event event) {
        if(event.getParticles().size()>1) {
            for(int i=1; i< event.getParticles().size(); i++) {
                Track t = event.getParticles().get(i);
                if(t.getDetector()==4 && 
                   t.charge()!=0 &&
                   t.p()>PMIN &&
                   Math.toDegrees(t.theta())>THETAMIN && Math.toDegrees(t.theta())<THETAMAX)
                    return t;
            }
        }
        return null;
    }
    
    @Override
    public void fillHistos(Event event) {
        
        Track electron = this.getElectron(event);
        Track forward  = this.getFDTrack(event);
        Track central  = this.getCDTrack(event);
        
        if(electron==null || forward==null || central==null) return;
        if(forward.charge()!=central.charge()) return;
        
        String charge = charges[(-forward.charge()+1)/2];
        this.getHistos().get("Selection").getH2F("hi_dtheta_dphi_"+charge).fill(Math.toDegrees(forward.theta()-central.theta()), Math.toDegrees(forward.phi()-central.phi()));
        this.getHistos().get("Selection").getH2F("hi_dtheta_dp_"+charge).fill(Math.toDegrees(forward.theta()-central.theta()), forward.p()-central.p());
        this.getHistos().get("Selection").getH2F("hi_theta_dtheta_"+charge).fill(Math.toDegrees(forward.theta()), Math.toDegrees(forward.theta()-central.theta()));
        this.getHistos().get("Selection").getH2F("hi_theta_dphi_"+charge).fill(Math.toDegrees(forward.theta()), Math.toDegrees(forward.phi()-central.phi()));
        if(central.p()>PCUT && Math.toDegrees(central.theta())<THETACUT) {
            this.getHistos().get("Resolution").getH1F("hi_dphi_"+charge).fill(Math.toDegrees(forward.phi()-central.phi()));
            this.getHistos().get("Resolution").getH1F("hi_dtheta_"+charge).fill(Math.toDegrees(forward.theta()-central.theta()));
            this.getHistos().get("Resolution").getH1F("hi_dz_"+charge).fill(forward.vz()-central.vz());
            this.getHistos().get("Resolution").getH1F("hi_dp_"+charge).fill((forward.p()-central.p())/central.p());
        }
        this.getHistos().get("P").getH2F("hi_p_dphi_"+charge).fill(central.p(), Math.toDegrees(forward.phi()-central.phi()));
        this.getHistos().get("P").getH2F("hi_p_dtheta_"+charge).fill(central.p(), Math.toDegrees(forward.theta()-central.theta()));
        this.getHistos().get("P").getH2F("hi_p_dz_"+charge).fill(central.p(), forward.vz()-central.vz());
        this.getHistos().get("P").getH2F("hi_p_dp_"+charge).fill(central.p(), (forward.p()-central.p())/central.p());
        this.getHistos().get("Phi").getH2F("hi_phi_dphi_"+charge).fill(Math.toDegrees(central.phi()), Math.toDegrees(forward.phi()-central.phi()));
        this.getHistos().get("Phi").getH2F("hi_phi_dtheta_"+charge).fill(Math.toDegrees(central.phi()), Math.toDegrees(forward.theta()-central.theta()));
        this.getHistos().get("Phi").getH2F("hi_phi_dz_"+charge).fill(Math.toDegrees(central.phi()), forward.vz()-central.vz());
        this.getHistos().get("Phi").getH2F("hi_phi_dp_"+charge).fill(Math.toDegrees(central.phi()), (forward.p()-central.p())/central.p());
        this.getHistos().get("Theta").getH2F("hi_theta_dphi_"+charge).fill(Math.toDegrees(central.theta()), Math.toDegrees(forward.phi()-central.phi()));
        this.getHistos().get("Theta").getH2F("hi_theta_dtheta_"+charge).fill(Math.toDegrees(central.theta()), Math.toDegrees(forward.theta()-central.theta()));
        this.getHistos().get("Theta").getH2F("hi_theta_dz_"+charge).fill(Math.toDegrees(central.theta()), forward.vz()-central.vz());
        this.getHistos().get("Theta").getH2F("hi_theta_dp_"+charge).fill(Math.toDegrees(central.theta()), (forward.p()-central.p())/central.p());
        
        this.setPhysicsEvent(this.writeToLund(electron));     
    }
    
    @Override
    public void analyzeHistos() {
        for(int i=0; i<2; i++) {
            fitGauss(this.getHistos().get("Resolution").getH1F("hi_dphi_"+charges[i]));
            fitGauss(this.getHistos().get("Resolution").getH1F("hi_dtheta_"+charges[i]));
            fitGauss(this.getHistos().get("Resolution").getH1F("hi_dz_"+charges[i]));
            fitGauss(this.getHistos().get("Resolution").getH1F("hi_dp_"+charges[i]));       
        }
    }
    
    @Override
    public void setPlottingOptions(String name) {
        for(EmbeddedPad p : this.getCanvas().getCanvas(name).getCanvasPads()) {
            p.getAxisZ().setLog(true);
            for(IDataSetPlotter dsp: p.getDatasetPlotters()) {
                IDataSet ds = dsp.getDataSet();
                if(ds instanceof H1F) {
                    H1F h1 = (H1F) ds;
                    h1.setLineWidth(2);
                }
            }
        }
    }

    void fitW(H1F hiw) {

        // get histogram maximum in the rane 0.8-1.2
        int i1=hiw.getXaxis().getBin(0.8);
        int i2=hiw.getXaxis().getBin(1.1);
        double hiMax=0;
        int    imax=i1;
        for(int i=i1; i<=i2; i++) {
            if(hiMax<hiw.getBinContent(i)) {
                imax=i;
                hiMax=hiw.getBinContent(i);
            }
        }           
        double mean = hiw.getDataX(imax); //hiw.getDataX(hiw.getMaximumBin());
        double amp  = hiMax;//hiw.getBinContent(hiw.getMaximumBin());
        double sigma = 0.05;
        F1D f1w = new F1D("f1","[amp]*gaus(x,[mean],[sigma])",0.8, 1.1);
        f1w.setLineColor(2);
        f1w.setLineWidth(2);
        f1w.setOptStat("1111");
        f1w.setParameter(0, amp);
        f1w.setParameter(1, mean);
        f1w.setParameter(2, sigma);
        double rmax = mean + 1.5 * Math.abs(sigma);
        double rmin = mean - 2.0 * Math.abs(sigma);
        f1w.setRange(rmin, rmax);
        DataFitter.fit(f1w, hiw, "Q"); //No options uses error for sigma 
    }


    @Override
    public void fitGauss(H1F hi) {
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
        double rmax = mean + 2.0 * Math.abs(sigma);
        double rmin = mean - 2.0 * Math.abs(sigma);
        f1.setRange(rmin, rmax);
        DataFitter.fit(f1, hi, "Q"); //No options uses error for sigma 
        rmax = f1.getParameter(1) + 2.0 * Math.abs(f1.getParameter(2));
        rmin = f1.getParameter(1) - 2.0 * Math.abs(f1.getParameter(2));
        f1.setRange(rmin, rmax);
        DataFitter.fit(f1, hi, "Q"); //No options uses error for sigma 
    }
    
    PhysicsEvent writeToLund(Particle electron) {
        PhysicsEvent ev = new PhysicsEvent();
        double m = PhysicsConstants.massProton();
        double e = this.getBeamEnergy()+m-electron.e();
        double p = Math.sqrt(e*e-m*m);
        double theta = Math.acos((this.getBeamEnergy()*this.getBeamEnergy()+p*p-Math.pow(this.getBeamEnergy()+m-e,2))/2/this.getBeamEnergy()/p);
        double phi = Math.PI+electron.phi();
        if(Double.isNaN(p) || Double.isNaN(theta)) return null;
        Particle part = new Particle(2212, 
                                     p*Math.sin(theta)*Math.cos(phi), 
                                     p*Math.sin(theta)*Math.sin(phi), 
                                     p*Math.cos(theta), 
                                     electron.vx(), 
                                     electron.vy(), 
                                     electron.vz());
        ev.addParticle(electron);	    
        ev.addParticle(part);	    
        ev.setBeam("e-", this.getBeamEnergy());    
        return ev;
    }

}
