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
public class ElasticModule extends Module {
    
    private final double MINW = 0.6;
    private final double MAXW = 4.5;
    private final double MINQ = 0.5;
    private final double MAXQ = 5.0;
    private final double PHIMIN = -180.0;
    private final double PHIMAX = 180.0;
    private final double DPHI = 10.0;
    private final double THETAMIN = 30.0;
    private final double THETAMAX = 70.0;
    private final double DTHETA   = 30.0;
    private final double VZMIN = -10; 
    private final double VZMAX =  10;

    public ElasticModule(double ebeam, boolean lund) {
        super("Elastic", false, ebeam, lund);
    }
    
    private DataGroup generalGroup() {
        System.out.println(this.getBeamEnergy());
        H2F hi_q2w   = histo2D("hi_q2w", "W (GeV)", "Q^2 (GeV^2)", 100, MINW, MAXW, 100, MINQ, MAXQ); 
        H1F hi_w     = histo1D("hi_w", "W (GeV)", "Counts", 100, MINW, MAXW, 0);
        H2F hi_wphi  = histo2D("hi_wphi", "#phi (deg)", "W (GeV)", 100, PHIMIN, PHIMAX, 100, MINW, MAXW); 
        H2F hi_el    = histo2D("hi_el", "p_e (GeV)", "#theta_e (deg)", 100, 0.5, this.getBeamEnergy()+0.5, 100, 0.0, THETAMIN); 
        F1D f1_el    = new F1D("f1_el", "2*(180/3.14)*atan(sqrt(0.93832*([e0]-x)/2/[e0]/x))", this.getBeamEnergy()*0.75, this.getBeamEnergy()*0.99);
        f1_el.setParameter(0, this.getBeamEnergy());
        H2F hi_pr    = histo2D("hi_pr", "p_p (GeV)", "#theta_p (deg)", 100, 0, this.getBeamEnergy()/2-0.1, 100, THETAMIN, THETAMAX); 
        F1D f1_pr = new F1D("f1_pr", "(180/3.14)*acos(([e0]*[e0]+x*x-pow(([e0]+0.93832-sqrt(x*x+0.9382*0.9382)),2))/2/[e0]/x)", this.getBeamEnergy()*0.08, this.getBeamEnergy()*0.35);
        f1_pr.setParameter(0, this.getBeamEnergy());
        H2F hi_dphi  = histo2D("hi_dphi", "#phi_p (deg)", "#Delta#phi (deg)", 100, PHIMIN, PHIMAX, 100, PHIMAX-DPHI, PHIMAX+DPHI); 

        DataGroup dg = new DataGroup(3,2);
        dg.addDataSet(hi_q2w,  0);
        dg.addDataSet(hi_w,    1);
        dg.addDataSet(hi_el,   2);
        dg.addDataSet(f1_el,   2);
        dg.addDataSet(hi_wphi, 3);
        dg.addDataSet(hi_pr,   4);
        dg.addDataSet(f1_pr,   4);
        dg.addDataSet(hi_dphi, 5);
        return dg;
    }

    private DataGroup protonGroup() {
        DataGroup dg = new DataGroup(4,2);
        H2F hi_p_dp    = histo2D("hi_p_dp", "p (GeV)", "Ebeam (GeV)", 100, 1, 3, 100, this.getBeamEnergy()*0.75, this.getBeamEnergy()*1.2);
        H2F hi_p_theta = histo2D("hi_p_theta", "p (GeV)", "#theta (deg)", 100, 0, this.getBeamEnergy()/2-0.1, 100, THETAMIN, THETAMAX); 
        F1D f1_pr = new F1D("f1_pr", "(180/3.14)*acos(([e0]*[e0]+x*x-pow(([e0]+0.93832-sqrt(x*x+0.9382*0.9382)),2))/2/[e0]/x)", this.getBeamEnergy()*0.08, this.getBeamEnergy()*0.35);
        f1_pr.setParameter(0, this.getBeamEnergy());
        H2F hi_p_dphi  = histo2D("hi_phi_dphi", "#phi (deg)", "#Delta#phi (deg)", 100, -180, 180, 100, PHIMAX-DPHI, PHIMAX+DPHI); 
        H2F hi_p_dz    = histo2D("hi_phi_dz", "#phi (deg)", "#Deltaz (cm)", 100, -180, 180, 100, VZMIN, VZMAX);
        H1F hi_dp      = histo1D("hi_dp", "Ebeam (GeV)", "Counts", 100, this.getBeamEnergy()*0.75, this.getBeamEnergy()*1.2, 0);
        H1F hi_dtheta  = histo1D("hi_dtheta", "#Delta#theta (deg)", "Counts", 100, -DTHETA, DTHETA, 0);
        H1F hi_dphi    = histo1D("hi_dphi", "#Delta#phi (deg)", "Counts", 100, PHIMAX-DPHI, PHIMAX+DPHI, 0); 
        H1F hi_dz      = histo1D("hi_dz", "#Deltaz (cm)", "Counts", 100, VZMIN, VZMAX, 0);
        dg.addDataSet(hi_p_dp,    0);
        dg.addDataSet(hi_p_theta, 1);
        dg.addDataSet(f1_pr,      1);
        dg.addDataSet(hi_p_dphi,  2);
        dg.addDataSet(hi_p_dz,    3);
        dg.addDataSet(hi_dp,      4);
        dg.addDataSet(hi_dtheta,  5);
        dg.addDataSet(hi_dphi,    6);
        dg.addDataSet(hi_dz,      7);
        return dg;
    }

    private DataGroup wGroup() {
        DataGroup dg = new DataGroup(3,2);
        for(int i=0; i<6; i++) {
            int sector = i+1;
            H1F hi_w     = histo1D("hi_w"+sector, "W - Sector " + sector + " (GeV)", "Counts", 100, MINW, MAXW, 0);
            dg.addDataSet(hi_w, i);
        }
        return dg;
    }

    private DataGroup phiGroup() {
        DataGroup dg = new DataGroup(3,2);
        for(int i=0; i<6; i++) {
            int sector = i+1;
            H1F hi_dphi     = histo1D("hi_dphi"+sector, "#Delta#phi - Sector" + sector + " (deg)", "Counts", 100, 180-DPHI, 180+DPHI, 0);
            dg.addDataSet(hi_dphi, i);
        }
        return dg;
    }

    private DataGroup sectorGroup() {
        DataGroup dg = new DataGroup(3,2);
        for(int i=0; i<6; i++) {
            int sector = i+1;
            H2F hi_p_theta = histo2D("hi_p_theta"+sector, "p (GeV)", "#theta (deg)", 100, 0, this.getBeamEnergy()/2-0.1, 100, THETAMIN, THETAMAX); 
            F1D f1_p_theta = new F1D("f1_pr", "(180/3.14)*acos(([e0]*[e0]+x*x-pow(([e0]+0.93832-sqrt(x*x+0.9382*0.9382)),2))/2/[e0]/x)", this.getBeamEnergy()*0.08, this.getBeamEnergy()*0.35);
            f1_p_theta.setParameter(0, this.getBeamEnergy());
            dg.addDataSet(hi_p_theta, i);
            dg.addDataSet(f1_p_theta, i);
        }
        return dg;
    }

    private DataGroup thetaGroup() {
        DataGroup dg = new DataGroup(3,2);
        for(int i=0; i<6; i++) {
            int sector = i+1;
            H1F hi_dtheta = histo1D("hi_dtheta"+sector, "#Delta#theta (deg)", "Counts", 100, -DTHETA, DTHETA, 0); 
            dg.addDataSet(hi_dtheta, i);
        }
        return dg;
    }

    private DataGroup ebeamGroup() {
        DataGroup dg = new DataGroup(3,2);
        for(int i=0; i<6; i++) {
            int sector = i+1;
            H1F hi_dp      = histo1D("hi_dp"+sector, "Ebeam (GeV)", "Counts", 100, this.getBeamEnergy()*0.75, this.getBeamEnergy()*1.2, 0);
            dg.addDataSet(hi_dp, i);
        }
        return dg;
    }

    private DataGroup efficiencyGroup() {
        DataGroup dg = new DataGroup(1,3);
        for(int i=0; i<3; i++) {
            H1F hi_eff      = histo1D("hi_eff"+i, "#phi (deg)", "Counts", 200, PHIMIN, PHIMAX, 0);
            dg.addDataSet(hi_eff, i);
        }
        return dg;
    }
    
    @Override
    public void createHistos() {
        this.getHistos().put("Efficiency", this.efficiencyGroup());
        this.getHistos().put("Ebeam",      this.ebeamGroup());
        this.getHistos().put("Theta",      this.thetaGroup());
        this.getHistos().put("Sector",     this.sectorGroup());
        this.getHistos().put("Phi",        this.phiGroup());
        this.getHistos().put("W",          this.wGroup());
        this.getHistos().put("Proton",     this.protonGroup());
        this.getHistos().put("General",    this.generalGroup());
    }
    
    @Override
    public void fillHistos(Event event) {
        if(!event.getParticles().isEmpty() 
         && event.getParticles().get(0).pid()==11
         && event.getParticles().get(0).getDetector()==2) {
            
            Track electron = event.getParticles().get(0);
            LorentzVector virtualPhoton = new LorentzVector(0.0, 0.0, this.getBeamEnergy(), this.getBeamEnergy());
            virtualPhoton.sub(electron.vector());
            LorentzVector hadronSystem = new LorentzVector(0.0, 0.0, this.getBeamEnergy(), PhysicsConstants.massProton()+this.getBeamEnergy());
            hadronSystem.sub(electron.vector());
            this.getHistos().get("General").getH2F("hi_q2w").fill(hadronSystem.mass(),-virtualPhoton.mass2());
            this.getHistos().get("General").getH1F("hi_w").fill(hadronSystem.mass());
            this.getHistos().get("General").getH2F("hi_el").fill(electron.p(),Math.toDegrees(electron.theta()));
            this.getHistos().get("General").getH2F("hi_wphi").fill(Math.toDegrees(electron.phi()), hadronSystem.mass());
            this.getHistos().get("W").getH1F("hi_w"+electron.getSector()).fill(hadronSystem.mass());
            
            if(hadronSystem.mass()>1.1) return;
            this.getHistos().get("Efficiency").getH1F("hi_eff0").fill(Math.toDegrees(electron.phi()));

            for(Track proton : event.getParticles()) {
                if(proton.charge()>0 && proton.getDetector()==4 && proton.getNDF()>0) {
                    double dphi = Math.abs(Math.toDegrees(proton.deltaPhi(electron)));
                    this.getHistos().get("General").getH2F("hi_dphi").fill(Math.toDegrees(proton.phi()),dphi);
                    this.getHistos().get("Proton").getH2F("hi_phi_dphi").fill(Math.toDegrees(proton.phi()), dphi);
                    this.getHistos().get("Proton").getH2F("hi_phi_dz").fill(Math.toDegrees(proton.phi()),proton.vz()-electron.vz());
                    this.getHistos().get("Proton").getH1F("hi_dphi").fill(dphi);
                    this.getHistos().get("Proton").getH1F("hi_dz").fill(proton.vz()-electron.vz());
                    this.getHistos().get("Phi").getH1F("hi_dphi"+electron.getSector()).fill(dphi);
                    if(Math.abs(dphi-180)<DPHI && Math.abs(proton.vz()-electron.vz())<VZMAX) {
                        this.getHistos().get("General").getH2F("hi_pr").fill(proton.p(),Math.toDegrees(proton.theta()));
                        this.getHistos().get("Proton").getH2F("hi_p_dp").fill(proton.p(),-PhysicsConstants.massProton()+proton.e()+electron.p());
                        this.getHistos().get("Proton").getH2F("hi_p_theta").fill(proton.p(),Math.toDegrees(proton.theta()));
                        this.getHistos().get("Proton").getH1F("hi_dp").fill(-PhysicsConstants.massProton()+proton.e()+electron.p());
                        this.getHistos().get("Proton").getH1F("hi_dtheta").fill(Math.toDegrees(proton.theta())-this.getHistos().get("Proton").getF1D("f1_pr").evaluate(proton.p()));
                        this.getHistos().get("Sector").getH2F("hi_p_theta"+electron.getSector()).fill(proton.p(),Math.toDegrees(proton.theta()));
                        this.getHistos().get("Theta").getH1F("hi_dtheta"+electron.getSector()).fill(Math.toDegrees(proton.theta())-this.getHistos().get("Proton").getF1D("f1_pr").evaluate(proton.p()));
                        this.getHistos().get("Ebeam").getH1F("hi_dp"+electron.getSector()).fill(-PhysicsConstants.massProton()+proton.e()+electron.p());
                        this.getHistos().get("Efficiency").getH1F("hi_eff1").fill(Math.toDegrees(electron.phi()));
                    }
                }
            }
            this.setPhysicsEvent(this.writeToLund(electron)); 
        }
    }
    
    @Override
    public void analyzeHistos() {
        fitW(this.getHistos().get("General").getH1F("hi_w"));
        for(int sector=1; sector <= 6; sector++) {
            fitW(this.getHistos().get("W").getH1F("hi_w" + sector));
            fitGauss(this.getHistos().get("Phi").getH1F("hi_dphi" + sector));
            fitGauss(this.getHistos().get("Theta").getH1F("hi_dtheta" + sector));
            fitGauss(this.getHistos().get("Ebeam").getH1F("hi_dp" + sector));
        }
        fitGauss(this.getHistos().get("Proton").getH1F("hi_dp"));
        fitGauss(this.getHistos().get("Proton").getH1F("hi_dtheta"));
        fitGauss(this.getHistos().get("Proton").getH1F("hi_dphi"));
        fitGauss(this.getHistos().get("Proton").getH1F("hi_dz"));
        
        H1F he1  = this.getHistos().get("Efficiency").getH1F("hi_eff0");
        H1F he2  = this.getHistos().get("Efficiency").getH1F("hi_eff1");
        H1F he3  = this.getHistos().get("Efficiency").getH1F("hi_eff2");
        for(int i=0; i<he1.getDataSize(0); i++) {
            double e1 = he1.getBinContent(i);
            double e2 = he2.getBinContent(i);
            double e3 = e2/e1;
            if(e1>10) {
                he3.setBinContent(i,e3);
                he3.setBinError(i,e3*Math.sqrt(e3*(1-e3)/e1));
            }
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
        double rmax = mean + 2.5 * Math.abs(sigma);
        double rmin = mean - 2.5 * Math.abs(sigma);
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
