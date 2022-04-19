package modules;

import analysis.Constants;
import objects.Track;
import objects.Event;
import analysis.Module;
import org.jlab.clas.pdg.PDGDatabase;
import org.jlab.clas.pdg.PhysicsConstants;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;

/**
 *
 * @author devita
 */
public class PIDModule extends Module {
    
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
    private final double DVZ = 5.0;
    
    private String[] charges; 

    public PIDModule() {
        super("PID");
    }  

    private DataGroup pidGroup() {
        DataGroup dg = new DataGroup(2,2);
        for(int i=0; i<charges.length; i++) {
            H2F hi_beta = histo2D("hi_beta_"+charges[i], "p (GeV)", "#beta", 200, 0, 2, 200, 0.1, 1.2);
            H1F hi_mass = histo1D("hi_mass_"+charges[i], "M^2 (GeV^2)", "Counts", 200, -1, 6, 42+i*2); 
            dg.addDataSet(hi_beta, i*2 + 0);
            dg.addDataSet(hi_mass, i*2 + 1);
        }
        return dg;
    }
    
    
    @Override
    public void createHistos() {
        charges = new String[]{"pos", "neg"};
        this.getHistos().put("PID", this.pidGroup());
    }
    
    @Override
    public void fillHistos(Event event) {
        if(!event.getParticles().isEmpty() 
         && event.getParticles().get(0).pid()==11
         && event.getParticles().get(0).getDetector()==2) {
            
            Track electron = event.getParticles().get(0);
            
            for(Track hadron : event.getParticles()) {
                if(hadron.charge()!=0 && hadron.getDetector()==4 && hadron.getNDF()>0) {
                    if(Math.abs(hadron.vz()-electron.vz())<VZMAX) {
                        int icharge = 0;
                        if(hadron.charge()<0) icharge = 1;
                        this.getHistos().get("PID").getH2F("hi_beta_"+charges[icharge]).fill(hadron.p(),hadron.getBeta());
                        this.getHistos().get("PID").getH1F("hi_mass_"+charges[icharge]).fill(hadron.p()*hadron.p()*(1/hadron.getBeta()/hadron.getBeta()-1));
                    }
                }
            }
        }
    }
    
    @Override
    public void drawHistos() {
        this.addCanvas("PID");
        this.getCanvas("PID").draw(this.getHistos().get("PID"));
        this.getCanvas("PID").cd(0);
        this.getCanvas("PID").getPad(0).getAxisZ().setLog(true);
        this.getCanvas("PID").draw(this.getF1Beta(PhysicsConstants.massPionCharged()), "same");
        this.getCanvas("PID").draw(this.getF1Beta(PhysicsConstants.massProton()), "same");
        this.getCanvas("PID").draw(this.getF1Beta(PhysicsConstants.massKaonCharged()), "same");
        this.getCanvas("PID").draw(this.getF1Beta(PDGDatabase.getParticleMass(45)), "same");
        this.getCanvas("PID").getPad(1).getAxisY().setLog(true);
        this.getCanvas("PID").cd(2);
        this.getCanvas("PID").getPad(2).getAxisZ().setLog(true);
        this.getCanvas("PID").draw(this.getF1Beta(PhysicsConstants.massPionCharged()), "same");
        this.getCanvas("PID").draw(this.getF1Beta(PhysicsConstants.massKaonCharged()), "same");
        this.getCanvas("PID").getPad(3).getAxisY().setLog(true);
    }
    

    private F1D getF1Beta(double mass) {
        F1D f1 = new F1D("fbeta","x/sqrt(x*x+[m]*[m])", 0.2, 2.0);
        f1.setParameter(0, mass);
        return f1;
    }

}
