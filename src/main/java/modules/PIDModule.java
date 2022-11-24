package modules;

import analysis.Constants;
import objects.Track;
import objects.Event;
import analysis.Module;
import java.util.List;
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
private int ntot, nacc, nchange;
    public PIDModule() {
        super("PID");
    }  

    private DataGroup pidGroup() {
        DataGroup dg = new DataGroup(3,2);
        for(int i=0; i<charges.length; i++) {
            H2F hi_elos = histo2D("hi_elos_"+charges[i], "p (GeV)", "#Deltap", 200, 0, 2, 200, -0.05, 0.25);
            H2F hi_beta = histo2D("hi_beta_"+charges[i], "p (GeV)", "#beta", 200, 0, 2, 200, 0.1, 1.2);
            H1F hi_mass = histo1D("hi_mass_"+charges[i], "M^2 (GeV^2)", "Counts", 200, -1, 6, 42+i*2); 
            dg.addDataSet(hi_elos, i*3 + 0);
            dg.addDataSet(hi_beta, i*3 + 1);
            dg.addDataSet(hi_mass, i*3 + 2);
        }
        return dg;
    }
    
    
    @Override
    public void createHistos() {
        charges = new String[]{"pos", "neg"};
        this.getHistos().put("Tracks", this.pidGroup());
        this.getHistos().put("FPTracks", this.pidGroup());
        this.getHistos().put("Seeds", this.pidGroup());
        this.getHistos().put("PID change", this.pidGroup());
        this.getHistos().put("Accidentals", this.pidGroup());
    }
    
    @Override
    public void fillHistos(Event event) {
        if(!event.getParticles().isEmpty() 
         && event.getParticles().get(0).pid()==11
         && event.getParticles().get(0).getDetector()==2) {
            
            Track electron = event.getParticles().get(0);
            
            for(Track hadron : event.getTracks()) {
                if(hadron.charge()!=0 && hadron.getNDF()>0) {
                    if(Math.abs(hadron.vz()-electron.vz())<3) {
                        int icharge = 0;
                        if(hadron.charge()<0) icharge = 1;
                        int tid = hadron.getId();
                        int sid = hadron.getSeedId();
                        Track fpass = null;
                        if(event.getFPTrackMap().containsKey(tid))
                            fpass = event.getFPTracks().get(event.getFPTrackMap().get(tid));
                        
                        this.fillGroup(hadron, fpass, icharge, "Tracks");
                        
                        if(hadron.getElossPid()!=hadron.getEBPid() && hadron.getEBPid()!=0) 
                            this.fillGroup(hadron, fpass, icharge, "PID change");
                        
                        if(hadron.getEBPid()==45 && hadron.getChi2pid()>5) 
                            this.fillGroup(hadron, fpass, icharge, "Accidentals");
                        
                        if(event.getSeedMap().containsKey(sid)) {
                            int si  = event.getSeedMap().get(sid);
                            Track seed = event.getSeeds().get(si);
                            this.fillGroup(seed, null, icharge, "Seeds");
                        }
                    }
                }
            }
            for(Track hadron : event.getFPTracks()) {
                if(hadron.charge()!=0 && hadron.getNDF()>0) {
                    if(Math.abs(hadron.vz()-electron.vz())<3) {
                        int icharge = 0;
                        if(hadron.charge()<0) icharge = 1;
                        this.fillGroup(hadron, null, icharge, "FPTracks");
                    }     
                }
            }
        }
    }
    
    public void fillGroup(Track spass, Track fpass, int icharge, String name) {
        if(fpass!=null) 
                this.getHistos().get(name).getH2F("hi_elos_"+charges[icharge]).fill(fpass.p(),spass.p()-fpass.p());
            this.getHistos().get(name).getH2F("hi_beta_"+charges[icharge]).fill(spass.p(),spass.getBeta());
            this.getHistos().get(name).getH1F("hi_mass_"+charges[icharge]).fill(spass.p()*spass.p()*(1/spass.getBeta()/spass.getBeta()-1));                   
    }
    
    @Override
    public void drawHistos() {
        this.drawGroup("Accidentals");
        this.drawGroup("PID change");
        this.drawGroup("Seeds");
        this.drawGroup("FPTracks");
        this.drawGroup("Tracks");
    }
    
    private void drawGroup(String name) {
        this.addCanvas(name);
        this.getCanvas(name).draw(this.getHistos().get(name));
        this.getCanvas(name).cd(0);
        this.getCanvas(name).getPad(0).getAxisZ().setLog(true);
        this.getCanvas(name).cd(1);
        this.getCanvas(name).getPad(1).getAxisZ().setLog(true);
        this.getCanvas(name).draw(this.getF1Beta(PhysicsConstants.massPionCharged()), "same");
        this.getCanvas(name).draw(this.getF1Beta(PhysicsConstants.massProton()), "same");
        this.getCanvas(name).draw(this.getF1Beta(PhysicsConstants.massKaonCharged()), "same");
        this.getCanvas(name).draw(this.getF1Beta(PDGDatabase.getParticleMass(45)), "same");
        this.getCanvas(name).getPad(2).getAxisY().setLog(true);
        this.getCanvas(name).cd(3);
        this.getCanvas(name).getPad(3).getAxisZ().setLog(true);
        this.getCanvas(name).cd(4);
        this.getCanvas(name).getPad(4).getAxisZ().setLog(true);
        this.getCanvas(name).draw(this.getF1Beta(PhysicsConstants.massPionCharged()), "same");
        this.getCanvas(name).draw(this.getF1Beta(PhysicsConstants.massKaonCharged()), "same");
        this.getCanvas(name).getPad(5).getAxisY().setLog(true);        
    }
    
    private F1D getF1Beta(double mass) {
        F1D f1 = new F1D("fbeta","x/sqrt(x*x+[m]*[m])", 0.2, 2.0);
        f1.setParameter(0, mass);
        return f1;
    }

}
