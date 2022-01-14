package modules;

import java.util.ArrayList;
import java.util.List;
import objects.Track;
import objects.Event;
import analysis.Module;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;

/**
 *
 * @author devita
 */
public class TrackModule extends Module {
    
    private final double PMIN = 0.0;
    private final double PMAX = 2.0;
    private final double PHIMIN = -180.0;
    private final double PHIMAX = 180.0;
    private final double THETAMIN = 20.0;
    private final double THETAMAX = 140.0;
    private final double VXYMIN = -0.5;//-10;
    private final double VXYMAX = 0.5;//10;
    private final double VZMIN = -10; //-26;
    private final double VZMAX = 4;//26;
    
    private final double CHI2PIDCUT = 10;
    
    
    
    public TrackModule() {
        super("Tracks");
    }
    
    public DataGroup createGroup(int col) {
        H1F hi_mult  = histo1D("hi_mult", "Track multiplicity", "Counts", 10, 0, 10, 48);
        H1F hi_chi2  = histo1D("hi_chi2", "Normalized #chi^2", "Counts", 100, 0, 50, 48);
        H1F hi_ndf   = histo1D("hi_ndf", "NDF", "Counts", 15, 0, 15.0, 48);
        H1F hi_q     = histo1D("hi_q", "Charge", "Counts", 3, -1.5, 1.5, 48);
        H1F hi_p     = histo1D("hi_p", "p (GeV)", "Counts", 100, PMIN, PMAX, col);
        H1F hi_pt    = histo1D("hi_pt", "pt (GeV)", "Counts", 100, PMIN, PMAX, col);
        H1F hi_theta = histo1D("hi_theta", "#theta (deg)", "Counts", 100, THETAMIN, THETAMAX, col);
        H1F hi_phi   = histo1D("hi_phi", "#phi (deg)", "Counts", 100, PHIMIN, PHIMAX, col);
        H1F hi_vx    = histo1D("hi_vx", "vx (cm)", "Counts", 100, VXYMIN, VXYMAX, col);
        H1F hi_vy    = histo1D("hi_vy", "vy (cm)", "Counts", 100, VXYMIN, VXYMAX, col);
        H1F hi_vz    = histo1D("hi_vz", "vz (cm)", "Counts", 100, VZMIN, VZMAX, col);
        H1F hi_type  = histo1D("hi_type", "Seed type", "Counts", 4, -0.5, 3.5, col);

        DataGroup dgTrack = new DataGroup(4,3);
        dgTrack.addDataSet(hi_mult, 0);
        dgTrack.addDataSet(hi_chi2, 1);
        dgTrack.addDataSet(hi_ndf,  2);
        dgTrack.addDataSet(hi_q,    3);
        dgTrack.addDataSet(hi_p,    4);
        dgTrack.addDataSet(hi_pt,   5);
        dgTrack.addDataSet(hi_theta,6);
        dgTrack.addDataSet(hi_phi,  7);
        dgTrack.addDataSet(hi_vx,   8);
        dgTrack.addDataSet(hi_vy,   9);
        dgTrack.addDataSet(hi_vz,   10);
        dgTrack.addDataSet(hi_type, 11);
        return dgTrack;
    }

    public DataGroup createGroup2D(int col) {
        H2F hi_ptheta   = histo2D("hi_ptheta", "#theta (deg)", "p (GeV)", 100, THETAMIN, THETAMAX, 100, PMIN, PMAX);
        H2F hi_pphi     = histo2D("hi_pphi", "#phi (deg)", "p (GeV)", 100, PHIMIN, PHIMAX, 100, PMIN, PMAX);
        H2F hi_thetaphi = histo2D("hi_thetaphi", "#phi (deg)", "#theta (deg)", 100, PHIMIN, PHIMAX, 100, THETAMIN, THETAMAX);
        H2F hi_vxy      = histo2D("hi_vxy", "vx (cm)", "vy (cm)", 100, VXYMIN, VXYMAX, 100, VXYMIN, VXYMAX);
        H2F hi_d0phi    = histo2D("hi_d0phi", "#phi (deg)", "d0 (cm)", 100, PHIMIN, PHIMAX, 100, VXYMIN, VXYMAX);
        H2F hi_vzphi    = histo2D("hi_vzphi", "#phi (deg)", "vz (cm)", 100, PHIMIN, PHIMAX, 100, VZMIN, VZMAX);

        DataGroup dgTrack = new DataGroup(3,2);
        dgTrack.addDataSet(hi_ptheta,   0);
        dgTrack.addDataSet(hi_pphi,     1);
        dgTrack.addDataSet(hi_thetaphi, 2);
        dgTrack.addDataSet(hi_vxy,      3);
        dgTrack.addDataSet(hi_d0phi,    4);
        dgTrack.addDataSet(hi_vzphi,   5);
        return dgTrack;
    }

    @Override
    public void createHistos() {
        this.getHistos().put("Tracks", this.createGroup(46));
        this.getHistos().put("Tracks2D", this.createGroup2D(46));
        this.getHistos().put("Seeds", this.createGroup(46));
        this.getHistos().put("TrackSeeds", this.createGroup(46));
        this.getHistos().put("TrackChi2pid", this.createGroup(49));
    }
    
    @Override
    public void fillHistos(Event event) {
        List<Track> trackSeeds = new ArrayList<>();
        List<Track> trackC2pid = new ArrayList<>();
        for(Track track : event.getTracks()) {
            int sid = track.getSeedId();
            if(event.getSeedMap().containsKey(sid)) {
                int si = event.getSeedMap().get(sid);           
                trackSeeds.add(event.getSeeds().get(si));
            }
            if(Math.abs(track.getChi2pid())<CHI2PIDCUT) trackC2pid.add(track);
        }
        this.fillGroup(this.getHistos().get("Tracks"),event.getTracks());
        this.fillGroup2D(this.getHistos().get("Tracks2D"),event.getTracks());
        this.fillGroup(this.getHistos().get("Seeds"),event.getSeeds());
        this.fillGroup(this.getHistos().get("TrackSeeds"),trackSeeds);
        this.fillGroup(this.getHistos().get("TrackChi2pid"),trackC2pid);
    }
    
    public void fillGroup(DataGroup group, List<Track> tracks) {
        group.getH1F("hi_mult").fill(tracks.size());
        for(Track track : tracks) {
            group.getH1F("hi_chi2").fill(track.getChi2()/track.getNDF());
            group.getH1F("hi_ndf").fill(track.getNDF());
            group.getH1F("hi_q").fill(track.charge());
            group.getH1F("hi_p").fill(track.p());
            group.getH1F("hi_pt").fill(track.pt());
            group.getH1F("hi_theta").fill(Math.toDegrees(track.theta()));
            group.getH1F("hi_phi").fill(Math.toDegrees(track.phi()));
            group.getH1F("hi_vx").fill(track.vx());
            group.getH1F("hi_vy").fill(track.vy());
            group.getH1F("hi_vz").fill(track.vz());
            group.getH1F("hi_type").fill(track.getSeedType());
        }
    }
    
    public void fillGroup2D(DataGroup group, List<Track> tracks) {
        for(Track track : tracks) {
            group.getH2F("hi_ptheta").fill(Math.toDegrees(track.theta()),track.p());
            group.getH2F("hi_pphi").fill(Math.toDegrees(track.phi()),track.p());
            group.getH2F("hi_thetaphi").fill(Math.toDegrees(track.phi()),Math.toDegrees(track.theta()));
            group.getH2F("hi_vxy").fill(track.vx(),track.vy());
            group.getH2F("hi_d0phi").fill(Math.toDegrees(track.phi()),track.d0());
            group.getH2F("hi_vzphi").fill(Math.toDegrees(track.phi()),track.vz());
        }
    }
    
    @Override
    public EmbeddedCanvasTabbed plotHistos() {
        EmbeddedCanvasTabbed canvas = new EmbeddedCanvasTabbed("Tracks", "Tracks2D", "Seeds", "EBTracks");
        canvas.getCanvas("Tracks").draw(this.getHistos().get("Tracks"));
        canvas.getCanvas("Tracks2D").draw(this.getHistos().get("Tracks2D"));
        canvas.getCanvas("Seeds").draw(this.getHistos().get("Seeds"));
        canvas.getCanvas("Seeds").draw(this.getHistos().get("TrackSeeds"));
        canvas.getCanvas("EBTracks").draw(this.getHistos().get("Tracks"));
        canvas.getCanvas("EBTracks").draw(this.getHistos().get("TrackChi2pid"));
        this.setPlottingOptions(canvas.getCanvas("Tracks"));
        this.setPlottingOptions(canvas.getCanvas("Seeds"));
        this.setPlottingOptions(canvas.getCanvas("EBTracks"));
        return canvas;
    }
       
    @Override
    public void setPlottingOptions(EmbeddedCanvas canvas) {
        canvas.setGridX(false);
        canvas.setGridY(false);
        canvas.getPad(1).getAxisY().setLog(true);
        canvas.getPad(8).getAxisY().setLog(true);
        canvas.getPad(9).getAxisY().setLog(true);
    }

}
