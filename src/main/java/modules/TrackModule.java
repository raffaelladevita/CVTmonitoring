package modules;

import java.util.ArrayList;
import java.util.List;
import objects.Track;
import objects.Event;
import analysis.Module;
import org.jlab.groot.data.H1F;
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
    private final double VXYMIN = -2;//-10;
    private final double VXYMAX = 2;//10;
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

    @Override
    public void createHistos() {
        this.getHistos().put("Tracks", this.createGroup(46));
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
            int si  = event.getSeedMap().get(sid);
            trackSeeds.add(event.getSeeds().get(si));
            if(Math.abs(track.getChi2pid())<CHI2PIDCUT) trackC2pid.add(track);
        }
        this.fillGroup(this.getHistos().get("Tracks"),event.getTracks());
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
    
    @Override
    public EmbeddedCanvasTabbed plotHistos() {
        EmbeddedCanvasTabbed canvas = new EmbeddedCanvasTabbed("Tracks", "Seeds", "EBTracks");
        canvas.getCanvas("Tracks").draw(this.getHistos().get("Tracks"));
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
