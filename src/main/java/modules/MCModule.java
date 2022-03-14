package modules;

import analysis.Constants;
import objects.Track;
import objects.Event;
import analysis.Module;
import org.jlab.clas.pdg.PhysicsConstants;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;

/**
 *
 * @author devita
 */
public class MCModule extends Module {
    
    private final double PMIN = 0.0;
    private final double PMAX = 2.0;
    private final double DP   = 0.2;
    private final double PHIMIN = -180.0;
    private final double PHIMAX = 180.0;
    private final double DPHI   = 2.0;
    private final double THETAMIN = 20.0;
    private final double THETAMAX = 140.0;
    private final double DTHETA   = 2.0;
    private final double VXYMIN = -0.2;//-10;
    private final double VXYMAX = 0.2;//10;
    private final double DVXY   = 0.1;//10;
    private final double VZMIN = -10; //-26;
    private final double VZMAX =  10;//26;
    private final double DVZ   = 0.5;//26;
    private final double DTX   = 0.03;;
    private final double DTZ   = 0.03;;


    public MCModule(boolean cosmics) {
        super("MC", cosmics);
    }
    
    private DataGroup trackGroup(int col) {
        H1F hi_q     = histo1D("hi_q", "Charge", "Counts", 3, -1.5, 1.5, col);
        H1F hi_p     = histo1D("hi_p", "p (GeV)", "Counts", 100, PMIN, PMAX, col);
        H1F hi_pt    = histo1D("hi_pt", "pt (GeV)", "Counts", 100, PMIN, PMAX, col);
        H1F hi_theta = histo1D("hi_theta", "#theta (deg)", "Counts", 100, THETAMIN, THETAMAX, col);
        H1F hi_phi   = histo1D("hi_phi", "#phi (deg)", "Counts", 100, PHIMIN, PHIMAX, col);
        H1F hi_d0    = histo1D("hi_d0", "d0 (cm)", "Counts", 100, VXYMIN, VXYMAX, col);
        H1F hi_vx    = histo1D("hi_vx", "vx (cm)", "Counts", 100, VXYMIN, VXYMAX, col);
        H1F hi_vy    = histo1D("hi_vy", "vy (cm)", "Counts", 100, VXYMIN, VXYMAX, col);
        H1F hi_vz    = histo1D("hi_vz", "vz (cm)", "Counts", 100, VZMIN, VZMAX, col);

        DataGroup dgTrack = new DataGroup(3,3);
        dgTrack.addDataSet(hi_q,    0);
        dgTrack.addDataSet(hi_p,    1);
        dgTrack.addDataSet(hi_pt,   2);
        dgTrack.addDataSet(hi_theta,3);
        dgTrack.addDataSet(hi_phi,  4);
        dgTrack.addDataSet(hi_d0,   5);
        dgTrack.addDataSet(hi_vx,   6);
        dgTrack.addDataSet(hi_vy,   7);
        dgTrack.addDataSet(hi_vz,   8);
        return dgTrack;
    }

    private DataGroup trackResolutionGroup(int icol) {
        H1F hi_p     = histo1D("hi_p", "#Deltap/p", "Counts", 100, -DP, DP, icol);
        H1F hi_theta = histo1D("hi_theta", "#Delta#theta (deg)", "Counts", 100, -DTHETA, DTHETA, icol);
        H1F hi_phi   = histo1D("hi_phi", "#Delta#phi (deg)", "Counts", 100, -DPHI, DPHI, icol);
        H1F hi_vx    = histo1D("hi_vx", "#Deltavx (cm)", "Counts", 100, -DVXY, DVXY, icol);
        H1F hi_vy    = histo1D("hi_vy", "#Deltavy (cm)", "Counts", 100, -DVXY, DVXY, icol);
        if(this.isCosmics()) {
            hi_vx.set(100, -DVZ, DVZ);
        }
        H1F hi_vz    = histo1D("hi_vz", "#Deltavz (cm)", "Counts", 100, -DVZ,  DVZ,  icol);

        DataGroup dg = new DataGroup(3,2);
        dg.addDataSet(hi_p,    0);
        dg.addDataSet(hi_theta,1);
        dg.addDataSet(hi_phi,  2);
        dg.addDataSet(hi_vx,   3);
        dg.addDataSet(hi_vy,   4);
        dg.addDataSet(hi_vz,   5);
        return dg;
    }

    private DataGroup trackResolution2DGroup() {
        H2F hi_p_dp     = histo2D("hi_p_dp",     "p (GeV)", "#Deltap/p",          50, PMIN, PMAX, 50, -DP, DP);
        H2F hi_p_dtheta = histo2D("hi_p_dtheta", "p (GeV)", "#Delta#theta (deg)", 50, PMIN, PMAX, 50, -DTHETA, DTHETA);
        H2F hi_p_dphi   = histo2D("hi_p_dphi",   "p (GeV)", "#Delta#phi (deg)",   50, PMIN, PMAX, 50, -DPHI, DPHI);
        H2F hi_p_dvz    = histo2D("hi_p_dvz",    "p (GeV)", "#Deltavz (cm)",      50, PMIN, PMAX, 50, -DVZ,  DVZ);
        H2F hi_theta_dp     = histo2D("hi_theta_dp",     "#theta (deg)", "#Deltap/p",          50, THETAMIN, THETAMAX, 50, -DP, DP);
        H2F hi_theta_dtheta = histo2D("hi_theta_dtheta", "#theta (deg)", "#Delta#theta (deg)", 50, THETAMIN, THETAMAX, 50, -DTHETA, DTHETA);
        H2F hi_theta_dphi   = histo2D("hi_theta_dphi",   "#theta (deg)", "#Delta#phi (deg)",   50, THETAMIN, THETAMAX, 50, -DPHI, DPHI);
        H2F hi_theta_dvz    = histo2D("hi_theta_dvz",    "#theta (deg)", "#Deltavz (cm)",      50, THETAMIN, THETAMAX, 50, -DVZ,  DVZ);
        H2F hi_phi_dp     = histo2D("hi_phi_dp",     "#phi (deg)", "#Deltap/p",          50, PHIMIN, PHIMAX, 50, -DP, DP);
        H2F hi_phi_dtheta = histo2D("hi_phi_dtheta", "#phi (deg)", "#Delta#theta (deg)", 50, PHIMIN, PHIMAX, 50, -DTHETA, DTHETA);
        H2F hi_phi_dphi   = histo2D("hi_phi_dphi",   "#phi (deg)", "#Delta#phi (deg)",   50, PHIMIN, PHIMAX, 50, -DPHI, DPHI);
        H2F hi_phi_dvz    = histo2D("hi_phi_dvz",    "#phi (deg)", "#Deltavz (cm)",      50, PHIMIN, PHIMAX, 50, -DVZ,  DVZ);

        DataGroup dg = new DataGroup(4,3);
        dg.addDataSet(hi_p_dp,        0);
        dg.addDataSet(hi_p_dtheta,    1);
        dg.addDataSet(hi_p_dphi,      2);
        dg.addDataSet(hi_p_dvz,       3);
        dg.addDataSet(hi_theta_dp,    4);
        dg.addDataSet(hi_theta_dtheta,5);
        dg.addDataSet(hi_theta_dphi,  6);
        dg.addDataSet(hi_theta_dvz,   7);
        dg.addDataSet(hi_phi_dp,      8);
        dg.addDataSet(hi_phi_dtheta,  9);
        dg.addDataSet(hi_phi_dphi,    10);
        dg.addDataSet(hi_phi_dvz,     11);
        return dg;
    }

    private DataGroup helixResolutionGroup(int icol) {
        H1F hi_chi2   = histo1D("hi_chi2", "#chi^2", "Counts", 100, 0, 5,  icol);
        H1F hi_d0     = histo1D("hi_d0", "#Deltad0 (cm)", "Counts", 100, -DVXY, DVXY, icol);
        H1F hi_phi0   = histo1D("hi_phi0", "#Delta#phi (deg)", "Counts", 100, -DPHI, DPHI, icol);
        H1F hi_rho    = histo1D("hi_rho", "#Delta#rho/#rho", "Counts", 100, -DP, DP, icol);
        H1F hi_z0     = histo1D("hi_z0", "#Deltav0 (cm)", "Counts", 100, -DVZ, DVZ, icol);
        H1F hi_tandip = histo1D("hi_tandip", "#DeltaTanDip", "Counts", 100, -DVXY, DVXY, icol);

        DataGroup dg = new DataGroup(3,2);
        dg.addDataSet(hi_chi2,   0);
        dg.addDataSet(hi_d0,     1);
        dg.addDataSet(hi_phi0,   2);
        dg.addDataSet(hi_rho,    3);
        dg.addDataSet(hi_z0,     4);
        dg.addDataSet(hi_tandip, 5);
        return dg;
    }

    private DataGroup helixResolution2DGroup() {
        DataGroup dg = new DataGroup(5,3);
        for(int i=0; i<3; i++) {
            String type = "";
            String titl = "#Delta";
            if(i<2) {
                if(i==1) {
                    type = "cov2D";
                    titl = "cov";
                }
                H2F hi_d0     = histo2D("hi_"+type+"d0", titl+"d0 (cm)", "d0 (cm)", 100, -DVXY, DVXY,100, 0, DVXY);
                H2F hi_phi0   = histo2D("hi_"+type+"phi0", titl+"#phi (deg)", "#phi0", 100, -DPHI, DPHI, 100, PHIMIN, PHIMAX);
                H2F hi_rho    = histo2D("hi_"+type+"rho",  titl+"#rho/#rho", "#rho", 100, -DP, DP,100, PhysicsConstants.speedOfLight()*Constants.B/PMAX/1E4, 
                                                                                                       PhysicsConstants.speedOfLight()*Constants.B/0.2/1E4);
                H2F hi_z0     = histo2D("hi_"+type+"z0", titl+"z0 (cm)", "z0 (cm)", 100, -DVZ, DVZ, 100, VZMIN, VZMAX);
                H2F hi_tandip = histo2D("hi_"+type+"tandip", titl+"tandip", "tandip", 100, -DVXY, DVXY, 100, -PMAX, PMAX);
                dg.addDataSet(hi_d0,     0 + i*5);
                dg.addDataSet(hi_phi0,   1 + i*5);
                dg.addDataSet(hi_rho,    2 + i*5);
                dg.addDataSet(hi_z0,     3 + i*5);
                dg.addDataSet(hi_tandip, 4 + i*5);
            }
            else {
                type = "cov1D";
                titl = "cov";
                H1F hi_d0     = histo1D("hi_"+type+"d0", titl+"d0 (cm)", "Counts", 100, -DVXY, DVXY, 46);
                H1F hi_phi0   = histo1D("hi_"+type+"phi0", titl+"#phi (deg)", "Counts", 100, -DPHI, DPHI, 46);
                H1F hi_rho    = histo1D("hi_"+type+"rho", titl+"#rho/#rho", "Counts", 100, -DP, DP, 46);
                H1F hi_z0     = histo1D("hi_"+type+"z0", titl+"z0 (cm)", "Counts", 100, -DVZ, DVZ, 46);
                H1F hi_tandip = histo1D("hi_"+type+"tandip", titl+"tandip", "Counts", 100, -DVXY, DVXY, 46);
                dg.addDataSet(hi_d0,     0 + i*5);
                dg.addDataSet(hi_phi0,   1 + i*5);
                dg.addDataSet(hi_rho,    2 + i*5);
                dg.addDataSet(hi_z0,     3 + i*5);
                dg.addDataSet(hi_tandip, 4 + i*5);                
            }
        }
        return dg;
    }

    private DataGroup rayResolutionGroup(int icol) {
        H1F hi_chi2   = histo1D("hi_chi2", "#chi^2", "Counts", 100, 0, 5,  icol);
        H1F hi_vx     = histo1D("hi_vx", "#Deltavx (cm)", "Counts", 100, -DVZ, DVZ, icol);
        H1F hi_vz     = histo1D("hi_vz", "#Deltavz (cm)", "Counts", 100, -DVZ, DVZ, icol);
        H1F hi_tx     = histo1D("hi_tx", "#Deltatx", "Counts", 100, -DTX, DTX, icol);
        H1F hi_tz     = histo1D("hi_tz", "#Deltatz", "Counts", 100, -DTZ, DTZ, icol);
        DataGroup dg = new DataGroup(3,2);
        dg.addDataSet(hi_chi2,   0);
        dg.addDataSet(hi_vx,     1);
        dg.addDataSet(hi_vz,     2);
        dg.addDataSet(hi_tx,     4);
        dg.addDataSet(hi_tz,     5);
        return dg;
    }

    private DataGroup rayResolution2DGroup() {
        DataGroup dg = new DataGroup(5,3);
        for(int i=0; i<3; i++) {
            String type = "";
            String titl = "#Delta";
            if(i<2) {
                if(i==1) {
                    type = "cov2D";
                    titl = "cov";
                }
                H2F hi_vx     = histo2D("hi_"+type+"vx", titl+"vx (cm)", "vx (cm)", 100, -DVZ, DVZ,100, VZMIN, VZMAX);
                H2F hi_vz     = histo2D("hi_"+type+"vz", titl+"vz (cm)", "vz (cm)", 100, -DVZ, DVZ,100, VZMIN, VZMAX);
                H2F hi_tx     = histo2D("hi_"+type+"tx", titl+"tx", "tx", 100, -DTX, DTX, 100, VZMIN, VZMAX);
                H2F hi_tz     = histo2D("hi_"+type+"tz", titl+"tz", "tz", 100, -DTZ, DTZ, 100, VZMIN, VZMAX);
                dg.addDataSet(hi_vx,     0 + i*5);
                dg.addDataSet(hi_vz,     1 + i*5);
                dg.addDataSet(hi_tx,     2 + i*5);
                dg.addDataSet(hi_tz,     3 + i*5);
            }
            else {
                type = "cov1D";
                titl = "cov";
                H1F hi_vx     = histo1D("hi_"+type+"vx", titl+"vx (cm)", "Counts", 100, -DVZ, DVZ, 46);
                H1F hi_vz     = histo1D("hi_"+type+"vz", titl+"vz (cm)", "Counts", 100, -DVZ, DVZ, 46);
                H1F hi_tx     = histo1D("hi_"+type+"tx", titl+"tx", "Counts", 100, -DTX, DTX, 46);
                H1F hi_tz     = histo1D("hi_"+type+"tz", titl+"tz", "Counts", 100, -DTZ, DTZ, 46);
                dg.addDataSet(hi_vx,     0 + i*5);
                dg.addDataSet(hi_vz,     1 + i*5);
                dg.addDataSet(hi_tx,     2 + i*5);
                dg.addDataSet(hi_tz,     3 + i*5);
            }
        }
        return dg;
    }
    
    private DataGroup repResolutionGroup(int icol) {
        if(!this.isCosmics())
            return this.helixResolutionGroup(icol);
        else 
            return this.rayResolutionGroup(icol);
    }

    private DataGroup repResolution2DGroup() {
        if(!this.isCosmics())
            return this.helixResolution2DGroup();
        else 
            return this.rayResolution2DGroup();
    }

    private DataGroup helixPullsGroup(int icol) {
        H1F hi_chi2   = histo1D("hi_chi2", "#chi^2", "Counts", 100, 0, 5,  icol);
        H1F hi_d0     = histo1D("hi_d0", "d0", "Counts", 100, -5, 5, icol);
        H1F hi_phi0   = histo1D("hi_phi0", "phi0", "Counts", 100, -5, 5, icol);
        H1F hi_rho    = histo1D("hi_rho", "rho", "Counts", 100, -5, 5, icol);
        H1F hi_z0     = histo1D("hi_z0", "z0", "Counts", 100, -5, 5, icol);
        H1F hi_tandip = histo1D("hi_tandip", "tandip", "Counts", 100, -5, 5, icol);

        DataGroup dg = new DataGroup(3,2);
        dg.addDataSet(hi_chi2,   0);
        dg.addDataSet(hi_d0,     1);
        dg.addDataSet(hi_phi0,   2);
        dg.addDataSet(hi_rho,    3);
        dg.addDataSet(hi_z0,     4);
        dg.addDataSet(hi_tandip, 5);
        return dg;
    }

    private DataGroup rayPullsGroup(int icol) {
        H1F hi_chi2   = histo1D("hi_chi2", "#chi^2", "Counts", 100, 0, 5,  icol);
        H1F hi_vx     = histo1D("hi_vx", "vx", "Counts", 100, -5, 5, icol);
        H1F hi_vz     = histo1D("hi_vz", "vz", "Counts", 100, -5, 5, icol);
        H1F hi_tx     = histo1D("hi_tx", "tx", "Counts", 100, -5, 5, icol);
        H1F hi_tz     = histo1D("hi_tz", "tz", "Counts", 100, -5, 5, icol);

        DataGroup dg = new DataGroup(3,2);
        dg.addDataSet(hi_chi2,   0);
        dg.addDataSet(hi_vx,     1);
        dg.addDataSet(hi_vz,     2);
        dg.addDataSet(hi_tx,     4);
        dg.addDataSet(hi_tz,     5);
        return dg;
    }

    private DataGroup pullsGroup(int icol) {
        if(!this.isCosmics())
            return this.helixPullsGroup(icol);
        else
            return this.rayPullsGroup(icol);
    }

    private DataGroup efficiencyGroup() {

        String[] type = {"MC", "Seed", "Rec", "Eff"};
        int[]    cols = { 35 ,    44 ,   46 ,   33 };
        
        DataGroup dg = new DataGroup(type.length,5);

        for(int it=0; it<type.length; it++) {
            String t = type[it];
            H1F hi_p     = histo1D("hi_p_"+t,     t + " p (GeV)",      "Counts", 100, PMIN,     PMAX,     cols[it]);
            H1F hi_theta = histo1D("hi_theta_"+t, t + " #theta (deg)", "Counts", 100, THETAMIN, THETAMAX, cols[it]);
            H1F hi_phi   = histo1D("hi_phi_"+t,   t + " #phi (deg)",   "Counts", 100, PHIMIN,   PHIMAX,   cols[it]);
            H1F hi_pt    = histo1D("hi_pt_"+t,    t + " pt (GeV)",     "Counts", 100, PMIN,     PMAX,     cols[it]);
            H1F hi_vz    = histo1D("hi_vz_"+t,    t + " vz (cm)",      "Counts", 100, VZMIN,    VZMAX,    cols[it]);

            dg.addDataSet(hi_p,     0*type.length + it);
            dg.addDataSet(hi_theta, 1*type.length + it);
            dg.addDataSet(hi_phi,   2*type.length + it);
            dg.addDataSet(hi_pt,    3*type.length + it);
            dg.addDataSet(hi_vz,    4*type.length + it);
        }
        return dg;
    }

    @Override
    public void createHistos() {
        this.getHistos().put("MC", this.trackGroup(45));
        this.getHistos().put("AllSeeds", this.trackGroup(44));
        this.getHistos().put("Seed", this.trackGroup(44));
        this.getHistos().put("SeedResolution1", this.trackResolutionGroup(44));
        this.getHistos().put("SeedResolution2", this.repResolutionGroup(44));
        this.getHistos().put("SeedResolution3", this.repResolution2DGroup());
        this.getHistos().put("SeedPulls", this.pullsGroup(44));
        this.getHistos().put("AllTracks", this.trackGroup(46));
        this.getHistos().put("Track", this.trackGroup(46));
        this.getHistos().put("Resolution1", this.trackResolutionGroup(46));
        this.getHistos().put("Resolution2", this.repResolutionGroup(46));
        this.getHistos().put("Resolution3", this.repResolution2DGroup());
        this.getHistos().put("Resolution4", this.trackResolution2DGroup());
        this.getHistos().put("Pulls", this.pullsGroup(46));
        this.getHistos().put("Efficiency", this.efficiencyGroup());
        this.getHistos().put("Efficiency2", this.efficiencyGroup());
        this.getHistos().put("EfficiencyG", this.efficiencyGroup());
    }
    
    @Override
    public void fillHistos(Event event) {
        Track mcTrack = event.getMCTrack(this.isCosmics());
        if(mcTrack!=null) {
            Track matchedTrack = null;
            for(Track track : event.getTracks()) {
                if(track.match(mcTrack)) {
                    matchedTrack = track;
                    break;
                }
            }            
            this.fillTrackGroup(this.getHistos().get("MC"),mcTrack);
            for(Track seed : event.getSeeds()) this.fillTrackGroup(this.getHistos().get("AllSeeds"),seed);
            for(Track track : event.getTracks()) this.fillTrackGroup(this.getHistos().get("AllTracks"),track);
            this.fillEfficiencyGroup(this.getHistos().get("Efficiency"), mcTrack, "MC");
            if(matchedTrack!=null) {
                int sid = matchedTrack.getSeedId();
                Track seedTrack = matchedTrack;
                if(event.getSeedMap().containsKey(sid)) {
                    int si  = event.getSeedMap().get(sid);
                    seedTrack = event.getSeeds().get(si);
               }
                this.fillTrackGroup(this.getHistos().get("Seed"),seedTrack);
                this.fillTrackGroup(this.getHistos().get("Track"),matchedTrack);
                this.fillEfficiencyGroup(this.getHistos().get("Efficiency"), seedTrack, "Seed");
                this.fillEfficiencyGroup(this.getHistos().get("Efficiency"), matchedTrack, "Rec");
                this.fillEfficiencyGroup(this.getHistos().get("EfficiencyG"), mcTrack, "Rec");
                if(seedTrack.getSeedType()==2) {
                    this.fillEfficiencyGroup(this.getHistos().get("Efficiency2"), seedTrack, "Seed");
                    this.fillEfficiencyGroup(this.getHistos().get("Efficiency2"), matchedTrack, "Rec");
                }
                this.fillTrackResolutionGroup(this.getHistos().get("SeedResolution1"), mcTrack, seedTrack);
                this.fillRepResolutionGroup(this.getHistos().get("SeedResolution2"), mcTrack, seedTrack);
                this.fillRepResolution2DGroup(this.getHistos().get("SeedResolution3"), mcTrack, seedTrack);
                this.fillPullsGroup(this.getHistos().get("SeedPulls"), mcTrack, seedTrack);
                this.fillTrackResolutionGroup(this.getHistos().get("Resolution1"), mcTrack, matchedTrack);
                this.fillRepResolutionGroup(this.getHistos().get("Resolution2"), mcTrack, matchedTrack);
                this.fillRepResolution2DGroup(this.getHistos().get("Resolution3"), mcTrack, matchedTrack);
                this.fillTrackResolution2DGroup(this.getHistos().get("Resolution4"), mcTrack, matchedTrack);
                this.fillPullsGroup(this.getHistos().get("Pulls"), mcTrack, matchedTrack);
            }
        }
    }
    
    private void fillEfficiencyGroup(DataGroup group, Track track, String type) {
        group.getH1F("hi_p_" + type).fill(track.p());
        group.getH1F("hi_pt_" + type).fill(track.pt());
        group.getH1F("hi_theta_" + type).fill(Math.toDegrees(track.theta()));
        group.getH1F("hi_phi_" + type).fill(Math.toDegrees(track.phi()));
        group.getH1F("hi_vz_" + type).fill(track.vz());
    }
    
    private void fillTrackGroup(DataGroup group, Track track) {
        group.getH1F("hi_q").fill(track.charge());
        group.getH1F("hi_p").fill(track.p());
        group.getH1F("hi_pt").fill(track.pt());
        group.getH1F("hi_theta").fill(Math.toDegrees(track.theta()));
        group.getH1F("hi_phi").fill(Math.toDegrees(track.phi()));
        group.getH1F("hi_d0").fill(track.d0());
        group.getH1F("hi_vx").fill(track.vx());
        group.getH1F("hi_vy").fill(track.vy());
        group.getH1F("hi_vz").fill(track.vz());
    }
    
    private void fillTrackResolutionGroup(DataGroup group, Track mc, Track track) {
        group.getH1F("hi_p").fill((mc.p()-track.p())/mc.p());
        group.getH1F("hi_theta").fill(Math.toDegrees(mc.theta()-track.theta()));
        group.getH1F("hi_phi").fill(Math.toDegrees(mc.deltaPhi(track)));
        group.getH1F("hi_vx").fill(mc.vx()-track.vx());
        group.getH1F("hi_vy").fill(mc.vy()-track.vy());
        group.getH1F("hi_vz").fill(mc.vz()-track.vz());
    }
    
    private void fillTrackResolution2DGroup(DataGroup group, Track mc, Track track) {
        group.getH2F("hi_p_dp").fill(mc.p(),(mc.p()-track.p())/mc.p());
        group.getH2F("hi_p_dtheta").fill(mc.p(),Math.toDegrees(mc.theta()-track.theta()));
        group.getH2F("hi_p_dphi").fill(mc.p(),Math.toDegrees(mc.deltaPhi(track)));
        group.getH2F("hi_p_dvz").fill(mc.p(),mc.vz()-track.vz());
        group.getH2F("hi_theta_dp").fill(Math.toDegrees(mc.theta()),(mc.p()-track.p())/mc.p());
        group.getH2F("hi_theta_dtheta").fill(Math.toDegrees(mc.theta()),Math.toDegrees(mc.theta()-track.theta()));
        group.getH2F("hi_theta_dphi").fill(Math.toDegrees(mc.theta()),Math.toDegrees(mc.deltaPhi(track)));
        group.getH2F("hi_theta_dvz").fill(Math.toDegrees(mc.theta()),mc.vz()-track.vz());
        group.getH2F("hi_phi_dp").fill(Math.toDegrees(mc.phi()),(mc.p()-track.p())/mc.p());
        group.getH2F("hi_phi_dtheta").fill(Math.toDegrees(mc.phi()),Math.toDegrees(mc.theta()-track.theta()));
        group.getH2F("hi_phi_dphi").fill(Math.toDegrees(mc.phi()),Math.toDegrees(mc.deltaPhi(track)));
        group.getH2F("hi_phi_dvz").fill(Math.toDegrees(mc.phi()),mc.vz()-track.vz());
    }
    
    private void fillHelixResolutionGroup(DataGroup group, Track mc, Track track) {
        group.getH1F("hi_chi2").fill(track.getChi2()/track.getNDF());
        group.getH1F("hi_d0").fill(mc.d0()-track.d0());
        group.getH1F("hi_phi0").fill(Math.toDegrees(mc.deltaPhi(track)));
        group.getH1F("hi_rho").fill((mc.rho()-track.rho())/mc.rho());
        group.getH1F("hi_z0").fill(mc.vz()-track.vz());
        group.getH1F("hi_tandip").fill(mc.tandip()-track.tandip());
    }
    
    private void fillHelixResolution2DGroup(DataGroup group, Track mc, Track track) {
        group.getH2F("hi_d0").fill(mc.d0()-track.d0(), mc.d0());
        group.getH2F("hi_phi0").fill(Math.toDegrees(mc.deltaPhi(track)), Math.toDegrees(mc.phi()));
        group.getH2F("hi_rho").fill((mc.rho()-track.rho())/mc.rho(), mc.rho());
        group.getH2F("hi_z0").fill(mc.vz()-track.vz(), mc.vz());
        group.getH2F("hi_tandip").fill(mc.tandip()-track.tandip(), mc.tandip());
        group.getH2F("hi_cov2Dd0").fill(track.getD0Err(), mc.d0());
        group.getH2F("hi_cov2Dphi0").fill(Math.toDegrees(track.getPhi0Err()), Math.toDegrees(mc.phi()));
        group.getH2F("hi_cov2Drho").fill(track.getRhoErr()/track.rho(), mc.rho());
        group.getH2F("hi_cov2Dz0").fill(track.getZ0Err(), mc.vz());
        group.getH2F("hi_cov2Dtandip").fill(track.getTanDipErr(), mc.tandip());
        group.getH1F("hi_cov1Dd0").fill(track.getD0Err());
        group.getH1F("hi_cov1Dphi0").fill(Math.toDegrees(track.getPhi0Err()));
        group.getH1F("hi_cov1Drho").fill(track.getRhoErr()/track.rho());
        group.getH1F("hi_cov1Dz0").fill(track.getZ0Err());
        group.getH1F("hi_cov1Dtandip").fill(track.getTanDipErr());
    }
    
    private void fillRayResolutionGroup(DataGroup group, Track mc, Track track) {
        group.getH1F("hi_chi2").fill(track.getChi2()/track.getNDF());
        group.getH1F("hi_vx").fill(mc.vx()-track.vx());
        group.getH1F("hi_vz").fill(mc.vz()-track.vz());
        group.getH1F("hi_tx").fill(mc.tx()-track.tx());
        group.getH1F("hi_tz").fill(mc.tz()-track.tz());
    }
    
    private void fillRayResolution2DGroup(DataGroup group, Track mc, Track track) {
        group.getH2F("hi_vx").fill(mc.vx()-track.vx(), mc.vx());
        group.getH2F("hi_vz").fill(mc.vz()-track.vz(), mc.vz());
        group.getH2F("hi_tx").fill(mc.tx()-track.tx(), mc.tx());
        group.getH2F("hi_tz").fill(mc.tz()-track.tz(), mc.tz());
        group.getH2F("hi_cov2Dvx").fill(track.getVxErr(), mc.vx());
        group.getH2F("hi_cov2Dvz").fill(track.getVzErr(), mc.vz());
        group.getH2F("hi_cov2Dtx").fill(track.getTxErr(), mc.tx());
        group.getH2F("hi_cov2Dtz").fill(track.getTzErr(), mc.tz());
        group.getH1F("hi_cov1Dvx").fill(track.getVxErr());
        group.getH1F("hi_cov1Dvz").fill(track.getVzErr());
        group.getH1F("hi_cov1Dtx").fill(track.getTxErr());
        group.getH1F("hi_cov1Dtz").fill(track.getTzErr());
    }
    
    private void fillRepResolutionGroup(DataGroup group, Track mc, Track track) {
        if(!this.isCosmics())
            this.fillHelixResolutionGroup(group, mc, track);
        else
            this.fillRayResolutionGroup(group, mc, track);
    }

    
    private void fillRepResolution2DGroup(DataGroup group, Track mc, Track track) {
        if(!this.isCosmics())
            this.fillHelixResolution2DGroup(group, mc, track);
        else
            this.fillRayResolution2DGroup(group, mc, track);
    }

    private void fillPullsGroup(DataGroup group, Track mc, Track track) {
        if(!this.isCosmics()) {
            group.getH1F("hi_chi2").fill(track.getChi2()/track.getNDF());
            group.getH1F("hi_d0").fill((mc.d0()-track.d0())/track.getD0Err());
            group.getH1F("hi_phi0").fill(mc.deltaPhi(track)/track.getPhi0Err());
            group.getH1F("hi_rho").fill((mc.rho()-track.rho())/track.getRhoErr());
            group.getH1F("hi_z0").fill((mc.vz()-track.vz())/track.getZ0Err());
            group.getH1F("hi_tandip").fill((mc.tandip()-track.tandip())/track.getTanDipErr());
        }
        else {
            group.getH1F("hi_chi2").fill(track.getChi2()/track.getNDF());
            group.getH1F("hi_vx").fill((mc.vx()-track.vx())/track.getVxErr());
            group.getH1F("hi_vz").fill((mc.vz()-track.vz())/track.getVzErr());
            group.getH1F("hi_tx").fill((mc.tx()-track.tx())/track.getTxErr());
            group.getH1F("hi_tz").fill((mc.tz()-track.tz())/track.getTzErr());
        }
    }
    
    @Override
    public void analyzeHistos() {
        this.fitDataGroup(this.getHistos().get("SeedResolution1"));
        this.fitDataGroup(this.getHistos().get("SeedResolution2"));
        this.fitDataGroup(this.getHistos().get("SeedPulls"));
        this.fitDataGroup(this.getHistos().get("Resolution1"));
        this.fitDataGroup(this.getHistos().get("Resolution2"));
        this.fitDataGroup(this.getHistos().get("Pulls"));
        this.getHistos().get("SeedResolution2").getH1F("hi_chi2").setFunction(null);
        this.getHistos().get("SeedPulls").getH1F("hi_chi2").setFunction(null);
        this.getHistos().get("Resolution2").getH1F("hi_chi2").setFunction(null);
        this.getHistos().get("Pulls").getH1F("hi_chi2").setFunction(null);
        
        DataGroup efficiency  = this.getHistos().get("Efficiency");
        DataGroup efficiencyG = this.getHistos().get("EfficiencyG");
        String[] type = {"p", "theta", "phi", "pt", "vz"};
        for(String t : type) {
            H1F h1 = efficiencyG.getH1F("hi_" + t + "_Rec");
            H1F h2 = efficiency.getH1F("hi_" + t + "_MC");
            H1F h3 = efficiency.getH1F("hi_" + t + "_Eff");
            h3.add(h1);
            h3.divide(h2);
        }
        double eff = 0;
        double evMC  = efficiency.getH1F("hi_phi_MC").getEntries();
        double evRec = efficiency.getH1F("hi_phi_Rec").getEntries();
        if(evMC>0) eff = 100*evRec/evMC;
        System.out.println();
        System.out.println(">>>>> MC/Rec tracks found: " + evMC + "/" + evRec);
        System.out.println(">>>>> Efficiency: " + String.format("%.1f%%", eff));        
    }
    
    @Override
    public void drawHistos() {
        this.addCanvas("MC", "Seed", "SeedResolution1", "SeedResolution2", "SeedResolution3", "SeedPulls", 
                       "Track", "Resolution1", "Resolution2", "Resolution3", "Resolution4", "Pulls", "Efficiency");
        this.getCanvas("MC").draw(this.getHistos().get("MC"));
        this.getCanvas("Seed").draw(this.getHistos().get("AllSeeds"));
        this.getCanvas("Seed").draw(this.getHistos().get("Seed"));
        this.getCanvas("SeedResolution1").draw(this.getHistos().get("SeedResolution1"));
        this.getCanvas("SeedResolution2").draw(this.getHistos().get("SeedResolution2"));
        this.getCanvas("SeedResolution3").draw(this.getHistos().get("SeedResolution3"));
        this.getCanvas("SeedPulls").draw(this.getHistos().get("SeedPulls"));
        this.getCanvas("Track").draw(this.getHistos().get("AllTracks"));
        this.getCanvas("Track").draw(this.getHistos().get("Track"));
        this.getCanvas("Resolution1").draw(this.getHistos().get("Resolution1"));
        this.getCanvas("Resolution2").draw(this.getHistos().get("Resolution2"));
        this.getCanvas("Resolution3").draw(this.getHistos().get("Resolution3"));
        this.getCanvas("Resolution4").draw(this.getHistos().get("Resolution4"));
        this.getCanvas("Pulls").draw(this.getHistos().get("Pulls"));
        this.getCanvas("Efficiency").draw(this.getHistos().get("Efficiency"));
        this.getCanvas("Efficiency").draw(this.getHistos().get("Efficiency2"));
        this.setPlottingOptions("MC");
        this.setPlottingOptions("Seed");
        this.setPlottingOptions("SeedResolution1");
        this.setPlottingOptions("SeedResolution2");
        this.setPlottingOptions("SeedResolution3");
        this.setPlottingOptions("SeedPulls");
        this.setPlottingOptions("Track");
        this.setPlottingOptions("Resolution1");
        this.setPlottingOptions("Resolution2");
        this.setPlottingOptions("Resolution3");
        this.setPlottingOptions("Resolution4");
        this.setPlottingOptions("Pulls");
        this.setPlottingOptions("Efficiency");
    }

    @Override
    public void fitGauss(H1F hi) {
        double mean = hi.getDataX(hi.getMaximumBin());
        double amp  = hi.getBinContent(hi.getMaximumBin());
        double rms  = hi.getRMS();
        double sigma = rms/2;
        String name = hi.getName().split("hi")[1];
        F1D f1 = new F1D("f1" + name,"[amp]*gaus(x,[mean],[sigma])",-0.3, 0.3);
        f1.setLineColor(4);
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
}
