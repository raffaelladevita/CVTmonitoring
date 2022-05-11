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
    
    private final double PMIN = 0.3;
    private final double PMAX = 1.7;
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

    private DataGroup trackTrack2DGroup() {
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

    private DataGroup helixTrack2DGroup(int icol) {
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
                H1F hi_d0     = histo1D("hi_"+type+"d0", titl+"d0 (cm)", "Counts", 100, -DVXY, DVXY, icol);
                H1F hi_phi0   = histo1D("hi_"+type+"phi0", titl+"#phi (deg)", "Counts", 100, -DPHI, DPHI, icol);
                H1F hi_rho    = histo1D("hi_"+type+"rho", titl+"#rho/#rho", "Counts", 100, -DP, DP, icol);
                H1F hi_z0     = histo1D("hi_"+type+"z0", titl+"z0 (cm)", "Counts", 100, -DVZ, DVZ, icol);
                H1F hi_tandip = histo1D("hi_"+type+"tandip", titl+"tandip", "Counts", 100, -DVXY, DVXY, icol);
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

    private DataGroup rayTrack2DGroup(int icol) {
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
                H1F hi_vx     = histo1D("hi_"+type+"vx", titl+"vx (cm)", "Counts", 100, -DVZ, DVZ, icol);
                H1F hi_vz     = histo1D("hi_"+type+"vz", titl+"vz (cm)", "Counts", 100, -DVZ, DVZ, icol);
                H1F hi_tx     = histo1D("hi_"+type+"tx", titl+"tx", "Counts", 100, -DTX, DTX, icol);
                H1F hi_tz     = histo1D("hi_"+type+"tz", titl+"tz", "Counts", 100, -DTZ, DTZ, icol);
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

    private DataGroup repTrack2DGroup(int icol) {
        if(!this.isCosmics())
            return this.helixTrack2DGroup(icol);
        else 
            return this.rayTrack2DGroup(icol);
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

    private DataGroup trackPullsGroup(int icol) {
        H1F hi_x  = histo1D("hi_x", "x", "Counts", 100, -5, 5,  icol);
        H1F hi_y  = histo1D("hi_y", "y", "Counts", 100, -5, 5, icol);
        H1F hi_z  = histo1D("hi_z", "z", "Counts", 100, -5, 5, icol);
        H1F hi_px = histo1D("hi_px", "px", "Counts", 100, -5, 5, icol);
        H1F hi_py = histo1D("hi_py", "py", "Counts", 100, -5, 5, icol);
        H1F hi_pz = histo1D("hi_pz", "pz", "Counts", 100, -5, 5, icol);

        DataGroup dg = new DataGroup(3,2);
        dg.addDataSet(hi_x,  0);
        dg.addDataSet(hi_y,  1);
        dg.addDataSet(hi_z,  2);
        dg.addDataSet(hi_px, 3);
        dg.addDataSet(hi_py, 4);
        dg.addDataSet(hi_pz, 5);
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
        this.getHistos().put("Seed",  this.trackGroup(44));
        this.getHistos().put("Seed1", this.trackResolutionGroup(44));
        this.getHistos().put("Seed2", this.repResolutionGroup(44));
        this.getHistos().put("Seed3", this.repTrack2DGroup(44));
        this.getHistos().put("Seed4", this.trackTrack2DGroup());
        this.getHistos().put("SeedPulls", this.pullsGroup(44));
        this.getHistos().put("AllFPTracks", this.trackGroup(49));
        this.getHistos().put("FPass",  this.trackGroup(49));
        this.getHistos().put("FPass1", this.trackResolutionGroup(49));
        this.getHistos().put("FPass2", this.repResolutionGroup(49));
        this.getHistos().put("FPass3", this.repTrack2DGroup(49));
        this.getHistos().put("FPass4", this.trackTrack2DGroup());
        this.getHistos().put("FPPulls",  this.pullsGroup(49));
        this.getHistos().put("AllUTracks", this.trackGroup(47));
        this.getHistos().put("UTrack",  this.trackGroup(47));
        this.getHistos().put("UTrack1", this.trackResolutionGroup(47));
        this.getHistos().put("UTrack2", this.repResolutionGroup(47));
        this.getHistos().put("UTrack3", this.repTrack2DGroup(47));
        this.getHistos().put("UTrack4", this.trackTrack2DGroup());
        this.getHistos().put("UPulls",  this.pullsGroup(47));
        this.getHistos().put("AllTracks", this.trackGroup(46));
        this.getHistos().put("Track",  this.trackGroup(46));
        this.getHistos().put("Track1", this.trackResolutionGroup(46));
        this.getHistos().put("Track2", this.repResolutionGroup(46));
        this.getHistos().put("Track3", this.repTrack2DGroup(46));
        this.getHistos().put("Track4", this.trackTrack2DGroup());
        this.getHistos().put("Pulls",  this.pullsGroup(46));
        this.getHistos().put("PullsXYZ",  this.trackPullsGroup(46));
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
            for(Track track : event.getFPTracks()) this.fillTrackGroup(this.getHistos().get("AllFPTracks"),track);
            for(Track track : event.getUTracks()) this.fillTrackGroup(this.getHistos().get("AllUTracks"),track);
            for(Track track : event.getTracks()) this.fillTrackGroup(this.getHistos().get("AllTracks"),track);
            this.fillEfficiencyGroup(this.getHistos().get("Efficiency"), mcTrack, "MC");
            if(matchedTrack!=null) {
                int tid = matchedTrack.getId();
                int sid = matchedTrack.getSeedId();
                Track uTrack    = matchedTrack;
                Track fpTrack   = matchedTrack;
                Track seedTrack = matchedTrack;
                if(event.getUTrackMap().containsKey(tid)) {
                    int ti  = event.getUTrackMap().get(tid);
                    uTrack = event.getUTracks().get(ti);
                }
                if(event.getSeedMap().containsKey(sid)) {
                    int si  = event.getSeedMap().get(sid);
                    seedTrack = event.getSeeds().get(si);
                }
                if(event.getFPTrackMap().containsKey(tid)) {
                    int fi = event.getFPTrackMap().get(tid);
                    fpTrack = event.getFPTracks().get(fi);
                }
                this.fillTrackGroup(this.getHistos().get("Seed"),seedTrack);
                this.fillTrackGroup(this.getHistos().get("FPass"),fpTrack);
                this.fillTrackGroup(this.getHistos().get("UTrack"),uTrack);
                this.fillTrackGroup(this.getHistos().get("Track"),matchedTrack);
                this.fillEfficiencyGroup(this.getHistos().get("Efficiency"), seedTrack, "Seed");
                this.fillEfficiencyGroup(this.getHistos().get("Efficiency"), matchedTrack, "Rec");
                this.fillEfficiencyGroup(this.getHistos().get("EfficiencyG"), mcTrack, "Rec");
                if(seedTrack.getSeedType()==2) {
                    this.fillEfficiencyGroup(this.getHistos().get("Efficiency2"), seedTrack, "Seed");
                    this.fillEfficiencyGroup(this.getHistos().get("Efficiency2"), matchedTrack, "Rec");
                }
                this.fillTrackResolutionGroup(this.getHistos().get("Seed1"), mcTrack, seedTrack);
                this.fillRepResolutionGroup(this.getHistos().get("Seed2"), mcTrack, seedTrack);
                this.fillRepTrack2DGroup(this.getHistos().get("Seed3"), mcTrack, seedTrack);
                this.fillTrackTrack2DGroup(this.getHistos().get("Seed4"), mcTrack, seedTrack);
                this.fillPullsGroup(this.getHistos().get("SeedPulls"), mcTrack, seedTrack);
                this.fillTrackResolutionGroup(this.getHistos().get("FPass1"), mcTrack, fpTrack);
                this.fillRepResolutionGroup(this.getHistos().get("FPass2"), mcTrack, fpTrack);
                this.fillRepTrack2DGroup(this.getHistos().get("FPass3"), mcTrack, fpTrack);
                this.fillTrackTrack2DGroup(this.getHistos().get("FPass4"), mcTrack, fpTrack);
                this.fillPullsGroup(this.getHistos().get("FPPulls"), mcTrack, fpTrack);
                this.fillTrackResolutionGroup(this.getHistos().get("UTrack1"), mcTrack, uTrack);
                this.fillRepResolutionGroup(this.getHistos().get("UTrack2"), mcTrack, uTrack);
                this.fillRepTrack2DGroup(this.getHistos().get("UTrack3"), mcTrack, uTrack);
                this.fillTrackTrack2DGroup(this.getHistos().get("UTrack4"), mcTrack, uTrack);
                this.fillPullsGroup(this.getHistos().get("UPulls"), mcTrack, uTrack);
                this.fillTrackResolutionGroup(this.getHistos().get("Track1"), mcTrack, matchedTrack);
                this.fillRepResolutionGroup(this.getHistos().get("Track2"), mcTrack, matchedTrack);
                this.fillRepTrack2DGroup(this.getHistos().get("Track3"), mcTrack, matchedTrack);
                this.fillTrackTrack2DGroup(this.getHistos().get("Track4"), mcTrack, matchedTrack);
                this.fillPullsGroup(this.getHistos().get("Pulls"), mcTrack, matchedTrack);
                this.fillPullsXYZGroup(this.getHistos().get("PullsXYZ"), mcTrack, matchedTrack);
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
    
    private void fillTrackTrack2DGroup(DataGroup group, Track mc, Track track) {
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
    
    private void fillHelixTrack2DGroup(DataGroup group, Track mc, Track track) {
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
    
    private void fillRayTrack2DGroup(DataGroup group, Track mc, Track track) {
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

    
    private void fillRepTrack2DGroup(DataGroup group, Track mc, Track track) {
        if(!this.isCosmics())
            this.fillHelixTrack2DGroup(group, mc, track);
        else
            this.fillRayTrack2DGroup(group, mc, track);
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

    private void fillPullsXYZGroup(DataGroup group, Track mc, Track track) {
        if(!this.isCosmics()) {
            group.getH1F("hi_x").fill((mc.vx()-track.vx())/track.getXErr());
            group.getH1F("hi_y").fill((mc.vy()-track.vy())/track.getYErr());
            group.getH1F("hi_z").fill((mc.vz()-track.vz())/track.getZErr());
            group.getH1F("hi_px").fill((mc.px()-track.px())/track.getPxErr());
            group.getH1F("hi_py").fill((mc.py()-track.py())/track.getPyErr());
            group.getH1F("hi_pz").fill((mc.pz()-track.pz())/track.getPzErr());            
        }
    }
    
    @Override
    public void analyzeHistos() {
        this.fitDataGroup(this.getHistos().get("Seed1"));
        this.fitDataGroup(this.getHistos().get("Seed2"));
        this.fitDataGroup(this.getHistos().get("SeedPulls"));
        this.fitDataGroup(this.getHistos().get("FPass1"));
        this.fitDataGroup(this.getHistos().get("FPass2"));
        this.fitDataGroup(this.getHistos().get("FPPulls"));
        this.fitDataGroup(this.getHistos().get("UTrack1"));
        this.fitDataGroup(this.getHistos().get("UTrack2"));
        this.fitDataGroup(this.getHistos().get("UPulls"));
        this.fitDataGroup(this.getHistos().get("Track1"));
        this.fitDataGroup(this.getHistos().get("Track2"));
        this.fitDataGroup(this.getHistos().get("Pulls"));
        this.fitDataGroup(this.getHistos().get("PullsXYZ"));
        this.getHistos().get("Seed2").getH1F("hi_chi2").setFunction(null);
        this.getHistos().get("SeedPulls").getH1F("hi_chi2").setFunction(null);
        this.getHistos().get("FPass2").getH1F("hi_chi2").setFunction(null);
        this.getHistos().get("FPPulls").getH1F("hi_chi2").setFunction(null);
        this.getHistos().get("UTrack2").getH1F("hi_chi2").setFunction(null);
        this.getHistos().get("UPulls").getH1F("hi_chi2").setFunction(null);
        this.getHistos().get("Track2").getH1F("hi_chi2").setFunction(null);
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
        this.addCanvas("MC", "S", "S1", "S2", "S3", "S4", "SPulls", 
                       "FP", "FP1", "FP2", "FP3", "FP4", "FPPulls",
                       "UT", "UT1", "UT2", "UT3", "UT4", "UPulls",
                       "T", "T1", "T2", "T3", "T4", "Pulls", "PullsXYZ", "Efficiency");
        this.getCanvas("MC").draw(this.getHistos().get("MC"));
        this.getCanvas("S").draw(this.getHistos().get("AllSeeds"));
        this.getCanvas("S").draw(this.getHistos().get("Seed"));
        this.getCanvas("S1").draw(this.getHistos().get("Seed1"));
        this.getCanvas("S2").draw(this.getHistos().get("Seed2"));
        this.getCanvas("S3").draw(this.getHistos().get("Seed3"));
        this.getCanvas("S4").draw(this.getHistos().get("Seed4"));
        this.getCanvas("SPulls").draw(this.getHistos().get("SeedPulls"));
        this.getCanvas("FP").draw(this.getHistos().get("AllFPTracks"));
        this.getCanvas("FP").draw(this.getHistos().get("FPass"));
        this.getCanvas("FP1").draw(this.getHistos().get("FPass1"));
        this.getCanvas("FP2").draw(this.getHistos().get("FPass2"));
        this.getCanvas("FP3").draw(this.getHistos().get("FPass3"));
        this.getCanvas("FP4").draw(this.getHistos().get("FPass4"));
        this.getCanvas("FPPulls").draw(this.getHistos().get("FPPulls"));
        this.getCanvas("UT").draw(this.getHistos().get("AllUTracks"));
        this.getCanvas("UT").draw(this.getHistos().get("UTrack"));
        this.getCanvas("UT1").draw(this.getHistos().get("UTrack1"));
        this.getCanvas("UT2").draw(this.getHistos().get("UTrack2"));
        this.getCanvas("UT3").draw(this.getHistos().get("UTrack3"));
        this.getCanvas("UT4").draw(this.getHistos().get("UTrack4"));
        this.getCanvas("UPulls").draw(this.getHistos().get("UPulls"));
        this.getCanvas("T").draw(this.getHistos().get("AllTracks"));
        this.getCanvas("T").draw(this.getHistos().get("Track"));
        this.getCanvas("T1").draw(this.getHistos().get("Track1"));
        this.getCanvas("T2").draw(this.getHistos().get("Track2"));
        this.getCanvas("T3").draw(this.getHistos().get("Track3"));
        this.getCanvas("T4").draw(this.getHistos().get("Track4"));
        this.getCanvas("Pulls").draw(this.getHistos().get("Pulls"));
        this.getCanvas("PullsXYZ").draw(this.getHistos().get("PullsXYZ"));
        this.getCanvas("Efficiency").draw(this.getHistos().get("Efficiency"));
        this.getCanvas("Efficiency").draw(this.getHistos().get("Efficiency2"));
        this.setPlottingOptions("MC");
        this.setPlottingOptions("S");
        this.setPlottingOptions("S1");
        this.setPlottingOptions("S2");
        this.setPlottingOptions("S3");
        this.setPlottingOptions("SPulls");
        this.setPlottingOptions("FP");
        this.setPlottingOptions("FP1");
        this.setPlottingOptions("FP2");
        this.setPlottingOptions("FP3");
        this.setPlottingOptions("FP4");
        this.setPlottingOptions("FPPulls");
        this.setPlottingOptions("UT");
        this.setPlottingOptions("UT1");
        this.setPlottingOptions("UT2");
        this.setPlottingOptions("UT3");
        this.setPlottingOptions("UT4");
        this.setPlottingOptions("UPulls");
        this.setPlottingOptions("T");
        this.setPlottingOptions("T1");
        this.setPlottingOptions("T2");
        this.setPlottingOptions("T3");
        this.setPlottingOptions("T4");
        this.setPlottingOptions("Pulls");
        this.setPlottingOptions("PullsXYZ");
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
