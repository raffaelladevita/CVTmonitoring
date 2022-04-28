package objects;

import analysis.Constants;
import java.util.ArrayList;
import java.util.List;
import org.jlab.clas.pdg.PhysicsConstants;
import org.jlab.clas.physics.Particle;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.io.base.DataBank;

/**
 *
 * @author devita
 */
public class Track extends Particle {

    private int id;
    private int seedId;
    private int seedType;
    private int NDF;
    private double chi2;
    private double xb;
    private double yb;
    private int index = -1;
    private int pindex = -1;
    private double solenoid = -1;
    private double beta = -1;
    private double chi2pid = Double.POSITIVE_INFINITY;
    private int recStatus=0;
    private int sector=0;
    private int status=0;
    private int type=0;
    private final double[][] covMatrix = new double[5][5];
    
    Track(int pid, double px, double py, double pz, double vx, double vy, double vz) {
        super(pid, px, py, pz, vx, vy, vz);
    }

    Track(int id, int pid, double px, double py, double pz,  double vx, double vy, double vz, 
          int type, int NDF, double chi2, int status) {
        super(pid, px, py, pz, vx, vy, vz);
        this.id = id;
        this.seedType = type;
        this.NDF  = NDF;
        this.chi2 = chi2;
        this.status = status;
    }

    Track(int id, double theta, double phi, int NDF, double chi2, int status) {
        super(13, 
              Math.cos(phi)*Math.sin(theta),
              Math.sin(phi)*Math.sin(theta),
              Math.cos(theta),
              0,0,0);
        this.id = id;
        this.seedType = 1;
        this.NDF  = NDF;
        this.chi2 = chi2;
        this.status = status;
        this.type = 1;
    }

    Track(int id, double x0, double z0, double tx, double tz, int NDF, double chi2, int status) {
        super(13, 
             -tx/Math.sqrt(1+tx*tx+tz*tz),
              -1/Math.sqrt(1+tx*tx+tz*tz),
             -tz/Math.sqrt(1+tx*tx+tz*tz),
              x0, 0, z0);
        this.id = id;
        this.seedType = 1;
        this.NDF  = NDF;
        this.chi2 = chi2;
        this.status = status;
        this.type = 1;
    }

    Track(int id, int charge, double pt, double tandip, double phi0, double d0, 
          double x0, double y0, double z0, int type, int NDF, double chi2, int status) {
        super(211*charge, 
              pt*Math.cos(phi0),
              pt*Math.sin(phi0),
              pt*tandip,
             -d0 * Math.sin(phi0) + x0,
              d0 * Math.cos(phi0) + y0,
              z0);
        this.id = id;
        this.seedType = type;
        this.NDF  = NDF;
        this.chi2 = chi2;
        this.xb = x0;
        this.yb = y0;
        this.status = status;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public int getPindex() {
        return pindex;
    }

    public void setPindex(int pindex) {
        this.pindex = pindex;
    }

    public int getSeedId() {
        return seedId;
    }

    public void setSeedId(int seedId) {
        this.seedId = seedId;
    }

    public int getSeedType() {
        return seedType;
    }

    public void setSeedType(int seedType) {
        this.seedType = seedType;
    }

    public int getNDF() {
        return NDF;
    }

    public void setNDF(int NDF) {
        this.NDF = NDF;
    }

    public double getChi2() {
        return chi2;
    }

    public void setChi2(double chi2) {
        this.chi2 = chi2;
    }

    public int getStatus() {
        return status;
    }

    public void setStatus(int status) {
        this.status = status;
    }

    public int getKFIterations() {
        return (int) status/1000;
    }

    public int getType() {
        return type;
    }

    public void setType(int type) {
        this.type = type;
    }

    public double getBeta() {
        return beta;
    }

    public void setBeta(double beta) {
        this.beta = beta;
    }

    public double getChi2pid() {
        return chi2pid;
    }

    public void setChi2pid(double chi2pid) {
        this.chi2pid = chi2pid;
    }

    public int getDetector() {
        return (int) Math.abs(this.recStatus)/1000;
    }

    public int getRECStatus() {
        return recStatus;
    }

    public void setRECStatus(int status) {
        this.recStatus = status;
    }

    public int getSector() {
        return sector;
    }

    public void setSector(int sector) {
        this.sector = sector;
    }

    public double pt() {
        return this.p()*Math.sin(this.theta());
    }
    
    public double rho() {
        return PhysicsConstants.speedOfLight()*Constants.B/this.pt()/1E5*10;
    }
   
    public double tandip() {
        return 1/Math.tan(this.theta());
    }
    
    public double d0() {
        double kappa = -this.charge()/this.pt();
        double xcen = this.vx() + Math.signum(kappa) * Constants.ALPHA * this.py();
        double ycen = this.vy() - Math.signum(kappa) * Constants.ALPHA * this.px();
        double phi0 = Math.atan2(ycen, xcen);
        if (Math.signum(kappa) < 0) {
            phi0 = Math.atan2(-ycen, -xcen);
        }
        double drh0 = (xcen-Constants.getXBeam()-xb)*Math.cos(phi0) + (ycen-Constants.getYBeam()-yb)*Math.sin(phi0) - Constants.ALPHA/ kappa;
        return -drh0;
//        return Math.signum(this.vy()/Math.cos(this.phi()))*Math.sqrt(this.vx()*this.vx()+this.vy()*this.vy());
    }
    
    public double tx() {
        return this.px()/this.py();
    }
    
    public double tz() {
        return this.pz()/this.py();
    }
    
    public void setCovMatrix(double c00, double c01, double c02, double c03, double c04, 
                             double c11, double c12, double c13, double c14, 
                             double c22, double c23, double c24,
                             double c33, double c34, 
                             double c44) {
        this.covMatrix[0][0] = c00;
        this.covMatrix[0][1] = c01;
        this.covMatrix[0][2] = c02;
        this.covMatrix[0][3] = c03;
        this.covMatrix[0][4] = c04;
        this.covMatrix[1][1] = c11;
        this.covMatrix[1][2] = c12;
        this.covMatrix[1][3] = c13;
        this.covMatrix[1][4] = c14;
        this.covMatrix[2][2] = c22;
        this.covMatrix[2][3] = c23;
        this.covMatrix[2][4] = c24;
        this.covMatrix[3][3] = c33;
        this.covMatrix[3][4] = c34;
        this.covMatrix[4][4] = c44;
    }

    public double[][] getCovMatrix() {
        return covMatrix;
    }
    
    public double getD0Err() {
        return Math.sqrt(this.covMatrix[0][0]);
    }
    
    public double getPhi0Err() {
        return Math.sqrt(this.covMatrix[1][1]);
    }
    
    public double getRhoErr() {
        return Math.sqrt(this.covMatrix[2][2]);
    }
    
    public double getZ0Err() {
        return Math.sqrt(this.covMatrix[3][3]);
    }
    
    public double getTanDipErr() {
        return Math.sqrt(this.covMatrix[4][4]);
    }
    
    public double getVxErr() {
        return Math.sqrt(this.covMatrix[0][0]);
    }
    
    public double getVzErr() {
        return Math.sqrt(this.covMatrix[1][1]);
    }
    
    public double getTxErr() {
        return Math.sqrt(this.covMatrix[2][2]);
    }
    
    public double getTzErr() {
        return Math.sqrt(this.covMatrix[3][3]);
    }
    
    public double deltaPhi(Track o) {
        double dphi = this.phi()-o.phi();
        if(Math.abs(this.phi()-o.phi())>2*Math.PI) dphi -= Math.signum(this.phi()-o.phi())*2*Math.PI;
        return dphi;
    }
    
    public boolean match(Particle p) {
        double dp = (this.p()-p.p())/p.p();
        double dth = Math.toDegrees(this.theta()-p.theta());
        double dph = Math.toDegrees(this.phi()-p.phi());
        if(Math.abs(dph)>Math.PI) dph -= Math.signum(dph)*2*Math.PI;
        
        if(this.solenoid!=0 && this.charge()!=p.charge()) return false;
        else if(this.solenoid!=0 && Math.abs(dp)>Constants.NSIGMA*Constants.SIGMA_P) return false;
        else if(Math.abs(dth)>Constants.NSIGMA*Constants.SIGMA_THETA) return false;
        else if(Math.abs(dph)>Constants.NSIGMA*Constants.SIGMA_PHI) return false;
        else return true;
    }
    
    public void toCosmic() {
        Vector3D dir = new Vector3D(this.px(), this.py(), this.pz()).asUnit();
        Point3D  vtx = new Point3D(this.vx(), this.vy(), this.vz());
        Line3D  line = new Line3D(vtx, dir);
        Point3D vtxn = line.lerpPoint(-vtx.y()/dir.y());
        this.setVector(this.pid(), this.px(), this.py(), this.pz(), vtxn.x(), vtxn.y(), vtxn.z()); 
        this.setType(1);
    }
    
    public static Track readTrack(DataBank bank, int row) {
        Track t = new Track(bank.getShort("ID", row),
                            bank.getByte("q", row),
                            bank.getFloat("pt", row),
                            bank.getFloat("tandip", row),
                            bank.getFloat("phi0", row),
                            bank.getFloat("d0", row),
                            bank.getFloat("xb", row),
                            bank.getFloat("yb", row),
                            bank.getFloat("z0", row),
                            bank.getByte("fittingMethod", row),
                            bank.getShort("ndf", row),
                            bank.getFloat("chi2", row),
                            bank.getShort("status", row));
        t.setCovMatrix(bank.getFloat("cov_d02", row),
                       bank.getFloat("cov_d0phi0", row),
                       bank.getFloat("cov_d0rho", row),
                       0, 0,
                       bank.getFloat("cov_phi02", row),
                       bank.getFloat("cov_phi0rho", row),
                       0, 0,
                       bank.getFloat("cov_rho2", row),
                       0, 0,
                       bank.getFloat("cov_z02", row),
                       0,
                       bank.getFloat("cov_tandip2", row));
        t.setIndex(row);
        t.setSeedId(bank.getShort("seedID", row));
        return t;
    }
    
    public static Track readRay(DataBank bank, int row) {
        Track t = new Track(bank.getShort("ID", row),
                            bank.getFloat("trkline_yx_interc", row),
                            bank.getFloat("trkline_yz_interc", row),
                            bank.getFloat("trkline_yx_slope", row),
                            bank.getFloat("trkline_yz_slope", row),
                            bank.getInt("ndf", row),
                            bank.getFloat("chi2", row),
                            bank.getShort("status", row));
        t.setCovMatrix(bank.getFloat("cov_x02",  row),
                       bank.getFloat("cov_x0z0", row),
                       bank.getFloat("cov_x0tx", row),
                       bank.getFloat("cov_x0tz", row),
                       0,
                       bank.getFloat("cov_z02",  row),
                       bank.getFloat("cov_z0tx", row),
                       bank.getFloat("cov_z0tz", row),
                       0,
                       bank.getFloat("cov_tx2",  row),
                       bank.getFloat("cov_txtz", row),
                       0,
                       bank.getFloat("cov_tz2",  row),
                       0, 0);
        t.setIndex(row);
        t.setSeedId(t.getId());
        return t;
    }
    
    public void addScale(DataBank config) {
        this.solenoid = config.getFloat("solenoid", 0);
    }
    
    public void addEBinfo(DataBank part, DataBank track) {
        for(int i=0; i<track.rows(); i++) {
            if(this.getIndex() == track.getShort("index", i) && 
               track.getByte("detector", i)==DetectorType.CVT.getDetectorId()) {
                this.setPindex(track.getShort("pindex", i));
                break;
            }
        }
        this.setChi2pid(part.getFloat("chi2pid", this.getPindex()));
        this.setRECStatus(part.getShort("status", this.getPindex()));
    }
    
    public static List<Track> readTracks(DataBank bank) {
        List<Track> tracks = new ArrayList<>();
        for(int i=0; i<bank.rows(); i++) {
            tracks.add(Track.readTrack(bank, i));
        }
        return tracks;
    }

    public static Track readParticle(DataBank recPart, DataBank recTrack, int row) {
        int pid    = recPart.getInt("pid", row);
        int charge = recPart.getByte("charge", row);
        if(pid==0) {
            pid = charge==0 ? 22 : charge*211;
        }
        Track t = new Track(pid,
                    recPart.getFloat("px", row),
                    recPart.getFloat("py", row),
                    recPart.getFloat("pz", row),
                    recPart.getFloat("vx", row)+0*Constants.getXBeam(),
                    recPart.getFloat("vy", row)+0*Constants.getYBeam(),
                    recPart.getFloat("vz", row));
        t.setBeta(recPart.getFloat("beta", row));
        t.setChi2pid(recPart.getFloat("chi2pid", row));
        t.setRECStatus(recPart.getShort("status", row));
        if(recTrack!=null) {
            for(int j=0; j<recTrack.rows(); j++) {
                if(recTrack.getShort("pindex", j)==row) {
                    t.setIndex(recTrack.getShort("index", j));
                    t.setSector(recTrack.getByte("sector", j));
                    t.setNDF(recTrack.getShort("NDF", j));
                    t.setChi2(recTrack.getFloat("chi2", j));
                    t.setStatus(recTrack.getShort("status", j));
                    break;
                }
            }
        }       
        return t;
    }
    
    public static Track readSeed(DataBank bank, int row) {
        if(bank.getByte("q", row)==0) bank.show();
        Track t = new Track(bank.getShort("ID", row),
                            bank.getByte("q", row),
                            bank.getFloat("pt", row),
                            bank.getFloat("tandip", row),
                            bank.getFloat("phi0", row),
                            bank.getFloat("d0", row),
                            bank.getFloat("xb", row),
                            bank.getFloat("yb", row),
                            bank.getFloat("z0", row),
                            bank.getByte("fittingMethod", row),
                            bank.getInt("ndf", row),
                            bank.getFloat("chi2", row), 
                            0);
        t.setCovMatrix(bank.getFloat("cov_d02", row),
                       bank.getFloat("cov_d0phi0", row),
                       bank.getFloat("cov_d0rho", row),
                       0, 0,
                       bank.getFloat("cov_phi02", row),
                       bank.getFloat("cov_phi0rho", row),
                       0, 0,
                       bank.getFloat("cov_rho2", row),
                       0, 0,
                       bank.getFloat("cov_z02", row),
                       0,
                       bank.getFloat("cov_tandip2", row));
        return t;
    }
}
