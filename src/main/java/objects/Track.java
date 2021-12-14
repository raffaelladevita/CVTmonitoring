package objects;

import analysis.Constants;
import java.util.ArrayList;
import java.util.List;
import org.jlab.clas.pdg.PhysicsConstants;
import org.jlab.clas.physics.Particle;
import org.jlab.detector.base.DetectorType;
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
    private int index = -1;
    private int pindex = -1;
    private double chi2pid = Double.POSITIVE_INFINITY;
    private final double[][] covMatrix = new double[5][5];
    private final double B = 5; //field in Tesla
    
    Track(int pid, double px, double py, double pz, double vx, double vy, double vz) {
        super(pid, px, py, pz, vx, vy, vz);
    }

    Track(int id, int pid, double px, double py, double pz, 
          double vx, double vy, double vz, int type, int NDF, double chi2) {
        super(pid, px, py, pz, vx, vy, vz);
        this.id = id;
        this.seedType = type;
        this.NDF  = NDF;
        this.chi2 = chi2;
    }

    Track(int id, int charge, double pt, double tandip, double phi0, double d0, 
          double x0, double y0, double z0, int type, int NDF, double chi2) {
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

    public double getChi2pid() {
        return chi2pid;
    }

    public void setChi2pid(double chi2pid) {
        this.chi2pid = chi2pid;
    }

    public double pt() {
        return this.p()*Math.cos(this.theta());
    }
    
    public double rho() {
        return PhysicsConstants.speedOfLight()*B/this.pt()/1E5;
    }
   
    public double tandip() {
        return 1/Math.tan(this.theta());
    }
    
    public double d0() {
        return Math.sqrt(this.vx()*this.vx()+this.vy()*this.vy());
    }
    
    public void setCovMatrix(double d02, double d0phi0, double d0rho, double phi02,
                             double phi0rho, double rho2, double z02, double tandip2) {
        this.covMatrix[0][0] = d02;
        this.covMatrix[1][1] = phi02;
        this.covMatrix[2][2] = rho2;
        this.covMatrix[3][3] = z02;
        this.covMatrix[4][4] = tandip2;
        this.covMatrix[0][1] = d0phi0;
        this.covMatrix[0][2] = d0rho;
        this.covMatrix[1][2] = phi0rho;
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
    
    public double deltaPhi(Track o) {
        double dphi = this.phi()-o.phi();
        if(Math.abs(this.phi()-o.phi())>2*Math.PI) dphi -= Math.signum(this.phi()-o.phi())*2*Math.PI;
        return dphi;
    }
    
    public boolean match(Particle p) {
        double dp = (this.p()-p.p())/p.p();
        double dth = Math.toDegrees(this.theta()-p.theta());
        double dph = Math.toDegrees(this.phi()-p.phi());
        if(Math.abs(dph)>2*Math.PI) dph -= Math.signum(dph)*2*Math.PI;
        
        if(this.charge()!=p.charge()) return false;
        else if(Math.abs(dp)>Constants.NSIGMA*Constants.SIGMA_P) return false;
        else if(Math.abs(dth)>Constants.NSIGMA*Constants.SIGMA_THETA) return false;
        else if(Math.abs(dph)>Constants.NSIGMA*Constants.SIGMA_PHI) return false;
        else return true;
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
                            bank.getInt("ndf", row),
                            bank.getFloat("chi2", row));
        t.setCovMatrix(bank.getFloat("cov_d02", row),
                       bank.getFloat("cov_d0phi0", row),
                       bank.getFloat("cov_d0rho", row),
                       bank.getFloat("cov_phi02", row),
                       bank.getFloat("cov_phi0rho", row),
                       bank.getFloat("cov_rho2", row),
                       bank.getFloat("cov_z02", row),
                       bank.getFloat("cov_tandip2", row));
        t.setIndex(row);
        t.setSeedId(bank.getShort("seedID", row));
        return t;
    }
    
    public void addEBinfo(DataBank part, DataBank track) {
        for(int i=0; i<track.rows(); i++) {
            if(this.getIndex() == track.getShort("index", i) && 
               track.getByte("detector", i)==DetectorType.CVT.getDetectorId()) {
                this.setPindex(track.getShort("pindex", i));
                break;
            }
        }
        double chi2pid = part.getFloat("chi2pid", this.getPindex());
        int    status  = part.getShort(("status"), this.getPindex());
        int det = status/1000;
        int sci = (status-det*1000)/100;
        if(det==DetectorType.CTOF.getDetectorId() && sci>0) this.setChi2pid(chi2pid);
    }
    
    public static List<Track> readTracks(DataBank bank) {
        List<Track> tracks = new ArrayList<>();
        for(int i=0; i<bank.rows(); i++) {
            tracks.add(Track.readTrack(bank, i));
        }
        return tracks;
    }

    public static Track readSeed(DataBank bank, int row) {
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
                            bank.getFloat("chi2", row));
        t.setCovMatrix(bank.getFloat("cov_d02", row),
                       bank.getFloat("cov_d0phi0", row),
                       bank.getFloat("cov_d0rho", row),
                       bank.getFloat("cov_phi02", row),
                       bank.getFloat("cov_phi0rho", row),
                       bank.getFloat("cov_rho2", row),
                       bank.getFloat("cov_z02", row),
                       bank.getFloat("cov_tandip2", row));
        return t;
    }
    
    public static List<Track> readSeeds(DataBank bank) {
        List<Track> tracks = new ArrayList<>();
        for(int i=0; i<bank.rows(); i++) {
            tracks.add(Track.readSeed(bank, i));
        }
        return tracks;
    }
}
