package objects;

import org.jlab.detector.base.DetectorType;
import org.jlab.io.base.DataBank;

/**
 *
 * @author devita
 */
public class Hit {
 
    private final int id;
    private final int sector;
    private final int layer;
    private final int strip;
    private final DetectorType type;
    
    private double energy;
    private double time;
    private double residual;
    
    private int clusterId;
    private int trackId;
    
    public Hit(int id, int sector, int layer, int strip, DetectorType type) {
        this.id = id;
        this.sector = sector;
        this.layer = layer;
        this.strip = strip;
        this.type = type;
    }

    public int getSector() {
        return sector;
    }

    public int getLayer() {
        return layer;
    }

    public int getStrip() {
        return strip;
    }

    public DetectorType getType() {
        return type;
    }

    public double getEnergy() {
        return energy;
    }

    public void setEnergy(double energy) {
        this.energy = energy;
    }

    public double getTime() {
        return time;
    }

    public void setTime(double time) {
        this.time = time;
    }

    public double getResidual() {
        return residual;
    }

    public void setResidual(double residual) {
        this.residual = residual;
    }

    public int getClusterId() {
        return clusterId;
    }

    public void setClusterId(int clusterId) {
        this.clusterId = clusterId;
    }
    
    public int getTrackId() {
        return trackId;
    }

    public void setTrackId(int trackId) {
        this.trackId = trackId;
    }
    
    public String getName() {
        if(this.type==DetectorType.BST)
            return "SVT";
        else {
            if(layer==1 || layer==4 || layer==6)
                return "BMTC";
            else
                return "BMTZ";
        }
    }
    public static Hit readHit(DataBank bank, int row, DetectorType type) {
        Hit hit = new Hit(bank.getShort("ID", row),
                          bank.getByte("sector", row),
                          bank.getByte("layer", row),
                          bank.getInt("strip", row),
                          type);
        hit.setEnergy(bank.getFloat("energy", row));
        hit.setTime(bank.getFloat("time", row));
        hit.setResidual(bank.getFloat("fitResidual", row));
        hit.setClusterId(bank.getShort("clusterID", row));
        hit.setTrackId(bank.getShort("trkID", row));
        return hit;
    }
    
}
