package objects;

import analysis.Constants;
import org.jlab.detector.base.DetectorType;
import org.jlab.io.base.DataBank;

/**
 *
 * @author devita
 */
public class Cluster {
 
    private final int id;
    private final int sector;
    private final int layer;
    private final DetectorType type;
    
    private int    size;
    private double energy;
    private double time;
    private double centroid;
    private double centroidValue;
    private double centroidError;
    private double centroidResidual;
    
    private int    seedStrip;
    private double seedEnergy;
    private double seedResidual;
    
    private int trackId;
    
    public Cluster(int id, int sector, int layer, DetectorType type) {
        this.id = id;
        this.sector = sector;
        this.layer = layer;
        this.type = type;
    }

    public int getSector() {
        return sector;
    }

    public int getLayer() {
        return layer;
    }

    public int getSize() {
        return size;
    }

    public void setSize(int size) {
        this.size = size;
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

    public double getCentroid() {
        return centroid;
    }

    public void setCentroid(double centroid) {
        this.centroid = centroid;
    }

    public double getCentroidValue() {
        return centroidValue;
    }

    public void setCentroidValue(double centroidValue) {
        this.centroidValue = centroidValue;
    }

    public double getCentroidError() {
        return centroidError;
    }

    public void setCentroidError(double centroidError) {
        this.centroidError = centroidError;
    }

    public double getCentroidResidual() {
        return centroidResidual;
    }

    public void setCentroidResidual(double centroidResidual) {
        this.centroidResidual = centroidResidual;
    }

    public double getCentroidPhi(boolean twoPI) {
        double phi = 0;
        if(this.getType()==CVTType.BMTZ) {
            int il = this.getLayer()-1;
            int ir = il/2;
            int is = this.getSector();
            phi = Math.toDegrees((this.getCentroid()-Constants.BMTZSTRIPS[ir]/2)*Constants.BMTZPITCH[ir]*1E-4/Constants.BMTRADIUS[il]);
            if(twoPI) phi += Constants.BMTMEANPHI[is];
        }
        return Math.toRadians(phi);
    }
    
    public int getSeedStrip() {
        return seedStrip;
    }

    public void setSeedStrip(int seedStrip) {
        this.seedStrip = seedStrip;
    }

    public double getSeedEnergy() {
        return seedEnergy;
    }

    public void setSeedEnergy(double seedEnergy) {
        this.seedEnergy = seedEnergy;
    }

    public double getSeedResidual() {
        return seedResidual;
    }

    public void setSeedResidual(double seedResidual) {
        this.seedResidual = seedResidual;
    }

    public int getTrackId() {
        return trackId;
    }

    public void setTrackId(int trackId) {
        this.trackId = trackId;
    }
    
    public CVTType getType() {
        if(this.type==DetectorType.BST)
            return CVTType.SVT;
        else {
            if(layer==1 || layer==4 || layer==6)
                return CVTType.BMTC;
            else
                return CVTType.BMTZ;
        }
    }
    public static Cluster readCluster(DataBank bank, int row, DetectorType type) {
        Cluster cluster = new Cluster(bank.getShort("ID", row),
                                      bank.getByte("sector", row),
                                      bank.getByte("layer", row),
                                      type);
        cluster.setSize(bank.getShort("size", row));
        cluster.setEnergy(bank.getFloat("ETot", row));
        cluster.setTime(bank.getFloat("time", row));
        cluster.setCentroid(bank.getFloat("centroid", row));
        cluster.setCentroidResidual(bank.getFloat("centroidResidual", row));
        if(type == DetectorType.BMT) {
            cluster.setCentroidValue(bank.getFloat("centroidValue", row));
        }
        cluster.setCentroidError(bank.getFloat("centroidError", row));
        cluster.setSeedStrip(bank.getInt("seedStrip", row));
        cluster.setSeedEnergy(bank.getFloat("seedE", row));
        cluster.setSeedResidual(bank.getFloat("seedResidual", row));
        
        cluster.setTrackId(bank.getShort("trkID", row));
        return cluster;
    }

}
