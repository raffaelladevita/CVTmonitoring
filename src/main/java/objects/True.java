package objects;

import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.io.base.DataBank;

/**
 *
 * @author devita
 */
public class True implements Comparable {
 
    private final int id;
    private final int sector;
    private final int layer;
    private final int strip;
    private final DetectorType type;
    
    private double time;
    private double energy;
    private Vector3D momentum;
    private Point3D  position;

    
    public True(int id, int sector, int layer, int strip, DetectorType type) {
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

    public Vector3D getMomentum() {
        return momentum;
    }

    public void setMomentum(Vector3D momentum) {
        this.momentum = momentum;
    }

    public Point3D getPosition() {
        return position;
    }

    public void setPosition(Point3D position) {
        this.position = position;
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
    public static True readTruth(DataBank adc, DataBank mc, int row, int offset, DetectorType type) {
        True hit = new True(row,
                          adc.getByte("sector", row),
                          adc.getByte("layer", row),
                          adc.getInt("component", row),
                          type);
        hit.setEnergy(mc.getFloat("trackE", row+offset));
        hit.setTime(mc.getFloat("avgT", row+offset));
        hit.setMomentum(new Vector3D(mc.getFloat("px", row+offset),
                                     mc.getFloat("py", row+offset),
                                     mc.getFloat("pz", row+offset)).divide(1000));
        hit.setPosition(new Point3D(mc.getFloat("avgX", row+offset),
                                    mc.getFloat("avgY", row+offset),
                                    mc.getFloat("avgZ", row+offset)));
        return hit;
    }

    @Override
    public int compareTo(Object o) {
        True ot = (True) o;
        if(this.getTime()<ot.getTime()) return 1;
        else if(this.getTime()==ot.getTime()) return 0;
        else return -1;
    }
    
}
