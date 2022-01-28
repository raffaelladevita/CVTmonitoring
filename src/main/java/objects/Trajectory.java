package objects;

import org.jlab.detector.base.DetectorDescriptor;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Point3D;
import org.jlab.io.base.DataBank;

/**
 *
 * @author devita
 */
public class Trajectory {

    private final int trackId;
    private final int sector;
    private final int layer;
    private final DetectorType detector;
    private double x;
    private double y;
    private double z;
    private double theta;
    private double phi;
    private double localAngle;
    private double path;

    
    Trajectory(int id, int type, int sector, int layer) {
        this.trackId  = id;
        this.detector = DetectorType.getType(type);
        this.sector   = sector;
        this.layer    = layer;
    }

    public DetectorDescriptor getDescriptor() {
       DetectorDescriptor desc = new DetectorDescriptor(this.detector); 
       desc.setSectorLayerComponent(sector, layer, 0);
       return desc;
    }
    
    public int getTrackId() {
        return trackId;
    }

    public int getSector() {
        return sector;
    }

    public int getLayer() {
        return layer;
    }

    public DetectorType getDetector() {
        return detector;
    }

    public double x() {
        return x;
    }

    public double y() {
        return y;
    }

    public double z() {
        return z;
    }

    public void setPosition(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
    
    public Point3D getPoint() {
        return new Point3D(x, y, z);
    }
    
    public double theta() {
        return theta;
    }

    public void setTheta(double theta) {
        this.theta = theta;
    }

    public double phi() {
        return phi;
    }

    public void setPhi(double phi) {
        this.phi = phi;
    }

    public double localAngle() {
        return localAngle;
    }

    public void setLocalAngle(double localAngle) {
        this.localAngle = localAngle;
    }

    public double path() {
        return path;
    }

    public void setPath(double path) {
        this.path = path;
    }
    
    public static Trajectory readTrajectory(DataBank bank, int row) {
        Trajectory t = new Trajectory(bank.getShort("id", row),
                            bank.getByte("detector", row),
                            bank.getByte("sector", row),
                            bank.getByte("layer", row));
        t.setPosition(bank.getFloat("x", row), bank.getFloat("y", row), bank.getFloat("z", row));
        t.setPhi(bank.getFloat("phi", row));
        t.setTheta(bank.getFloat("theta", row));
        t.setLocalAngle(bank.getFloat("langle", row));
        t.setPath(bank.getFloat("path", row));
        return t;
    }
    
}
