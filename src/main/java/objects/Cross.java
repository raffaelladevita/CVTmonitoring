package objects;

import analysis.Constants;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Point3D;
import org.jlab.io.base.DataBank;

/**
 *
 * @author devita
 */
public class Cross {
 
    private final int id;
    private final int sector;
    private final int region;
    private final int layer;
    private final DetectorType type;
    
    private Point3D point;
    private Point3D point0;
    
    private int trackId;
    
    public Cross(int id, int sector, int region, int layer, DetectorType type) {
        this.id = id;
        this.sector = sector;
        this.region = region;
        this.layer = layer;
        this.type = type;
    }

    public int getSector() {
        return sector;
    }

    public int getRegion() {
        return region;
    }

    public int getLayer() {
        return layer;
    }

    public Point3D getPoint() {
        return point;
    }

    public void setPoint(double x, double y, double z) {
        this.point = new Point3D(x, y, z);
    }

    public Point3D getPoint0() {
        return point0;
    }

    public void setPoint0(double x, double y, double z) {
        this.point0 = new Point3D(x, y, z);
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
    public static Cross readCross(DataBank bank, int row, DetectorType type) {
        Cross cross = new Cross(bank.getShort("ID", row),
                                bank.getByte("sector", row),
                                bank.getByte("region", row),
                                0,
                                type);
        cross.setPoint(bank.getFloat("x", row), bank.getFloat("y", row), bank.getFloat("z", row));
        cross.setPoint0(bank.getFloat("x0", row), bank.getFloat("y0", row), bank.getFloat("z0", row));
        
        cross.setTrackId(bank.getShort("trkID", row));
        return cross;
    }
    
}
