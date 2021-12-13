package analysis;

/**
 *
 * @author devita
 */
public class Constants {
    
    public static final int BMTSECTORS = 3;
    public static final int BMTREGIONS = 3;
    public static final int SVTLAYERS = 6;
    public static final int[] SVTSECTORS = {10, 10, 14, 14, 18, 18};
    public static final int[] BMTCLAYERS = {1,4,6};
    public static final int[] BMTZLAYERS = {2,3,5};
    public static final double[] BMTRADIUS = {147.746, 162.646, 177.646, 192.646, 207.646, 222.646}; // radius at half drift
    
    
    public static int PID = 2212;
    
    public static final int NSIGMA = 3;
    public static final double SIGMA_P = 5; // 3;
    public static final double SIGMA_PHI = 0.29;//0.57; //0.3; //0.12;
    public static final double SIGMA_THETA = 1.146;//0.57; //0.3; //0.2;

    
    public static void setPID(int PID) {
        Constants.PID = PID;
    }
    
    
}
