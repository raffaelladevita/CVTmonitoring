package modules;

import java.util.ArrayList;
import java.util.List;
import objects.Event;
import analysis.Module;
import objects.CVTType;
import objects.Cross;
import org.jlab.groot.data.H1F;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.group.DataGroup;

/**
 *
 * @author devita
 */
public class CrossModule extends Module {
    
    
    public CrossModule() {
        super("Crosses", false);
    }
    
    public DataGroup crossGroup(int col) {
        String[] names = new String[]{"SVT", "BMTC", "BMTZ"};
        DataGroup dg = new DataGroup(1,3);
        for(int i=0; i<1; i++) {
            String name = names[i];
            for(int ir=0; ir<3; ir++) {
                H1F hi_z0  = histo1D("hi_z0_" + name + "_R" + (ir+1), name + " z0 (cm)",   "Counts", 200, -50, 50, col);
                dg.addDataSet(hi_z0, ir);
            }
        }
        return dg;
    }

    @Override
    public void createHistos() {
        this.getHistos().put("Crosses", this.crossGroup(44));
        this.getHistos().put("CrossesNotOnTrack", this.crossGroup(44));
        this.getHistos().put("CrossesOnTrack", this.crossGroup(3));
    }
    
    @Override
    public void fillHistos(Event event) {
        List<Cross> crossesOnTrack    = new ArrayList<>();
        List<Cross> crossesNotOnTrack = new ArrayList<>();
        for(Cross cross : event.getCrosses()) {
            if(cross.getTrackId()>0)
                crossesOnTrack.add(cross);
            else
                crossesNotOnTrack.add(cross);
        }
        this.fillGroup(this.getHistos().get("Crosses"),event.getCrosses());
        this.fillGroup(this.getHistos().get("CrossesOnTrack"),crossesOnTrack);
        this.fillGroup(this.getHistos().get("CrossesNotOnTrack"),crossesNotOnTrack);
    }
    
    public void fillGroup(DataGroup group, List<Cross> crosses) {
        for(Cross cross : crosses) {
            if(cross.getType()==CVTType.SVT) group.getH1F("hi_z0_" + cross.getType().getName() + "_R" + cross.getRegion()).fill(cross.getPoint0().z());    
        }
    }
    
    @Override
    public void drawHistos() {
        this.addCanvas("Crosses");
        this.getCanvas().getCanvas("Crosses").draw(this.getHistos().get("Crosses"));
        this.getCanvas().getCanvas("Crosses").draw(this.getHistos().get("CrossesNotOnTrack"));
        this.getCanvas().getCanvas("Crosses").draw(this.getHistos().get("CrossesOnTrack"));
        this.setPlottingOptions("Crosses");
    }
        
    @Override
    public void setPlottingOptions(String name) {
        this.getCanvas().getCanvas(name).setGridX(false);
        this.getCanvas().getCanvas(name).setGridY(false);
        for(EmbeddedPad pad : this.getCanvas().getCanvas(name).getCanvasPads()) {
            pad.getAxisY().setLog(true);
        }
    }

}
