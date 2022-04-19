package modules;

import analysis.Constants;
import java.util.List;
import objects.Event;
import analysis.Module;
import java.util.HashMap;
import java.util.Map;
import objects.True;
import org.jlab.detector.base.DetectorType;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.graphics.EmbeddedPad;

/**
 *
 * @author devita
 */
public class TruthModule extends Module {
    
    public TruthModule() {
        super("Truth", false);
    }
    
    public DataGroup elossGroup(int col) {
        DetectorType[] types = new DetectorType[]{DetectorType.BST, DetectorType.BMT};
        DataGroup dg = new DataGroup(2,2);
        for(int i=0; i<types.length; i++) {
            String name = types[i].getName();
            H1F hi_eloss   = histo1D("hi_eloss_"  + name, name + " ELoss (MeV)", "Counts"  , 100, 0, 5, col);
            H2F hi_peloss  = histo2D("hi_peloss_" + name, "p (GeV)", name + " ELoss (MeV)", 100, 0.2, 1.2, 100, 0, 20);
            dg.addDataSet(hi_eloss,  i + 0);
            dg.addDataSet(hi_peloss, i + 2);
        }
        return dg;
    }


    @Override
    public void createHistos() {
        this.getHistos().put("ELoss", this.elossGroup(44));
    }
    
    @Override
    public void fillHistos(Event event) {
        List<True> trues = event.getTrues();
        Map<Integer,True> trueMap = new HashMap<>();
        for(True t : trues) {
            int layer = t.getLayer();
            if(t.getType()==DetectorType.BMT)
                layer += Constants.SVTLAYERS;
            if(!trueMap.containsKey(layer))
                trueMap.put(layer, t);
        }
        for(int i=0; i<2; i++) {
            int layer = 2*i+1;
            this.fillEloss(layer, layer+2, DetectorType.BST, trueMap);
        }
        for(int i=0; i<Constants.BMTREGIONS*2; i++) {
            int layer = Constants.SVTLAYERS + i +1;
            this.fillEloss(layer, layer+1, DetectorType.BMT, trueMap);
        }         
    }
    
    private boolean fillEloss(int layer1, int layer2, DetectorType type, Map<Integer,True> trueMap) {
        boolean value = false;
        if(trueMap.containsKey(layer1) && trueMap.containsKey(layer2)) {
            double edep = (trueMap.get(layer1).getEnergy()-trueMap.get(layer2).getEnergy());
            double pmom = trueMap.get(layer1).getMomentum().mag();
            if(trueMap.get(layer1).getTime()>trueMap.get(layer2).getTime()) {
                edep = -edep;
                pmom = trueMap.get(layer2).getMomentum().mag();
            }
//            if(type==DetectorType.BST && Math.abs(pmom-0.4)<0.01 && edep<4) {
//                System.out.println(trueMap.get(layer1).getEnergy() + " " + trueMap.get(layer2).getEnergy() + " " + edep);
//                value=true;
//            }
            this.getHistos().get("ELoss").getH1F("hi_eloss_"+type.getName()).fill(edep);
            this.getHistos().get("ELoss").getH2F("hi_peloss_"+type.getName()).fill(pmom, edep);                
        }
        return value;
    }

       
    @Override
    public void setPlottingOptions(String name) {
        super.setPlottingOptions(name);
        for(EmbeddedPad pad : this.getCanvas(name).getCanvasPads())
            pad.getAxisZ().setLog(true);
    }

}
