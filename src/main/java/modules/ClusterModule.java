package modules;

import java.util.ArrayList;
import java.util.List;
import objects.Cluster;
import objects.Event;
import analysis.Module;
import org.jlab.groot.data.H1F;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;

/**
 *
 * @author devita
 */
public class ClusterModule extends Module {
    
    public ClusterModule() {
        super("Clusters");
    }
    
    public DataGroup clusterGroup(int col) {
        String[] names = new String[]{"SVT", "BMTC", "BMTZ"};
        double[] EMAX = {2000, 10000, 10000};
        DataGroup dg = new DataGroup(3,3);
        for(int i=0; i<names.length; i++) {
            String name = names[i];
            H1F hi_size    = histo1D("hi_size_" + name,   name + " Cluster Size",   "Counts",  20, 0, 20, col);
            H1F hi_energy  = histo1D("hi_energy_" + name, name + " Cluster Energy", "Counts", 100, 0, EMAX[i], col);
            H1F hi_time    = histo1D("hi_time_" + name,   name + " Cluster Time",   "Counts", 440, 0, 440, col);

            dg.addDataSet(hi_size,   0 + i*3);
            dg.addDataSet(hi_energy, 1 + i*3);
            dg.addDataSet(hi_time,   2 + i*3);
        }
        return dg;
    }

    @Override
    public void createHistos() {
        this.getHistos().put("Clusters",        this.clusterGroup(44));
        this.getHistos().put("ClustersOnTrack", this.clusterGroup(3));
        this.getHistos().put("ClustersNotOnTrack", this.clusterGroup(44));
    }
    
    @Override
    public void fillHistos(Event event) {
        List<Cluster> clustersOnTrack    = new ArrayList<>();
        List<Cluster> clustersNotOnTrack = new ArrayList<>();
        for(Cluster cluster : event.getClusters()) {
            if(cluster.getTrackId()>0)
                clustersOnTrack.add(cluster);
            else
                clustersNotOnTrack.add(cluster);
        }
        this.fillGroup(this.getHistos().get("Clusters"),event.getClusters());
        this.fillGroup(this.getHistos().get("ClustersOnTrack"),clustersOnTrack);
        this.fillGroup(this.getHistos().get("ClustersNotOnTrack"),clustersNotOnTrack);
    }
    
    public void fillGroup(DataGroup group, List<Cluster> clusters) {
        for(Cluster cluster : clusters) {
            group.getH1F("hi_size_" + cluster.getName()).fill(cluster.getSize());
            group.getH1F("hi_energy_" + cluster.getName()).fill(cluster.getEnergy());
            group.getH1F("hi_time_" + cluster.getName()).fill(cluster.getTime());            
        }
    }
    
    @Override
    public EmbeddedCanvasTabbed plotHistos() {
        EmbeddedCanvasTabbed canvas = new EmbeddedCanvasTabbed("Clusters");
        canvas.getCanvas("Clusters").draw(this.getHistos().get("Clusters"));
        canvas.getCanvas("Clusters").draw(this.getHistos().get("ClustersNotOnTrack"));
        canvas.getCanvas("Clusters").draw(this.getHistos().get("ClustersOnTrack"));
        this.setPlottingOptions(canvas.getCanvas("Clusters"));
        return canvas;
    }
       
    @Override
    public void setPlottingOptions(EmbeddedCanvas canvas) {
        canvas.setGridX(false);
        canvas.setGridY(false);
        for(int i=0; i<9; i++)
            canvas.getPad(i).getAxisY().setLog(true);
    }

}
