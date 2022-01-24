package modules;

import java.util.ArrayList;
import java.util.List;
import objects.Event;
import objects.Hit;
import analysis.Module;
import org.jlab.groot.data.H1F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.math.F1D;

/**
 *
 * @author devita
 */
public class HitModule extends Module {
    
    
    public HitModule() {
        super("Hits", false);
    }
    
    public DataGroup hitGroup(int col) {
        String[] names = new String[]{"SVT", "BMTC", "BMTZ"};
        double[] EMAX = { 1000, 5000, 5000};
        double[] TMAX = {  600, 440, 440};
        double[] RMAX = { 1.5, 3, 3};
        DataGroup dg = new DataGroup(3,3);
        for(int i=0; i<names.length; i++) {
            String name = names[i];
            H1F hi_energy  = histo1D("hi_energy_" + name, name + " Hit Energy",   "Counts", 100,   0, EMAX[i], col);
            H1F hi_time    = histo1D("hi_time_" + name,   name + " Hit Time",     "Counts", 10+(int) TMAX[i], -10, TMAX[i], col);
            H1F hi_resi    = histo1D("hi_resi_" + name,   name + " Hit Residual", "Counts", 100, -RMAX[i], RMAX[i], col);

            dg.addDataSet(hi_energy, 0 + i*3);
            dg.addDataSet(hi_time,   1 + i*3);
            dg.addDataSet(hi_resi,   2 + i*3);
        }
        return dg;
    }

    @Override
    public void createHistos() {
        this.getHistos().put("Hits",           this.hitGroup(44));
        this.getHistos().put("HitsOnTrack",    this.hitGroup(3));
        this.getHistos().put("HitsNotOnTrack", this.hitGroup(44));
    }
    
    @Override
    public void fillHistos(Event event) {
        List<Hit> hitsOnTrack    = new ArrayList<>();
        List<Hit> hitsNotOnTrack = new ArrayList<>();
        for(Hit hit : event.getHits()) {
            if(hit.getTrackId()>0)
                hitsOnTrack.add(hit);
            else
                hitsNotOnTrack.add(hit);
        }
        this.fillGroup(this.getHistos().get("Hits"),event.getHits());
        this.fillGroup(this.getHistos().get("HitsOnTrack"),hitsOnTrack);
        this.fillGroup(this.getHistos().get("HitsNotOnTrack"),hitsNotOnTrack);
    }
    
    public void fillGroup(DataGroup group, List<Hit> hits) {
        for(Hit hit : hits) {
            group.getH1F("hi_energy_" + hit.getName()).fill(hit.getEnergy());
            group.getH1F("hi_time_" + hit.getName()).fill(hit.getTime());    
            if(hit.getTrackId()>0)
                group.getH1F("hi_resi_" + hit.getName()).fill(hit.getResidual());    
        }
    }
    
    @Override
    public void analyzeHistos() {
        this.fitTime(this.getHistos().get("HitsOnTrack").getH1F("hi_time_BMTC"));
        this.fitTime(this.getHistos().get("HitsOnTrack").getH1F("hi_time_BMTZ"));
    }
    
    @Override
    public EmbeddedCanvasTabbed plotHistos() {
        EmbeddedCanvasTabbed canvas = new EmbeddedCanvasTabbed("Hits");
        canvas.getCanvas("Hits").draw(this.getHistos().get("Hits"));
        canvas.getCanvas("Hits").draw(this.getHistos().get("HitsNotOnTrack"));
        canvas.getCanvas("Hits").draw(this.getHistos().get("HitsOnTrack"));
        this.setPlottingOptions(canvas.getCanvas("Hits"));
        return canvas;
    }
       
    @Override
    public void setPlottingOptions(EmbeddedCanvas canvas) {
        canvas.setGridX(false);
        canvas.setGridY(false);
        for(int i=0; i<3; i++) {
            canvas.getPad(i*3+0).getAxisY().setLog(true);
            canvas.getPad(i*3+1).getAxisY().setLog(true);
        }
    }

    private void fitTime(H1F hi) {
        F1D ftime = new F1D("ftime", "[amp]*gaus(x,[mean],[sigma])",100, 300);
        double amp   = hi.getBinContent(hi.getDataSize(0)/2);
        double mean  = hi.getDataX(hi.getDataSize(0)/2);
        double sigma = hi.getRMS();
        ftime.setParameter(0, amp);
        ftime.setParameter(1, mean);
        ftime.setParameter(2, sigma);
        ftime.setRange(mean-sigma*2, mean+sigma*2);
        ftime.setOptStat("1110");
        ftime.setLineColor(2);
        ftime.setLineWidth(2);
        DataFitter.fit(ftime, hi, "Q");
    }
}
