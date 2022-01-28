package modules;

import analysis.Constants;
import java.util.ArrayList;
import java.util.List;
import objects.Event;
import objects.Hit;
import analysis.Module;
import org.jlab.detector.base.DetectorType;
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
        double[] RMAX = { 1500, 3000, 3000};
        DataGroup dg = new DataGroup(3,3);
        for(int i=0; i<names.length; i++) {
            String name = names[i];
            H1F hi_energy  = histo1D("hi_energy_" + name, name + " Hit Energy",   "Counts", 100,   0, EMAX[i], col);
            H1F hi_time    = histo1D("hi_time_" + name,   name + " Hit Time",     "Counts", 10+(int) TMAX[i], -10, TMAX[i], col);
            H1F hi_resi    = histo1D("hi_resi_" + name,   name + " Hit Residual (um)", "Counts", 100, -RMAX[i], RMAX[i], col);

            dg.addDataSet(hi_energy, 0 + i*3);
            dg.addDataSet(hi_time,   1 + i*3);
            dg.addDataSet(hi_resi,   2 + i*3);
        }
        return dg;
    }

    public DataGroup occSVTGroup(int col) {
        DataGroup dg = new DataGroup(3,2);
        for(int il=0; il<Constants.SVTLAYERS; il++) {
            int layer = il+1;
            String name = "layer"+layer;
            int nbins = Constants.SVTSTRIPS*Constants.SVTSECTORS[il];
            H1F hi_occ  = histo1D("hi_occ_" + name, "L" + layer + "Strip",   "Counts", nbins,   0, nbins, col);           
            dg.addDataSet(hi_occ, 0 + il);
        }
        return dg;
    }

    public DataGroup occBMTGroup(String name, int col) {
        DataGroup dg = new DataGroup(3,3);
        for(int ir=0; ir<Constants.BMTREGIONS; ir++) {
            for(int is=0; is<Constants.BMTSECTORS; is++) {
                int sector  = is+1;
                int layer   = Constants.BMTCLAYERS[ir];
                int nstrips = Constants.BMTCSTRIPS[ir];
                if(name.equals("BMTZ")) {
                    layer = Constants.BMTZLAYERS[ir];
                    nstrips = Constants.BMTZSTRIPS[ir];
                }
                H1F hi_occ  = histo1D("hi_occ_" +layer + sector, "L" + layer + "S" + sector + " Strip",   "Counts", nstrips+1,   0, nstrips+1, col);
                dg.addDataSet(hi_occ, is + ir*Constants.BMTSECTORS);
            }
        }
        return dg;
    }

    @Override
    public void createHistos() {
        this.getHistos().put("Hits",           this.hitGroup(44));
        this.getHistos().put("HitsOnTrack",    this.hitGroup(3));
        this.getHistos().put("HitsNotOnTrack", this.hitGroup(44));
        this.getHistos().put("SVT",            this.occSVTGroup(44));
        this.getHistos().put("SVTOnTrack",     this.occSVTGroup(3));
        this.getHistos().put("BMTC",           this.occBMTGroup("BMTC", 44));
        this.getHistos().put("BMTCOnTrack",    this.occBMTGroup("BMTC", 3));
        this.getHistos().put("BMTZ",           this.occBMTGroup("BMTZ", 44));
        this.getHistos().put("BMTZOnTrack",    this.occBMTGroup("BMTZ", 3));
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
        this.fillSVTGroup(this.getHistos().get("SVT"),event.getHits());
        this.fillSVTGroup(this.getHistos().get("SVTOnTrack"),hitsOnTrack);
        this.fillBMTGroup(this.getHistos().get("BMTC"), "BMTC", event.getHits());
        this.fillBMTGroup(this.getHistos().get("BMTCOnTrack"), "BMTC", hitsOnTrack);
        this.fillBMTGroup(this.getHistos().get("BMTZ"), "BMTZ", event.getHits());
        this.fillBMTGroup(this.getHistos().get("BMTZOnTrack"), "BMTZ", hitsOnTrack);
    }
    
    public void fillGroup(DataGroup group, List<Hit> hits) {
        for(Hit hit : hits) {
            group.getH1F("hi_energy_" + hit.getName()).fill(hit.getEnergy());
            group.getH1F("hi_time_" + hit.getName()).fill(hit.getTime());    
            if(hit.getTrackId()>0)
                group.getH1F("hi_resi_" + hit.getName()).fill(hit.getResidual());    
        }
    }

    public void fillSVTGroup(DataGroup group, List<Hit> hits) {
        for(Hit hit : hits) {
            if(hit.getType()==DetectorType.BST) 
                group.getH1F("hi_occ_layer"+hit.getLayer()).fill(hit.getStrip()+Constants.SVTSTRIPS*(hit.getSector()-1));  
        }
    }
    
    public void fillBMTGroup(DataGroup group, String name, List<Hit> hits) {
        for(Hit hit : hits) {
            if(hit.getName().equals(name)) {
                group.getH1F("hi_occ_"+hit.getLayer()+hit.getSector()).fill(hit.getStrip()); 
            }
        }
    }
    
    @Override
    public void analyzeHistos() {
        this.fitTime(this.getHistos().get("HitsOnTrack").getH1F("hi_time_BMTC"));
        this.fitTime(this.getHistos().get("HitsOnTrack").getH1F("hi_time_BMTZ"));
    }
    
    @Override
    public EmbeddedCanvasTabbed plotHistos() {
        EmbeddedCanvasTabbed canvas = new EmbeddedCanvasTabbed("Hits", "SVT", "BMTC", "BMTZ");
        canvas.getCanvas("Hits").draw(this.getHistos().get("Hits"));
        canvas.getCanvas("Hits").draw(this.getHistos().get("HitsNotOnTrack"));
        canvas.getCanvas("Hits").draw(this.getHistos().get("HitsOnTrack"));
        canvas.getCanvas("SVT").draw(this.getHistos().get("SVT"));
        canvas.getCanvas("SVT").draw(this.getHistos().get("SVTOnTrack"));
        canvas.getCanvas("BMTC").draw(this.getHistos().get("BMTC"));
        canvas.getCanvas("BMTC").draw(this.getHistos().get("BMTCOnTrack"));
        canvas.getCanvas("BMTZ").draw(this.getHistos().get("BMTZ"));
        canvas.getCanvas("BMTZ").draw(this.getHistos().get("BMTZOnTrack"));
        this.setPlottingOptions(canvas.getCanvas("Hits"));
        super.setPlottingOptions(canvas.getCanvas("SVT"));
        super.setPlottingOptions(canvas.getCanvas("BMTC"));
        super.setPlottingOptions(canvas.getCanvas("BMTZ"));
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
