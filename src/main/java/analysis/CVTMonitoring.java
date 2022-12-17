package analysis;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import objects.Event;
import java.util.ArrayList;
import java.util.List;
import javax.swing.JFrame;
import javax.swing.JTabbedPane;
import modules.ClusterModule;
import modules.CrossModule;
import modules.DoubletsModule;
import modules.ElasticModule;
import modules.HitModule;
import modules.LorentzModule;
import modules.MCModule;
import modules.PIDModule;
import modules.TruthModule;
import modules.PullsModule;
import modules.ResidualModule;
import modules.TrackModule;
import modules.VertexModule;
import org.jlab.groot.base.GStyle;

import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;


import org.jlab.groot.data.TDirectory;
import org.jlab.jnp.utils.benchmark.ProgressPrintout;
import org.jlab.utils.options.OptionParser;

/**
 *
 * @author devita
 */

public class CVTMonitoring {
    
    private final boolean debug = false;
    private boolean fastmode = false;
    ByteArrayOutputStream pipeOut = new ByteArrayOutputStream();
    private static PrintStream outStream = System.out;
    private static PrintStream errStream = System.err;
        
    ArrayList<Module>    modules = new ArrayList<>();

    public CVTMonitoring(String active, boolean mode, int pid, double ebeam, double[] beamSpot, 
                         boolean cosmics, int residualScale, String opts, boolean lund) {
        this.init(active, mode, pid, ebeam, beamSpot, cosmics, residualScale, opts, lund);
    }
    
    

    private void init(String active, boolean mode, int pid, double ebeam, double[] beamSpot, 
                      boolean cosmics, int residualScale, String opts, boolean lund) {
        
        GStyle.getH1FAttributes().setOptStat(opts);
        GStyle.getAxisAttributesX().setTitleFontSize(24);
        GStyle.getAxisAttributesX().setLabelFontSize(18);
        GStyle.getAxisAttributesY().setTitleFontSize(24);
        GStyle.getAxisAttributesY().setLabelFontSize(18);
        GStyle.getAxisAttributesZ().setLabelFontSize(14);
//        GStyle.getAxisAttributesX().setLabelFontName("Arial");
//        GStyle.getAxisAttributesY().setLabelFontName("Arial");
//        GStyle.getAxisAttributesZ().setLabelFontName("Arial");
//        GStyle.getAxisAttributesX().setTitleFontName("Arial");
//        GStyle.getAxisAttributesY().setTitleFontName("Arial");
//        GStyle.getAxisAttributesZ().setTitleFontName("Arial");
        GStyle.setGraphicsFrameLineWidth(2);
        GStyle.getH1FAttributes().setLineWidth(1);

        Constants.setFASTMODE(mode);        
        Constants.setPID(pid);
        Constants.setCharge(pid);
        
        if(beamSpot.length==2) Constants.setBEAMSPOT(beamSpot);
        
        this.addModule(active, new TrackModule(cosmics));
        this.addModule(active, new VertexModule());
        this.addModule(active, new HitModule(residualScale));
        this.addModule(active, new CrossModule());
        this.addModule(active, new ClusterModule());
        this.addModule(active, new LorentzModule());
        this.addModule(active, new ResidualModule(residualScale));
        this.addModule(active, new PullsModule(residualScale));
        this.addModule(active, new MCModule(cosmics));
        this.addModule(active, new TruthModule());
        this.addModule(active, new PIDModule());
        this.addModule(active, new ElasticModule(ebeam, lund));
        this.addModule(active, new DoubletsModule(ebeam, lund));
    }

    private void addModule(String active, Module module) {
        boolean flag = true;
        if(active!=null && !active.isEmpty()) {
            flag = false;
            String[] mods = active.split(":");
            for(String m : mods) {
                if(m.trim().equalsIgnoreCase(module.getName())) {
                    flag = true;
                    break;
                }
            }
        }
        if(flag) {
            System.out.println("Adding module " + module.getName());
            this.modules.add(module);
        }
    }
    
    private void processEvent(DataEvent de) {
        Event event = new Event(de);
        
        for(Module m : modules) m.processEvent(event);
    }

    private void analyzeHistos() {
        for(Module m : modules) m.analyzeHistos();
    }

    public JTabbedPane plotHistos() {
        JTabbedPane panel = new JTabbedPane();
        for(Module m : modules) {
            panel.add(m.getName(), m.plotHistos());
        }
        return panel;
    }
    
    public void readHistos(String fileName) {
        System.out.println("Opening file: " + fileName);
        PrintStream pipeStream = new PrintStream(pipeOut);
        System.setOut(pipeStream);
        System.setErr(pipeStream);
        TDirectory dir = new TDirectory();
        dir.readFile(fileName);
        System.out.println(dir.getDirectoryList());
        dir.cd();
        dir.pwd();
        for(Module m : modules) {
            m.readDataGroup(dir);
        }
        System.setOut(outStream);
        System.setErr(errStream);
    }

    public void saveHistos(String fileName) {
        System.out.println("\n>>>>> Saving histograms to file " + fileName);
        PrintStream pipeStream = new PrintStream(pipeOut);
        System.setOut(pipeStream);
        System.setErr(pipeStream);
        TDirectory dir = new TDirectory();
        for(Module m : modules) {
            m.writeDataGroup(dir);
        }
        dir.writeFile(fileName);
        System.setOut(outStream);
        System.setErr(errStream);
    }

    private void printHistos() {
        System.out.println("\n>>>>> Printing canvases to directory plots");
        for(Module m : modules) {
            m.printHistos("plots");
        }
    }
    
    private void testHistos() {
        for(Module m : modules) {
            m.testHistos();
        }
    }
    
    public static void main(String[] args) {
        
        OptionParser parser = new OptionParser("cvtMonitoring [options] file1 file2 ... fileN");
        parser.setRequiresInputList(false);
        // valid options for event-base analysis
        parser.addOption("-o"          ,"",     "histogram file name prefix");
        parser.addOption("-n"          ,"-1",   "maximum number of events to process");
        // histogram based analysis
        parser.addOption("-histo"      ,"0",    "read histogram file (0/1)");
        parser.addOption("-fast"       ,"0",    "run in fast mode, i.e. only tracks (0/1)");
        parser.addOption("-plot"       ,"1",    "display histograms (0/1)");
        parser.addOption("-print"      ,"0",    "print histograms (0/1)");
        parser.addOption("-stats"      ,"",     "histogram stat option (e.g. \"10\" will display entries)");
        parser.addOption("-pid"        ,"0",    "MC particle PID (default: use first particle in the bank)");
        parser.addOption("-beam"       ,"10.6", "Beam energy (GeV)");
//        parser.addOption("-spot"       ,"0:0" , "Beam spot coordinates (cm)");
        parser.addOption("-cosmics"    ,"0",    "analyze as cosmics (0=false, 1=true)");
        parser.addOption("-residual"   ,"1",    "residual scale (1=cm, 10=mm)");
        parser.addOption("-lund"       ,"0",    "save events to lund (0=false, 1=true");
        parser.addOption("-modules"    ,"",     "comma-separated list of modules to be activated");
        
        parser.parse(args);
        
        String namePrefix  = parser.getOption("-o").stringValue();        
        String histoName   = "histo.hipo";
        if(!namePrefix.isEmpty()) {
            histoName  = namePrefix + "_" + histoName; 
        }
        int     maxEvents     = parser.getOption("-n").intValue();        
        boolean readHistos    = (parser.getOption("-histo").intValue()!=0);            
        boolean openWindow    = (parser.getOption("-plot").intValue()!=0);
        boolean printHistos   = (parser.getOption("-print").intValue()!=0);
        String  optStats      = parser.getOption("-stats").stringValue(); 
        int     pid           = parser.getOption("-pid").intValue();  
        double  ebeam         = parser.getOption("-beam").doubleValue();
        double[] beamSpot = {0, 0};
//        String[]  spot    = parser.getOption("-spot").stringValue().split(":");
//        if(spot.length==2) for(int i=0; i<2; i++) beamSpot[i] = Double.parseDouble(spot[i]);
        boolean cosmics       = parser.getOption("-cosmics").intValue()!=0;
        int     residualScale = parser.getOption("-residual").intValue();  
        boolean lund          = parser.getOption("-lund").intValue()!=0;
        String  modules       = parser.getOption("-modules").stringValue();
        boolean fast          = parser.getOption("-fast").intValue()!=0;
        
        if(!openWindow) System.setProperty("java.awt.headless", "true");

        CVTMonitoring cvtMon = new CVTMonitoring(modules, fast, pid, ebeam, beamSpot, cosmics, residualScale, optStats, lund);
        
        List<String> inputList = parser.getInputList();
        if(inputList.isEmpty()==true){
            parser.printUsage();
            System.out.println("\n >>>> error: no input file is specified....\n");
            System.exit(0);
        }

        if(readHistos) {
            cvtMon.readHistos(inputList.get(0));
            cvtMon.analyzeHistos();
            cvtMon.testHistos();
        }
        else{

            ProgressPrintout progress = new ProgressPrintout();

            int counter = -1;
            for(String inputFile : inputList){
                HipoDataSource reader = new HipoDataSource();
                reader.open(inputFile);

                
                while (reader.hasEvent()) {

                    counter++;

                    DataEvent event = reader.getNextEvent();
                    cvtMon.processEvent(event);
                    
                    progress.updateStatus();
                    if(maxEvents>0){
                        if(counter>=maxEvents) break;
                    }
                }
                progress.showStatus();
                reader.close();
            }    
            cvtMon.analyzeHistos();
            cvtMon.testHistos();
            cvtMon.saveHistos(histoName);
        }

        if(openWindow) {
            JFrame frame = new JFrame("CVTMonitoring");
            frame.setSize(1400, 900);
            frame.add(cvtMon.plotHistos());
            frame.setLocationRelativeTo(null);
            frame.setVisible(true);
            if(printHistos) cvtMon.printHistos();
        }
    }

}