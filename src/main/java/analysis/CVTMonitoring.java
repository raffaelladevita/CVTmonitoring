package analysis;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import objects.Event;
import java.util.ArrayList;
import java.util.List;
import javax.swing.JFrame;
import javax.swing.JTabbedPane;
import modules.ClusterModule;
import modules.HitModule;
import modules.MCModule;
import modules.PullsModule;
import modules.ResidualModule;
import modules.TrackModule;
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
    ByteArrayOutputStream pipeOut = new ByteArrayOutputStream();
    private static PrintStream outStream = System.out;
    private static PrintStream errStream = System.err;
        
    ArrayList<Module>    modules = new ArrayList<>();

    public CVTMonitoring(int pid, String opts) {
        this.init(pid, opts);
    }
    
    

    private void init(int pid, String opts) {
        
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

        Constants.setPID(pid);
        this.modules.add(new TrackModule());
        this.modules.add(new ClusterModule());
        this.modules.add(new HitModule());
        this.modules.add(new ResidualModule());
        this.modules.add(new PullsModule());
        this.modules.add(new MCModule());
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
        TDirectory dir = new TDirectory();
        dir.readFile(fileName);
        System.out.println(dir.getDirectoryList());
        dir.cd();
        dir.pwd();
        for(Module m : modules) {
            m.readDataGroup(dir);
        }
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
        parser.addOption("-plot"       ,"1",    "display histograms (0/1)");
        parser.addOption("-stats"      ,"",     "histogram stat option (e.g. \"10\" will display entries)");
        parser.addOption("-pid"        ,"2212", "MC particle PID");
        
        parser.parse(args);
        
        String namePrefix  = parser.getOption("-o").stringValue();        
        String histoName   = "histo.hipo";
        if(!namePrefix.isEmpty()) {
            histoName  = namePrefix + "_" + histoName; 
        }
        int     maxEvents    = parser.getOption("-n").intValue();        
        boolean readHistos   = (parser.getOption("-histo").intValue()!=0);            
        boolean openWindow   = (parser.getOption("-plot").intValue()!=0);
        String  optStats     = parser.getOption("-stats").stringValue(); 
        int     pid          = parser.getOption("-pid").intValue();        

        if(!openWindow) System.setProperty("java.awt.headless", "true");

        CVTMonitoring cvtMon = new CVTMonitoring(pid, optStats);
        
        List<String> inputList = parser.getInputList();
        if(inputList.isEmpty()==true){
            parser.printUsage();
            System.out.println("\n >>>> error: no input file is specified....\n");
            System.exit(0);
        }

        if(readHistos) {
            cvtMon.readHistos(inputList.get(0));
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
        }
    }

}