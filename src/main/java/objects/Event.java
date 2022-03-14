package objects;

import analysis.Constants;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.detector.base.DetectorType;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

/**
 *
 * @author devita
 */
public class Event {
    private final boolean debug = false;
    
    private int run;
    private int event;
    private double startTime;
    private Track mcParticle;
    private final List<Track> particles  = new ArrayList<>();
    private final List<Track> tracks     = new ArrayList<>();
    private final List<Track> seeds      = new ArrayList<>();
    private final List<Track> mcTracks   = new ArrayList<>();
    private final List<Cluster> clusters = new ArrayList<>();
    private final List<Hit> hits         = new ArrayList<>();
    private final List<Trajectory> trajs = new ArrayList<>();
    private final List<True> trues       = new ArrayList<>();
    private final Map<Integer, Integer> trackMap = new HashMap<>();
    private final Map<Integer, Integer> seedMap = new HashMap<>();
    private final Map<Integer, List<Integer>> clusterMap = new HashMap<>();
    private final Map<Integer, List<Integer>> hitMap = new HashMap<>();
    private final Map<Integer, List<Integer>> trajMap = new HashMap<>();
    private DataEvent hipoEvent;

    public Event(DataEvent event) {
        this.hipoEvent = event;
        this.readEvent(event);
        if(debug) System.out.println("Read event with " + tracks.size() + " particles");
    }
    


    private DataBank getBank(DataEvent de, String bankName) {
        DataBank bank = null;
        if (de.hasBank(bankName)) {
            bank = de.getBank(bankName);
        }
        return bank;
    }

    private void readHeader(DataEvent event) {
        DataBank head = this.getBank(event, "RUN::config");
        if(head!=null) {
            this.run   = head.getInt("run", 0);
            this.event = head.getInt("event", 0);
        }
    }
    
    private void readMCParticle(DataEvent event) {
        DataBank mc  = this.getBank(event, "MC::Particle");
        if(mc!=null) {
            int index = 0;
            int pid = 0;
            if(Constants.PID!=0) {
                pid = Constants.PID;
                for (int loop = 0; loop < mc.rows(); loop++) {
                    double px = mc.getFloat("px", loop);
                    double py = mc.getFloat("py", loop);
                    double pz = mc.getFloat("pz", loop);
                    double theta = Math.toDegrees(Math.acos(pz/Math.sqrt(px*px+py*py+pz*pz)));
                    if(mc.getInt("pid", loop)==Constants.PID  && theta>Constants.THMIN) {
                        index = loop;
                        break;
                    }
                }
            }
            else {
                if(mc.getInt("pid", 0)!=0) 
                    pid = mc.getInt("pid", 0);
                else
                    pid = -13;
            }
            mcParticle = new Track(pid,
                            mc.getFloat("px", index),
                            mc.getFloat("py", index),
                            mc.getFloat("pz", index),
                            mc.getFloat("vx", index),
                            mc.getFloat("vy", index),
                            mc.getFloat("vz", index));
        }
    }
    
    private void readParticles(DataEvent event) {
        DataBank recPart   = this.getBank(event, "REC::Particle");
        DataBank recTrack  = this.getBank(event, "REC::Track");
        if(recPart!=null) {
            for (int loop = 0; loop < recPart.rows(); loop++) {    
                int pid    = recPart.getInt("pid", loop);
                int charge = recPart.getByte("charge", loop);
                if(pid==0) {
                    pid = charge==0 ? 22 : charge*211;
                }
                Track t = new Track(pid,
                            recPart.getFloat("px", loop),
                            recPart.getFloat("py", loop),
                            recPart.getFloat("pz", loop),
                            recPart.getFloat("vx", loop),
                            recPart.getFloat("vy", loop),
                            recPart.getFloat("vz", loop));
                t.setRECStatus(recPart.getShort("status", loop));
                t.setChi2pid(recPart.getFloat("chi2pid", loop));
                if(recTrack!=null) {
                    for(int j=0; j<recTrack.rows(); j++) {
                        if(recTrack.getShort("pindex", j)==loop) {
                            t.setSector(recTrack.getByte("sector", j));
                            t.setNDF(recTrack.getShort("NDF", j));
                            t.setChi2(recTrack.getFloat("chi2", j));
                            t.setStatus(recTrack.getShort("status", j));
                            break;
                        }
                    }
                }                    
                particles.add(t);
            }
        }
    }
    
    private void readStartTime(DataEvent event) {
        DataBank recEvent = this.getBank(event, "REC::Event");
        if(recEvent!=null) {
            startTime = recEvent.getFloat("startTime",0);
        }
    }
    
    private void readTracks(DataEvent event) {
        DataBank cvtBank   = this.getBank(event, "CVTRec::Tracks");
        DataBank cosBank   = this.getBank(event, "CVTRec::Cosmics");
        DataBank recPart   = this.getBank(event, "REC::Particle");
        DataBank recTrack  = this.getBank(event, "REC::Track");
        DataBank runConfig = this.getBank(event, "RUN::config");
        if(cvtBank!=null) {
            for(int i=0; i<cvtBank.rows(); i++) {
                Track track = Track.readTrack(cvtBank, i);
                if(runConfig!=null) track.addScale(runConfig);
                if(recPart!=null && recTrack!=null) track.addEBinfo(recPart, recTrack);
                tracks.add(track);
                trackMap.put(track.getId(), i);
            }
        }
        else if(cosBank!=null) {
            for(int i=0; i<cosBank.rows(); i++) {
                Track track = Track.readRay(cosBank, i);
                if(runConfig!=null) track.addScale(runConfig);
                tracks.add(track);
                trackMap.put(track.getId(), i);
            }                
        }
    }
        
    private void readSeeds(DataEvent event) {
        DataBank trackBank = this.getBank(event, "CVTRec::Seeds");
        DataBank cosmiBank = this.getBank(event, "CVTRec::CosmicSeeds");
        if(trackBank!=null) {
            for(int i=0; i<trackBank.rows(); i++) {
                Track track = Track.readSeed(trackBank, i);
                seeds.add(track);
                seedMap.put(track.getId(), i);
            }
        }
        else if(cosmiBank!=null) {
            for(int i=0; i<cosmiBank.rows(); i++) {
                Track track = Track.readRay(cosmiBank, i);
                seeds.add(track);
                seedMap.put(track.getId(), i);
            }                
        }
    }
        
    private void readTrajectory(DataEvent event) {
        DataBank trajBank = this.getBank(event, "CVTRec::Trajectory");
        if(trajBank!=null) {
            for(int i=0; i<trajBank.rows(); i++) {
                Trajectory traj = Trajectory.readTrajectory(trajBank, i);
                trajs.add(traj);
                if(!this.trajMap.containsKey(traj.getTrackId()))
                    this.trajMap.put(traj.getTrackId(), new ArrayList<>());
                this.trajMap.get(traj.getTrackId()).add(i);
            }
        }
    }
        
    private void readClusters(DataEvent event) {
        DataBank bank = this.getBank(event, "BSTRec::Clusters");
        if(bank==null) return;
        for(int i=0; i<bank.rows(); i++) {
            Cluster cluster = Cluster.readCluster(bank, i, DetectorType.BST);
            clusters.add(cluster);
        }
        bank = this.getBank(event, "BMTRec::Clusters");
        if(bank==null) return;
        for(int i=0; i<bank.rows(); i++) {
            Cluster cluster = Cluster.readCluster(bank, i, DetectorType.BMT);
            clusters.add(cluster);
        }
        for(int i=0; i<clusters.size(); i++) {
            Cluster cluster = clusters.get(i);
            if(cluster.getTrackId()>0) {
                if(!this.clusterMap.containsKey(cluster.getTrackId()))
                    this.clusterMap.put(cluster.getTrackId(), new ArrayList<>());
                this.clusterMap.get(cluster.getTrackId()).add(i);
            }
        }
    }
                
    private void readHits(DataEvent event) {
        DataBank bank = this.getBank(event, "BSTRec::Hits");
        if(bank==null) return;
        for(int i=0; i<bank.rows(); i++) {
            Hit hit = Hit.readHit(bank, i, DetectorType.BST);
            hits.add(hit);
        }
        bank = this.getBank(event, "BMTRec::Hits");
        if(bank==null) return;
        for(int i=0; i<bank.rows(); i++) {
            Hit hit = Hit.readHit(bank, i, DetectorType.BMT);
            hits.add(hit);
        }
        for(int i=0; i<hits.size(); i++) {
            Hit hit = hits.get(i);
            if(hit.getTrackId()>0) {
                if(!this.hitMap.containsKey(hit.getTrackId()))
                    this.hitMap.put(hit.getTrackId(), new ArrayList<>());
                this.hitMap.get(hit.getTrackId()).add(i);
            }
        }
    }

    private void readTrues(DataEvent event) {
        DataBank mc  = this.getBank(event, "MC::True");
        DataBank bmt = this.getBank(event, "BMT::adc");
        DataBank bst = this.getBank(event, "BST::adc");
        if(mc==null || (bmt==null && bst==null)) return;
        int offset = 0;
        if(bmt!=null) {
            for(int i=0; i<bmt.rows(); i++) {
                True t = True.readTruth(bmt, mc, i, offset, DetectorType.BMT);
                trues.add(t);
            }
            offset += bmt.rows();
        }
        if(bst!=null) {
            for(int i=0; i<bst.rows(); i++) {
                True t = True.readTruth(bst, mc, i, offset, DetectorType.BST);
                trues.add(t);
            }
        }
        Collections.sort(trues);
    }
    
    private void loadMap(Map<Integer, List<DetectorResponse>> map, DetectorResponse response) {
        final int iTo = response.getAssociation();
        if (map.containsKey(iTo)) {
            map.get(iTo).add(response);
        } else {
            List<DetectorResponse> iFrom = new ArrayList<>();
            map.put(iTo, iFrom);
            map.get(iTo).add(response);
        }
    }

    private void readEvent(DataEvent de) {
        this.readHeader(de);
        this.readMCParticle(de);
        this.readParticles(de);
        this.readStartTime(de);
        this.readTracks(de);
        this.readTrajectory(de);
        this.readSeeds(de);
        this.readClusters(de);
        this.readHits(de);
        this.readTrues(de);
    }

    public List<Track> getTracks() {
        return tracks;
    }
 
    public List<Track> getSeeds() {
        return seeds;
    }
 
    public Map<Integer, Integer> getTrackMap() {
        return trackMap;
    }
    
    public Map<Integer, Integer> getSeedMap() {
        return seedMap;
    }

    public List<Cluster> getClusters() {
        return clusters;
    }

    public Map<Integer, List<Integer>> getClusterMap() {
        return clusterMap;
    }

    public List<Hit> getHits() {
        return hits;
    }

    public Map<Integer, List<Integer>> getHitMap() {
        return hitMap;
    }      

    public List<Trajectory> getTrajectories(int trackId) {
        List<Trajectory> trackTrajs = new ArrayList<>();
        for(int i : this.getTrajectoryMap().get(trackId))
            trackTrajs.add(this.trajs.get(i));
        return trackTrajs;
    }

    public Map<Integer, List<Integer>> getTrajectoryMap() {
        return trajMap;
    }
    
    public Track getMCTrack(boolean cosmics) {
        if(cosmics && mcParticle!=null) mcParticle.toCosmic();
        return mcParticle;
    }

    public List<Track> getParticles() {
        return particles;
    }

    public List<True> getTrues() {
        return trues;
    }

    public int getRun() {
        return run;
    }

    public int getEvent() {
        return event;
    }

    
    public double getStartTime() {
        return startTime;
    }

    public DataEvent getHipoEvent() {
        return hipoEvent;
    }
    
    
}
