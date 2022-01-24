package objects;

import analysis.Constants;
import java.util.ArrayList;
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
    
    private Track mcParticle;
    private double startTime;
    private final List<Track> tracks   = new ArrayList<>();
    private final List<Track> seeds    = new ArrayList<>();
    private final List<Track> mcTracks = new ArrayList<>();
    private final List<Cluster> clusters = new ArrayList<>();
    private final List<Hit> hits = new ArrayList<>();
    private final Map<Integer, Integer> trackMap = new HashMap<>();
    private final Map<Integer, Integer> seedMap = new HashMap<>();
    private final Map<Integer, List<Integer>> clusterMap = new HashMap<>();
    private final Map<Integer, List<Integer>> hitMap = new HashMap<>();


    public Event(DataEvent event) {
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


    private void readMCParticle(DataEvent event) {
        DataBank mc  = this.getBank(event, "MC::Particle");
        if(mc!=null) {
            int index = 0;
            if(Constants.PID!=0) {
                for (int loop = 0; loop < mc.rows(); loop++) {
                    if(mc.getInt("pid", loop)==Constants.PID) {
                        index = loop;
                        break;
                    }
                }
            }
            else {
                Constants.PID = mc.getInt("pid", 0);
            }
            mcParticle = new Track(Constants.PID,
                            mc.getFloat("px", index),
                            mc.getFloat("py", index),
                            mc.getFloat("pz", index),
                            mc.getFloat("vx", index),
                            mc.getFloat("vy", index),
                            mc.getFloat("vz", index));
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
        this.readMCParticle(de);
        this.readStartTime(de);
        this.readTracks(de);
        this.readSeeds(de);
        this.readClusters(de);
        this.readHits(de);
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
    
    public Track getMCTrack(boolean cosmics) {
        if(cosmics && mcParticle!=null) mcParticle.toCosmic();
        return mcParticle;
    }

    public double getStartTime() {
        return startTime;
    }
    
    
}
