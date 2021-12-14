# CVTmonitoring

Displays CVT tracking information for beam data and MC

## To compile:

      mvn install
  
## To run:

      ./bin/cvtMonitoring  [options] file1 file2 ... fileN 
      
It requires the following banks to be present:
  - BMTRec::Hits, BSTRec::Hits
  - BMTRec::Clusters, BSTRec::Clusters
  - CVTRec::Seeds, CVTRec::Tracks
  - (for matching with CTOF) REC::Particle, REC::Track
  - (for MC only) MC::Particle
  
## Usage:

      Usage : cvtMonitoring  [options] file1 file2 ... fileN 

      Options :
       -histo : read histogram file (0/1) (default = 0)
           -n : maximum number of events to process (default = -1)
           -o : histogram file name prefix (default = )
         -pid : MC particle PID (default = 2212)
        -plot : display histograms (0/1) (default = 1)
       -stats : histogram stat option (e.g. "10" will display entries, default = )  

## Plots:
Plots are organized by tabs and subtabs.
- Tracks tab:
  - Tracks: parameters distribution for all tracks
  - Seeds: parameter distributions for all seeds and for track seeds (darker histograms) 
  - EBtracks: parameter distributions for all tracks and for tracks matched to CTOF (purple histograms)
- Cluster tab:
  - Size, energy and time distribution for all clusters (light blue), on-track clusters (green) and off-track clusters (darker blue)
- Hits tab:
  - Energy, time and residuals distribution for all hits (light blue), on-track hits (green) and off-track hits (darker blue)
- Residual tab:
  - Centroid residuals for SVT, BMT-C and BMT-Z clusters
- Pulls tab:
  - Error-normalized residuals
- MC:
  - MC and MC-Rec comparison. These are filled only for MC events where the MC:Particle bank is present. Only the MC particle with the selected PID is used. MC-Rec differences are computed only if the reconstructed particle matches the generated within 3 sigma of the resolution.
  - MC, Seed, Track: parameters distributions for the selected track type
  - (Seed)Resolution1,2,3: seed and track resolutions as a function of different track parameters.
  - (Seed)Pulls: helical track parameter pulls, calculated as difference between generated and reconstructed parameters, normalized to the corresponding uncertainty from the tracking covariance matrix.
  - Efficiency: generated (yellow), reconstructed seed (blue), reconstructed track (pink) distributions and corresponding efficiencies. The darker histograms show seed and tracks from the second seeding algorithm.
