# CVTmonitoring

Displays CVT tracking information for beam data and MC

## To compile:

      mvn install
  
## To run:

      ./bin/cvtMonitoring  [options] file1 file2 ... fileN 
      
It requires the following banks to be present:
  - BMTRec::Hits, BSTRec::Hits
  - BMTRec::Clusters, BSTRec::Clusters
  - BMTRec::Crosses, BSTRec::Crosses
  - CVT::Seeds, CVTRec::Seeds, CVT::Tracks, CVTRec::UTracks, CVTRec::Tracks, CVTRec::Trajectory
  - (for matching with CTOF, elastic and doublets analysis) REC::Particle, REC::Track
  - (for MC only) MC::Particle
  
## Usage:

      Usage : cvtMonitoring  [options] file1 file2 ... fileN 

      Options :
        -beam : Beam energy (GeV) (default = 10.6)
     -cosmics : analyze as cosmics (0=false, 1=true) (default = 0)
       -histo : read histogram file (0/1) (default = 0)
        -lund : save events to lund (default = 0)
     -modules : comma-separated list of modules to be activated (default = )
           -n : maximum number of events to process (default = -1)
           -o : histogram file name prefix (default = )
         -pid : MC particle PID (default: use first particle in the bank) (default = 0)
        -plot : display histograms (0/1) (default = 1)
       -print : print histograms (0/1) (default = 0)
    -residual : residual scale (1=cm, 10=mm) (default = 1)
       -stats : histogram stat option (e.g. "10" will display entries) (default = )

## Plots (description to be updated):
Plots are organized by tabs and subtabs.
- Tracks tab:
  - Tracks: parameters distribution for all tracks
  - Seeds: parameter distributions for all seeds and for track seeds (darker histograms) 
  - EBtracks: parameter distributions for all tracks and for tracks matched to CTOF (purple histograms)
- Vertex tab:
  - Vertex analysis for positives and negative tracks. For straight tracks, only the positive subtab will be populated. The top row plots show the distance of closest approach to the beamline, d0, defined by the beam offset read during reconstruction from the CCDB table /geometry/beam/position, d0 vs. phi and vz. The bottom row show the vertex x vs. y, x vs. phi, y vs. phi and vz vs. phi. If no ```CVTRec::Tracks``` bank is present, but ```REC::Particle``` is, d0 will be defined as the distance of closets approach to the ideal beamline wth x=y=0.
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

## Beam Spot analysis

The plots in the Vertex tab and corresponding log can be used to check the global alignment of CVT with respect to the beam spot used in reconstruction. For this purpose, the CVT reconstruction services, both first and second pass, should be run turning off the beam spot constraint with the yaml setting:

     beamSpotConst: "0"
     
The presence of a phi modulation in the d0 vs. phi plots indicates an offset of the detector with respect to the beam spot read from the reconstructed banks. This modulation is fitted and analyzed to determine the offset and calculate the corresponding correction to be applied to the beam offset or, alternatively, to the detector position. 

The plot below shows examples of the results from a Spring 2019 run after internal CVT alignment was completed. The corresponfing code output from the log is:
```
Analyzing Vertex group: Positives
d0(phi) = p0 sin(p1 x + p2):
	 p0 = (-0.0948 +/- 0.0325)
	 p1 = (0.0174 +/- 0.0003)
	 p2 = (-10.6555 +/- 0.3069)
x_offset: (-0.316 +/- 0.295) mm, y_offset: (-0.893 +/- 0.321) mm
  with respect to beam spot read from banks: (0.000, 0.000) mm
Update the beam (x,y) position to: (-0.316, -0.893) mm
or shift the detector position by: (0.316, 0.893) mm

Analyzing Vertex group: Negatives
d0(phi) = p0 sin(p1 x + p2):
	 p0 = (0.1163 +/- 0.0270)
	 p1 = (0.0175 +/- 0.0002)
	 p2 = (-1.1563 +/- 0.2215)
x_offset: (-0.469 +/- 0.260) mm, y_offset: (-1.065 +/- 0.268) mm
  with respect to beam spot read from banks: (0.000, 0.000) mm
Update the beam (x,y) position to: (-0.469, -1.065) mm
or shift the detector position by: (0.469, 1.065) mm
```
Results from positive and negative tracks are expected to be compatible within the errors. The results for negative tracks have typiccaly better accuracy because of the Lorentz angle having a smaller impact on the cluster size and therefore the resolution.
#### Positives
![Plot_07-24-2022_10 25 25_PM](https://user-images.githubusercontent.com/7524926/180664539-4d34b854-1c59-4bb7-ac97-588a3b148cb6.png)

#### Negatives
![Plot_07-24-2022_10 22 20_PM](https://user-images.githubusercontent.com/7524926/180664526-f5852115-91ff-49bc-81cc-627783a1a526.png)
