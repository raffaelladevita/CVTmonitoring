# CVTmonitoring

Displays CVT tracking information for beam data and MC

- To compile:

      mvn install
  
- To run:

      ./bin/cvtMonitoring
  
- Usage:

      Usage : cvtMonitoring 

      Options :
       -histo : read histogram file (0/1) (default = 0)
           -n : maximum number of events to process (default = -1)
           -o : histogram file name prefix (default = )
         -pid : MC particle PID (default = 2212)
        -plot : display histograms (0/1) (default = 1)
       -stats : histogram stat option (e.g. "10" will display entries, default = )  
