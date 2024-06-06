# MANGO at INAF simulation

## GEANT

- Working only with geant4 11.1.0 (no idea why)
  - following: <https://geant4.web.cern.ch/download/release-notes/notes-v11.2.0.html> There is a line of code in the main to add for v11.2 and next
    - `G4HadronicParameters::Instance()->SetTimeThresholdForRadioactiveDecay( 1.0e+60*CLHEP::year );`
  - comment it if you are using <=11.1.0
- Using `myMac.mac`

## ANALYSIS

- root -l 'RecoTrack.C("path/to/output.root")'
- g++ -o Study Study.cpp `root-config --cflags --libs` -lVc
