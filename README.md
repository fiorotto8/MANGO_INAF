# MANGO at INAF simulation

## GEANT

- Working with geant4 11.1.0
  - following: <https://geant4.web.cern.ch/download/release-notes/notes-v11.2.0.html> There is a line of code in the main to add for v11.2 and next
    - `G4HadronicParameters::Instance()->SetTimeThresholdForRadioactiveDecay( 1.0e+60*CLHEP::year );`
  - comment it if you are using <=11.1.0
- Using `myMac.mac`
- The new version has the possibility to simulate calibration dataset with fixed energy gammas shot from the cathode
  - Comment the standard block in `PrimaryGenerator.cc` and uncomment the block for gamma calibration

## ANALYSIS

- root -l 'RecoTrack.C("path/to/output.root")'
- g++ -o study Study.cpp `root-config --cflags --libs` -lVc

If you use it as it is you can run `RunAnal.sh` (be sude to run `chmod 777 RunAnal.sh` before) when ROOT has finished you should `.q` it and wait for the study to occur

## Convert for Digitization

- Move the `output_t0.root` in the analysis folder
- g++ -o convert ConvertForDigi.cpp `root-config --cflags --libs`