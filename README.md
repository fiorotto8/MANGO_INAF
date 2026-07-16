# MANGO at INAF simulation

Geant4 simulation of radioactive sources next to the MANGO gas volume.
The current supported source profiles are:

- `Am241.mac`: elemental Am-241 source, primary decay only, alpha hits.
- `Sr90.mac`: elemental Sr-90 source, full Sr-90 -> Y-90 -> Zr-90 chain,
  electron hits for the angular reconstruction.

`myMac.mac` runs the Am-241 profile by default.

## Build and run

```bash
cmake -S . -B build
cmake --build build -j
cd build
./rdecay01 Am241.mac
# or
./rdecay01 Sr90.mac
```

The simulation is intentionally restricted to one thread for now so that the
ROOT hit metadata and output remain consistent. Output files are written under
`build/output_files/`, for example `Am241_t0.root` and `Sr90_t0.root`.

The source macros configure:

- isotope atomic and mass number;
- elemental source material;
- whether daughter nuclei continue to decay;
- which particle is stored in the hit ntuple.

## Sr-90 electron reconstruction

Build the angular study program once:

```bash
cd analysis
g++ -O2 -o study Study.cpp $(root-config --cflags --libs)
```

Then run the complete Sr-90 reduction and angular study:

```bash
./RunAnal.sh ../build/output_files/Sr90_t0.root e-
```

The reducer creates one `elabHits` entry per direct radioactive-decay electron
track. Energies are stored in MeV, track lengths in mm, and the ion-pair count
uses a W-value of 38 eV.

For Am-241 track reduction without the electron angular study:

```bash
root -l -b -q 'RecoTrack.C("../build/output_files/Am241_t0.root","alpha",true)'
```

## Geant4 compatibility

The active physics list is the `QGSP_BIC_EMZ` reference list with radioactive
decay physics added in `rdecay01.cc`. The radioactive-decay time-threshold API
is guarded so the code can also build with Geant4 versions older than 11.2.
