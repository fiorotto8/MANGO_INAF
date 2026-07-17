# Sr-90 polar-angle acceptance

## Provenance

- Repository tag: `MANGO_INAF_v1.0`
- Repository commit: `1e8cc4b05db3314365d03662f6397e0f00987398`
- Geant4: 11.4.0 (`geant4-11-04`, 5 December 2025)
- ROOT: 6.36.08
- Physics: the repository's `QGSP_BIC_EMZ` reference list plus
  `G4RadioactiveDecayPhysics`
- Production macro: the unchanged `myMac.mac` from the checkout
- Production events: 3,000,000

The production macro contains
`/process/eLoss/StepFunction 0.001 0.01 mm` and
`/run/beamOn 3000000`. Its SHA-256 in both the working tree and copied build
directory was
`922146c78bcee754f5d681a4f71157312e6c4c9ad974e6add36f3ced496eee00`.
`git diff --ignore-space-at-eol --exit-code -- myMac.mac` returned zero. The
different SHA-256 of the LF-rendered Git blob,
`e2a829565c1ca99af3cbaa291dcb3c56d0d04ac14eb5e77bb6253c8b11fe2977`,
is solely due to the checkout's pre-existing CRLF line endings.

No detector, source, collimator, gas, physics-list, radioactive-decay,
step-function, containment, or retained simulation-selection definition was
changed. The gas remains He/CF4 60/40 at 1 atm.

## Definitions

Geant4's existing geometry has its cylinder/drift normal along global `y`.
To implement the requested manuscript convention without changing geometry,
all saved directions are mapped as

```text
(x, y, z)_analysis = (y, x, z)_Geant4
```

Thus `x` is the drift/GEM normal, `y-z` is the GEM/camera plane,
`theta = acos(ux)`, and `phi = atan2(uy,uz)`. Angles in the `Acceptance` tree
are radians; summary angles are degrees.

A direct beta electron is an electron created by `RadioactiveDecay` whose
parent ion is Sr-90 or Y-90. Each simulated decay-chain event has one of each,
so EventID alone is not a unique electron key; use `(EventID,TrackID)`.
Generated quantities are recorded at creation. Gas-entry energy and direction
are taken from the pre-step point of that electron's first step in the
sensitive gas, after any source/collimator scattering. Entry fields are `NaN`
when `GasEntered=0`.

`FullyContained` is a direct port of
`analysis/RecoTrack.C::areAllPointsInsideCylinder`: for every gas hit, with
`dx = G4 x`, `dz = G4 z - 51.4 mm`, the hit must have radius at most
`36.9-5.0 = 31.9 mm`, or `atan2(dx,dz)` normalized to `[0,360)` must lie in
the inclusive 150--210 degree sector. As in the retained code, `G4 y` is not
used by this cut.

`RecoAvailable` requires at least four non-degenerate sensitive-gas pre-step
positions. `RecoU*` is an orthogonal straight-line/PCA fit to the first four
positions projected into the requested `y-z` camera plane. The axis is
oriented from hit 1 toward hit 4. Consequently `RecoUx=0`,
`RecoTheta=90 degrees`, and the reconstructed azimuth is directional rather
than axial.

The retained repository simulation analysis first requires
`FullyContained` and then needs a reconstructed direction. Therefore

```text
Selected = FullyContained && RecoAvailable
```

No additional simulation-side event cut is present in this repository.
Accordingly, no untraceable selection was invented. Any selection beyond
containment plus availability of a direction cannot be reproduced from this
repository and is not represented by `Selected`.

## Energy bins and statistics

The manuscript bins were recovered from the retained directionality analysis
and are strict open intervals:

```text
(10,15), (15,20), (20,30), (30,40), (40,60) keV
```

The generated stage is binned in generated kinetic energy and uses the
generated direction. Gas-entry, contained, and selected stages are binned in
deposited energy, as in the retained manuscript analysis, and use the
direction at gas entry. This distinction is explicit in the CSV
`energy_definition` column and ROOT `Summary.EnergyDefinition` branch.

Percentiles use sorted-sample linear interpolation at
`p*(N-1)` (R/NumPy type 7). `rms_sin_theta` is the population RMS spread about
the mean, `sqrt(sum((sin(theta)-mean)^2)/N)`, matching the usual ROOT
distribution-RMS convention. All requested fractions use strict `<` tests.
Counts are given both as beta-electron records and distinct EventIDs.

## Production counts

Overall:

```text
Geant4 events                         3,000,000
Direct Sr-90/Y-90 beta records       6,000,000
GasEntered                              20,413
FullyContained                           6,293
RecoAvailable                           20,409
Selected                                 6,289
GasEntered without RecoAvailable             4
FullyContained without RecoAvailable          4
```

Counts in the manuscript bins (`electron records / distinct EventIDs`):

| Energy (keV) | Generated | Gas entry | Contained | Selected |
|---:|---:|---:|---:|---:|
| 10--15 | 45,296 / 45,260 | 7,005 / 6,925 | 1,480 / 1,473 | 1,480 / 1,473 |
| 15--20 | 45,014 / 44,974 | 3,189 / 3,177 | 713 / 712 | 713 / 712 |
| 20--30 | 93,574 / 93,331 | 2,262 / 2,255 | 720 / 718 | 720 / 718 |
| 30--40 | 94,567 / 94,247 | 944 / 943 | 465 / 465 | 465 / 465 |
| 40--60 | 189,870 / 188,638 | 978 / 977 | 585 / 584 | 585 / 584 |

The complete theta quantiles, `sin(theta)` mean/RMS, and all six requested
fractions for all four stages are in `sr90_acceptance.csv` and the ROOT
`Summary` tree.

## Output contents

- `sr90_acceptance.root`: `Acceptance` tree with 6,000,000 direct-beta rows;
  20-row `Summary` tree; definitions; count histograms; selected-stage
  quantile, RMS, and fraction graphs; and the plotted canvas.
- `sr90_acceptance.csv`: 20 summary rows (four stages times five bins).
- `sr90_acceptance.png`: plotted count and selected-angle summaries.
- `run.log`: exact commands and validation/production/analysis summaries,
  followed by the complete final production stdout/stderr and resource report.

Final artifact SHA-256 values:

```text
75efc513ed19aa772abb9f6f52f2f17e330301bf73edfc7d954c0d12952beb1d  sr90_acceptance.root
cf288cd621c2f343898f56ad6c154b08d6be295c2c0e1fc29045d62150502e94  sr90_acceptance.csv
137e4696bd7fcc847b4a2add9cd681bcef0fb87ea159b5df1f98cfccc7868455  sr90_acceptance.png
```

The `Acceptance` branches are:

```text
EventID TrackID ParentID
GeneratedEnergy_keV GeneratedUx GeneratedUy GeneratedUz
GeneratedTheta_rad GeneratedPhi_rad GeneratedSinTheta GeneratedAbsUx
GasEntered EntryEnergy_keV Ux Uy Uz Theta_rad Phi_rad SinTheta AbsUx
DepositedEnergy_keV HitCount FullyContained Selected RecoAvailable
RecoUx RecoUy RecoUz RecoTheta_rad RecoPhi_rad RecoSinTheta RecoAbsUx
```

## Exact commands

The repository's existing `build/` cache named a different source directory,
so a clean out-of-tree build was used.

```bash
export ROOTSYS=/home/fior/root_v6.36.08
export PATH="$ROOTSYS/bin:$PATH"
export LD_LIBRARY_PATH="$ROOTSYS/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
source /home/fior/geant4_v11.4.0/bin/geant4.sh

cmake -S . -B /tmp/mango_sr90_round2_1e8cc4b_build \
  -DGeant4_DIR=/home/fior/geant4_v11.4.0/lib/cmake/Geant4 \
  -DWITH_GEANT4_UIVIS=ON
cmake --build /tmp/mango_sr90_round2_1e8cc4b_build -j2

g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
  analysis/Sr90Acceptance.cpp \
  $(/home/fior/root_v6.36.08/bin/root-config --cflags --libs) \
  -o /tmp/mango_sr90_round2_1e8cc4b_build/sr90_acceptance

/tmp/mango_sr90_round2_1e8cc4b_build/sr90_acceptance --self-test

cd /tmp/mango_sr90_round2_1e8cc4b_build
mkdir -p outfiles
/usr/bin/time -v ./rdecay01 \
  /mnt/c/Users/david/MyDrive/WORK/POLARIMETRY/MANGO_INAF_Sr90+Am241/round2_validation.mac

/usr/bin/time -v ./rdecay01 myMac.mac

cd /mnt/c/Users/david/MyDrive/WORK/POLARIMETRY/MANGO_INAF_Sr90+Am241
/usr/bin/time -v \
  /tmp/mango_sr90_round2_1e8cc4b_build/sr90_acceptance \
  /tmp/mango_sr90_round2_1e8cc4b_build/outfiles/output_t0.root \
  round2_outputs
```

The angle self-test passed for `+x`, `-x`, `+y`, `-y`, `+z`, `-z`, and
`(1,1,1)/sqrt(3)`. The 1,000-event validation produced 2,000 direct-beta
records, 10 gas entries, 3 fully contained records, and 3 selected records;
its angle consistency warning count was zero. The final production took
13:58.77 wall time, used 178,332 kB maximum RSS, and exited with status zero.
The analysis took 42.71 seconds and used 750,288 kB maximum RSS.

Independent ROOT checks of all 20,413 gas-entry rows found zero failures at
`1e-10` for unit-vector magnitude, `theta=acos(ux)`,
`phi=atan2(uy,uz)`, `sin(theta)`, and `abs(ux)`. The analyzer also reported
zero angle-consistency warnings. The final CSV and PNG were byte-identical to
an earlier full-run analysis that filtered the same direct beta records.

## Warnings and repository state

- Geant4 warns that the pre-existing zero-density `Ar_gas` component is
  replaced by its minimum density, `1e-25 g/cm3`. The configured Ar fraction
  is zero; the detector mixture remains He/CF4 60/40 at 1 atm.
- The startup physics dump appears before macro execution and therefore shows
  the Geant4 default electron step function `(0.2,0.01 mm)`. The unchanged
  production macro then applies `0.001 0.01 mm` and Geant4 executes
  `/run/physicsModified` before the events.
- The pre-existing run-summary fields label the gun as `geantino`, report
  infinite geantino activity, and show an invalid minimum visible-energy
  value. The application swaps the gun to Sr-90 during event generation and
  its visible-energy accumulator is disabled. None of these summary fields is
  used here.
- Sourcing this ROOT installation's `thisroot.sh` under this shell emitted a
  null-byte/root-config warning. The exact ROOT paths were exported directly
  in the commands above; build, self-test, and analysis all succeeded.
- The checked-out worktree was already dirty before this work because many
  tracked files had LF-to-CRLF-only differences, and
  `analysis/convert`/`analysis/desktop.ini` were already untracked. Those
  unrelated changes were preserved.

Files modified for this work:

```text
include/EventAction.hh
src/EventAction.cc
src/RunAction.cc
src/SensitiveDetector.cc
src/TrackingAction.cc
analysis/Sr90Acceptance.cpp                  (new)
round2_validation.mac                        (new)
round2_outputs/README.md                     (new)
round2_outputs/run.log                       (new)
round2_outputs/sr90_acceptance.root           (new)
round2_outputs/sr90_acceptance.csv            (new)
round2_outputs/sr90_acceptance.png            (new)
```
