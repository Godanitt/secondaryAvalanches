# secondaryAvalanches

Microscopic simulation of **primary and photon-induced secondary electron avalanches** in a uniform electric field using [Garfield++](https://garfieldpp.docs.cern.ch/), Magboltz and ROOT.

The project starts from microscopic electron avalanches, generates photons from the recorded collision populations, transports the relevant photons through the gas, applies gas photoabsorption/photoionisation and cathode quantum efficiency, and iterates the resulting feedback electrons over successive avalanche generations.

The project provides:

- a C++17 simulator, `secondaryAvalanches.cxx`;
- microscopic electron avalanches with `Garfield::AvalancheMicroscopic`;
- photon generation from the project kinetic/spectral model and CSV parameters;
- photon transport through Garfield++ with gas absorption and photoionisation;
- cathode photoelectron production using a wavelength-dependent QE model;
- iterative secondary-avalanche chains from cathode photoelectrons and gas-photoionisation electrons;
- fixed-field and target-gain campaigns driven by YAML;
- parallel campaign execution and a PySide6 campaign interface;
- ROOT outputs containing the primary avalanche, total feedback gain, avalanche genealogy, photon spectra, interaction positions and time structure;
- effectively infinite parallel electrodes in the transverse directions;
- a prompt-time view for the primary avalanche and a long-time view for delayed feedback.

GIF generation is not part of this project.

---

## Contents

1. [Physics workflow](#physics-workflow)
2. [Gain definitions](#gain-definitions)
3. [Geometry](#geometry)
4. [Repository structure](#repository-structure)
5. [Requirements](#requirements)
6. [Required Garfield++ interface](#required-garfield-interface)
7. [Installation](#installation)
8. [Quick start](#quick-start)
9. [Campaign configuration](#campaign-configuration)
10. [Photon and feedback options](#photon-and-feedback-options)
11. [Campaign modes](#campaign-modes)
12. [Output files](#output-files)
13. [ROOT contents](#root-contents)
14. [Understanding `mc_samples`](#understanding-mc_samples)
15. [Performance recommendations](#performance-recommendations)
16. [Troubleshooting](#troubleshooting)
17. [Scope and limitations](#scope-and-limitations)

---

## Physics workflow

For every primary electron, the simulation performs the following chain.

### 1. Primary microscopic avalanche

A primary electron is launched at the cathode plane and transported towards the anode with `Garfield::AvalancheMicroscopic`.

The simulation records:

- final electron and ion counts;
- electron endpoints and times;
- Magboltz collision levels;
- optional excitation positions and `z-t` correlations;
- the collision sites required by the photon-emission model.

The primary avalanche is assigned feedback generation

```text
0
```

### 2. Photon generation

`PhotonModel` converts the avalanche collision populations into kinetic emission components using the files under:

```text
data/parameters/
```

For each component, the code samples:

- emission position from the relevant microscopic collision sites;
- wavelength from the configured spectral model;
- emission delay from the component lifetime;
- the expected number of emitted photons.

The complete generated spectrum is always stored in `hSpectra`, including photons that are not subsequently transported.

### 3. Automatic photoelectric threshold

With

```yaml
propagate_only_above_phit: true
```

the effective photon-transport threshold is

```text
effective cut = max(photon_transport_cut_ev, PhiT)
```

where `PhiT` is obtained from the same QE material/model used to generate cathode photoelectrons.

Therefore, the simulation retains the full emitted spectrum but transports only photons satisfying

```text
E_gamma > PhiT
```

when the automatic threshold is enabled.

### 4. Photon propagation

With

```yaml
photo_absorption: true
```

photons are transported using Garfield++ optical cross sections through `AvalancheMicroscopic::TransportPhotonExternal`.

A transported photon may:

- be absorbed without ionising the gas;
- photoionise the gas and create one or more electron seeds;
- reach the cathode;
- reach the anode;
- leave the numerical optical volume;
- fail transport because no valid optical result is available.

With

```yaml
photo_absorption: false
```

the code uses isotropic ballistic propagation to the parallel-plane boundaries and does not simulate gas absorption.

### 5. Cathode photoelectrons

For a photon reaching the cathode, `QuantumEfficiency` evaluates the selected material/model at the photon wavelength.

A cathode electron is produced according to:

- the wavelength-dependent quantum efficiency;
- the configured electron-extraction efficiency;
- the photoelectron energy model.

Accepted cathode electrons are inserted into the avalanche queue as generation `g + 1`.

### 6. Gas-photoionisation electrons

When Garfield++ classifies an optical interaction as photoionisation, the generated gas electrons can be transported as new avalanche seeds when

```yaml
propagate_photoionisation_electrons: true
```

These electrons are also assigned generation `g + 1`.

### 7. Iterative feedback chain

The queue is processed until no seeds remain or a safety limit is reached:

```yaml
max_feedback_generations: 5
max_avalanches_per_primary: 10000
max_mc_photons_per_primary: 100000000
```

Seed types stored in `dataPerAvalanche` are:

```text
0 = primary electron
1 = cathode photoelectron
2 = gas-photoionisation electron
```

---

## Gain definitions

The simulation distinguishes the primary avalanche from the complete feedback chain.

For each initial electron:

```text
nePrimaryAvalanche
```

contains only the generation-0 avalanche, while

```text
ne
neTotalWithFeedback
```

contain the sum of all avalanches initiated by that primary electron.

The corresponding mean gains are therefore

```math
G_{\mathrm{primary}} = \left\langle n_{e,\,g=0} \right\rangle
```

and

```math
G_{\mathrm{total}} = \left\langle \sum_g n_{e,g} \right\rangle.
```

The campaign runner uses the `ne` branch of `dataPerPrimaryElectron`; consequently, the gain displayed by the campaign is the **total gain including feedback**. When no secondary avalanche occurs, primary and total gains are identical.

---

## Geometry

The uniform-field convention is:

```text
anode:   z = 0
cathode: z = gap
```

Primary electrons start at the cathode and drift towards the anode.

### Infinite parallel electrodes

For the uniform-field studies, the recommended configuration is:

```yaml
infinite_electrodes: true
optical_half_width_gaps: 100.0
```

In this mode:

- the cathode and anode are treated as physically infinite in `x/y`;
- there are no physical lateral walls in the analytic propagation;
- the finite Garfield++ transverse area is only a numerical safety volume;
- the numerical half-width is `optical_half_width_gaps × gap`.

For a `0.150 mm` gap and `optical_half_width_gaps: 100`, the optical volume has a half-width of `15 mm`.

A photon intersecting `z = gap` is classified as a cathode impact, and a photon intersecting `z = 0` is classified as an anode impact. Lateral losses should therefore be negligible except for numerical boundary cases in Garfield++ transport.

### Finite-electrode mode

When

```yaml
infinite_electrodes: false
```

`transport_half_size_cm` and `cathode_half_size_cm` define finite transverse limits. This mode should only be used when a finite optical aperture is physically intended.

---

## Repository structure

```text
.
├── CMakeLists.txt
├── secondaryAvalanches.cxx       # integrated avalanche/photon/feedback simulator
├── run_campaign.py               # YAML campaign runner
├── gui.py                        # PySide6 campaign interface
├── run_gui.sh                    # stable GUI launcher
├── campaign.yaml                 # active campaign configuration
├── alpha_model.py                # target-gain predictor
├── fit_diagnostics.py            # PDF fit diagnostics
├── requirements.txt
├── src/
│   └── aux/
│       ├── PhotonModel.hh
│       ├── PhotonModel.cxx
│       ├── QuantumEfficiency.hh
│       └── QuantumEfficiency.cxx
├── data/
│   ├── levels/                   # Magboltz-level mappings
│   ├── parameters/               # kinetic and spectral photon parameters
│   └── qe/
│       └── qe_materials.csv      # material-dependent QE data
├── build/                        # regenerated CMake build
├── fits/                         # target-gain fit PDFs
└── outputs/
    ├── roots/                    # production ROOT files
    ├── alpha/                    # machine-readable alpha fits
    ├── metadata/                 # campaign-attempt log
    └── legacy_roots/             # incompatible previous ROOT schemas
```

`uniformE.cxx` may remain in the repository as historical/reference code, but the production executable is built from `secondaryAvalanches.cxx`.

---

## Requirements

### C++ and scientific software

- Linux as the primary supported platform;
- a C++17 compiler;
- a Fortran compiler for Garfield++/Magboltz;
- CMake 3.16 or newer;
- ROOT with `Core`, `RIO`, `Tree`, `Hist` and `MathCore`;
- Garfield++ with Magboltz;
- GSL and the usual Garfield++ build dependencies.

### Python

Recommended: Python 3.10 or newer in a virtual environment.

Required packages:

```text
PySide6
numpy
uproot
PyYAML
matplotlib
```

`matplotlib` is only required for PDF fit diagnostics. The physical simulations can run without generating those PDFs.

---

## Required Garfield++ interface

The integrated photon study calls a public wrapper named

```cpp
TransportPhotonExternal(...)
```

which forwards to Garfield++'s internal `TransportPhoton` implementation.

Add the following method to:

```text
<garfield-source>/Include/Garfield/AvalancheMicroscopic.hh
```

inside the public section, immediately before `private:`:

```cpp
/// Transport a photon supplied externally.
void TransportPhotonExternal(const double x, const double y,
                             const double z, const double t,
                             const double e, const std::size_t w,
                             std::vector<Seed>& stack) {
  TransportPhoton(x, y, z, t, e, w, stack);
}
```

Keep the original private declaration unchanged:

```cpp
void TransportPhoton(const double x, const double y, const double z,
                     const double t, const double e, const std::size_t w,
                     std::vector<Seed>& stack);
```

Reinstall the modified header:

```bash
cd ~/garfield
cmake --install build --prefix ~/garfield/install
```

Verify both source and installed copies:

```bash
grep -n "TransportPhotonExternal" \
  ~/garfield/Include/Garfield/AvalancheMicroscopic.hh \
  ~/garfield/install/include/Garfield/AvalancheMicroscopic.hh
```

The production code must call:

```cpp
photon_transport.TransportPhotonExternal(...);
```

---

## Installation

### 1. Install basic tools

Example for Ubuntu/Debian:

```bash
sudo apt update
sudo apt install -y \
  build-essential \
  gfortran \
  cmake \
  libgsl-dev \
  python3 \
  python3-venv \
  python3-pip \
  git
```

Install ROOT and Garfield++ following their official documentation.

### 2. Load ROOT

```bash
source /path/to/root/bin/thisroot.sh
root-config --version
```

`run_gui.sh` also checks common locations such as:

```text
$HOME/root-install/bin/thisroot.sh
/opt/root/bin/thisroot.sh
```

### 3. Make Garfield++ discoverable

```bash
export CMAKE_PREFIX_PATH="$HOME/garfield/install:${CMAKE_PREFIX_PATH:-}"
export LD_LIBRARY_PATH="$HOME/garfield/install/lib:${LD_LIBRARY_PATH:-}"
```

### 4. Create the Python environment

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

### 5. Build

```bash
rm -rf build
cmake -S . -B build \
  -DCMAKE_PREFIX_PATH="$HOME/garfield/install"
cmake --build build -j2
```

The executable is:

```text
build/secondaryAvalanches
```

---

## Quick start

### Graphical interface

Launch from the repository root:

```bash
chmod +x run_gui.sh
./run_gui.sh
```

The interface provides:

- an editable YAML campaign;
- campaign open/save controls;
- excitation-position and Magboltz-transport switches;
- automatic CMake configuration and compilation;
- one table row per simulation;
- requested and measured electric field;
- measured total gain;
- primary-electron progress;
- status and error details;
- campaign stop control;
- access to the fit directory.

There is no GIF tab or GIF-generation workflow.

### Command line

```bash
source .venv/bin/activate
python3 run_campaign.py campaign.yaml
```

The runner rebuilds the executable unless explicitly instructed otherwise by its command-line options.

---

## Campaign configuration

The following is the project baseline for the current Ar–CF4 secondary-avalanche study:

```yaml
scan_mode: field
measure_gas_transport: false
record_excitation_positions: true
space_charge: false

# Integrated photon generation, Garfield++ transport and feedback chain.
parameters_dir: data/parameters
mc_samples: 10
random_seed: 12345
photo_absorption: true
# Keep the full spectrum, but only transport photons with E_gamma > PhiT.
# PhiT is read automatically from the selected QE material/model.
propagate_only_above_phit: true
photon_transport_cut_ev: 0.0

# Parallel plates are treated as effectively infinite in x/y. The finite
# Garfield++ area is only a numerical safety box, with a half-width of 100 gaps.
infinite_electrodes: true
optical_half_width_gaps: 100.0

# Keep the long optical window for delayed feedback, but use this prompt
# window for the main electron-endpoint histogram so the primary peak is visible.
prompt_time_max_ns: 10.0
qe_material: stainless_steel # Ti
qe_csv: data/qe/qe_materials.csv
qe_model: measured_extended
electron_extraction_efficiency: 1.0
propagate_photoionisation_electrons: true
max_feedback_generations: 5
max_avalanches_per_primary: 10000
max_mc_photons_per_primary: 100000000

field_sets:
  - components: {ar: 90, cf4: 10}
    pressures_bar: [1]
    gap_mm: 0.150
    reference_field_kv_cm: 39.8
    field_scales: [1.1, 1.0, 0.75]

  - components: {ar: 85, cf4: 15}
    pressures_bar: [1]
    gap_mm: 0.150
    reference_field_kv_cm: 43
    field_scales: [1.1, 1.0, 0.75]

  - components: {ar: 80, cf4: 20}
    pressures_bar: [1]
    gap_mm: 0.150
    reference_field_kv_cm: 44
    field_scales: [1.1, 1.0, 0.75]

npe: 5
workers: 2
```

For statistically meaningful secondary-avalanche distributions, increase `npe`. Five primaries are useful only as a fast functional test and can easily produce zero photoelectrons when the QE is small.

---

## Photon and feedback options

### `parameters_dir`

Directory containing the kinetic and spectral parameter files used by `PhotonModel`.

```yaml
parameters_dir: data/parameters
```

### `photo_absorption`

```yaml
photo_absorption: true
```

Uses Garfield++ optical transport and allows gas absorption and gas photoionisation.

```yaml
photo_absorption: false
```

Uses isotropic geometric propagation between the planes without gas optical interactions.

### `propagate_only_above_phit`

```yaml
propagate_only_above_phit: true
```

Uses the selected QE material/model to determine `PhiT` automatically. All generated photons remain in `hSpectra`, but only photons above the effective threshold are transported.

### `photon_transport_cut_ev`

A technical lower cut applied in addition to `PhiT`:

```yaml
photon_transport_cut_ev: 0.0
```

With the automatic `PhiT` option enabled, the effective cut is the larger of these two values.

### QE model

```yaml
qe_material: stainless_steel
qe_csv: data/qe/qe_materials.csv
qe_model: measured_extended
```

Available model modes are:

```text
measured_extended
measured_table
constant_threshold
```

`measured_extended` uses the measured table together with the project's extension outside the tabulated interval. `measured_table` restricts the response to the tabulated model. `constant_threshold` uses the configured fallback threshold and QE.

### Electron extraction efficiency

```yaml
electron_extraction_efficiency: 1.0
```

This multiplies the probability that a QE-accepted electron is extracted into the gas and available to start a secondary avalanche.

### Feedback generations

```yaml
max_feedback_generations: 5
```

Generation `0` is the primary avalanche. A seed produced by generation `g` creates generation `g + 1`.

### Gas-photoionisation feedback

```yaml
propagate_photoionisation_electrons: true
```

When enabled, electrons created by optical ionisation in the gas are inserted into the feedback queue.

### Safety limits

```yaml
max_avalanches_per_primary: 10000
max_mc_photons_per_primary: 100000000
```

These limits prevent runaway feedback or accidental memory/time explosions. They are safety constraints, not physical parameters.

---

## Campaign modes

### Fixed-field mode

```yaml
scan_mode: field
```

The runner simulates exactly the requested fields. It does not modify them to obtain a target gain.

A `field_set` accepts either explicit fields:

```yaml
fields_kv_cm: [30, 35, 40]
```

or a reference field and scale factors:

```yaml
reference_field_kv_cm: 39.8
field_scales: [1.1, 1.0, 0.75]
```

Fixed-field mode is the recommended mode for studying photon feedback at selected experimental operating points.

### Target-gain mode

```yaml
scan_mode: gain
```

The inherited campaign controller can use existing ROOT points and the alpha-model hierarchy to propose fields for requested gains. The measured gain is read from `dataPerPrimaryElectron/ne`, so it is the total gain including feedback.

Use target-gain mode only after validating the feedback configuration and collecting stable fixed-field statistics.

---

## Output files

Production ROOT files are stored under:

```text
outputs/roots/<mixture>/gap_<gap>mm/
```

Example:

```text
outputs/roots/ArCF4/gap_0.150mm/
```

Typical file name:

```text
ar_90.0_cf4_10.0_39.8kVcm_1.000bar_0.1500mm_100npe.root
```

The name contains:

- gas fractions;
- electric field in `kV/cm`;
- pressure in `bar`;
- gap in `mm`;
- actual number of primary electrons.

Incompatible ROOT schemas may be moved to:

```text
outputs/legacy_roots/
```

Target-gain fit state and diagnostics are stored under:

```text
outputs/alpha/
fits/
```

---

## ROOT contents

### `gasData`

One entry per ROOT file containing only the gas and run identity:

| Branch | Unit | Meaning |
|---|---:|---|
| `gas1` | — | first Magboltz gas identifier |
| `composition1_pct` | `%` | first-gas fraction |
| `gas2` | — | second Magboltz gas identifier |
| `composition2_pct` | `%` | second-gas fraction |
| `gas3` | — | optional third gas identifier |
| `composition3_pct` | `%` | optional third-gas fraction |
| `pressure_bar` | `bar` | gas pressure |
| `temperature_K` | `K` | gas temperature |
| `electricField_V_cm` | `V/cm` | uniform field magnitude |
| `gap_mm` | `mm` | multiplication gap |
| `npe` | — | actual simulated primary electrons |
| `randomSeed` | — | run random seed |

Physics results are intentionally not duplicated in `gasData`.

### `photonTransportData`

One summary entry containing only:

| Branch | Meaning |
|---|---|
| `nPhotoabsorptions` | total weighted gas photoabsorptions |
| `nPhotoionisations` | weighted gas-photoionisation subset |

Both branches are stored as floating-point values so weighted `mc_samples > 1` runs remain valid.

### `dataPerPrimaryElectron`

One entry per initial primary electron:

| Branch | Meaning |
|---|---|
| `ne` | total final electrons including every feedback generation |
| `ni` | total ions including every feedback generation |
| `neTotalWithFeedback` | explicit alias of total `ne` |
| `niTotalWithFeedback` | explicit alias of total `ni` |
| `nePrimaryAvalanche` | generation-0 electrons only |
| `niPrimaryAvalanche` | generation-0 ions only |
| `npe` | represented primaries, currently `1` per entry |
| `nAvalanches` | total avalanches in this primary's chain |
| `nPhotoelectrons` | cathode-photoelectron seeds that entered the chain |
| `nPhotoionisationSeeds` | gas-photoionisation seeds that entered the chain |
| `maxGeneration` | highest simulated feedback generation |

### `dataPerAvalanche`

One entry for every primary or feedback avalanche:

| Branch | Meaning |
|---|---|
| `primaryId` | parent initial-electron index |
| `avalancheId` | avalanche index within the primary chain |
| `generation` | feedback generation |
| `seedType` | primary, cathode photoelectron or gas photoionisation |
| `seedXcm`, `seedYcm`, `seedZcm` | seed position |
| `seedTimeNs` | seed time |
| `seedEnergyEv` | seed energy |
| `ne`, `ni` | avalanche electron and ion counts |
| `nPhotons` | weighted photons generated by this avalanche |
| `nCathodeImpacts` | weighted cathode impacts |
| `nPhotoelectrons` | weighted QE-accepted cathode electrons |
| `nPhotoionisationElectrons` | weighted electrons created by gas photoionisation |

This is the principal tree for reconstructing the secondary-avalanche chain.

### `dataPerElectron`

One entry per final electron endpoint:

| Branch | Meaning |
|---|---|
| `status` | Garfield++ endpoint status |
| `primaryId` | parent initial electron |
| `avalancheId` | parent avalanche |
| `generation` | feedback generation |

This tree can become large at high gain.

### Primary-avalanche histograms

```text
hElectronEnergyDistribution
hLevels
hExcXYZ                 # only when record_excitation_positions is true
hExcZT                  # only when record_excitation_positions is true
```

`hLevels` retains the historical project convention: valid non-elastic Magboltz levels from collision types 1–5, rather than only pure excitation-type collisions.

### Time structure

```text
hElectronsVsTime
hElectronsVsTimeFull
hAvalancheElectronsVsTime
hAvalancheElectronsVsTimePrompt
hAvalancheElectronsVsTimeGeneration
hFeedbackGeneration
```

- `hElectronsVsTime` contains electron endpoints in the fixed prompt window, normally `0–10 ns`.
- `hElectronsVsTimeFull` retains the complete delayed chain.
- `hAvalancheElectronsVsTime` places each avalanche at its seed time with weight equal to its electron yield.
- `hAvalancheElectronsVsTimeGeneration` separates avalanche strength by time and generation.

### Photon production and propagation

```text
hSpectra
hPhotonXYZ
hPhotonWavelengthTime
hCosTheta
hPhi
hQE
```

### Cathode response

```text
hImpactsXY
hPhotoElectronXY
hPhotoElectronTimeEnergy
```

### Gas optical interactions

```text
hPhotoAbsorptionXYZ
hPhotoAbsorptionZT
hPhotoIonisationXYZ
hPhotoIonisationZT
```

### ROOT drawing style

All one-dimensional ROOT histograms are stored with default option:

```text
HIST
```

They therefore open as stepped histograms in ROOT/JSROOT. Bin errors and `Sumw2` information remain available and may be drawn explicitly with:

```cpp
hist->Draw("E1");
```

Two- and three-dimensional histograms retain their normal ROOT drawing behaviour.

---

## Understanding `mc_samples`

`mc_samples` controls photon Monte Carlo oversampling.

For

```yaml
mc_samples: 10
```

approximately ten Monte Carlo photons are sampled per expected physical photon, and each sample carries weight

```text
1 / mc_samples = 0.1
```

Therefore:

- ROOT histogram `Entries` counts raw Monte Carlo samples;
- the histogram integral represents the weighted physical-equivalent count;
- `4000` impact entries at `mc_samples: 10` correspond to approximately `400` physical-equivalent impacts.

For the discrete feedback chain, a QE-accepted Monte Carlo photoelectron is retained as an actual avalanche seed with probability `1 / mc_samples`. The same unbiased downsampling is applied to gas-photoionisation electron seeds.

This prevents photon oversampling from artificially multiplying the number of secondary avalanches.

### Recommended use

For smooth optical distributions and geometric efficiencies:

```yaml
mc_samples: 10
```

or larger may be useful.

For event-by-event secondary-avalanche studies, where one impact should correspond directly to one physical trial:

```yaml
mc_samples: 1
```

is easier to interpret:

```text
one sampled photon = one physical photon
one impact entry   = one impact
one accepted seed  = one secondary avalanche seed
```

Do not remove the YAML key; set it to `1`.

---

## Performance recommendations

### Start with `mc_samples: 1` for feedback-chain debugging

This gives the cleanest event-by-event interpretation and avoids spending time oversampling distributions before the secondary-avalanche logic is validated.

### Increase `npe` for rare photoemission

With low-QE materials, `npe: 5` frequently produces zero photoelectrons even when thousands of weighted Monte Carlo impact entries are visible. Use at least tens or hundreds of primaries when studying secondary-avalanche statistics.

### Worker count

Garfield++ and Magboltz jobs can use substantial memory. Begin with:

```yaml
workers: 2
```

Reduce to `1` if the system swaps or becomes unresponsive.

### Keep Magboltz transport disabled unless needed

```yaml
measure_gas_transport: false
```

The microscopic avalanche already uses Magboltz collision data. This switch controls the additional macroscopic transport-table measurement, which can dominate the runtime of small tests.

### Excitation-position output

```yaml
record_excitation_positions: true
```

is required for `hExcXYZ` and `hExcZT`. Disable it only when those spatial diagnostics are not required.

### Feedback safety

A rapidly increasing chain can generate enormous numbers of photons and avalanches. Keep the generation, avalanche and photon safety limits enabled during development.

---

## Troubleshooting

### `AvalancheMicroscopic has no member TransportPhotonExternal`

The project is compiling against a Garfield++ header that does not contain the public wrapper.

Verify:

```bash
grep -n "TransportPhotonExternal" \
  ~/garfield/Include/Garfield/AvalancheMicroscopic.hh \
  ~/garfield/install/include/Garfield/AvalancheMicroscopic.hh
```

Then rebuild this project from a clean directory:

```bash
rm -rf build
cmake -S . -B build \
  -DCMAKE_PREFIX_PATH="$HOME/garfield/install"
cmake --build build -j2
```

### CMake finds the wrong Garfield++ installation

Inspect:

```bash
cmake --build build --verbose -j2
```

and confirm that include and library paths point to the intended installation prefix.

### Zero photoelectrons despite many impact entries

Check:

1. `mc_samples`: `Entries` are not physical counts when `mc_samples > 1`;
2. the integral of `hImpactsXY`, not only the displayed `Entries`;
3. the selected material QE in `hQE`;
4. the number of primaries;
5. `electron_extraction_efficiency`;
6. whether `E_gamma > PhiT` leaves sufficient transported photons.

For low-QE stainless steel and very small `npe`, zero photoelectrons can be the statistically expected result.

### Too many lateral photon losses

Use:

```yaml
infinite_electrodes: true
optical_half_width_gaps: 100.0
```

and regenerate the ROOT. Old files retain the previous finite geometry.

### Prompt-time histogram appears empty

Confirm that the ROOT was generated with the current prompt/full separation. `hElectronsVsTime` should cover only `prompt_time_max_ns`, while delayed events belong to `hElectronsVsTimeFull`.

### Old ROOT schema is reused

Delete or move the specific old ROOT and rerun the point. Existing ROOT files are not rewritten merely because the C++ schema has changed.

### PySide6 / Qt library conflict

Launch with:

```bash
./run_gui.sh
```

rather than directly running `python3 gui.py`. The launcher selects the Qt libraries bundled with PySide6 while preserving ROOT/Garfield++ library paths.

---

## Scope and limitations

The current project assumes:

- a uniform electric field;
- parallel anode and cathode planes;
- no imported Micromegas/GEM field map;
- no microscopic ion-drift simulation;
- optional charged-ring space charge rather than a complete plasma treatment;
- photon emission determined by the supplied kinetic model and parameter tables;
- optical transport determined by the Garfield++/Magboltz optical data available for the gas;
- cathode response determined by the selected QE table/model.

The effectively infinite transverse geometry is appropriate for studying intrinsic photon feedback in a uniform gap. A realistic detector mesh, hole structure or finite photocathode requires a non-uniform field and explicit detector geometry.

---

## Citation and software documentation

This project depends on:

- Garfield++: https://garfieldpp.docs.cern.ch/
- ROOT: https://root.cern/
- CMake: https://cmake.org/
- Qt for Python / PySide6: https://doc.qt.io/qtforpython-6/

When publishing results, document at least:

- Garfield++ and ROOT versions;
- gas composition and pressure;
- temperature;
- electric field and gap;
- QE material/model and extraction efficiency;
- `mc_samples`;
- photon-absorption mode;
- `PhiT` propagation setting;
- feedback-generation and safety limits;
- number of primary electrons;
- random seed.
