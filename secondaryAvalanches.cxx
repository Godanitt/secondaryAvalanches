#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2I.h>
#include <TH3F.h>
#include <TH3I.h>
#include <TMath.h>
#include <TObject.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TTree.h>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentChargedRing.hh"
#include "Garfield/ComponentUser.hh"
#include "Garfield/Medium.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/RandomEngineSTL.hh"
#include "Garfield/Sensor.hh"

#include "PhotonModel.hh"
#include "QuantumEfficiency.hh"

using namespace Garfield;

namespace {

constexpr double kTorrPerBar = 750.061683;
constexpr double kLightSpeedCmPerNs = 29.9792458;
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr int kElasticType = 0;
constexpr int kIonisationType = 1;
constexpr int kAttachmentType = 2;
constexpr int kInelasticType = 3;
constexpr int kExcitationType = 4;
constexpr int kSuperelasticType = 5;
constexpr int kPhotonAbsorbedStatus = -2;
constexpr int kEnergyBins = 1000;
constexpr int kExcitationSpatialBins = 128;
constexpr int kExcitationTimeBins = 256;
constexpr Long64_t kMaxEnergySamples = 200000;
constexpr double kInitialTimeRangeNs = 10.0;
constexpr double kSpaceChargeToleranceCm = 1.0e-5;

// ============================================================================
//                              Configuration
// ============================================================================

int parse_int(const char* text, const char* name) {
  try {
    std::size_t used = 0;
    const int value = std::stoi(text, &used);
    if (used != std::string(text).size()) throw std::invalid_argument("suffix");
    return value;
  } catch (...) {
    throw std::runtime_error(std::string("Invalid integer for ") + name +
                             ": " + text);
  }
}

void use_hist_draw_option(TH1& hist) {
  // Store the ROOT-browser default as a stepped histogram. Sumw2/bin errors
  // remain available and can still be drawn explicitly with "E1".
  hist.SetOption("HIST");
}

double parse_double(const char* text, const char* name) {
  try {
    std::size_t used = 0;
    const double value = std::stod(text, &used);
    if (used != std::string(text).size()) throw std::invalid_argument("suffix");
    return value;
  } catch (...) {
    throw std::runtime_error(std::string("Invalid number for ") + name +
                             ": " + text);
  }
}

bool parse_on_off(const std::string& value, const std::string& option) {
  if (value == "on") return true;
  if (value == "off") return false;
  throw std::runtime_error(option + " must be on or off");
}

std::string value_after(int& index, const int argc, char** argv,
                        const std::string& option) {
  if (index + 1 >= argc) {
    throw std::runtime_error("Missing value for " + option);
  }
  return argv[++index];
}

struct Config {
  std::string root_file;
  std::string mixture_name;
  double field_v_cm = 0.0;
  double gap_mm = 0.0;
  double pressure_bar = 0.0;
  int min_npe = 10;
  int max_npe = 100;
  double target_relative_error = 0.03;
  std::string gas1;
  double composition1 = 0.0;
  std::string gas2;
  double composition2 = 0.0;
  std::string gas3;
  double composition3 = 0.0;
  double height_factor = 1.5;
  bool space_charge = false;
  bool record_excitation_positions = true;
  bool measure_gas_transport = false;
  int magboltz_collisions = 1;
  int job_id = 0;

  double temperature_k = 293.15;
  double initial_energy_ev = 0.1;
  double max_electron_energy_ev = 400.0;

  // Photon generation and transport.
  std::string parameters_dir = "data/parameters";
  long long mc_samples = 10;
  unsigned int seed = 12345;
  bool photo_absorption = true;
  double photon_transport_cut_ev = 0.0;
  bool propagate_only_above_phit = true;
  double wavelength_min_nm = 100.0;
  double wavelength_max_nm = 900.0;
  double wavelength_bin_nm = 1.0;
  int position_bins = 400;
  int angular_bins = 200;
  int time_bins = 1000;
  double time_max_ns = 0.0;
  // Dedicated prompt window. The full optical chain can extend to tens of us,
  // while microscopic electron avalanches typically finish within a few ns.
  double prompt_time_max_ns = 10.0;
  int electron_energy_bins = 300;
  double electron_energy_max_ev = 20.0;
  double cathode_half_size_cm = -1.0;
  double transport_half_size_cm = -1.0;
  // Uniform-field studies represent parallel plates. In this mode the x/y
  // limits are only a numerical Garfield++ boundary, not physical walls.
  bool infinite_electrodes = true;
  double optical_half_width_gaps = 100.0;
  photonfeedback::QeConfig qe;
  double electron_extraction_efficiency = 1.0;

  // Iterative avalanche chain.
  bool propagate_photoionisation_electrons = true;
  int max_feedback_generations = 5;
  int max_avalanches_per_primary = 10000;
  long long max_mc_photons_per_primary = 100000000;
  bool quiet = false;

  double pressure_torr() const { return pressure_bar * kTorrPerBar; }
  double gap_cm() const { return 0.1 * gap_mm; }
  double launch_z_cm() const { return gap_cm(); }
  double z_max_cm() const { return height_factor * gap_cm(); }
  double xy_half_width_cm() const { return 2.0 * gap_cm(); }
  double transport_half_width_cm() const {
    if (transport_half_size_cm > 0.0) return transport_half_size_cm;
    if (infinite_electrodes) return optical_half_width_gaps * gap_cm();
    const double cathode =
        cathode_half_size_cm > 0.0 ? cathode_half_size_cm : gap_cm();
    return std::max(xy_half_width_cm(), cathode);
  }
  double cathode_half_width_cm() const {
    if (cathode_half_size_cm > 0.0) return cathode_half_size_cm;
    return infinite_electrodes ? transport_half_width_cm() : gap_cm();
  }
};

void print_help(const char* executable) {
  std::cout
      << "Usage:\n  " << executable
      << " output.root mixture field(V/cm) gap(mm) pressure(bar) "
         "minNpe maxNpe targetRelativeError gas1 comp1 gas2 comp2 "
         "heightFactor spaceCharge legacyMakeGif legacyGifTmax "
         "legacyGifFrames jobId [legacyGifFile] [legacyGifMoveIons] "
         "[legacyGifIonSpeed] [recordExcitationPositions] "
         "[measureGasTransport] [gas3] [comp3] [magboltzCollisions] "
         "[photon options]\n\n"
      << "The three legacy GIF arguments are accepted only to keep old campaign "
         "commands compatible; GIF generation has been removed.\n\n"
      << "Photon options:\n"
      << "  --mc-samples N\n"
      << "  --seed N\n"
      << "  --parameters-dir DIR\n"
      << "  --photo-absorption on|off\n"
      << "  --photon-transport-cut-eV VALUE\n"
      << "  --propagate-only-above-PhiT on|off\n"
      << "  --qe-material NAME\n"
      << "  --qe-csv FILE\n"
      << "  --qe-model measured_extended|measured_table|constant_threshold\n"
      << "  --qe VALUE\n"
      << "  --qe-threshold-eV VALUE\n"
      << "  --work-function-eV VALUE\n"
      << "  --eee VALUE\n"
      << "  --cathode-half-size-cm VALUE\n"
      << "  --transport-half-size-cm VALUE\n"
      << "  --infinite-electrodes on|off\n"
      << "  --optical-half-width-gaps VALUE\n"
      << "  --prompt-time-max-ns VALUE\n"
      << "  --max-feedback-generations N\n"
      << "  --max-avalanches-per-primary N\n"
      << "  --propagate-photoionisation-electrons on|off\n";
}

Config read_config(int argc, char* argv[]) {
  if (argc < 19) {
    print_help(argv[0]);
    throw std::runtime_error("Not enough positional arguments");
  }

  Config c;
  c.root_file = argv[1];
  c.mixture_name = argv[2];
  c.field_v_cm = parse_double(argv[3], "field");
  c.gap_mm = parse_double(argv[4], "gap");
  c.pressure_bar = parse_double(argv[5], "pressure");
  c.min_npe = parse_int(argv[6], "minNpe");
  c.max_npe = parse_int(argv[7], "maxNpe");
  c.target_relative_error = parse_double(argv[8], "targetRelativeError");
  c.gas1 = argv[9];
  c.composition1 = parse_double(argv[10], "composition1");
  c.gas2 = argv[11];
  c.composition2 = parse_double(argv[12], "composition2");
  c.height_factor = parse_double(argv[13], "heightFactor");
  c.space_charge = parse_int(argv[14], "spaceCharge") != 0;

  const bool legacy_make_gif = parse_int(argv[15], "legacyMakeGif") != 0;
  if (legacy_make_gif) {
    throw std::runtime_error(
        "GIF generation was deliberately removed from secondaryAvalanches");
  }
  // argv[16] and argv[17] are deliberately ignored legacy GIF values.
  c.job_id = parse_int(argv[18], "jobId");

  int option_start = argc;
  for (int i = 19; i < argc; ++i) {
    if (std::string(argv[i]).rfind("--", 0) == 0) {
      option_start = i;
      break;
    }
  }

  // Keep the old run_campaign.py positional layout usable.
  if (option_start > 22) {
    c.record_excitation_positions =
        parse_int(argv[22], "recordExcitationPositions") != 0;
  }
  if (option_start > 23) {
    c.measure_gas_transport =
        parse_int(argv[23], "measureGasTransport") != 0;
  }
  if (option_start > 24) c.gas3 = argv[24];
  if (option_start > 25) c.composition3 = parse_double(argv[25], "composition3");
  if (option_start > 26) {
    c.magboltz_collisions = parse_int(argv[26], "magboltzCollisions");
  }

  for (int i = option_start; i < argc; ++i) {
    const std::string option = argv[i];
    if (option == "--help" || option == "-h") {
      print_help(argv[0]);
      std::exit(EXIT_SUCCESS);
    } else if (option == "--mc-samples") {
      c.mc_samples = std::stoll(value_after(i, argc, argv, option));
    } else if (option == "--seed") {
      c.seed = static_cast<unsigned int>(
          std::stoul(value_after(i, argc, argv, option)));
    } else if (option == "--parameters-dir" || option == "--parameter-dir") {
      c.parameters_dir = value_after(i, argc, argv, option);
    } else if (option == "--photo-absorption") {
      c.photo_absorption =
          parse_on_off(value_after(i, argc, argv, option), option);
    } else if (option == "--photon-transport-cut-eV" ||
               option == "--photon-transport-cut-ev") {
      c.photon_transport_cut_ev =
          std::stod(value_after(i, argc, argv, option));
    } else if (option == "--propagate-only-above-PhiT" ||
               option == "--propagate-only-above-phit") {
      c.propagate_only_above_phit =
          parse_on_off(value_after(i, argc, argv, option), option);
    } else if (option == "--qe-material") {
      c.qe.material = value_after(i, argc, argv, option);
    } else if (option == "--qe-csv") {
      c.qe.csv_path = value_after(i, argc, argv, option);
    } else if (option == "--qe-model") {
      c.qe.model_mode = value_after(i, argc, argv, option);
    } else if (option == "--qe") {
      c.qe.fallback_qe = std::stod(value_after(i, argc, argv, option));
    } else if (option == "--qe-threshold-eV" ||
               option == "--qe-threshold-ev") {
      c.qe.fallback_threshold_ev =
          std::stod(value_after(i, argc, argv, option));
    } else if (option == "--work-function-eV" ||
               option == "--work-function-ev") {
      c.qe.work_function_override_ev =
          std::stod(value_after(i, argc, argv, option));
    } else if (option == "--eee") {
      c.electron_extraction_efficiency =
          std::stod(value_after(i, argc, argv, option));
    } else if (option == "--wavelength-min") {
      c.wavelength_min_nm = std::stod(value_after(i, argc, argv, option));
    } else if (option == "--wavelength-max") {
      c.wavelength_max_nm = std::stod(value_after(i, argc, argv, option));
    } else if (option == "--wavelength-bin") {
      c.wavelength_bin_nm = std::stod(value_after(i, argc, argv, option));
    } else if (option == "--position-bins") {
      c.position_bins = std::stoi(value_after(i, argc, argv, option));
    } else if (option == "--angular-bins") {
      c.angular_bins = std::stoi(value_after(i, argc, argv, option));
    } else if (option == "--time-bins") {
      c.time_bins = std::stoi(value_after(i, argc, argv, option));
    } else if (option == "--time-max-ns") {
      c.time_max_ns = std::stod(value_after(i, argc, argv, option));
    } else if (option == "--prompt-time-max-ns") {
      c.prompt_time_max_ns = std::stod(value_after(i, argc, argv, option));
    } else if (option == "--electron-energy-bins") {
      c.electron_energy_bins = std::stoi(value_after(i, argc, argv, option));
    } else if (option == "--electron-energy-max-eV" ||
               option == "--electron-energy-max-ev") {
      c.electron_energy_max_ev =
          std::stod(value_after(i, argc, argv, option));
    } else if (option == "--cathode-half-size-cm") {
      c.cathode_half_size_cm =
          std::stod(value_after(i, argc, argv, option));
    } else if (option == "--transport-half-size-cm") {
      c.transport_half_size_cm =
          std::stod(value_after(i, argc, argv, option));
    } else if (option == "--infinite-electrodes") {
      c.infinite_electrodes =
          parse_on_off(value_after(i, argc, argv, option), option);
    } else if (option == "--optical-half-width-gaps") {
      c.optical_half_width_gaps =
          std::stod(value_after(i, argc, argv, option));
    } else if (option == "--max-feedback-generations") {
      c.max_feedback_generations =
          std::stoi(value_after(i, argc, argv, option));
    } else if (option == "--max-avalanches-per-primary") {
      c.max_avalanches_per_primary =
          std::stoi(value_after(i, argc, argv, option));
    } else if (option == "--max-mc-photons-per-primary") {
      c.max_mc_photons_per_primary =
          std::stoll(value_after(i, argc, argv, option));
    } else if (option == "--propagate-photoionisation-electrons") {
      c.propagate_photoionisation_electrons =
          parse_on_off(value_after(i, argc, argv, option), option);
    } else if (option == "--quiet") {
      c.quiet = true;
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }

  if (c.field_v_cm <= 0.0) throw std::runtime_error("field must be positive");
  if (c.gap_mm <= 0.0) throw std::runtime_error("gap must be positive");
  if (c.pressure_bar <= 0.0) throw std::runtime_error("pressure must be positive");
  if (c.min_npe <= 0 || c.max_npe < c.min_npe) {
    throw std::runtime_error("Require 0 < minNpe <= maxNpe");
  }
  if (c.target_relative_error < 0.0) {
    throw std::runtime_error("targetRelativeError cannot be negative");
  }
  if (c.height_factor < 1.0) {
    throw std::runtime_error("heightFactor must be at least 1");
  }
  if (c.magboltz_collisions < 1) {
    throw std::runtime_error("magboltzCollisions must be at least 1");
  }
  if (c.composition1 < 0.0 || c.composition2 < 0.0 ||
      c.composition3 < 0.0) {
    throw std::runtime_error("Gas compositions cannot be negative");
  }
  const double composition_sum =
      c.composition1 + c.composition2 + c.composition3;
  if (std::abs(composition_sum - 100.0) > 1.0e-6) {
    throw std::runtime_error("Gas compositions must add to 100 percent");
  }
  if (c.mc_samples <= 0) throw std::runtime_error("mcSamples must be positive");
  if (c.max_feedback_generations < 0 || c.max_avalanches_per_primary < 1) {
    throw std::runtime_error("Invalid feedback-chain limits");
  }
  if (c.max_mc_photons_per_primary < 1) {
    throw std::runtime_error("max MC photons per primary must be positive");
  }
  if (!std::isfinite(c.electron_extraction_efficiency) ||
      c.electron_extraction_efficiency < 0.0 ||
      c.electron_extraction_efficiency > 1.0) {
    throw std::runtime_error("--eee must be between 0 and 1");
  }
  if (!std::isfinite(c.photon_transport_cut_ev) ||
      c.photon_transport_cut_ev < 0.0) {
    throw std::runtime_error("Photon transport cut must be non-negative");
  }
  if (c.wavelength_max_nm <= c.wavelength_min_nm ||
      c.wavelength_bin_nm <= 0.0) {
    throw std::runtime_error("Invalid wavelength range");
  }
  if (!std::isfinite(c.prompt_time_max_ns) || c.prompt_time_max_ns <= 0.0) {
    throw std::runtime_error("Prompt time maximum must be positive");
  }
  if (!std::isfinite(c.optical_half_width_gaps) ||
      c.optical_half_width_gaps < 1.0) {
    throw std::runtime_error("Optical half width must be at least one gap");
  }
  if (c.transport_half_size_cm > 0.0 &&
      c.transport_half_size_cm + 1.0e-12 < c.cathode_half_width_cm()) {
    throw std::runtime_error(
        "transport half size cannot be smaller than the cathode");
  }
  return c;
}

photonemission::GasData make_model_gas(const Config& config) {
  photonemission::GasData gas;
  gas.gas1 = config.gas1;
  gas.gas2 = config.gas2;
  gas.gas3 = config.gas3;
  gas.composition1_percent = config.composition1;
  gas.composition2_percent = config.composition2;
  gas.composition3_percent = config.composition3;
  gas.pressure_bar = config.pressure_bar;
  gas.temperature_k = config.temperature_k;
  gas.electric_field_v_cm = config.field_v_cm;
  gas.gap_mm = config.gap_mm;
  gas.start_z_fraction_from_anode = 1.0;
  gas.start_z_mm = config.gap_mm;
  gas.avalanche_distance_mm = config.gap_mm;
  gas.space_charge_enabled = config.space_charge;
  return gas;
}

std::vector<photonemission::LevelInfo> read_level_info(
    const Config& config, MediumMagboltz& gas) {
  const int n_levels = std::max(0, static_cast<int>(gas.GetNumberOfLevels()));
  std::vector<int> active_slots;
  if (config.composition1 > 0.0) active_slots.push_back(0);
  if (config.composition2 > 0.0) active_slots.push_back(1);
  if (config.composition3 > 0.0) active_slots.push_back(2);

  std::vector<photonemission::LevelInfo> levels;
  levels.reserve(static_cast<std::size_t>(n_levels));
  for (int level = 0; level < n_levels; ++level) {
    int gas_index = -1;
    int process_type = -1;
    std::string description;
    double energy_ev = std::numeric_limits<double>::quiet_NaN();
    gas.GetLevel(static_cast<unsigned int>(level), gas_index, process_type,
                 description, energy_ev);

    photonemission::LevelInfo info;
    info.level = level;
    info.state_name = description;
    info.energy_ev = energy_ev;
    if (process_type == kElasticType) info.type = "elastic";
    if (process_type == kIonisationType) info.type = "ionisation";
    if (process_type == kAttachmentType) info.type = "attachment";
    if (process_type == kInelasticType) info.type = "inelastic";
    if (process_type == kExcitationType) info.type = "excitation";
    if (process_type == kSuperelasticType) info.type = "superelastic";

    if (gas_index >= 0 && gas_index < static_cast<int>(active_slots.size())) {
      const int slot = active_slots[static_cast<std::size_t>(gas_index)];
      info.gas = slot == 0 ? config.gas1 : (slot == 1 ? config.gas2 : config.gas3);
    }
    levels.push_back(std::move(info));
  }
  return levels;
}

// ============================================================================
//                    Electron-energy and collision callbacks
// ============================================================================

struct EnergyReservoir {
  Long64_t seen = 0;
  std::vector<float> values;

  void reset() {
    seen = 0;
    values.clear();
    values.reserve(static_cast<std::size_t>(kMaxEnergySamples));
  }

  void add(const double energy_ev, const bool hole) {
    if (hole || !std::isfinite(energy_ev) || energy_ev < 0.0) return;
    ++seen;
    if (static_cast<Long64_t>(values.size()) < kMaxEnergySamples) {
      values.push_back(static_cast<float>(energy_ev));
      return;
    }
    const Long64_t index = static_cast<Long64_t>(
        std::floor(gRandom->Uniform(0.0, static_cast<double>(seen))));
    if (index < kMaxEnergySamples) {
      values[static_cast<std::size_t>(index)] = static_cast<float>(energy_ev);
    }
  }

  double random_energy(const double fallback) const {
    if (values.empty()) return fallback;
    return values[static_cast<std::size_t>(gRandom->Integer(values.size()))];
  }
};

EnergyReservoir g_energy;

void handle_step(double, double, double, double, double energy,
                 double, double, double, bool hole) {
  g_energy.add(energy, hole);
}

struct CollisionRecorder {
  TH1D* h_levels = nullptr;
  TH3I* h_xyz = nullptr;
  TH2I* h_zt = nullptr;
  bool keep_ions = false;

  Long64_t n_non_elastic = 0;
  Long64_t n_inelastic_type3 = 0;
  Long64_t n_excitation_type4 = 0;
  Long64_t n_excitation_like = 0;
  Long64_t n_ionisations = 0;
  Long64_t n_attachments = 0;
  Long64_t n_superelastic = 0;

  std::vector<photonemission::PhotonSourceSite> sites_this_avalanche;
  std::vector<std::array<double, 4>> ions_this_avalanche;

  void reset(TH1D& levels, TH3I* xyz, TH2I* zt, const bool store_ions) {
    h_levels = &levels;
    h_xyz = xyz;
    h_zt = zt;
    keep_ions = store_ions;
    n_non_elastic = 0;
    n_inelastic_type3 = 0;
    n_excitation_type4 = 0;
    n_excitation_like = 0;
    n_ionisations = 0;
    n_attachments = 0;
    n_superelastic = 0;
    sites_this_avalanche.clear();
    ions_this_avalanche.clear();
  }

  void start_avalanche() {
    sites_this_avalanche.clear();
    ions_this_avalanche.clear();
  }

  void add(double x, double y, double z, double t, int type, int level) {
    if (type == kIonisationType) {
      ++n_ionisations;
      if (keep_ions && std::isfinite(x) && std::isfinite(y) &&
          std::isfinite(z)) {
        ions_this_avalanche.push_back({x, y, z, t});
      }
    } else if (type == kAttachmentType) {
      ++n_attachments;
    } else if (type == kSuperelasticType) {
      ++n_superelastic;
    }

    if (type != kElasticType && level >= 0) {
      ++n_non_elastic;
      if (h_levels != nullptr) h_levels->Fill(level);

      // Photon-producing channels are not limited to excitation collisions.
      // CF4+* and the Ar third continuum are sourced by ionisation levels, so
      // every resolved non-elastic collision keeps its own x-y-z-t site for
      // the kinetic model. Diagnostic hExcXYZ/hExcZT remain excitation-only.
      sites_this_avalanche.push_back({level, x, y, z, t});
    }
    if (type == kInelasticType && level >= 0) ++n_inelastic_type3;
    if (type == kExcitationType && level >= 0) ++n_excitation_type4;

    if ((type == kInelasticType || type == kExcitationType) && level >= 0) {
      ++n_excitation_like;
      if (h_xyz != nullptr && std::isfinite(x) && std::isfinite(y) &&
          std::isfinite(z)) {
        h_xyz->Fill(x, y, z);
      }
      if (h_zt != nullptr && std::isfinite(z) && std::isfinite(t)) {
        h_zt->Fill(z, t);
      }
    }
  }
};

CollisionRecorder g_collisions;

void handle_collision(double x, double y, double z, double t,
                      int type, int level, Medium*, double, double, double,
                      double, double, double, double, double) {
  g_collisions.add(x, y, z, t, type, level);
}

// ============================================================================
//                        Statistics and gas transport
// ============================================================================

struct RunningStatistics {
  int n = 0;
  double mean = 0.0;
  double m2 = 0.0;

  void add(const double value) {
    ++n;
    const double delta = value - mean;
    mean += delta / n;
    m2 += delta * (value - mean);
  }
  double standard_deviation() const {
    return n > 1 ? std::sqrt(m2 / (n - 1)) : 0.0;
  }
  double error_on_mean() const {
    return n > 1 ? standard_deviation() / std::sqrt(static_cast<double>(n))
                 : std::numeric_limits<double>::infinity();
  }
  double relative_error() const {
    return mean > 0.0 ? error_on_mean() / mean
                      : std::numeric_limits<double>::infinity();
  }
};

struct GasTransport {
  double vz_cm_ns = std::numeric_limits<double>::quiet_NaN();
  double longitudinal_diffusion = std::numeric_limits<double>::quiet_NaN();
  double transverse_diffusion = std::numeric_limits<double>::quiet_NaN();
  double townsend_cm_inv = std::numeric_limits<double>::quiet_NaN();
  double attachment_cm_inv = std::numeric_limits<double>::quiet_NaN();
  double runtime_seconds = std::numeric_limits<double>::quiet_NaN();

  void measure(const Config& config, MediumMagboltz& gas) {
    if (!config.measure_gas_transport) return;
    gas.SetFieldGrid(config.field_v_cm, config.field_v_cm, 1, false);
    std::cout << "MAGBOLTZ_START " << config.job_id << " "
              << config.field_v_cm << std::endl;
    const auto start = std::chrono::steady_clock::now();
    gas.GenerateGasTable(config.magboltz_collisions, false);
    runtime_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    std::cout << "MAGBOLTZ_DONE " << config.job_id << " "
              << runtime_seconds << std::endl;

    double vx = 0.0, vy = 0.0;
    const double ex = 0.0, ey = 0.0, ez = config.field_v_cm;
    const double bx = 0.0, by = 0.0, bz = 0.0;
    if (!gas.ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz_cm_ns) ||
        !gas.ElectronDiffusion(ex, ey, ez, bx, by, bz,
                               longitudinal_diffusion,
                               transverse_diffusion) ||
        !gas.ElectronTownsend(ex, ey, ez, bx, by, bz, townsend_cm_inv) ||
        !gas.ElectronAttachment(ex, ey, ez, bx, by, bz,
                                attachment_cm_inv)) {
      throw std::runtime_error(
          "Could not read the requested Magboltz transport values");
    }
  }
};

struct SpaceCharge {
  std::unique_ptr<ComponentChargedRing> rings;
  Long64_t n_ions = 0;

  void initialise(const Config& config, MediumMagboltz& gas) {
    if (!config.space_charge) return;
    rings = std::make_unique<ComponentChargedRing>();
    rings->SetArea(-config.xy_half_width_cm(), -config.xy_half_width_cm(), 0.0,
                   config.xy_half_width_cm(), config.xy_half_width_cm(),
                   config.z_max_cm());
    rings->SetSpacingTolerance(kSpaceChargeToleranceCm);
    rings->SetMedium(&gas);
    rings->ClearActiveRings();
    rings->UpdateCentre(0.0, 0.0);
  }

  void add_new_ions(const Config& config) {
    if (rings == nullptr) return;
    for (const auto& ion : g_collisions.ions_this_avalanche) {
      if (std::abs(ion[0]) > config.xy_half_width_cm() ||
          std::abs(ion[1]) > config.xy_half_width_cm() || ion[2] < 0.0 ||
          ion[2] > config.z_max_cm()) {
        continue;
      }
      rings->AddChargedRing(ion[0], ion[1], ion[2], +1.0);
      ++n_ions;
    }
  }
};

// ============================================================================
//                           Photon propagation
// ============================================================================

enum class PhotonTransportProcess {
  Invalid,
  Transparent,
  Photoabsorption,
  Photoionisation,
  CathodeImpact,
  AnodeImpact,
  LateralEscape
};

struct PhotoionisationElectron {
  double x_cm = 0.0;
  double y_cm = 0.0;
  double z_cm = 0.0;
  double time_ns = 0.0;
  double energy_ev = 0.0;
  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;
  std::size_t multiplicity = 1;
};

struct PhotonResult {
  PhotonTransportProcess process = PhotonTransportProcess::Invalid;
  double energy_ev = 0.0;
  double wavelength_nm = 0.0;
  double x0_cm = 0.0;
  double y0_cm = 0.0;
  double z0_cm = 0.0;
  double t0_ns = 0.0;
  double x1_cm = 0.0;
  double y1_cm = 0.0;
  double z1_cm = 0.0;
  double t1_ns = 0.0;
  int status = 0;
  bool used_garfield = false;
  std::vector<PhotoionisationElectron> electrons;
};

struct PhotonGeometry {
  double xmin_cm = -0.2;
  double xmax_cm = 0.2;
  double ymin_cm = -0.2;
  double ymax_cm = 0.2;
  double gap_cm = 0.1;
  double cathode_half_size_cm = 0.1;
  bool infinite_electrodes = false;
};

struct PhotonDirection {
  double ux = 0.0;
  double uy = 0.0;
  double uz = 1.0;
  double cos_theta = 1.0;
  double phi = 0.0;
};

enum class BoundarySurface { Invalid, Cathode, Anode, Lateral };

struct BoundaryExit {
  bool valid = false;
  BoundarySurface surface = BoundarySurface::Invalid;
  double x_cm = 0.0;
  double y_cm = 0.0;
  double z_cm = 0.0;
  double distance_cm = 0.0;
};

PhotonDirection sample_isotropic_direction(TRandom3& random) {
  PhotonDirection direction;
  direction.cos_theta = random.Uniform(-1.0, 1.0);
  direction.phi = random.Uniform(-kPi, kPi);
  const double sin_theta = std::sqrt(
      std::max(0.0, 1.0 - direction.cos_theta * direction.cos_theta));
  direction.ux = sin_theta * std::cos(direction.phi);
  direction.uy = sin_theta * std::sin(direction.phi);
  direction.uz = direction.cos_theta;
  return direction;
}

BoundaryExit propagate_to_box_boundary(
    const double x_cm, const double y_cm, const double z_cm,
    const PhotonDirection& direction, const PhotonGeometry& geometry) {
  BoundaryExit exit;
  constexpr double direction_tolerance = 1.0e-15;
  constexpr double distance_tolerance = 1.0e-12;

  struct Candidate {
    double distance = std::numeric_limits<double>::infinity();
    BoundarySurface surface = BoundarySurface::Invalid;
  };
  std::array<Candidate, 6> candidates{};
  std::size_t count = 0;
  auto add = [&](const double distance, const BoundarySurface surface) {
    if (std::isfinite(distance) && distance >= -distance_tolerance) {
      candidates[count++] = {std::max(0.0, distance), surface};
    }
  };

  if (!geometry.infinite_electrodes) {
    if (direction.ux > direction_tolerance) {
      add((geometry.xmax_cm - x_cm) / direction.ux,
          BoundarySurface::Lateral);
    } else if (direction.ux < -direction_tolerance) {
      add((geometry.xmin_cm - x_cm) / direction.ux,
          BoundarySurface::Lateral);
    }
    if (direction.uy > direction_tolerance) {
      add((geometry.ymax_cm - y_cm) / direction.uy,
          BoundarySurface::Lateral);
    } else if (direction.uy < -direction_tolerance) {
      add((geometry.ymin_cm - y_cm) / direction.uy,
          BoundarySurface::Lateral);
    }
  }
  if (direction.uz > direction_tolerance) {
    add((geometry.gap_cm - z_cm) / direction.uz,
        BoundarySurface::Cathode);
  } else if (direction.uz < -direction_tolerance) {
    add((0.0 - z_cm) / direction.uz, BoundarySurface::Anode);
  }
  if (count == 0) return exit;

  const auto best = std::min_element(
      candidates.begin(), candidates.begin() + static_cast<std::ptrdiff_t>(count),
      [](const Candidate& left, const Candidate& right) {
        return left.distance < right.distance;
      });
  if (!std::isfinite(best->distance)) return exit;

  exit.valid = true;
  exit.distance_cm = best->distance;
  exit.x_cm = x_cm + best->distance * direction.ux;
  exit.y_cm = y_cm + best->distance * direction.uy;
  exit.z_cm = z_cm + best->distance * direction.uz;
  exit.surface = best->surface;
  if (!geometry.infinite_electrodes &&
      exit.surface == BoundarySurface::Cathode &&
      (std::abs(exit.x_cm) > geometry.cathode_half_size_cm ||
       std::abs(exit.y_cm) > geometry.cathode_half_size_cm)) {
    exit.surface = BoundarySurface::Lateral;
  }
  return exit;
}

class PhotonPropagation {
 public:
  PhotonPropagation(const PhotonGeometry& geometry, MediumMagboltz& medium,
                    const double transport_cut_ev)
      : geometry_(geometry), cut_ev_(transport_cut_ev), medium_(&medium) {
    component_.SetElectricField(
        [](double, double, double, double& ex, double& ey, double& ez) {
          ex = 0.0;
          ey = 0.0;
          ez = 0.0;
        });
    component_.SetArea(geometry_.xmin_cm, geometry_.ymin_cm, 0.0,
                       geometry_.xmax_cm, geometry_.ymax_cm,
                       geometry_.gap_cm);
    component_.SetMedium(medium_);
    sensor_.AddComponent(&component_);
    sensor_.SetArea(geometry_.xmin_cm, geometry_.ymin_cm, 0.0,
                    geometry_.xmax_cm, geometry_.ymax_cm,
                    geometry_.gap_cm);
  }

  PhotonResult transport(double x_cm, double y_cm, double z_cm,
                         double time_ns, double wavelength_nm) {
    PhotonResult result;
    result.wavelength_nm = wavelength_nm;
    result.energy_ev = photonemission::photon_energy_ev(wavelength_nm);
    result.x0_cm = result.x1_cm = x_cm;
    result.y0_cm = result.y1_cm = y_cm;
    result.z0_cm = result.z1_cm = z_cm;
    result.t0_ns = result.t1_ns = time_ns;

    if (!std::isfinite(result.energy_ev) || result.energy_ev <= 0.0) {
      return result;
    }
    if (result.energy_ev <= cut_ev_) {
      result.process = PhotonTransportProcess::Transparent;
      return result;
    }
    const double rate = medium_->GetPhotonCollisionRate(result.energy_ev);
    if (!std::isfinite(rate) || rate <= 0.0) {
      result.process = PhotonTransportProcess::Transparent;
      return result;
    }

    // TransportPhoton is an internal Garfield++ routine and keeps its photon
    // tracks in the AvalancheMicroscopic instance. Use a fresh instance for
    // every generated photon so no previous track can leak into this result.
    AvalancheMicroscopic photon_transport(&sensor_);
    photon_transport.SetPhotonTransportCut(cut_ev_);

    std::vector<Garfield::Seed> electron_stack;
    photon_transport.TransportPhotonExternal(x_cm, y_cm, z_cm, time_ns,
                                     result.energy_ev, 1, electron_stack);
    result.used_garfield = true;
    const std::size_t n_tracks = photon_transport.GetNumberOfPhotons();
    if (n_tracks == 0) return result;

    double stored_energy = result.energy_ev;
    photon_transport.GetPhoton(n_tracks - 1, stored_energy, result.x0_cm,
                               result.y0_cm, result.z0_cm, result.t0_ns,
                               result.x1_cm, result.y1_cm, result.z1_cm,
                               result.t1_ns, result.status);
    result.energy_ev = stored_energy;

    for (const auto& seed : electron_stack) {
      result.electrons.push_back(
          {seed.pt.x, seed.pt.y, seed.pt.z, seed.pt.t, seed.pt.energy,
           seed.pt.kx, seed.pt.ky, seed.pt.kz, seed.w});
    }

    if (result.status == kPhotonAbsorbedStatus) {
      result.process = result.electrons.empty()
                           ? PhotonTransportProcess::Photoabsorption
                           : PhotonTransportProcess::Photoionisation;
      return result;
    }

    const double tolerance = std::max(
        1.0e-8, 5.0e-6 * std::max({geometry_.xmax_cm - geometry_.xmin_cm,
                                   geometry_.ymax_cm - geometry_.ymin_cm,
                                   geometry_.gap_cm, 1.0e-6}));
    const bool on_cathode =
        std::abs(result.z1_cm - geometry_.gap_cm) <= tolerance;
    const bool inside_cathode =
        geometry_.infinite_electrodes ||
        (std::abs(result.x1_cm) <=
             geometry_.cathode_half_size_cm + tolerance &&
         std::abs(result.y1_cm) <=
             geometry_.cathode_half_size_cm + tolerance);
    if (on_cathode && inside_cathode) {
      result.process = PhotonTransportProcess::CathodeImpact;
    } else if (std::abs(result.z1_cm) <= tolerance) {
      result.process = PhotonTransportProcess::AnodeImpact;
    } else if (on_cathode ||
               std::abs(result.x1_cm - geometry_.xmin_cm) <= tolerance ||
               std::abs(result.x1_cm - geometry_.xmax_cm) <= tolerance ||
               std::abs(result.y1_cm - geometry_.ymin_cm) <= tolerance ||
               std::abs(result.y1_cm - geometry_.ymax_cm) <= tolerance) {
      result.process = PhotonTransportProcess::LateralEscape;
    }
    return result;
  }

 private:
  PhotonGeometry geometry_;
  double cut_ev_ = 0.0;
  MediumMagboltz* medium_ = nullptr;
  ComponentUser component_;
  Sensor sensor_;
};

// ============================================================================
//                              ROOT objects
// ============================================================================

struct OutputHistograms {
  // Time structure of the complete avalanche chain.
  // hElectronsVsTime counts individual electron endpoints.
  // hAvalancheElectronsVsTime places each avalanche at its seed time with
  // weight equal to the number of electrons produced by that avalanche.
  // Prompt endpoint view (default 0--10 ns) and full optical-chain view.
  TH1D electrons_vs_time;
  TH1D electrons_vs_time_full;
  TH1D avalanche_electrons_vs_time;
  TH1D avalanche_electrons_vs_time_prompt;
  TH2D avalanche_electrons_vs_time_generation;

  TH1D spectra;
  TH3F photon_xyz;
  TH2D photon_wavelength_time;
  TH1D cos_theta;
  TH1D phi;
  TH1D qe;
  TH2D impacts_xy;
  TH2D photoelectron_xy;
  TH2D photoelectron_time_energy;
  TH3F photoabsorption_xyz;
  TH2D photoabsorption_zt;
  TH3F photoionisation_xyz;
  TH2D photoionisation_zt;
  TH1D feedback_generation;

  OutputHistograms(const Config& c, const double time_max_ns,
                   const PhotonGeometry& geometry,
                   const std::string& qe_name)
      : electrons_vs_time(
            "hElectronsVsTime",
            "Prompt electron endpoints versus time;time [ns];electrons / bin",
            c.time_bins, 0.0, c.prompt_time_max_ns),
        electrons_vs_time_full(
            "hElectronsVsTimeFull",
            "Electron endpoints over the full feedback chain;time [ns];electrons / bin",
            c.time_bins, 0.0, time_max_ns),
        avalanche_electrons_vs_time(
            "hAvalancheElectronsVsTime",
            "Avalanche electron yield over the full chain;avalanche start time [ns];electrons / bin",
            c.time_bins, 0.0, time_max_ns),
        avalanche_electrons_vs_time_prompt(
            "hAvalancheElectronsVsTimePrompt",
            "Prompt avalanche electron yield;avalanche start time [ns];electrons / bin",
            c.time_bins, 0.0, c.prompt_time_max_ns),
        avalanche_electrons_vs_time_generation(
            "hAvalancheElectronsVsTimeGeneration",
            "Avalanche strength by feedback generation;avalanche start time [ns];feedback generation;number of electrons",
            c.time_bins, 0.0, time_max_ns,
            c.max_feedback_generations + 1, -0.5,
            c.max_feedback_generations + 0.5),
        spectra("hSpectra", "Photon spectrum;wavelength [nm];photons",
                std::max(1, static_cast<int>(std::llround(
                                (c.wavelength_max_nm - c.wavelength_min_nm) /
                                c.wavelength_bin_nm))),
                c.wavelength_min_nm, c.wavelength_max_nm),
        photon_xyz("hPhotonXYZ",
                   "Photon emission positions;x [cm];y [cm];z [cm]",
                   kExcitationSpatialBins, -c.xy_half_width_cm(),
                   c.xy_half_width_cm(), kExcitationSpatialBins,
                   -c.xy_half_width_cm(), c.xy_half_width_cm(),
                   kExcitationSpatialBins, 0.0, c.gap_cm()),
        photon_wavelength_time(
            "hPhotonWavelengthTime",
            "Photon wavelength and emission time;wavelength [nm];time [ns]",
            spectra.GetNbinsX(), c.wavelength_min_nm, c.wavelength_max_nm,
            c.time_bins, 0.0, time_max_ns),
        cos_theta("hCosTheta",
                  "Propagated isotropic direction;cos(#theta);photons",
                  c.angular_bins, -1.0, 1.0),
        phi("hPhi", "Propagated isotropic direction;#phi [rad];photons",
            c.angular_bins, -kPi, kPi),
        qe("hQE", ("Quantum efficiency: " + qe_name +
                    ";wavelength [nm];QE")
                       .c_str(),
           spectra.GetNbinsX(), c.wavelength_min_nm, c.wavelength_max_nm),
        impacts_xy("hImpactsXY",
                   "Propagated photon impacts on cathode;x [cm];y [cm]",
                   c.position_bins, -geometry.cathode_half_size_cm,
                   geometry.cathode_half_size_cm, c.position_bins,
                   -geometry.cathode_half_size_cm,
                   geometry.cathode_half_size_cm),
        photoelectron_xy("hPhotoElectronXY",
                         "Photoelectrons emitted from cathode;x [cm];y [cm]",
                         c.position_bins, -geometry.cathode_half_size_cm,
                         geometry.cathode_half_size_cm, c.position_bins,
                         -geometry.cathode_half_size_cm,
                         geometry.cathode_half_size_cm),
        photoelectron_time_energy(
            "hPhotoElectronTimeEnergy",
            "Photoelectron time and energy;time [ns];energy [eV]",
            c.time_bins, 0.0, time_max_ns, c.electron_energy_bins, 0.0,
            c.electron_energy_max_ev),
        photoabsorption_xyz(
            "hPhotoAbsorptionXYZ",
            "All gas photoabsorption positions;x [cm];y [cm];z [cm]",
            kExcitationSpatialBins, geometry.xmin_cm, geometry.xmax_cm,
            kExcitationSpatialBins, geometry.ymin_cm, geometry.ymax_cm,
            kExcitationSpatialBins, 0.0, geometry.gap_cm),
        photoabsorption_zt("hPhotoAbsorptionZT",
                           "All gas photoabsorption;z [cm];time [ns]",
                           kExcitationSpatialBins, 0.0, geometry.gap_cm,
                           c.time_bins, 0.0, time_max_ns),
        photoionisation_xyz(
            "hPhotoIonisationXYZ",
            "Gas photoionisation positions;x [cm];y [cm];z [cm]",
            kExcitationSpatialBins, geometry.xmin_cm, geometry.xmax_cm,
            kExcitationSpatialBins, geometry.ymin_cm, geometry.ymax_cm,
            kExcitationSpatialBins, 0.0, geometry.gap_cm),
        photoionisation_zt("hPhotoIonisationZT",
                           "Gas photoionisation;z [cm];time [ns]",
                           kExcitationSpatialBins, 0.0, geometry.gap_cm,
                           c.time_bins, 0.0, time_max_ns),
        feedback_generation(
            "hFeedbackGeneration",
            "Avalanche seeds by generation;generation;seeds",
            c.max_feedback_generations + 1, -0.5,
            c.max_feedback_generations + 0.5) {
    use_hist_draw_option(electrons_vs_time);
    use_hist_draw_option(electrons_vs_time_full);
    use_hist_draw_option(avalanche_electrons_vs_time);
    use_hist_draw_option(avalanche_electrons_vs_time_prompt);
    use_hist_draw_option(spectra);
    use_hist_draw_option(cos_theta);
    use_hist_draw_option(phi);
    use_hist_draw_option(qe);
    use_hist_draw_option(feedback_generation);

    electrons_vs_time.Sumw2();
    electrons_vs_time_full.Sumw2();
    avalanche_electrons_vs_time.Sumw2();
    avalanche_electrons_vs_time_prompt.Sumw2();
    avalanche_electrons_vs_time_generation.Sumw2();

    // Only the full-chain objects are allowed to grow. Extending the prompt
    // histogram would reproduce the old problem: a late optical event would
    // stretch the axis to tens of microseconds and hide the ns-scale peak.
    electrons_vs_time_full.SetCanExtend(TH1::kXaxis);
    avalanche_electrons_vs_time.SetCanExtend(TH1::kXaxis);
    avalanche_electrons_vs_time_generation.SetCanExtend(TH1::kXaxis);

    // The weighted avalanche histograms use one Fill call per avalanche and
    // weight it by ne. ROOT's "Entries" would therefore mean avalanches, not
    // electrons, so hide that potentially misleading statistics box.
    avalanche_electrons_vs_time.SetStats(false);
    avalanche_electrons_vs_time_prompt.SetStats(false);
    avalanche_electrons_vs_time_generation.SetStats(false);

    spectra.Sumw2();
    photon_wavelength_time.Sumw2();
    cos_theta.Sumw2();
    phi.Sumw2();
    impacts_xy.Sumw2();
    photoelectron_xy.Sumw2();
    photoelectron_time_energy.Sumw2();
  }
};

struct PhotonTransportSummary {
  bool photo_absorption_enabled = false;
  double generated_photons = 0.0;
  double transport_candidates = 0.0;
  double garfield_transport_attempts = 0.0;
  double garfield_tracks = 0.0;
  double transparent_fallbacks = 0.0;
  double geometric_rejections = 0.0;
  double photoabsorptions = 0.0;
  double nonionising_photoabsorptions = 0.0;
  double photoionisations = 0.0;
  double photoionisation_electrons = 0.0;
  double cathode_impacts = 0.0;
  double anode_impacts = 0.0;
  double lateral_escapes = 0.0;
  double transport_failures = 0.0;
  double accounted_transport_outcomes = 0.0;
  double unaccounted_transport_candidates = 0.0;
  Long64_t feedback_photoelectron_seeds = 0;
  Long64_t feedback_photoionisation_seeds = 0;
  Long64_t suppressed_generation_seeds = 0;
  Long64_t suppressed_avalanche_limit_seeds = 0;
};

// ============================================================================
//                         Iterative avalanche chain
// ============================================================================

enum class SeedType : int { Primary = 0, CathodePhotoelectron = 1,
                            GasPhotoionisation = 2 };

struct ElectronSeed {
  double x_cm = 0.0;
  double y_cm = 0.0;
  double z_cm = 0.0;
  double time_ns = 0.0;
  double energy_ev = 0.1;
  double dx = 0.0;
  double dy = 0.0;
  double dz = -1.0;
  int generation = 0;
  SeedType type = SeedType::Primary;
};

ElectronSeed primary_seed(const Config& config, const int event,
                          TRandom3& random) {
  ElectronSeed seed;
  seed.x_cm = random.Uniform(-config.gap_cm(), config.gap_cm());
  seed.y_cm = random.Uniform(-config.gap_cm(), config.gap_cm());
  seed.z_cm = config.launch_z_cm() - 1.0e-9;
  seed.energy_ev = event == 0
                       ? config.initial_energy_ev
                       : g_energy.random_energy(config.initial_energy_ev);
  const double phi = random.Uniform(0.0, 2.0 * kPi);
  const double theta = 0.75 * kPi;
  seed.dx = std::cos(phi) * std::sin(theta);
  seed.dy = std::sin(phi) * std::sin(theta);
  seed.dz = std::cos(theta);
  return seed;
}

ElectronSeed cathode_photoelectron_seed(
    const Config& config, double x_cm, double y_cm, double time_ns,
    double energy_ev, int generation, TRandom3& random) {
  ElectronSeed seed;
  seed.x_cm = x_cm;
  seed.y_cm = y_cm;
  seed.z_cm = config.gap_cm() - 1.0e-9;
  seed.time_ns = time_ns;
  seed.energy_ev = std::max(1.0e-6, energy_ev);
  const double cos_theta = -random.Uniform(0.0, 1.0);
  const double phi = random.Uniform(-kPi, kPi);
  const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
  seed.dx = sin_theta * std::cos(phi);
  seed.dy = sin_theta * std::sin(phi);
  seed.dz = cos_theta;
  seed.generation = generation;
  seed.type = SeedType::CathodePhotoelectron;
  return seed;
}

ElectronSeed gas_photoionisation_seed(const PhotoionisationElectron& electron,
                                      int generation, TRandom3& random) {
  ElectronSeed seed;
  seed.x_cm = electron.x_cm;
  seed.y_cm = electron.y_cm;
  seed.z_cm = electron.z_cm;
  seed.time_ns = electron.time_ns;
  seed.energy_ev = std::max(1.0e-6, electron.energy_ev);
  seed.dx = electron.dx;
  seed.dy = electron.dy;
  seed.dz = electron.dz;
  const double norm = std::sqrt(seed.dx * seed.dx + seed.dy * seed.dy +
                                seed.dz * seed.dz);
  if (!std::isfinite(norm) || norm <= 0.0) {
    const double cos_theta = random.Uniform(-1.0, 1.0);
    const double phi = random.Uniform(-kPi, kPi);
    const double sin_theta =
        std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    seed.dx = sin_theta * std::cos(phi);
    seed.dy = sin_theta * std::sin(phi);
    seed.dz = cos_theta;
  } else {
    seed.dx /= norm;
    seed.dy /= norm;
    seed.dz /= norm;
  }
  seed.generation = generation;
  seed.type = SeedType::GasPhotoionisation;
  return seed;
}

void fill_track_direction(const PhotonResult& result,
                          OutputHistograms& histograms, const double weight) {
  const double dx = result.x1_cm - result.x0_cm;
  const double dy = result.y1_cm - result.y0_cm;
  const double dz = result.z1_cm - result.z0_cm;
  const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (!std::isfinite(distance) || distance <= 0.0) return;
  histograms.cos_theta.Fill(std::clamp(dz / distance, -1.0, 1.0), weight);
  histograms.phi.Fill(std::atan2(dy, dx), weight);
}

}  // namespace

int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);
  try {
    const std::time_t wall_start = std::time(nullptr);
    const std::clock_t cpu_start = std::clock();
    Config config = read_config(argc, argv);

    if (!config.quiet) {
      std::cout << "[secondaryAvalanches] " << config.root_file << "\n"
                << "  mixture                  = " << config.gas1 << " "
                << config.composition1 << "% + " << config.gas2 << " "
                << config.composition2 << "%";
      if (!config.gas3.empty() && config.composition3 > 0.0) {
        std::cout << " + " << config.gas3 << " " << config.composition3
                  << "%";
      }
      std::cout << "\n"
                << "  field / pressure / gap    = " << config.field_v_cm
                << " V/cm / " << config.pressure_bar << " bar / "
                << config.gap_mm << " mm\n"
                << "  photon MC samples         = " << config.mc_samples
                << "\n"
                << "  gas photoabsorption       = "
                << (config.photo_absorption ? "on" : "off") << "\n"
                << "  propagate only E > PhiT   = "
                << (config.propagate_only_above_phit ? "on" : "off") << "\n"
                << "  requested photon cut      = "
                << config.photon_transport_cut_ev << " eV\n"
                << "  QE material / model       = " << config.qe.material
                << " / " << config.qe.model_mode << "\n"
                << "  maximum feedback depth    = "
                << config.max_feedback_generations << std::endl;
    }

    gRandom->SetSeed(config.seed);
    TRandom3 random(config.seed);
    Garfield::RandomEngineSTL garfield_random(config.seed);
    Garfield::Random::SetEngine(garfield_random);

    // ========================================================================
    //                              Setup gas
    // ========================================================================

    MediumMagboltz gas;
    std::vector<std::pair<std::string, double>> active_gases;
    if (config.composition1 > 0.0) {
      active_gases.push_back({config.gas1, config.composition1});
    }
    if (config.composition2 > 0.0) {
      active_gases.push_back({config.gas2, config.composition2});
    }
    if (config.composition3 > 0.0) {
      active_gases.push_back({config.gas3, config.composition3});
    }
    if (active_gases.size() == 1) {
      gas.SetComposition(active_gases[0].first, active_gases[0].second);
    } else if (active_gases.size() == 2) {
      gas.SetComposition(active_gases[0].first, active_gases[0].second,
                         active_gases[1].first, active_gases[1].second);
    } else if (active_gases.size() == 3) {
      gas.SetComposition(active_gases[0].first, active_gases[0].second,
                         active_gases[1].first, active_gases[1].second,
                         active_gases[2].first, active_gases[2].second);
    } else {
      throw std::runtime_error("Require between one and three active gases");
    }
    gas.SetTemperature(config.temperature_k);
    gas.SetPressure(config.pressure_torr());
    gas.SetMaxElectronEnergy(config.max_electron_energy_ev);
    gas.Initialise();

    GasTransport gas_transport;
    gas_transport.measure(config, gas);
    const auto level_info = read_level_info(config, gas);
    const int n_levels = std::max(1, static_cast<int>(level_info.size()));
    photonemission::GasData model_gas = make_model_gas(config);

    // ========================================================================
    //                    Uniform field and optical geometry
    // ========================================================================

    ComponentUser field;
    const double uniform_field = config.field_v_cm;
    field.SetElectricField(
        [uniform_field](double, double, double, double& ex, double& ey,
                        double& ez) {
          ex = 0.0;
          ey = 0.0;
          ez = uniform_field;
        });
    field.SetArea(-config.xy_half_width_cm(), -config.xy_half_width_cm(), 0.0,
                  config.xy_half_width_cm(), config.xy_half_width_cm(),
                  config.z_max_cm());
    field.SetMedium(&gas);

    SpaceCharge space_charge;
    space_charge.initialise(config, gas);

    Sensor sensor;
    sensor.AddComponent(&field);
    if (space_charge.rings != nullptr) {
      sensor.AddComponent(space_charge.rings.get());
    }
    sensor.SetArea(-config.xy_half_width_cm(), -config.xy_half_width_cm(), 0.0,
                   config.xy_half_width_cm(), config.xy_half_width_cm(),
                   config.z_max_cm());

    PhotonGeometry photon_geometry;
    photon_geometry.xmin_cm = -config.transport_half_width_cm();
    photon_geometry.xmax_cm = config.transport_half_width_cm();
    photon_geometry.ymin_cm = -config.transport_half_width_cm();
    photon_geometry.ymax_cm = config.transport_half_width_cm();
    photon_geometry.gap_cm = config.gap_cm();
    photon_geometry.cathode_half_size_cm = config.cathode_half_width_cm();
    photon_geometry.infinite_electrodes = config.infinite_electrodes;

    const auto qe_material =
        photonfeedback::load_photoelectron_material(config.qe);
    const double qe_threshold_ev =
        photonfeedback::photoelectron_threshold_ev(config.qe, qe_material);
    const double effective_photon_transport_cut_ev =
        config.propagate_only_above_phit
            ? std::max(config.photon_transport_cut_ev, qe_threshold_ev)
            : config.photon_transport_cut_ev;

    if (!config.quiet) {
      std::cout << "  PhiT / QE threshold       = " << qe_threshold_ev
                << " eV\n"
                << "  effective photon cut      = "
                << effective_photon_transport_cut_ev << " eV\n"
                << "  infinite electrodes       = "
                << (config.infinite_electrodes ? "on" : "off") << "\n"
                << "  optical half width        = "
                << config.transport_half_width_cm() << " cm ("
                << config.transport_half_width_cm() / config.gap_cm()
                << " gaps)\n"
                << "  prompt time window        = "
                << config.prompt_time_max_ns << " ns\n";
    }

    std::unique_ptr<PhotonPropagation> photon_propagation;
    if (config.photo_absorption) {
      photon_propagation = std::make_unique<PhotonPropagation>(
          photon_geometry, gas, effective_photon_transport_cut_ev);
    }

    const double time_max_ns = config.time_max_ns > 0.0
                                   ? config.time_max_ns
                                   : std::max(100.0, 20.0 * 3140.0);

    // ========================================================================
    //                              ROOT output
    // ========================================================================

    const std::filesystem::path output_path(config.root_file);
    if (!output_path.parent_path().empty()) {
      std::filesystem::create_directories(output_path.parent_path());
    }

    // Keep every histogram/tree in memory while the full simulation runs.
    // The ROOT file is opened only after a successful chain, so a failed run
    // cannot leave a valid-looking but incomplete output file behind.
    TH1D h_electron_energy(
        "hElectronEnergyDistribution",
        "Electron energy distribution from null-collision steps;E_{e} [eV];samples",
        kEnergyBins, 0.0, 50.0);
    TH1D h_levels("hLevels", "Excitation Distribution;hLevel;excitations",
                  n_levels, 0.0, static_cast<double>(n_levels));
    use_hist_draw_option(h_electron_energy);
    use_hist_draw_option(h_levels);

    std::unique_ptr<TH3I> h_exc_xyz;
    std::unique_ptr<TH2I> h_exc_zt;
    if (config.record_excitation_positions) {
      h_exc_xyz = std::make_unique<TH3I>(
          "hExcXYZ",
          "Inelastic/excitation spatial distribution;x [cm];y [cm];z [cm]",
          kExcitationSpatialBins, -config.xy_half_width_cm(),
          config.xy_half_width_cm(), kExcitationSpatialBins,
          -config.xy_half_width_cm(), config.xy_half_width_cm(),
          kExcitationSpatialBins, 0.0, config.gap_cm());
      h_exc_zt = std::make_unique<TH2I>(
          "hExcZT",
          "Inelastic/excitation longitudinal-time distribution;z [cm];t [ns]",
          kExcitationSpatialBins, 0.0, config.gap_cm(),
          kExcitationTimeBins, 0.0, kInitialTimeRangeNs);
      h_exc_zt->SetCanExtend(TH1::kYaxis);
    }

    OutputHistograms photon_hists(
        config, time_max_ns, photon_geometry,
        qe_material.loaded ? qe_material.material : std::string("constant"));
    for (int bin = 1; bin <= photon_hists.qe.GetNbinsX(); ++bin) {
      const double wavelength_nm = photon_hists.qe.GetBinCenter(bin);
      photon_hists.qe.SetBinContent(
          bin, photonfeedback::quantum_efficiency(
                   wavelength_nm, config.qe, qe_material));
    }

    TTree data_per_primary("dataPerPrimaryElectron",
                           "Data per primary electron including feedback");
    TTree data_per_avalanche("dataPerAvalanche",
                             "Every primary and feedback avalanche");
    TTree data_per_electron("dataPerElectron", "Data per electron endpoint");
    TTree gas_data("gasData", "Gas mixture and run conditions");

    Long64_t primary_ne = 0;
    Long64_t primary_ni = 0;
    Long64_t primary_avalanche_ne = 0;
    Long64_t primary_avalanche_ni = 0;
    Int_t one_primary = 1;
    Int_t primary_n_avalanches = 0;
    Int_t primary_n_photoelectrons = 0;
    Int_t primary_n_photoionisation_seeds = 0;
    Int_t primary_max_generation = 0;
    // ne/ni are the complete effective charge after all feedback generations.
    data_per_primary.Branch("ne", &primary_ne, "ne/L");
    data_per_primary.Branch("ni", &primary_ni, "ni/L");
    data_per_primary.Branch("neTotalWithFeedback", &primary_ne,
                            "neTotalWithFeedback/L");
    data_per_primary.Branch("niTotalWithFeedback", &primary_ni,
                            "niTotalWithFeedback/L");
    data_per_primary.Branch("nePrimaryAvalanche", &primary_avalanche_ne,
                            "nePrimaryAvalanche/L");
    data_per_primary.Branch("niPrimaryAvalanche", &primary_avalanche_ni,
                            "niPrimaryAvalanche/L");
    data_per_primary.Branch("npe", &one_primary, "npe/I");
    data_per_primary.Branch("nAvalanches", &primary_n_avalanches,
                            "nAvalanches/I");
    data_per_primary.Branch("nPhotoelectrons", &primary_n_photoelectrons,
                            "nPhotoelectrons/I");
    data_per_primary.Branch("nPhotoionisationSeeds",
                            &primary_n_photoionisation_seeds,
                            "nPhotoionisationSeeds/I");
    data_per_primary.Branch("maxGeneration", &primary_max_generation,
                            "maxGeneration/I");

    Int_t avalanche_primary_id = 0;
    Int_t avalanche_id = 0;
    Int_t avalanche_generation = 0;
    Int_t avalanche_seed_type = 0;
    Double_t avalanche_seed_x = 0.0;
    Double_t avalanche_seed_y = 0.0;
    Double_t avalanche_seed_z = 0.0;
    Double_t avalanche_seed_t = 0.0;
    Double_t avalanche_seed_energy = 0.0;
    Int_t avalanche_ne = 0;
    Int_t avalanche_ni = 0;
    Double_t avalanche_photons = 0.0;
    Double_t avalanche_cathode_impacts = 0.0;
    Double_t avalanche_photoelectrons = 0.0;
    Double_t avalanche_photoionisation_electrons = 0.0;
    data_per_avalanche.Branch("primaryId", &avalanche_primary_id,
                              "primaryId/I");
    data_per_avalanche.Branch("avalancheId", &avalanche_id,
                              "avalancheId/I");
    data_per_avalanche.Branch("generation", &avalanche_generation,
                              "generation/I");
    data_per_avalanche.Branch("seedType", &avalanche_seed_type,
                              "seedType/I");
    data_per_avalanche.Branch("seedXcm", &avalanche_seed_x, "seedXcm/D");
    data_per_avalanche.Branch("seedYcm", &avalanche_seed_y, "seedYcm/D");
    data_per_avalanche.Branch("seedZcm", &avalanche_seed_z, "seedZcm/D");
    data_per_avalanche.Branch("seedTimeNs", &avalanche_seed_t,
                              "seedTimeNs/D");
    data_per_avalanche.Branch("seedEnergyEv", &avalanche_seed_energy,
                              "seedEnergyEv/D");
    data_per_avalanche.Branch("ne", &avalanche_ne, "ne/I");
    data_per_avalanche.Branch("ni", &avalanche_ni, "ni/I");
    data_per_avalanche.Branch("nPhotons", &avalanche_photons,
                              "nPhotons/D");
    data_per_avalanche.Branch("nCathodeImpacts",
                              &avalanche_cathode_impacts,
                              "nCathodeImpacts/D");
    data_per_avalanche.Branch("nPhotoelectrons",
                              &avalanche_photoelectrons,
                              "nPhotoelectrons/D");
    data_per_avalanche.Branch("nPhotoionisationElectrons",
                              &avalanche_photoionisation_electrons,
                              "nPhotoionisationElectrons/D");

    Int_t endpoint_status = 0;
    Int_t endpoint_primary_id = 0;
    Int_t endpoint_avalanche_id = 0;
    Int_t endpoint_generation = 0;
    data_per_electron.Branch("status", &endpoint_status, "status/I");
    data_per_electron.Branch("primaryId", &endpoint_primary_id,
                             "primaryId/I");
    data_per_electron.Branch("avalancheId", &endpoint_avalanche_id,
                             "avalancheId/I");
    data_per_electron.Branch("generation", &endpoint_generation,
                             "generation/I");

    // ========================================================================
    //                         Microscopic avalanches
    // ========================================================================

    g_energy.reset();
    g_collisions.reset(h_levels, h_exc_xyz.get(), h_exc_zt.get(),
                       config.space_charge);

    AvalancheMicroscopic avalanche;
    avalanche.SetSensor(&sensor);
    avalanche.EnableSignalCalculation(false);
    avalanche.EnableNullCollisionSteps(true, 1);
    avalanche.SetUserHandleStep(handle_step);
    avalanche.SetUserHandleCollision(handle_collision);

    RunningStatistics electron_statistics;
    RunningStatistics ion_statistics;
    RunningStatistics primary_avalanche_electron_statistics;
    RunningStatistics primary_avalanche_ion_statistics;
    PhotonTransportSummary photon_summary;
    photon_summary.photo_absorption_enabled = config.photo_absorption;
    Long64_t ne_total = 0;
    Long64_t ni_total = 0;
    Long64_t total_avalanches = 0;
    Long64_t total_mc_photons = 0;
    int global_max_generation = 0;
    const double photon_weight = 1.0 / static_cast<double>(config.mc_samples);

    auto last_progress_update =
        std::chrono::steady_clock::now() - std::chrono::seconds(1);

    for (int event = 0; event < config.max_npe; ++event) {
      std::deque<ElectronSeed> queue;
      queue.push_back(primary_seed(config, event, random));
      primary_ne = 0;
      primary_ni = 0;
      primary_avalanche_ne = 0;
      primary_avalanche_ni = 0;
      primary_n_avalanches = 0;
      primary_n_photoelectrons = 0;
      primary_n_photoionisation_seeds = 0;
      primary_max_generation = 0;
      long long primary_mc_photons = 0;

      while (!queue.empty()) {
        if (primary_n_avalanches >= config.max_avalanches_per_primary) {
          photon_summary.suppressed_avalanche_limit_seeds += queue.size();
          queue.clear();
          break;
        }

        const ElectronSeed seed = queue.front();
        queue.pop_front();
        g_collisions.start_avalanche();

        avalanche.AvalancheElectron(seed.x_cm, seed.y_cm, seed.z_cm,
                                    seed.time_ns, seed.energy_ev,
                                    seed.dx, seed.dy, seed.dz);
        avalanche.GetAvalancheSize(avalanche_ne, avalanche_ni);
        if (seed.type == SeedType::Primary && seed.generation == 0) {
          primary_avalanche_ne = avalanche_ne;
          primary_avalanche_ni = avalanche_ni;
        }
        primary_ne += avalanche_ne;
        primary_ni += avalanche_ni;
        ++primary_n_avalanches;
        ++total_avalanches;
        primary_max_generation =
            std::max(primary_max_generation, seed.generation);
        global_max_generation =
            std::max(global_max_generation, seed.generation);
        photon_hists.feedback_generation.Fill(seed.generation);

        endpoint_primary_id = event;
        endpoint_avalanche_id = primary_n_avalanches - 1;
        endpoint_generation = seed.generation;
        for (int electron = 0; electron < avalanche_ne; ++electron) {
          double ex0 = 0.0, ey0 = 0.0, ez0 = 0.0, et0 = 0.0, ee0 = 0.0;
          double ex1 = 0.0, ey1 = 0.0, ez1 = 0.0, et1 = 0.0, ee1 = 0.0;
          avalanche.GetElectronEndpoint(
              electron, ex0, ey0, ez0, et0, ee0, ex1, ey1, ez1, et1,
              ee1, endpoint_status);
          data_per_electron.Fill();

          // One entry per final electron. Peaks in this histogram show when
          // the primary and feedback avalanches deliver their electrons.
          if (std::isfinite(et1) && et1 >= 0.0) {
            photon_hists.electrons_vs_time_full.Fill(et1);
            if (et1 <= config.prompt_time_max_ns) {
              photon_hists.electrons_vs_time.Fill(et1);
            }
          }
        }

        // One weighted entry per avalanche. This representation makes the
        // start time and strength of the primary, secondary, ... avalanches
        // explicit. The 2-D version keeps the feedback generation separated.
        if (avalanche_ne > 0 && std::isfinite(seed.time_ns) &&
            seed.time_ns >= 0.0) {
          photon_hists.avalanche_electrons_vs_time.Fill(
              seed.time_ns, static_cast<double>(avalanche_ne));
          if (seed.time_ns <= config.prompt_time_max_ns) {
            photon_hists.avalanche_electrons_vs_time_prompt.Fill(
                seed.time_ns, static_cast<double>(avalanche_ne));
          }
          photon_hists.avalanche_electrons_vs_time_generation.Fill(
              seed.time_ns, static_cast<double>(seed.generation),
              static_cast<double>(avalanche_ne));
        }
        space_charge.add_new_ions(config);

        avalanche_photons = 0.0;
        avalanche_cathode_impacts = 0.0;
        avalanche_photoelectrons = 0.0;
        avalanche_photoionisation_electrons = 0.0;

        const auto populations = photonemission::populations_from_sites(
            level_info, g_collisions.sites_this_avalanche);
        model_gas.ne_total = avalanche_ne;
        model_gas.ni_total = avalanche_ni;
        const auto components = photonemission::create_photon_emission(
            model_gas, populations, config.parameters_dir,
            static_cast<int>(std::min<long long>(
                config.mc_samples, std::numeric_limits<int>::max())));

        auto enqueue_child = [&](ElectronSeed child, bool is_photoelectron) {
          if (child.generation > config.max_feedback_generations) {
            ++photon_summary.suppressed_generation_seeds;
            return;
          }
          if (primary_n_avalanches + static_cast<int>(queue.size()) >=
              config.max_avalanches_per_primary) {
            ++photon_summary.suppressed_avalanche_limit_seeds;
            return;
          }
          queue.push_back(child);
          if (is_photoelectron) {
            ++primary_n_photoelectrons;
            ++photon_summary.feedback_photoelectron_seeds;
          } else {
            ++primary_n_photoionisation_seeds;
            ++photon_summary.feedback_photoionisation_seeds;
          }
        };

        auto record_cathode_impact =
            [&](double x_cm, double y_cm, double impact_time_ns,
                double wavelength_nm) {
              photon_hists.impacts_xy.Fill(x_cm, y_cm, photon_weight);
              photon_summary.cathode_impacts += photon_weight;
              avalanche_cathode_impacts += photon_weight;

              const auto photoelectron = photonfeedback::emit_photoelectron(
                  wavelength_nm, config.electron_extraction_efficiency,
                  config.qe, qe_material, random);
              if (!photoelectron.emitted) return;

              photon_hists.photoelectron_xy.Fill(x_cm, y_cm, photon_weight);
              photon_hists.photoelectron_time_energy.Fill(
                  impact_time_ns, photoelectron.energy_ev, photon_weight);
              avalanche_photoelectrons += photon_weight;

              // The spectra/histograms use all mcSamples with weight 1/N. For
              // the physical feedback chain, randomly keep one in N accepted
              // photoelectrons. This is unbiased and prevents mcSamples from
              // artificially multiplying the number of avalanches.
              if (random.Uniform() >= photon_weight) return;
              enqueue_child(cathode_photoelectron_seed(
                                config, x_cm, y_cm, impact_time_ns,
                                photoelectron.energy_ev, seed.generation + 1,
                                random),
                            true);
            };

        for (const auto& component : components) {
          const long long samples = photonemission::number_of_mc_photons(
              component.expected_photons, config.mc_samples, random);
          if (primary_mc_photons + samples >
              config.max_mc_photons_per_primary) {
            throw std::runtime_error(
                "Photon safety limit reached for one primary electron; reduce "
                "mcSamples or the feedback limits");
          }
          primary_mc_photons += samples;
          total_mc_photons += samples;

          for (long long sample = 0; sample < samples; ++sample) {
            const auto* site = photonemission::sample_source_site(
                component, level_info, g_collisions.sites_this_avalanche,
                random);
            if (site == nullptr) continue;

            const double x_cm = site->x_cm;
            const double y_cm = site->y_cm;
            const double z_cm = std::clamp(site->z_cm, 0.0, config.gap_cm());
            const double wavelength_nm = photonemission::sample_wavelength_nm(
                component, random, config.wavelength_min_nm,
                config.wavelength_max_nm);
            const double emission_time_ns =
                std::max(0.0, site->time_ns) +
                photonemission::sample_emission_delay_ns(component, random);
            const double photon_energy_ev =
                photonemission::photon_energy_ev(wavelength_nm);

            photon_hists.spectra.Fill(wavelength_nm, photon_weight);
            photon_hists.photon_xyz.Fill(x_cm, y_cm, z_cm, photon_weight);
            photon_hists.photon_wavelength_time.Fill(
                wavelength_nm, emission_time_ns, photon_weight);
            photon_summary.generated_photons += photon_weight;
            avalanche_photons += photon_weight;

            // Keep the complete emitted spectrum, but only transport photons
            // above the selected physical/technical threshold. When enabled,
            // PhiT comes from the same QE material/model used for emission.
            if (photon_energy_ev <= effective_photon_transport_cut_ev) {
              continue;
            }

            if (!config.photo_absorption) {
              photon_summary.transport_candidates += photon_weight;
              const auto direction = sample_isotropic_direction(random);
              photon_hists.cos_theta.Fill(direction.cos_theta,
                                          photon_weight);
              photon_hists.phi.Fill(direction.phi, photon_weight);
              const auto boundary = propagate_to_box_boundary(
                  x_cm, y_cm, z_cm, direction, photon_geometry);
              if (!boundary.valid) {
                photon_summary.geometric_rejections += photon_weight;
                continue;
              }
              const double boundary_time_ns =
                  emission_time_ns + boundary.distance_cm / kLightSpeedCmPerNs;
              if (boundary.surface == BoundarySurface::Cathode) {
                record_cathode_impact(boundary.x_cm, boundary.y_cm,
                                      boundary_time_ns, wavelength_nm);
              } else if (boundary.surface == BoundarySurface::Anode) {
                photon_summary.anode_impacts += photon_weight;
              } else {
                photon_summary.lateral_escapes += photon_weight;
              }
              continue;
            }

            photon_summary.transport_candidates += photon_weight;
            const auto result = photon_propagation->transport(
                x_cm, y_cm, z_cm, emission_time_ns, wavelength_nm);
            if (result.used_garfield) {
              photon_summary.garfield_transport_attempts += photon_weight;
            }

            if (result.process == PhotonTransportProcess::Transparent) {
              photon_summary.transparent_fallbacks += photon_weight;
              const auto direction = sample_isotropic_direction(random);
              photon_hists.cos_theta.Fill(direction.cos_theta,
                                          photon_weight);
              photon_hists.phi.Fill(direction.phi, photon_weight);
              const auto boundary = propagate_to_box_boundary(
                  x_cm, y_cm, z_cm, direction, photon_geometry);
              if (!boundary.valid) {
                photon_summary.transport_failures += photon_weight;
                continue;
              }
              const double boundary_time_ns =
                  emission_time_ns + boundary.distance_cm / kLightSpeedCmPerNs;
              if (boundary.surface == BoundarySurface::Cathode) {
                record_cathode_impact(boundary.x_cm, boundary.y_cm,
                                      boundary_time_ns, wavelength_nm);
              } else if (boundary.surface == BoundarySurface::Anode) {
                photon_summary.anode_impacts += photon_weight;
              } else {
                photon_summary.lateral_escapes += photon_weight;
              }
              continue;
            }

            if (result.process == PhotonTransportProcess::Invalid) {
              photon_summary.transport_failures += photon_weight;
              continue;
            }
            photon_summary.garfield_tracks += photon_weight;
            fill_track_direction(result, photon_hists, photon_weight);

            switch (result.process) {
              case PhotonTransportProcess::Photoabsorption:
                photon_hists.photoabsorption_xyz.Fill(
                    result.x1_cm, result.y1_cm, result.z1_cm, photon_weight);
                photon_hists.photoabsorption_zt.Fill(
                    result.z1_cm, result.t1_ns, photon_weight);
                photon_summary.photoabsorptions += photon_weight;
                photon_summary.nonionising_photoabsorptions += photon_weight;
                break;

              case PhotonTransportProcess::Photoionisation:
                photon_hists.photoabsorption_xyz.Fill(
                    result.x1_cm, result.y1_cm, result.z1_cm, photon_weight);
                photon_hists.photoabsorption_zt.Fill(
                    result.z1_cm, result.t1_ns, photon_weight);
                photon_hists.photoionisation_xyz.Fill(
                    result.x1_cm, result.y1_cm, result.z1_cm, photon_weight);
                photon_hists.photoionisation_zt.Fill(
                    result.z1_cm, result.t1_ns, photon_weight);
                photon_summary.photoabsorptions += photon_weight;
                photon_summary.photoionisations += photon_weight;
                for (const auto& electron : result.electrons) {
                  photon_summary.photoionisation_electrons +=
                      photon_weight * electron.multiplicity;
                  avalanche_photoionisation_electrons +=
                      photon_weight * electron.multiplicity;
                  if (!config.propagate_photoionisation_electrons) continue;
                  for (std::size_t copy = 0; copy < electron.multiplicity;
                       ++copy) {
                    if (random.Uniform() >= photon_weight) continue;
                    enqueue_child(gas_photoionisation_seed(
                                      electron, seed.generation + 1, random),
                                  false);
                  }
                }
                break;

              case PhotonTransportProcess::CathodeImpact:
                record_cathode_impact(result.x1_cm, result.y1_cm,
                                      result.t1_ns, wavelength_nm);
                break;
              case PhotonTransportProcess::AnodeImpact:
                photon_summary.anode_impacts += photon_weight;
                break;
              case PhotonTransportProcess::LateralEscape:
                photon_summary.lateral_escapes += photon_weight;
                break;
              case PhotonTransportProcess::Transparent:
              case PhotonTransportProcess::Invalid:
                break;
            }
          }
        }

        avalanche_primary_id = event;
        avalanche_id = primary_n_avalanches - 1;
        avalanche_generation = seed.generation;
        avalanche_seed_type = static_cast<int>(seed.type);
        avalanche_seed_x = seed.x_cm;
        avalanche_seed_y = seed.y_cm;
        avalanche_seed_z = seed.z_cm;
        avalanche_seed_t = seed.time_ns;
        avalanche_seed_energy = seed.energy_ev;
        data_per_avalanche.Fill();
      }

      ne_total += primary_ne;
      ni_total += primary_ni;
      electron_statistics.add(static_cast<double>(primary_ne));
      ion_statistics.add(static_cast<double>(primary_ni));
      primary_avalanche_electron_statistics.add(
          static_cast<double>(primary_avalanche_ne));
      primary_avalanche_ion_statistics.add(
          static_cast<double>(primary_avalanche_ni));
      data_per_primary.Fill();

      const int completed = event + 1;
      const auto now = std::chrono::steady_clock::now();
      const bool progress_due =
          completed == 1 || completed == config.max_npe ||
          std::chrono::duration_cast<std::chrono::milliseconds>(
              now - last_progress_update)
                  .count() >= 250;
      if (progress_due) {
        std::cout << "PROGRESS " << config.job_id << " " << completed << " "
                  << config.max_npe << " " << electron_statistics.mean << " "
                  << electron_statistics.relative_error() << std::endl;
        last_progress_update = now;
      }

      const bool enough_primaries = completed >= config.min_npe;
      const bool precision_reached =
          config.target_relative_error > 0.0 && enough_primaries &&
          electron_statistics.relative_error() <= config.target_relative_error;
      if (precision_reached) break;
    }

    // ========================================================================
    //                           Final summaries
    // ========================================================================

    const int actual_npe = electron_statistics.n;
    double gain = electron_statistics.mean;
    double gain_error = electron_statistics.error_on_mean();
    double ion_mean = ion_statistics.mean;
    double primary_avalanche_gain =
        primary_avalanche_electron_statistics.mean;
    double primary_avalanche_gain_error =
        primary_avalanche_electron_statistics.error_on_mean();
    double primary_avalanche_ion_mean =
        primary_avalanche_ion_statistics.mean;
    double alpha_effective =
        gain > 1.0 ? std::log(gain) / config.gap_cm()
                   : std::numeric_limits<double>::quiet_NaN();
    double alpha_error =
        gain > 1.0 && std::isfinite(gain_error)
            ? gain_error / (gain * config.gap_cm())
            : std::numeric_limits<double>::quiet_NaN();
    double relative_gain_error =
        gain > 0.0 ? gain_error / gain
                   : std::numeric_limits<double>::quiet_NaN();
    double primary_avalanche_alpha =
        primary_avalanche_gain > 1.0
            ? std::log(primary_avalanche_gain) / config.gap_cm()
            : std::numeric_limits<double>::quiet_NaN();
    double primary_avalanche_alpha_error =
        primary_avalanche_gain > 1.0 &&
                std::isfinite(primary_avalanche_gain_error)
            ? primary_avalanche_gain_error /
                  (primary_avalanche_gain * config.gap_cm())
            : std::numeric_limits<double>::quiet_NaN();

    for (const float energy : g_energy.values) h_electron_energy.Fill(energy);

    photon_summary.accounted_transport_outcomes =
        photon_summary.nonionising_photoabsorptions +
        photon_summary.photoionisations + photon_summary.cathode_impacts +
        photon_summary.anode_impacts + photon_summary.lateral_escapes +
        photon_summary.geometric_rejections +
        photon_summary.transport_failures;
    photon_summary.unaccounted_transport_candidates =
        photon_summary.transport_candidates -
        photon_summary.accounted_transport_outcomes;
    const double closure_tolerance = std::max(
        1.0e-8,
        1.0e-8 * std::max(1.0, photon_summary.transport_candidates));
    if (std::abs(photon_summary.unaccounted_transport_candidates) <=
        closure_tolerance) {
      photon_summary.unaccounted_transport_candidates = 0.0;
    }

    // Keep gasData deliberately minimal: only the information required to
    // identify and reproduce the gas conditions of this run. Physics results
    // are stored in the dedicated trees and histograms.
    Long64_t actual_npe_root = actual_npe;
    UInt_t random_seed_root = config.seed;

    gas_data.Branch("gas1", &config.gas1);
    gas_data.Branch("composition1_pct", &config.composition1,
                    "composition1_pct/D");
    gas_data.Branch("gas2", &config.gas2);
    gas_data.Branch("composition2_pct", &config.composition2,
                    "composition2_pct/D");
    gas_data.Branch("gas3", &config.gas3);
    gas_data.Branch("composition3_pct", &config.composition3,
                    "composition3_pct/D");
    gas_data.Branch("pressure_bar", &config.pressure_bar, "pressure_bar/D");
    gas_data.Branch("temperature_K", &config.temperature_k,
                    "temperature_K/D");
    gas_data.Branch("electricField_V_cm", &config.field_v_cm,
                    "electricField_V_cm/D");
    gas_data.Branch("gap_mm", &config.gap_mm, "gap_mm/D");
    gas_data.Branch("npe", &actual_npe_root, "npe/L");
    gas_data.Branch("randomSeed", &random_seed_root, "randomSeed/i");
    gas_data.Fill();

    // Only the two gas-interaction totals that are not represented as simple
    // cathode/photoelectron histograms are kept in photonTransportData.
    // Double_t is intentional because weighted MC runs can be non-integer.
    Double_t photoabsorptions = photon_summary.photoabsorptions;
    Double_t photoionisations = photon_summary.photoionisations;

    TTree photon_transport_data("photonTransportData",
                                "Photon gas-interaction totals");
    photon_transport_data.Branch("nPhotoabsorptions", &photoabsorptions,
                                 "nPhotoabsorptions/D");
    photon_transport_data.Branch("nPhotoionisations", &photoionisations,
                                 "nPhotoionisations/D");
    photon_transport_data.Fill();

    // ========================================================================
    //                              Save ROOT
    // ========================================================================

    TFile output(config.root_file.c_str(), "RECREATE");
    if (output.IsZombie()) {
      throw std::runtime_error("Could not create ROOT file: " +
                               config.root_file);
    }
    output.cd();
    h_electron_energy.Write();
    h_levels.Write();
    if (h_exc_xyz != nullptr) h_exc_xyz->Write();
    if (h_exc_zt != nullptr) h_exc_zt->Write();
    data_per_primary.Write();
    data_per_avalanche.Write();
    data_per_electron.Write();
    gas_data.Write();
    photon_transport_data.Write();
    photon_hists.electrons_vs_time.Write();
    photon_hists.electrons_vs_time_full.Write();
    photon_hists.avalanche_electrons_vs_time.Write();
    photon_hists.avalanche_electrons_vs_time_prompt.Write();
    photon_hists.avalanche_electrons_vs_time_generation.Write();
    photon_hists.spectra.Write();
    photon_hists.photon_xyz.Write();
    photon_hists.photon_wavelength_time.Write();
    photon_hists.cos_theta.Write();
    photon_hists.phi.Write();
    photon_hists.qe.Write();
    photon_hists.impacts_xy.Write();
    photon_hists.photoelectron_xy.Write();
    photon_hists.photoelectron_time_energy.Write();
    photon_hists.photoabsorption_xyz.Write();
    photon_hists.photoabsorption_zt.Write();
    photon_hists.photoionisation_xyz.Write();
    photon_hists.photoionisation_zt.Write();
    photon_hists.feedback_generation.Write();

    // All ROOT objects above are owned by C++ stack variables/unique_ptrs.
    // Detach them before TFile::Close so the file cannot delete them first and
    // leave their C++ owners with dangling pointers (the classic exit-139
    // double-destruction failure).
    h_electron_energy.SetDirectory(nullptr);
    h_levels.SetDirectory(nullptr);
    if (h_exc_xyz != nullptr) h_exc_xyz->SetDirectory(nullptr);
    if (h_exc_zt != nullptr) h_exc_zt->SetDirectory(nullptr);
    data_per_primary.SetDirectory(nullptr);
    data_per_avalanche.SetDirectory(nullptr);
    data_per_electron.SetDirectory(nullptr);
    gas_data.SetDirectory(nullptr);
    photon_transport_data.SetDirectory(nullptr);
    photon_hists.electrons_vs_time.SetDirectory(nullptr);
    photon_hists.electrons_vs_time_full.SetDirectory(nullptr);
    photon_hists.avalanche_electrons_vs_time.SetDirectory(nullptr);
    photon_hists.avalanche_electrons_vs_time_prompt.SetDirectory(nullptr);
    photon_hists.avalanche_electrons_vs_time_generation.SetDirectory(nullptr);
    photon_hists.spectra.SetDirectory(nullptr);
    photon_hists.photon_xyz.SetDirectory(nullptr);
    photon_hists.photon_wavelength_time.SetDirectory(nullptr);
    photon_hists.cos_theta.SetDirectory(nullptr);
    photon_hists.phi.SetDirectory(nullptr);
    photon_hists.qe.SetDirectory(nullptr);
    photon_hists.impacts_xy.SetDirectory(nullptr);
    photon_hists.photoelectron_xy.SetDirectory(nullptr);
    photon_hists.photoelectron_time_energy.SetDirectory(nullptr);
    photon_hists.photoabsorption_xyz.SetDirectory(nullptr);
    photon_hists.photoabsorption_zt.SetDirectory(nullptr);
    photon_hists.photoionisation_xyz.SetDirectory(nullptr);
    photon_hists.photoionisation_zt.SetDirectory(nullptr);
    photon_hists.feedback_generation.SetDirectory(nullptr);
    output.Close();

    const double wall_time = std::difftime(std::time(nullptr), wall_start);
    const double cpu_time =
        static_cast<double>(std::clock() - cpu_start) / CLOCKS_PER_SEC;

    std::cout << std::setprecision(8)
              << "RESULT gain=" << gain << " gainError=" << gain_error
              << " primaryGain=" << primary_avalanche_gain
              << " primaryGainError=" << primary_avalanche_gain_error
              << " alphaEffective=" << alpha_effective
              << " alphaError=" << alpha_error << " npe=" << actual_npe
              << "\n"
              << "Avalanches: " << total_avalanches << "\n"
              << "Photons: " << photon_summary.generated_photons << "\n"
              << "Gas absorptions: " << photon_summary.photoabsorptions << "\n"
              << "Photoionisations: " << photon_summary.photoionisations << "\n"
              << "Cathode impacts: " << photon_summary.cathode_impacts << "\n"
              << "Photoelectrons: "
              << photon_hists.photoelectron_xy.Integral() << "\n"
              << "Feedback PE seeds: "
              << photon_summary.feedback_photoelectron_seeds << "\n"
              << "Feedback gas-ionisation seeds: "
              << photon_summary.feedback_photoionisation_seeds << "\n"
              << "Wall time: " << wall_time << " s\n"
              << "CPU time: " << cpu_time << " s\n"
              << "DONE " << config.job_id << std::endl;

    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "[secondaryAvalanches:ERROR] " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
