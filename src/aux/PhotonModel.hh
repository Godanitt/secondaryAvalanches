#ifndef SECONDARY_AVALANCHES_PHOTON_MODEL_HH
#define SECONDARY_AVALANCHES_PHOTON_MODEL_HH

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <TRandom3.h>

namespace photonemission {

struct GasData {
  std::string gas1;
  std::string gas2;
  std::string gas3;
  double composition1_percent = 0.0;
  double composition2_percent = 0.0;
  double composition3_percent = 0.0;
  double pressure_bar = 0.0;
  double temperature_k = 293.15;
  double electric_field_v_cm = 0.0;
  double gap_mm = 0.0;
  double start_z_fraction_from_anode = 1.0;
  double start_z_mm = std::numeric_limits<double>::quiet_NaN();
  double avalanche_distance_mm = std::numeric_limits<double>::quiet_NaN();
  bool space_charge_enabled = false;
  int npe = 0;
  long long ne_total = 0;
  long long ni_total = 0;
};

struct PrimaryGainRow {
  int primary_id = 0;
  int ne = 0;
  int ni = 0;
  double gain = 0.0;
};

struct LevelInfo {
  int level = -1;
  std::string gas;
  std::string state_name;
  std::string type;
  double energy_ev = std::numeric_limits<double>::quiet_NaN();
};

struct LevelPopulation {
  int level = -1;
  double n_events = 0.0;
  bool known_level = false;
  LevelInfo level_info;
};

struct SourceSelector {
  std::string source_gas;
  std::vector<std::string> state_tokens;
  double min_energy_ev = -1.0e99;
  double max_energy_ev = 1.0e99;
};

struct SpectralPeak {
  double center_nm = 0.0;
  double sigma_nm = 1.0;
  double weight = 1.0;
};

enum class KineticStageKind { Excitation, Emission };

struct KineticStage {
  std::string name;
  KineticStageKind kind = KineticStageKind::Excitation;
  // Effective lifetime of this state after every competing channel has been
  // included: tau_tot = 1 / sum_j Gamma_j.  The branching probability is
  // already accounted for in KineticComponent::expected_photons.
  double tau_total_ns = 0.0;
};

struct KineticComponent {
  std::string name;
  std::string region;
  SourceSelector source;
  std::vector<SpectralPeak> peaks;
  double source_population = 0.0;
  double expected_photons = 0.0;
  double probability = 0.0;

  // A photon is emitted at
  //   t_gamma = t_Garfield_collision + tau_exc(sampled) + tau_em(sampled).
  // Intermediate cascade/transfer/formation stages contribute to tau_exc;
  // the final radiating state contributes to tau_em.  Every waiting time is
  // sampled independently from Exp(tau_total) because sequential stages do
  // not combine into one exponential distribution.
  std::vector<KineticStage> kinetic_stages;
  double mean_excitation_delay_ns = 0.0;
  double mean_emission_delay_ns = 0.0;
  double mean_total_delay_ns = 0.0;
};

struct EmissionDelaySample {
  double excitation_delay_ns = 0.0;
  double emission_delay_ns = 0.0;

  double total_delay_ns() const {
    return excitation_delay_ns + emission_delay_ns;
  }
};

inline std::string lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

inline std::string trim_copy(std::string value) {
  auto not_space = [](unsigned char c) { return !std::isspace(c); };
  value.erase(value.begin(), std::find_if(value.begin(), value.end(), not_space));
  value.erase(std::find_if(value.rbegin(), value.rend(), not_space).base(), value.end());
  return value;
}

inline std::string collapse_spaces(std::string value) {
  std::string out;
  bool in_space = false;
  for (unsigned char c : value) {
    if (std::isspace(c)) {
      if (!in_space && !out.empty()) out.push_back(' ');
      in_space = true;
    } else {
      out.push_back(static_cast<char>(c));
      in_space = false;
    }
  }
  return trim_copy(out);
}

inline std::string normalise_text_for_match(const std::string& value) {
  return collapse_spaces(lower(value));
}

inline std::string normalise_gas_name(const std::string& name) {
  std::string s = normalise_text_for_match(name);
  s.erase(std::remove_if(s.begin(), s.end(), [](unsigned char c) {
    return c == '-' || c == '_' || std::isspace(c);
  }), s.end());
  if (s == "argon") return "ar";
  if (s == "nitrogen") return "n2";
  if (s == "helium") return "he";
  if (s == "methane") return "ch4";
  if (s == "carbondioxide") return "co2";
  if (s == "isobutane" || s == "isobutano" || s == "iso" ||
      s == "c4h10" || s == "ic4h10") return "ic4h10";
  return s;
}

inline bool contains_gas(const GasData& gas, const std::string& token) {
  const std::string t = normalise_gas_name(token);
  return normalise_gas_name(gas.gas1) == t || normalise_gas_name(gas.gas2) == t ||
         normalise_gas_name(gas.gas3) == t;
}

inline std::string mixture_name(const GasData& gas) {
  if (contains_gas(gas, "ar") && contains_gas(gas, "cf4") &&
      contains_gas(gas, "ic4h10")) return "ArIsoCF4";
  if (contains_gas(gas, "ar") && contains_gas(gas, "cf4")) return "ArCF4";
  if (contains_gas(gas, "ar") && contains_gas(gas, "n2")) return "ArN2";
  if (contains_gas(gas, "ar") && contains_gas(gas, "ic4h10")) return "ArIso";
  if (contains_gas(gas, "ar") && contains_gas(gas, "ch4")) return "ArCH4";
  if (contains_gas(gas, "ar") && contains_gas(gas, "co2")) return "ArCO2";
  if (contains_gas(gas, "he") && contains_gas(gas, "cf4")) return "HeCF4";
  if (contains_gas(gas, "cf4")) return "CF4";
  if (contains_gas(gas, "ar")) return "Ar";
  if (contains_gas(gas, "n2")) return "N2";
  if (contains_gas(gas, "ic4h10")) return "Iso";
  return "unknown";
}

inline double fraction_of(const GasData& gas, const std::string& token) {
  const std::string t = normalise_gas_name(token);
  if (normalise_gas_name(gas.gas1) == t) return gas.composition1_percent / 100.0;
  if (normalise_gas_name(gas.gas2) == t) return gas.composition2_percent / 100.0;
  if (normalise_gas_name(gas.gas3) == t) return gas.composition3_percent / 100.0;
  return 0.0;
}

inline double gaussian_pdf(const double x, const double mu, const double sigma) {
  if (sigma <= 0.0) return 0.0;
  static constexpr double inv_sqrt_2pi = 0.39894228040143267794;
  const double u = (x - mu) / sigma;
  return inv_sqrt_2pi / sigma * std::exp(-0.5 * u * u);
}

inline double clamp_value(const double x, const double lo, const double hi) {
  return std::max(lo, std::min(hi, x));
}

inline double safe_ratio(const double numerator, const double denominator) {
  if (!std::isfinite(numerator) || !std::isfinite(denominator) || std::abs(denominator) <= 0.0) return 0.0;
  return numerator / denominator;
}

inline bool state_matches_tokens(const std::string& state_name,
                                const std::vector<std::string>& tokens) {
  if (tokens.empty()) return true;
  const std::string state = normalise_text_for_match(state_name);
  for (const auto& token : tokens) {
    const std::string t = normalise_text_for_match(token);
    if (!t.empty() && state.find(t) == std::string::npos) return false;
  }
  return true;
}

inline bool selector_matches_level(const SourceSelector& selector, const LevelInfo& info) {
  if (!selector.source_gas.empty() &&
      normalise_gas_name(info.gas) != normalise_gas_name(selector.source_gas)) return false;
  if (!state_matches_tokens(info.state_name, selector.state_tokens)) return false;
  const double e = info.energy_ev;
  if (std::isfinite(e) && (e <= selector.min_energy_ev || e >= selector.max_energy_ev)) return false;
  return true;
}

inline std::vector<SpectralPeak> peaks_single(const double center, const double sigma = 2.5) {
  return {{center, sigma, 1.0}};
}

inline std::vector<SpectralPeak> peaks_cf4_uv() {
  return {{235.0, 17.0, 0.55}, {290.0, 17.0, 0.75}, {364.0, 50.0, 0.35}};
}

inline std::vector<SpectralPeak> peaks_cf3_uv() {
  return {{260.0, 30, 1.0}};
}

inline double fwhm_to_sigma(const double fwhm) {
  return fwhm > 0.0 ? fwhm / 2.3548200450309493 : 1.0;
}

inline std::vector<SpectralPeak> peaks_cf4_vuv_branch(const double center_nm = 155.0,
                                                          const double fwhm_nm = 10.0) {
  return {{center_nm, fwhm_to_sigma(fwhm_nm), 1.0}};
}

inline std::vector<SpectralPeak> peaks_cf4_vuv_150() {
  return peaks_cf4_vuv_branch(150.0, 10.0);
}

inline std::vector<SpectralPeak> peaks_ar_second_continuum(const double center_nm = 128.0,
                                                           const double fwhm_nm = 10.0) {
  return {{center_nm, fwhm_to_sigma(fwhm_nm), 1.0}};
}

inline std::vector<SpectralPeak> peaks_cf3_visible() {
  // Broad CF3/CF4 visible structure; the exact detailed spectrum can be swapped
  // later by reading a spectral-shape CSV without changing the ROOT output format.
  return {{630.0, 70.0, 1.00}};
}

inline std::vector<SpectralPeak> peaks_ar_3rd_uv() {
  return {{176.0, 30.0, 1.0}, {188.0, 30.0, 1.0}, {199.0, 30.0, 1.0},
          {212.0, 30.0, 1.0}, {225.0, 30.0, 1.0}, {245.0, 30.0, 1.0}};
}

inline std::vector<SpectralPeak> peaks_n2_second_positive() {
  return {{336.5, 3.75, 0.42}, {359.6, 3.75, 0.30}, {379.0, 3.75, 0.10},
          {404.0, 3.75, 0.05}, {431.0, 3.75, 0.05}};
}


inline std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> out;
  std::string current;
  bool in_quotes = false;
  for (std::size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (c == '"') {
      if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
        current.push_back('"');
        ++i;
      } else {
        in_quotes = !in_quotes;
      }
    } else if (c == ',' && !in_quotes) {
      out.push_back(current);
      current.clear();
    } else {
      current.push_back(c);
    }
  }
  out.push_back(current);
  return out;
}

inline std::string trim(std::string value) {
  auto not_space = [](unsigned char c) { return !std::isspace(c); };
  value.erase(value.begin(), std::find_if(value.begin(), value.end(), not_space));
  value.erase(std::find_if(value.rbegin(), value.rend(), not_space).base(), value.end());
  return value;
}

inline double population_for_selector(const SourceSelector& selector,
                                      const std::vector<LevelPopulation>& populations) {
  double total = 0.0;
  for (const auto& pop : populations) {
    if (!pop.known_level || pop.n_events <= 0.0) continue;
    if (!selector_matches_level(selector, pop.level_info)) continue;
    total += pop.n_events;
  }
  return total;
}

inline bool csv_value_is_false(const std::string& value) {
  const std::string v = normalise_text_for_match(trim(value));
  return v == "0" || v == "false" || v == "no" || v == "off" ||
         v == "disabled" || v == "inactive";
}

inline std::string csv_additive_suffix(std::string additive) {
  additive = lower(trim(additive));
  additive.erase(std::remove_if(additive.begin(), additive.end(),
                                [](unsigned char c) {
                                  return c == '-' || c == '_' ||
                                         std::isspace(c);
                                }),
                 additive.end());
  if (additive == "cf4") return "CF4";
  if (additive == "n2" || additive == "nitrogen") return "N2";
  if (additive == "co2" || additive == "carbondioxide") return "CO2";
  if (additive == "ch4" || additive == "methane") return "CH4";
  if (additive == "ic4h10" || additive == "c4h10" ||
      additive == "isobutane" || additive == "isobutano" ||
      additive == "iso") {
    return "IC4H10";
  }
  std::transform(additive.begin(), additive.end(), additive.begin(),
                 [](unsigned char c) {
                   return static_cast<char>(std::toupper(c));
                 });
  return additive;
}

inline std::string scoped_additive_parameter_name(
    const std::string& name, const std::string& additive_suffix) {
  static constexpr const char* unit_suffix = "_m3_s";
  if (name.size() >= 5 &&
      name.compare(name.size() - 5, 5, unit_suffix) == 0) {
    return name.substr(0, name.size() - 5) + "_" + additive_suffix +
           unit_suffix;
  }
  return name + "_" + additive_suffix;
}

inline std::map<std::string, double> read_parameter_csv_values(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("No puedo abrir parameter CSV: " + path);

  std::string header;
  if (!std::getline(in, header)) throw std::runtime_error("Parameter CSV vacío: " + path);
  const auto columns = split_csv_line(header);
  std::map<std::string, std::size_t> index;
  for (std::size_t i = 0; i < columns.size(); ++i) index[trim(columns[i])] = i;

  auto optional_col = [&](const std::string& name) -> std::size_t {
    const auto it = index.find(name);
    return it == index.end() ? std::numeric_limits<std::size_t>::max() : it->second;
  };

  const std::size_t absent = std::numeric_limits<std::size_t>::max();
  const std::size_t i_name = optional_col("name") != absent
                               ? optional_col("name")
                               : optional_col("parameter");
  const std::size_t i_value = optional_col("value") != absent
                                ? optional_col("value")
                                : optional_col("secondary_optimal");
  const std::size_t i_scope = optional_col("scope");
  const std::size_t i_additive = optional_col("additive");
  const std::size_t i_enabled = optional_col("enabled");
  if (i_name == absent) {
    throw std::runtime_error("Columna ausente en parameter CSV: name/parameter");
  }
  if (i_value == absent) {
    throw std::runtime_error("Columna ausente en parameter CSV: value/secondary_optimal");
  }

  std::map<std::string, double> out;
  std::string line;
  while (std::getline(in, line)) {
    if (trim(line).empty()) continue;
    const auto fields = split_csv_line(line);
    if (fields.size() <= std::max(i_name, i_value)) continue;
    if (i_enabled != absent && i_enabled < fields.size() &&
        csv_value_is_false(fields[i_enabled])) {
      continue;
    }

    std::string name = trim(fields[i_name]);
    const std::string value_text = trim(fields[i_value]);
    if (name.empty() || value_text.empty()) continue;

    // Compatibility with the first resonant-line patch.
    if (name == "r_Ar_1s5_direct") name = "r";

    std::string storage_name = name;
    if (i_scope != absent && i_scope < fields.size() &&
        normalise_text_for_match(fields[i_scope]) == "additive") {
      if (i_additive == absent || i_additive >= fields.size()) {
        throw std::runtime_error(
            "Fila aditiva sin columna additive en parameter CSV: " + path);
      }
      const std::string suffix = csv_additive_suffix(fields[i_additive]);
      if (suffix.empty()) {
        throw std::runtime_error(
            "Fila aditiva sin gas en parameter CSV: " + path);
      }
      storage_name = scoped_additive_parameter_name(name, suffix);
    }

    try {
      out[storage_name] = std::stod(value_text);
    } catch (...) {
      // Empty/planned rows are intentionally ignored. Active numeric rows are
      // the only values consumed by the C++ kinetic model.
    }
  }
  return out;
}

inline double pget(const std::map<std::string, double>& p, const std::string& name,
                   const double fallback = 0.0) {
  const auto it = p.find(name);
  if (it == p.end() || !std::isfinite(it->second)) return fallback;
  return it->second;
}

inline bool file_exists_readable(const std::string& path) {
  std::ifstream in(path);
  return static_cast<bool>(in);
}

inline std::string first_existing_csv(const std::vector<std::string>& candidates) {
  for (const auto& path : candidates) {
    if (file_exists_readable(path)) return path;
  }
  return candidates.empty() ? std::string{} : candidates.back();
}

inline std::string parameter_csv_for_model(const std::string& parameter_dir,
                                           const std::string& mixture,
                                           const bool infrared) {
  const std::string base = parameter_dir.empty() ? "." : parameter_dir;

  // For avalanche/secondary scintillation prefer the secondary-fit tables when
  // they are present.  Old archives only shipped the primary-fit CSVs, so keep a
  // primary fallback rather than silently failing.
  if (mixture == "ArN2") {
    return first_existing_csv(infrared
        ? std::vector<std::string>{base + "/ArN2_IR_secondary.csv", base + "/ArN2_IR_primary_secondary.csv", base + "/ArN2_IR_primary.csv"}
        : std::vector<std::string>{base + "/ArN2_secondary.csv", base + "/ArN2_primary_secondary.csv", base + "/ArN2_primary.csv"});
  }
  return first_existing_csv(infrared
      ? std::vector<std::string>{base + "/ArCF4_IR_secondary.csv", base + "/ArCF4_IR_primary_secondary.csv", base + "/ArCF4_IR_primary.csv"}
      : std::vector<std::string>{base + "/ArCF4_secondary.csv", base + "/ArCF4_primary_secondary.csv", base + "/ArCF4_primary.csv"});
}

inline std::string ar2nd_parameter_csv(const std::string& parameter_dir) {
  const std::string base = parameter_dir.empty() ? "." : parameter_dir;
  // Prefer the full scintillation-model table.  The old minimal
  // Ar2nd_continuum_secondary.csv only had f_S/Y_anchor-like placeholders and
  // is deliberately kept last for backwards compatibility.
  return first_existing_csv({base + "/Ar2nd_continium.csv",
                             base + "/Ar2nd_continuum.csv",
                             base + "/Ar2nd_continuum_secondary.csv"});
}

inline std::map<std::string, double> read_optional_parameter_csv_values(const std::string& path) {
  if (path.empty() || !file_exists_readable(path)) return {};
  try {
    return read_parameter_csv_values(path);
  } catch (...) {
    return {};
  }
}

inline double secondary_normalization_reference(const std::string& parameter_dir) {
  // The secondary projection uses the primary-fit formulae as optical weights,
  // but the avalanche populations are already absolute. Keep the Nnorm factors
  // inside the original UV/VIS models, then divide every component by the same
  // ArCF4 Nnorm reference at the end. This also puts IR components, which do
  // not carry an explicit Nnorm in their parameterization, on the same scale.
  try {
    const auto p = read_parameter_csv_values(parameter_csv_for_model(parameter_dir, "ArCF4", false));
    const double nref = pget(p, "Nnorm", 1.0);
    if (std::isfinite(nref) && nref > 0.0) return nref;
  } catch (...) {
  }
  return 1.0;
}

inline void apply_secondary_normalization(std::vector<KineticComponent>& components,
                                          const std::string& parameter_dir) {
  const double nref = secondary_normalization_reference(parameter_dir);
  if (!std::isfinite(nref) || nref <= 0.0) return;
  for (auto& component : components) {
    // Ar2nd is a kinetic probability applied directly to Garfield hLevel counts,
    // not an optical branch fitted with the primary Nnorm.
    if (component.name.find("Ar2nd") != std::string::npos) continue;
    component.expected_photons /= nref;
  }
}

inline KineticStage excitation_stage(std::string name,
                                      const double tau_total_ns) {
  return {std::move(name), KineticStageKind::Excitation, tau_total_ns};
}

inline KineticStage emission_stage(std::string name,
                                    const double tau_total_ns) {
  return {std::move(name), KineticStageKind::Emission, tau_total_ns};
}

inline void push_component(std::vector<KineticComponent>& out,
                           std::string name,
                           std::string region,
                           SourceSelector source,
                           std::vector<SpectralPeak> peaks,
                           const double source_population,
                           const double expected_photons_one_sample,
                           const double probability,
                           const std::vector<KineticStage>& kinetic_stages,
                           const int /*mc_samples*/) {
  if (!std::isfinite(expected_photons_one_sample) ||
      expected_photons_one_sample <= 0.0) {
    return;
  }
  KineticComponent c;
  c.name = std::move(name);
  c.region = std::move(region);
  c.source = std::move(source);
  c.peaks = std::move(peaks);
  c.source_population = source_population;
  c.expected_photons = expected_photons_one_sample;
  c.probability = probability;

  for (const auto& stage : kinetic_stages) {
    if (!std::isfinite(stage.tau_total_ns) || stage.tau_total_ns <= 0.0) {
      continue;
    }
    c.kinetic_stages.push_back(stage);
    c.mean_total_delay_ns += stage.tau_total_ns;
    if (stage.kind == KineticStageKind::Emission) {
      c.mean_emission_delay_ns += stage.tau_total_ns;
    } else {
      c.mean_excitation_delay_ns += stage.tau_total_ns;
    }
  }

  // Every generated component must end in a real radiative state.  This is a
  // model-construction error, not a condition that should be silently fixed.
  const bool has_emission_stage = std::any_of(
      c.kinetic_stages.begin(), c.kinetic_stages.end(),
      [](const KineticStage& stage) {
        return stage.kind == KineticStageKind::Emission;
      });
  if (!has_emission_stage) {
    throw std::runtime_error("Photon component without final emission stage: " +
                             c.name);
  }
  out.push_back(std::move(c));
}

// Convenience overload for existing multi-stage channels.  The final waiting
// time belongs to the radiating state; every previous waiting time is an
// excitation/transfer/cascade stage.  Explicitly named stages are preferred
// for physically important channels such as the Ar second continuum.
inline void push_component(std::vector<KineticComponent>& out,
                           std::string name,
                           std::string region,
                           SourceSelector source,
                           std::vector<SpectralPeak> peaks,
                           const double source_population,
                           const double expected_photons_one_sample,
                           const double probability,
                           const std::vector<double>& lifetime_stages_ns,
                           const int mc_samples) {
  std::vector<KineticStage> stages;
  stages.reserve(lifetime_stages_ns.size());
  for (std::size_t i = 0; i < lifetime_stages_ns.size(); ++i) {
    const bool final_stage = i + 1 == lifetime_stages_ns.size();
    stages.push_back(final_stage
                         ? emission_stage("final emitting state",
                                          lifetime_stages_ns[i])
                         : excitation_stage("cascade / transfer",
                                            lifetime_stages_ns[i]));
  }
  push_component(out, std::move(name), std::move(region), std::move(source),
                 std::move(peaks), source_population,
                 expected_photons_one_sample, probability, stages, mc_samples);
}

inline void push_component(std::vector<KineticComponent>& out,
                           std::string name,
                           std::string region,
                           SourceSelector source,
                           std::vector<SpectralPeak> peaks,
                           const double source_population,
                           const double expected_photons_one_sample,
                           const double probability,
                           const double lifetime_ns,
                           const int mc_samples) {
  push_component(out, std::move(name), std::move(region), std::move(source),
                 std::move(peaks), source_population,
                 expected_photons_one_sample, probability,
                 std::vector<KineticStage>{
                     emission_stage("final emitting state", lifetime_ns)},
                 mc_samples);
}

inline double lifetime_from_total_rate_ns(const double total_rate_ns_inv) {
  if (!std::isfinite(total_rate_ns_inv) || total_rate_ns_inv <= 0.0) return 0.0;
  return 1.0 / total_rate_ns_inv;
}

inline double positive_rate(const double value) {
  return std::isfinite(value) ? std::max(0.0, value) : 0.0;
}

inline bool has_finite_parameter(const std::map<std::string, double>& p,
                                 const std::string& name) {
  const auto it = p.find(name);
  return it != p.end() && std::isfinite(it->second);
}

inline double ar2nd_two_body_rate_ns_inv(const std::map<std::string, double>& p2nd,
                                         const std::string& si_name,
                                         const double default_si_m3_s,
                                         const std::string& legacy_name = "") {
  if (!legacy_name.empty() && !has_finite_parameter(p2nd, si_name) &&
      has_finite_parameter(p2nd, legacy_name)) {
    return positive_rate(pget(p2nd, legacy_name));
  }
  constexpr double seconds_per_nanosecond = 1.0e-9;
  const double loschmidt_m3 = positive_rate(pget(p2nd, "loschmidt_m3", 2.6868e25));
  const double k_m3_s = positive_rate(pget(p2nd, si_name, default_si_m3_s));
  return k_m3_s * loschmidt_m3 * seconds_per_nanosecond;
}

inline double ar2nd_three_body_rate_ns_inv(const std::map<std::string, double>& p2nd,
                                           const std::string& si_name,
                                           const double default_si_m6_s,
                                           const std::string& legacy_name = "") {
  if (!legacy_name.empty() && !has_finite_parameter(p2nd, si_name) &&
      has_finite_parameter(p2nd, legacy_name)) {
    return positive_rate(pget(p2nd, legacy_name));
  }
  constexpr double seconds_per_nanosecond = 1.0e-9;
  const double loschmidt_m3 = positive_rate(pget(p2nd, "loschmidt_m3", 2.6868e25));
  const double k_m6_s = positive_rate(pget(p2nd, si_name, default_si_m6_s));
  return k_m6_s * loschmidt_m3 * loschmidt_m3 * seconds_per_nanosecond;
}

struct Ar2ndKineticRates {
  double upper_ar = 0.0;
  double four_s_ar = 0.0;
  double excimer_formation = 0.0;
  // These are already composition weighted: sum_i f_i K_i.
  double upper_additives = 0.0;
  double four_s_additives = 0.0;
  double excimer_additives = 0.0;
  double total_additive_fraction = 0.0;
  double argon_fraction = 0.0;
};

struct Ar2ndAdditiveDefaults {
  double upper_m3_s = 0.0;
  double four_s_m3_s = 0.0;
  double excimer_m3_s = 0.0;
};

inline std::string ar2nd_additive_suffix(const std::string& gas_name) {
  const std::string gas = normalise_gas_name(gas_name);
  if (gas == "cf4") return "CF4";
  if (gas == "n2") return "N2";
  if (gas == "ch4") return "CH4";
  if (gas == "ic4h10") return "IC4H10";
  return "";
}

inline Ar2ndAdditiveDefaults ar2nd_additive_defaults(const std::string& suffix) {
  if (suffix == "CF4") return {1.80e-16, 3.00e-17, 3.00e-17};
  if (suffix == "N2") return {1.58e-16, 3.60e-17, 2.24e-18};
  if (suffix == "CH4") return {6.53e-16, 3.30e-16, 3.30e-16};
  // Explicit effective approximation adopted in the current model: the
  // measured Ar(1s5)+iC4H10 coefficient is used for Ar**, Ar(4s), and Ar2*.
  if (suffix == "IC4H10") return {7.10e-16, 7.10e-16, 7.10e-16};
  throw std::runtime_error(
      "No hay cinética completa del segundo continuo para el aditivo " + suffix + ".");
}

inline std::map<std::string, double> ar2nd_additive_fractions(const GasData& gas) {
  std::map<std::string, double> fractions;
  const auto add = [&](const std::string& name, const double percent) {
    const std::string normalised = normalise_gas_name(name);
    const double fraction = clamp_value(percent / 100.0, 0.0, 1.0);
    if (normalised.empty() || normalised == "ar" || fraction <= 0.0) return;
    fractions[normalised] += fraction;
  };
  add(gas.gas1, gas.composition1_percent);
  add(gas.gas2, gas.composition2_percent);
  add(gas.gas3, gas.composition3_percent);
  return fractions;
}

inline double ar2nd_additive_rate_ns_inv(
    const std::map<std::string, double>& p2nd,
    const std::string& base_name,
    const std::string& suffix,
    const double fallback_si_m3_s) {
  return ar2nd_two_body_rate_ns_inv(
      p2nd, base_name + "_" + suffix + "_m3_s", fallback_si_m3_s);
}

inline Ar2ndKineticRates ar2nd_rates_for_mixture(
    const std::map<std::string, double>& p2nd,
    const GasData& gas) {
  Ar2ndKineticRates rates;
  rates.upper_ar = ar2nd_two_body_rate_ns_inv(
      p2nd, "k_Ar_dbleStar_Q_Ar_m3_s", 1.63e-17, "K_Ar_dbleStar_Q_Ar");
  rates.four_s_ar = ar2nd_two_body_rate_ns_inv(
      p2nd, "k_Ar_4s_Q_Ar_m3_s", 1.63e-17, "K_Ar_4s_Q_Ar");
  rates.excimer_formation = ar2nd_three_body_rate_ns_inv(
      p2nd, "k_Ar_4s_Q_2Ar_m6_s", 1.00e-44, "K_Ar_star_Q_2Ar");

  const auto fractions = ar2nd_additive_fractions(gas);
  for (const auto& [additive, fraction] : fractions) {
    const std::string suffix = ar2nd_additive_suffix(additive);
    if (suffix.empty()) {
      throw std::runtime_error(
          "Faltan coeficientes del segundo continuo para el aditivo " + additive +
          ". No se asumirá quenching cero.");
    }
    const auto defaults = ar2nd_additive_defaults(suffix);
    rates.upper_additives += fraction * ar2nd_additive_rate_ns_inv(
        p2nd, "k_Ar_dbleStar_Q", suffix, defaults.upper_m3_s);
    rates.four_s_additives += fraction * ar2nd_additive_rate_ns_inv(
        p2nd, "k_Ar_4s_Q", suffix, defaults.four_s_m3_s);
    rates.excimer_additives += fraction * ar2nd_additive_rate_ns_inv(
        p2nd, "k_Ar2star_Q", suffix, defaults.excimer_m3_s);
    rates.total_additive_fraction += fraction;
  }

  if (rates.total_additive_fraction > 1.0 + 1.0e-10) {
    throw std::runtime_error("La suma de fracciones aditivas supera 1.");
  }
  // Same convention as the current ScintillationModel: the productive argon
  // fraction is one minus the sum of every non-Ar component.
  rates.argon_fraction = clamp_value(1.0 - rates.total_additive_fraction, 0.0, 1.0);
  return rates;
}

inline void ar2nd_formation_fractions(const std::map<std::string, double>& p2nd,
                                      double& f_singlet,
                                      double& f_triplet) {
  f_singlet = clamp_value(pget(p2nd, "f_1Sigma", pget(p2nd, "f_S", 0.1)), 0.0, 1.0);
  f_triplet = clamp_value(pget(p2nd, "f_3Sigma", 1.0 - f_singlet), 0.0, 1.0);
  const double norm = f_singlet + f_triplet;
  if (norm <= 0.0) {
    f_singlet = 0.1;
    f_triplet = 0.9;
  } else {
    f_singlet /= norm;
    f_triplet /= norm;
  }
}

inline void push_ar2nd_components(std::vector<KineticComponent>& out,
                                  const GasData& gas,
                                  const std::vector<LevelPopulation>& populations,
                                  const std::map<std::string, double>& p2nd,
                                  const int mc_samples) {
  const double n = std::max(gas.pressure_bar, 0.0);
  if (n <= 0.0 || fraction_of(gas, "ar") <= 0.0) return;

  const Ar2ndKineticRates rates = ar2nd_rates_for_mixture(p2nd, gas);
  const double f_ar = rates.argon_fraction;
  if (f_ar <= 0.0) return;

  // The supplied model/CSV stores the 11.548 eV metastable population as
  // Ar_1s4 and the 11.624 eV resonant population as Ar_1s5.  Because the
  // Paschen labels are interchanged in some conventions, component names use
  // the unambiguous energy/state descriptions instead of relying on 1s4/1s5.
  const double energy_1s4_ev = pget(p2nd, "energy_Ar_1s4_eV", 11.548);
  const double energy_1s5_ev = pget(p2nd, "energy_Ar_1s5_eV", 11.624);
  const double lower_state_split_ev =
      (energy_1s4_ev < energy_1s5_ev)
          ? 0.5 * (energy_1s4_ev + energy_1s5_ev)
          : 11.586;

  // Disjoint source populations. A Garfield collision can enter exactly one
  // kinetic source family; this prevents double counting in the VUV model.
  const SourceSelector sel_ar_1s4{"Ar", {"EXC"}, 11.50,
                                  lower_state_split_ev};
  const SourceSelector sel_ar_1s5{"Ar", {"EXC"}, lower_state_split_ev,
                                  11.70};
  const SourceSelector sel_ar_high_4s{"Ar", {"EXC"}, 11.70, 12.00};
  const SourceSelector sel_ar_upper{"Ar", {"EXC"}, 12.00, 100.0};

  const double population_1s4 =
      population_for_selector(sel_ar_1s4, populations);
  const double population_1s5 =
      population_for_selector(sel_ar_1s5, populations);
  const double population_high_4s =
      population_for_selector(sel_ar_high_4s, populations);
  const double population_upper =
      population_for_selector(sel_ar_upper, populations);
  if (population_1s4 <= 0.0 && population_1s5 <= 0.0 &&
      population_high_4s <= 0.0 && population_upper <= 0.0) {
    return;
  }

  // Ar**: productive cascade/self-relaxation competes with every additive.
  const double tau_upper_rad_ns =
      std::max(pget(p2nd, "tau_Ar_dbleStar_ns", 30.0), 1.0e-30);
  const double gamma_upper_rad = 1.0 / tau_upper_rad_ns;
  const double gamma_upper_productive =
      gamma_upper_rad + n * f_ar * rates.upper_ar;
  const double gamma_upper_total =
      gamma_upper_productive + n * rates.upper_additives;
  const double p_upper_productive =
      safe_ratio(gamma_upper_productive, gamma_upper_total);
  const double tau_upper_tot_ns =
      lifetime_from_total_rate_ns(gamma_upper_total);

  // Ar(1s2,1s3): relaxation into the two lower 4s states competes with
  // additive quenching.
  const double gamma_high_productive = n * f_ar * rates.four_s_ar;
  const double gamma_high_total =
      gamma_high_productive + n * rates.four_s_additives;
  const double p_high_productive =
      safe_ratio(gamma_high_productive, gamma_high_total);
  const double tau_high_tot_ns =
      lifetime_from_total_rate_ns(gamma_high_total);

  // Lower 4s states are now resolved. Both can form Ar2*, but only the
  // resonant 11.624 eV state has the competing trapped direct line.
  const double gamma_excimer_formation =
      n * n * f_ar * f_ar * rates.excimer_formation;
  const double gamma_4s_quenching = n * rates.four_s_additives;

  const double gamma_1s4_total =
      gamma_excimer_formation + gamma_4s_quenching;
  const double p_1s4_to_excimer =
      safe_ratio(gamma_excimer_formation, gamma_1s4_total);
  const double tau_1s4_tot_ns =
      lifetime_from_total_rate_ns(gamma_1s4_total);

  // r is the mean number of resonant absorption/re-emission cycles. It acts on
  // the radiative lifetime, tau_direct_eff = r * tau_atomic. Therefore changing
  // only r in Ar2nd_continium.csv updates both the direct branching ratio and
  // the waiting-time distribution used by hPhotonWavelengthTime.
  const double tau_atomic_1s5_ns = std::max(
      pget(p2nd, "tau_Ar_1s5_ns",
           pget(p2nd, "tau_Ar_1s4_ns",
                pget(p2nd, "tau_Ar4s_ns", 7.0))),
      1.0e-30);
  const double resonant_reemissions =
      std::max(positive_rate(pget(p2nd, "r", 150.0)), 1.0e-30);
  const double tau_direct_eff_ns =
      resonant_reemissions * tau_atomic_1s5_ns;
  const double gamma_direct_1s5 = 1.0 / tau_direct_eff_ns;
  const double gamma_1s5_total =
      gamma_excimer_formation + gamma_4s_quenching + gamma_direct_1s5;
  const double p_1s5_to_excimer =
      safe_ratio(gamma_excimer_formation, gamma_1s5_total);
  const double p_1s5_direct =
      safe_ratio(gamma_direct_1s5, gamma_1s5_total);
  const double tau_1s5_tot_ns =
      lifetime_from_total_rate_ns(gamma_1s5_total);

  // Final excimer states: radiative emission competes with additive quenching.
  const double tau_s_rad_ns =
      std::max(pget(p2nd, "tau_S_ns", 4.5), 1.0e-30);
  const double tau_t_rad_ns =
      std::max(pget(p2nd, "tau_T_ns", 3140.0), 1.0e-30);
  const double gamma_s_rad = 1.0 / tau_s_rad_ns;
  const double gamma_t_rad = 1.0 / tau_t_rad_ns;
  const double gamma_s_total =
      gamma_s_rad + n * rates.excimer_additives;
  const double gamma_t_total =
      gamma_t_rad + n * rates.excimer_additives;
  const double p_s_radiative = safe_ratio(gamma_s_rad, gamma_s_total);
  const double p_t_radiative = safe_ratio(gamma_t_rad, gamma_t_total);
  const double tau_s_tot_ns = lifetime_from_total_rate_ns(gamma_s_total);
  const double tau_t_tot_ns = lifetime_from_total_rate_ns(gamma_t_total);

  double f_s = 0.1;
  double f_t = 0.9;
  ar2nd_formation_fractions(p2nd, f_s, f_t);

  // The supplied Python model divides every successful high/upper cascade
  // equally between the independently represented 1s4 and 1s5 populations.
  // Optional CSV keys are accepted for later sensitivity studies without a
  // source-code change; absent keys reproduce the current 0.5/0.5 model.
  double f_cascade_to_1s4 = clamp_value(
      pget(p2nd, "f_Ar_cascade_to_1s4", 0.5), 0.0, 1.0);
  double f_cascade_to_1s5 = clamp_value(
      pget(p2nd, "f_Ar_cascade_to_1s5", 0.5), 0.0, 1.0);
  const double cascade_norm = f_cascade_to_1s4 + f_cascade_to_1s5;
  if (cascade_norm <= 0.0) {
    f_cascade_to_1s4 = 0.5;
    f_cascade_to_1s5 = 0.5;
  } else {
    f_cascade_to_1s4 /= cascade_norm;
    f_cascade_to_1s5 /= cascade_norm;
  }

  const double scale = std::max(0.0, pget(p2nd, "scale_Ar2nd", 1.0));
  const double triplet_weight =
      clamp_value(pget(p2nd, "triplet_weight", 1.0), 0.0, 1.0);
  const double w_upper_to_1s =
      std::max(0.0, pget(p2nd, "W_Ar_dbleStar_to_1s", 1.0));
  const double singlet_factor = scale * f_s * p_s_radiative;
  const double triplet_factor =
      scale * triplet_weight * f_t * p_t_radiative;

  const double lambda_2nd_nm = pget(p2nd, "lambda_Ar2nd_nm", 128.0);
  const double fwhm_2nd_nm = pget(p2nd, "fwhm_Ar2nd_nm", 10.0);
  const auto second_continuum_peaks =
      peaks_ar_second_continuum(lambda_2nd_nm, fwhm_2nd_nm);

  static constexpr double hc_ev_nm = 1239.8419843320026;
  const double lambda_direct_nm =
      hc_ev_nm / std::max(energy_1s5_ev, 1.0e-30);
  const double sigma_direct_nm = std::max(
      pget(p2nd, "sigma_Ar_1s5_direct_nm",
           pget(p2nd, "sigma_Ar_1s4_direct_nm", 0.25)),
      1.0e-12);
  const auto direct_line_peaks = peaks_single(lambda_direct_nm,
                                               sigma_direct_nm);

  struct FeedPath {
    std::string tag;
    SourceSelector source;
    double source_population = 0.0;
    double feed_probability = 0.0;
    std::vector<KineticStage> prior_stages;
    bool feeds_1s5 = false;
  };

  std::vector<FeedPath> paths;
  paths.push_back({"metastable_4s_11p548", sel_ar_1s4, population_1s4, 1.0, {}, false});
  paths.push_back({"resonant_4s_11p624", sel_ar_1s5, population_1s5, 1.0, {}, true});
  paths.push_back({
      "1s2_1s3_to_metastable_4s_11p548", sel_ar_high_4s, population_high_4s,
      p_high_productive * f_cascade_to_1s4,
      {excitation_stage("Ar(1s2,1s3) relaxation", tau_high_tot_ns)}, false});
  paths.push_back({
      "1s2_1s3_to_resonant_4s_11p624", sel_ar_high_4s, population_high_4s,
      p_high_productive * f_cascade_to_1s5,
      {excitation_stage("Ar(1s2,1s3) relaxation", tau_high_tot_ns)}, true});
  paths.push_back({
      "upper_to_metastable_4s_11p548", sel_ar_upper, population_upper,
      w_upper_to_1s * p_upper_productive * f_cascade_to_1s4,
      {excitation_stage("Ar** productive cascade", tau_upper_tot_ns)}, false});
  paths.push_back({
      "upper_to_resonant_4s_11p624", sel_ar_upper, population_upper,
      w_upper_to_1s * p_upper_productive * f_cascade_to_1s5,
      {excitation_stage("Ar** productive cascade", tau_upper_tot_ns)}, true});

  const std::string prefix = mixture_name(gas) + "_VUV_Ar2nd";
  for (const auto& path : paths) {
    if (path.source_population <= 0.0 || path.feed_probability <= 0.0) {
      continue;
    }

    const double p_to_excimer = path.feeds_1s5
                                    ? p_1s5_to_excimer
                                    : p_1s4_to_excimer;
    const double tau_lower_tot_ns = path.feeds_1s5
                                        ? tau_1s5_tot_ns
                                        : tau_1s4_tot_ns;
    const std::string lower_state_name = path.feeds_1s5
        ? "Ar resonant 4s state (11.624 eV)"
        : "Ar metastable 4s state (11.548 eV)";

    if (p_to_excimer > 0.0 && tau_lower_tot_ns > 0.0) {
      const double path_to_excimer = path.feed_probability * p_to_excimer;

      auto singlet_stages = path.prior_stages;
      singlet_stages.push_back(excitation_stage(
          lower_state_name + " -> Ar2* formation", tau_lower_tot_ns));
      singlet_stages.push_back(emission_stage(
          "Ar2*(1Sigma) radiative emission", tau_s_tot_ns));
      push_component(
          out, prefix + "_" + path.tag + "_to_singlet", "VUV",
          path.source, second_continuum_peaks, path.source_population,
          path.source_population * path_to_excimer * singlet_factor,
          path_to_excimer * singlet_factor, singlet_stages, mc_samples);

      auto triplet_stages = path.prior_stages;
      triplet_stages.push_back(excitation_stage(
          lower_state_name + " -> Ar2* formation", tau_lower_tot_ns));
      triplet_stages.push_back(emission_stage(
          "Ar2*(3Sigma) radiative emission", tau_t_tot_ns));
      push_component(
          out, prefix + "_" + path.tag + "_to_triplet", "VUV",
          path.source, second_continuum_peaks, path.source_population,
          path.source_population * path_to_excimer * triplet_factor,
          path_to_excimer * triplet_factor, triplet_stages, mc_samples);
    }

    if (path.feeds_1s5 && p_1s5_direct > 0.0 && tau_1s5_tot_ns > 0.0) {
      const double direct_path_probability =
          path.feed_probability * p_1s5_direct;
      auto direct_stages = path.prior_stages;
      // The same lower-state tau_tot is used for both possible successful
      // exits. Conditioned on direct emission, the event time is still drawn
      // from Exp(sum Gamma) rather than from Exp(r * tau_atomic) alone.
      direct_stages.push_back(emission_stage(
          "Ar resonant 11.624 eV trapped direct emission",
          tau_1s5_tot_ns));
      push_component(
          out, prefix + "_" + path.tag + "_direct_resonant", "VUV",
          path.source, direct_line_peaks, path.source_population,
          path.source_population * direct_path_probability * scale,
          direct_path_probability * scale, direct_stages, mc_samples);
    }
  }
}

inline std::vector<KineticComponent> build_arcf4_kinetic_components(
    const GasData& gas,
    const std::vector<LevelPopulation>& populations,
    const std::string& parameter_dir,
    const int mc_samples) {
  const auto p = read_parameter_csv_values(parameter_csv_for_model(parameter_dir, "ArCF4", false));
  const auto pir = read_parameter_csv_values(parameter_csv_for_model(parameter_dir, "ArCF4", true));

  const double f = clamp_value(fraction_of(gas, "cf4"), 0.0, 1.0);
  const double n = std::max(gas.pressure_bar, 0.0);
  constexpr double tau_3rd = 5.02;
  constexpr double tercer_continuo = 0.4866;

  const SourceSelector sel_cf3{"CF4", {"NEUTRAL DISS"}, 12.9, 100.0};
  const SourceSelector sel_ar_dbl{"Ar", {"EXC"}, 12.9, 100.0};
  const SourceSelector sel_cf4{"CF4", {"ION CF3 +"}, 15.0, 100.0};
  const SourceSelector sel_ar3rd{"Ar", {"CHARGE STATE ="}, 40.0, 100.0};

  const double P_CF3 = population_for_selector(sel_cf3, populations);
  const double P_Ar_dbl = population_for_selector(sel_ar_dbl, populations);
  const double P_CF4 = population_for_selector(sel_cf4, populations);
  const double P_Ar3rd = population_for_selector(sel_ar3rd, populations);

  const double N = pget(p, "Nnorm");
  const double p_CF3_vis = pget(p, "P_CF3_vis_dir");
  const double p_DbleStar = pget(p, "P_Ar_dbleStar");
  const double K = pget(p, "K_Ar_dbleStar_Q_Ar");
  const double K2 = pget(p, "K_Ar_dbleStar_Q_CF4");
  const double K1 = pget(p, "inv_tau_dis_K_relax");
  const double K3 = pget(p, "tau_uv_K_CF4_Q_CF4");
  const double p_CF4_dir = pget(p, "P_CF4_dir");
  const double K4 = pget(p, "K_Ar3rd_Q_CF4");
  const double PAr_3rd = pget(p, "P_Ar3rd");
  const double p_CF3_uv = pget(p, "P_CF3_uv_dir");

  // Radiative lifetimes are converted into effective total lifetimes by
  // adding every available collisional loss rate. Missing optional CF3/CF4
  // quenching coefficients default to zero, so the explicit 15 ns lifetime is
  // retained until a measured coefficient is supplied in the CSV.
  const double tau_cf3 = std::max(pget(p, "tau_CF3_ns", 15.0), 1.0e-30);
  const double cf3_total_rate =
      1.0 / tau_cf3 +
      n * (fraction_of(gas, "ar") * pget(p, "K_CF3_Q_Ar", 0.0) +
           fraction_of(gas, "cf4") * pget(p, "K_CF3_Q_CF4", 0.0) +
           fraction_of(gas, "he") * pget(p, "K_CF3_Q_He", 0.0) +
           fraction_of(gas, "n2") * pget(p, "K_CF3_Q_N2", 0.0) +
           fraction_of(gas, "ic4h10") * pget(p, "K_CF3_Q_IC4H10", 0.0));
  const double tau_cf3_tot = lifetime_from_total_rate_ns(cf3_total_rate);

  const double tau_cf4_uv =
      std::max(pget(p, "tau_CF4_UV_ns", 15.0), 1.0e-30);
  const double cf4_uv_total_rate =
      (1.0 / tau_cf4_uv) * (1.0 + pget(p, "tau_uv_K_CF4_Q_CF4") * n * f);
  const double tau_cf4_uv_tot =
      lifetime_from_total_rate_ns(cf4_uv_total_rate);

  const auto p2nd = read_optional_parameter_csv_values(ar2nd_parameter_csv(parameter_dir));
  const double br_cf4_d_to_x = clamp_value(pget(p2nd, "Br_CF4_D_to_X", pget(p, "Br_CF4_D_to_X", 0.1)), 0.0, 1.0);
  const double lambda_cf4_d_to_x_nm = pget(p2nd, "lambda_CF4_D_to_X_nm", 155.0);
  const double fwhm_cf4_d_to_x_nm = pget(p2nd, "fwhm_CF4_D_to_X_nm", 10.0);

  const double denom_dbl = n * f * K2 + n * (1.0 - f) * K + 1.0 / 30.0;
  const double frac_dbl_to_cf4 = safe_ratio(K2 * n * f, denom_dbl);
  const double tau_dbl_tot = lifetime_from_total_rate_ns(denom_dbl);

  const double frac1 = safe_ratio(f * n, f * n + K1);
  const double frac2 = safe_ratio(1.0, 1.0 + K3 * n * f);
  const double denom3 = (1.0 / tau_3rd) + f * n * K4;
  const double frac3 = safe_ratio(f * n * K4, denom3);
  const double frac4 = safe_ratio(1.0 / tau_3rd, denom3);
  const double tau_3rd_tot = lifetime_from_total_rate_ns(denom3);

  std::vector<KineticComponent> out;

  const double vis_cf3 = N * p_CF3_vis * P_CF3;
  const double vis_ar = N * frac_dbl_to_cf4 * p_DbleStar * P_Ar_dbl;
  push_component(out, "ArCF4_VIS_CF3_direct", "VIS", sel_cf3,
                 peaks_cf3_visible(), P_CF3, vis_cf3, p_CF3_vis,
                 std::vector<KineticStage>{
                     emission_stage("CF3* visible emission", tau_cf3_tot)},
                 mc_samples);
  push_component(out, "ArCF4_VIS_ArDbleStar_transfer", "VIS", sel_ar_dbl,
                 peaks_cf3_visible(), P_Ar_dbl, vis_ar,
                 frac_dbl_to_cf4 * p_DbleStar,
                 std::vector<KineticStage>{
                     excitation_stage("Ar** -> CF3* transfer", tau_dbl_tot),
                     emission_stage("CF3* visible emission", tau_cf3_tot)},
                 mc_samples);

  // CF4+(D)->CF4+(X) branch from the current scintillation scheme:
  // Y_150 = Br_D_to_X * Y_CF4,UV.  It is an added VUV branch relative to the
  // fitted CF4 ionic UV channel; do not reduce the ordinary CF4 UV component.
  const double uv_cf3 = p_CF3_uv * N * p_CF3_vis * P_CF3;
  const double uv_ar_dbl = p_CF3_uv * N * frac_dbl_to_cf4 * p_DbleStar * P_Ar_dbl;
  const double uv_cf4_total = N * (frac1 * frac2) * (p_CF4_dir * P_CF4);
  const double uv_cf4_150 = br_cf4_d_to_x * uv_cf4_total;
  const double uv_ar3_transfer = N * (frac1 * frac2) * (frac3 * P_Ar3rd * PAr_3rd);
  const double uv_ar3_direct = N * (tercer_continuo * frac4 * P_Ar3rd);

  push_component(out, "ArCF4_UV_CF3_direct", "UV", sel_cf3,
                 peaks_cf3_uv(), P_CF3, uv_cf3,
                 p_CF3_uv * p_CF3_vis,
                 std::vector<KineticStage>{
                     emission_stage("CF3* UV emission", tau_cf3_tot)},
                 mc_samples);
  push_component(out, "ArCF4_UV_ArDbleStar_transfer", "UV", sel_ar_dbl,
                 peaks_cf3_uv(), P_Ar_dbl, uv_ar_dbl,
                 p_CF3_uv * frac_dbl_to_cf4 * p_DbleStar,
                 std::vector<KineticStage>{
                     excitation_stage("Ar** -> CF3* transfer", tau_dbl_tot),
                     emission_stage("CF3* UV emission", tau_cf3_tot)},
                 mc_samples);
  push_component(out, "ArCF4_UV_CF4_direct", "UV", sel_cf4,
                 peaks_cf4_uv(), P_CF4, uv_cf4_total,
                 frac1 * frac2 * p_CF4_dir,
                 std::vector<KineticStage>{
                     emission_stage("CF4+* UV emission", tau_cf4_uv_tot)},
                 mc_samples);
  push_component(out, "ArCF4_VUV155_CF4_direct", "VUV", sel_cf4,
                 peaks_cf4_vuv_branch(lambda_cf4_d_to_x_nm,
                                      fwhm_cf4_d_to_x_nm),
                 P_CF4, uv_cf4_150,
                 br_cf4_d_to_x * frac1 * frac2 * p_CF4_dir,
                 std::vector<KineticStage>{
                     emission_stage("CF4+*(D->X) VUV emission",
                                    tau_cf4_uv_tot)},
                 mc_samples);
  push_component(out, "ArCF4_UV_Ar3rd_transfer", "UV", sel_ar3rd,
                 peaks_ar_3rd_uv(), P_Ar3rd, uv_ar3_transfer,
                 frac1 * frac2 * frac3 * PAr_3rd,
                 std::vector<KineticStage>{
                     excitation_stage("Ar third-continuum -> CF4+* transfer",
                                      tau_3rd_tot),
                     emission_stage("CF4+* UV emission", tau_cf4_uv_tot)},
                 mc_samples);
  push_component(out, "ArCF4_UV_Ar3rd_third_continuum", "UV", sel_ar3rd,
                 peaks_ar_3rd_uv(), P_Ar3rd, uv_ar3_direct,
                 tercer_continuo * frac4,
                 std::vector<KineticStage>{
                     emission_stage("Ar third-continuum emission",
                                    tau_3rd_tot)},
                 mc_samples);

  push_ar2nd_components(out, gas, populations, p2nd, mc_samples);

  struct IrLine { const char* tag; double lambda; double e0; const char* tau_name; const char* k_cf4_name; };
  const IrLine lines[] = {
      {"696", 696.0, 13.32, "tau_CF4_696", "K_Ar_Q_CF4_696"},
      {"727", 727.0, 13.32, "tau_CF4_727", "K_Ar_Q_CF4_727"},
      {"750", 750.0, 13.47, "tau_CF4_750", "K_Ar_Q_CF4_750"},
      {"763", 763.0, 13.17, "tau_CF4_764", "K_Ar_Q_CF4_764"},
      {"772", 772.0, 13.32, "tau_CF4_772", "K_Ar_Q_CF4_772"},
      {"794", 794.0, 13.28, "tau_CF4_794", "K_Ar_Q_CF4_794"},
  };

  for (const auto& line : lines) {
    const SourceSelector sel{"Ar", {"EXC"}, line.e0, line.e0 + 10.0};
    const double P = population_for_selector(sel, populations);
    const std::string s = line.tag;
    const double pstar = pget(pir, "PAr_star_" + s);
    const double tau = pget(pir, line.tau_name, 0.0);
    const double k_ar = pget(pir, "K_Ar_Q_Ar_" + s);
    const double k_cf4 = pget(pir, line.k_cf4_name);
    const double radiative = tau > 0.0 ? 1.0 / tau : 0.0;
    const double denom = radiative + n * f * k_cf4 + n * (1.0 - f) * k_ar;
    const double prob = pstar * safe_ratio(radiative, denom);
    const double yield = prob * P;
    push_component(
        out, "ArCF4_IR_" + s, "IR", sel,
        peaks_single(line.lambda, 2.5), P, yield, prob,
        std::vector<KineticStage>{emission_stage(
            "Ar IR " + s + " nm emission", lifetime_from_total_rate_ns(denom))},
        mc_samples);
  }

  return out;
}

inline std::vector<KineticComponent> build_arn2_kinetic_components(
    const GasData& gas,
    const std::vector<LevelPopulation>& populations,
    const std::string& parameter_dir,
    const int mc_samples) {
  const auto p = read_parameter_csv_values(parameter_csv_for_model(parameter_dir, "ArN2", false));
  const auto pir = read_parameter_csv_values(parameter_csv_for_model(parameter_dir, "ArN2", true));

  const double f = clamp_value(fraction_of(gas, "n2"), 0.0, 1.0);
  const double n = std::max(gas.pressure_bar, 0.0);

  const auto p2nd = read_optional_parameter_csv_values(ar2nd_parameter_csv(parameter_dir));

  const SourceSelector sel_n2{"N2", {"C 3PI"}, 11.0, 15.5};
  const SourceSelector sel_ar_meta{"Ar", {"EXC"}, 0.0, 11.6};
  const SourceSelector sel_ar_res{"Ar", {"EXC"}, 11.6, 11.7};
  const SourceSelector sel_ar_dbl{"Ar", {"EXC"}, 11.7, 100.0};

  const double P_N2_pop = population_for_selector(sel_n2, populations);
  const double P_Ar_meta = population_for_selector(sel_ar_meta, populations);
  const double P_Ar_res = population_for_selector(sel_ar_res, populations);
  const double P_Ar_dbl = population_for_selector(sel_ar_dbl, populations);

  const double Nnorm = pget(p, "Nnorm");
  const double P_N2 = pget(p, "P_N2");
  const double tau_N2 = pget(p, "tau_N2");
  const double K_N2_Q_N2 = pget(p, "K_N2_Q_N2");
  const double K_N2_Q_Ar = pget(p, "K_N2_Q_Ar");
  const double K_ArMeta_Q_N2c = pget(p, "K_ArMeta_Q_N2c");
  const double K_ArMeta_Q_N2b = pget(p, "K_ArMeta_Q_N2b");
  const double K_ArMeta_Q_2Ar = pget(p, "K_ArMeta_Q_2Ar");
  const double K_ArRes_Q_N2c = pget(p, "K_ArRes_Q_N2c");
  const double K_ArRes_Q_N2b = pget(p, "K_ArRes_Q_N2b");
  const double K_ArRes_Q_2Ar = pget(p, "K_ArRes_Q_2Ar");
  const double P_Ar_dbleStar = pget(p, "P_Ar_dbleStar");
  const double frac_Ar_dbleStar = pget(p, "frac_Ar_dbleStar");

  const double radiative_n2 = tau_N2 > 0.0 ? 1.0 / tau_N2 : 0.0;
  const double n2_total_rate =
      radiative_n2 + n * f * K_N2_Q_N2 +
      n * (1.0 - f) * K_N2_Q_Ar;
  const double factor_N2 = safe_ratio(radiative_n2, n2_total_rate);
  const double tau_n2_tot = lifetime_from_total_rate_ns(n2_total_rate);

  const double meta_total_rate =
      (K_ArMeta_Q_N2b + K_ArMeta_Q_N2c) * f * n +
      K_ArMeta_Q_2Ar * std::pow(1.0 - f, 2) * n * n;
  const double res_total_rate =
      (K_ArRes_Q_N2b + K_ArRes_Q_N2c) * f * n +
      K_ArRes_Q_2Ar * std::pow(1.0 - f, 2) * n * n;
  const double factor_Ar_meta =
      safe_ratio(K_ArMeta_Q_N2c * f * n, meta_total_rate);
  const double factor_Ar_res =
      safe_ratio(K_ArRes_Q_N2c * f * n, res_total_rate);
  const double tau_meta_tot = lifetime_from_total_rate_ns(meta_total_rate);
  const double tau_res_tot = lifetime_from_total_rate_ns(res_total_rate);

  const Ar2ndKineticRates ar2nd_rates = ar2nd_rates_for_mixture(p2nd, gas);
  const double tau_upper =
      std::max(pget(p2nd, "tau_Ar_dbleStar_ns", 30.0), 1.0e-30);
  const double upper_productive_rate =
      1.0 / tau_upper + n * ar2nd_rates.argon_fraction * ar2nd_rates.upper_ar;
  const double upper_total_rate =
      upper_productive_rate + n * ar2nd_rates.upper_additives;
  const double tau_upper_tot = lifetime_from_total_rate_ns(upper_total_rate);

  std::vector<KineticComponent> out;
  const auto n2_peaks = peaks_n2_second_positive();

  push_component(out, "ArN2_UV_N2_C_direct", "UV", sel_n2, n2_peaks,
                 P_N2_pop, Nnorm * factor_N2 * (P_N2_pop * P_N2),
                 factor_N2 * P_N2,
                 std::vector<KineticStage>{
                     emission_stage("N2(C) second-positive emission",
                                    tau_n2_tot)},
                 mc_samples);
  push_component(out, "ArN2_UV_ArMeta_to_N2C", "UV", sel_ar_meta,
                 n2_peaks, P_Ar_meta,
                 Nnorm * factor_N2 * (P_Ar_meta * factor_Ar_meta),
                 factor_N2 * factor_Ar_meta,
                 std::vector<KineticStage>{
                     excitation_stage("Ar metastable -> N2(C) transfer",
                                      tau_meta_tot),
                     emission_stage("N2(C) second-positive emission",
                                    tau_n2_tot)},
                 mc_samples);
  push_component(out, "ArN2_UV_ArRes_to_N2C", "UV", sel_ar_res,
                 n2_peaks, P_Ar_res,
                 Nnorm * factor_N2 * (P_Ar_res * factor_Ar_res),
                 factor_N2 * factor_Ar_res,
                 std::vector<KineticStage>{
                     excitation_stage("Ar resonant -> N2(C) transfer",
                                      tau_res_tot),
                     emission_stage("N2(C) second-positive emission",
                                    tau_n2_tot)},
                 mc_samples);

  // Keep the metastable and resonant Ar** branches separate: a weighted
  // average lifetime would not reproduce the true mixture of exponentials.
  const double ar_dbl_meta_probability =
      factor_N2 * P_Ar_dbleStar * frac_Ar_dbleStar * factor_Ar_meta;
  const double ar_dbl_res_probability =
      factor_N2 * P_Ar_dbleStar * (1.0 - frac_Ar_dbleStar) * factor_Ar_res;
  push_component(out, "ArN2_UV_ArDbleStar_to_N2C_via_meta", "UV",
                 sel_ar_dbl, n2_peaks, P_Ar_dbl,
                 Nnorm * P_Ar_dbl * ar_dbl_meta_probability,
                 ar_dbl_meta_probability,
                 std::vector<KineticStage>{
                     excitation_stage("Ar** upper cascade", tau_upper_tot),
                     excitation_stage("Ar metastable -> N2(C) transfer",
                                      tau_meta_tot),
                     emission_stage("N2(C) second-positive emission",
                                    tau_n2_tot)},
                 mc_samples);
  push_component(out, "ArN2_UV_ArDbleStar_to_N2C_via_res", "UV",
                 sel_ar_dbl, n2_peaks, P_Ar_dbl,
                 Nnorm * P_Ar_dbl * ar_dbl_res_probability,
                 ar_dbl_res_probability,
                 std::vector<KineticStage>{
                     excitation_stage("Ar** upper cascade", tau_upper_tot),
                     excitation_stage("Ar resonant -> N2(C) transfer",
                                      tau_res_tot),
                     emission_stage("N2(C) second-positive emission",
                                    tau_n2_tot)},
                 mc_samples);

  push_ar2nd_components(out, gas, populations, p2nd, mc_samples);

  struct IrLine { const char* tag; double lambda; double e0; };
  const IrLine lines[] = {
      {"696", 696.0, 13.32}, {"727", 727.0, 13.32}, {"750", 750.0, 13.47},
      {"763", 763.0, 13.17}, {"772", 772.0, 13.32}, {"794", 794.0, 13.28},
  };
  for (const auto& line : lines) {
    const SourceSelector sel{"Ar", {"EXC"}, line.e0, line.e0 + 10.0};
    const double P = population_for_selector(sel, populations);
    const std::string s = line.tag;
    const double pstar = pget(pir, "PAr_star_" + s);
    const double tau = pget(pir, "tau_N2_" + s, 0.0);
    const double k_ar = pget(pir, "K_Ar_Q_Ar_" + s);
    const double k_n2 = pget(pir, "K_Ar_Q_N2_" + s);
    const double radiative = tau > 0.0 ? 1.0 / tau : 0.0;
    const double denom = radiative + n * f * k_n2 + n * (1.0 - f) * k_ar;
    const double prob = pstar * safe_ratio(radiative, denom);
    const double yield = prob * P;
    push_component(
        out, "ArN2_IR_" + s, "IR", sel,
        peaks_single(line.lambda, 2.8), P, yield, prob,
        std::vector<KineticStage>{emission_stage(
            "Ar IR " + s + " nm emission", lifetime_from_total_rate_ns(denom))},
        mc_samples);
  }

  return out;
}

inline std::vector<KineticComponent> build_hecf4_kinetic_components(
    const GasData& gas,
    const std::vector<LevelPopulation>& populations,
    const std::string& parameter_dir,
    const int mc_samples) {
  // No dedicated HeCF4 fit was supplied here. Use the CF4-driven part of the
  // ArCF4 model and ignore Ar-specific channels. This keeps HeCF4 usable while
  // making the approximation explicit in kineticData.
  const auto p = read_parameter_csv_values(parameter_csv_for_model(parameter_dir, "ArCF4", false));
  const double f = clamp_value(fraction_of(gas, "cf4"), 0.0, 1.0);
  const double n = std::max(gas.pressure_bar, 0.0);

  const SourceSelector sel_cf3{"CF4", {"NEUTRAL DISS"}, 12.9, 100.0};
  const SourceSelector sel_cf4{"CF4", {"ION CF3 +"}, 15.0, 100.0};
  const double P_CF3 = population_for_selector(sel_cf3, populations);
  const double P_CF4 = population_for_selector(sel_cf4, populations);

  const double N = pget(p, "Nnorm");
  const double p_CF3_vis = pget(p, "P_CF3_vis_dir");
  const double K1 = pget(p, "inv_tau_dis_K_relax");
  const double K3 = pget(p, "tau_uv_K_CF4_Q_CF4");
  const double p_CF4_dir = pget(p, "P_CF4_dir");
  const double p_CF3_uv = pget(p, "P_CF3_uv_dir");

  const double tau_cf3 = std::max(pget(p, "tau_CF3_ns", 15.0), 1.0e-30);
  const double cf3_total_rate =
      1.0 / tau_cf3 +
      n * (fraction_of(gas, "he") * pget(p, "K_CF3_Q_He", 0.0) +
           fraction_of(gas, "cf4") * pget(p, "K_CF3_Q_CF4", 0.0));
  const double tau_cf3_tot = lifetime_from_total_rate_ns(cf3_total_rate);
  const double tau_cf4_uv =
      std::max(pget(p, "tau_CF4_UV_ns", 15.0), 1.0e-30);
  const double cf4_uv_total_rate =
      (1.0 / tau_cf4_uv) * (1.0 + K3 * n * f);
  const double tau_cf4_uv_tot =
      lifetime_from_total_rate_ns(cf4_uv_total_rate);

  const double frac1 = safe_ratio(f * n, f * n + K1);
  const double frac2 = safe_ratio(1.0, 1.0 + K3 * n * f);
  const auto p2nd = read_optional_parameter_csv_values(ar2nd_parameter_csv(parameter_dir));
  const double br_cf4_d_to_x = clamp_value(pget(p2nd, "Br_CF4_D_to_X", pget(p, "Br_CF4_D_to_X", 0.1)), 0.0, 1.0);
  const double lambda_cf4_d_to_x_nm = pget(p2nd, "lambda_CF4_D_to_X_nm", 155.0);
  const double fwhm_cf4_d_to_x_nm = pget(p2nd, "fwhm_CF4_D_to_X_nm", 10.0);

  std::vector<KineticComponent> out;
  push_component(
      out, "HeCF4_VIS_CF3_direct_ArCF4Fit", "VIS", sel_cf3,
      peaks_cf3_visible(), P_CF3, N * p_CF3_vis * P_CF3, p_CF3_vis,
      std::vector<KineticStage>{
          emission_stage("CF3* visible emission", tau_cf3_tot)},
      mc_samples);
  push_component(
      out, "HeCF4_UV_CF3_direct_ArCF4Fit", "UV", sel_cf3,
      peaks_cf3_uv(), P_CF3, p_CF3_uv * N * p_CF3_vis * P_CF3,
      p_CF3_uv * p_CF3_vis,
      std::vector<KineticStage>{
          emission_stage("CF3* UV emission", tau_cf3_tot)},
      mc_samples);
  const double uv_cf4_total = N * frac1 * frac2 * p_CF4_dir * P_CF4;
  push_component(
      out, "HeCF4_UV_CF4_direct_ArCF4Fit", "UV", sel_cf4,
      peaks_cf4_uv(), P_CF4, uv_cf4_total, frac1 * frac2 * p_CF4_dir,
      std::vector<KineticStage>{
          emission_stage("CF4+* UV emission", tau_cf4_uv_tot)},
      mc_samples);
  push_component(
      out, "HeCF4_VUV155_CF4_direct_ArCF4Fit", "VUV", sel_cf4,
      peaks_cf4_vuv_branch(lambda_cf4_d_to_x_nm,
                           fwhm_cf4_d_to_x_nm),
      P_CF4, br_cf4_d_to_x * uv_cf4_total,
      br_cf4_d_to_x * frac1 * frac2 * p_CF4_dir,
      std::vector<KineticStage>{
          emission_stage("CF4+*(D->X) VUV emission", tau_cf4_uv_tot)},
      mc_samples);
  return out;
}

inline std::vector<KineticComponent> build_kinetic_components(
    const GasData& gas,
    const std::vector<LevelPopulation>& populations,
    const std::string& parameter_dir,
    const int mc_samples) {
  const std::string mix = mixture_name(gas);
  std::vector<KineticComponent> out;
  if (mix == "ArN2") {
    out = build_arn2_kinetic_components(gas, populations, parameter_dir, mc_samples);
  } else if (mix == "HeCF4") {
    out = build_hecf4_kinetic_components(gas, populations, parameter_dir, mc_samples);
  } else {
    out = build_arcf4_kinetic_components(gas, populations, parameter_dir, mc_samples);
  }
  apply_secondary_normalization(out, parameter_dir);
  return out;
}

inline void rescale_components_to_photons_per_electron(std::vector<KineticComponent>& components,
                                                       const GasData& gas,
                                                       const int /*mc_samples*/,
                                                       const double photons_per_electron) {
  if (!std::isfinite(photons_per_electron) || photons_per_electron < 0.0) return;
  double current = 0.0;
  for (const auto& c : components) current += c.expected_photons;
  const double target = static_cast<double>(std::max<long long>(0, gas.ne_total)) * photons_per_electron;
  if (current <= 0.0 || target < 0.0) return;
  const double scale = target / current;
  for (auto& c : components) c.expected_photons *= scale;
}

std::vector<KineticComponent> create_photon_emission(
    const GasData& gas, const std::vector<LevelPopulation>& populations,
    const std::string& parameters_dir, int mc_samples);

double slowest_emission_lifetime_ns(
    const std::vector<KineticComponent>& components);

// One resolved non-elastic collision recorded during a single avalanche.
// The photon model uses these sites directly, preserving the x-y-z-t
// correlation instead of resampling a global histogram after the simulation.
struct PhotonSourceSite {
  int level = -1;
  double x_cm = 0.0;
  double y_cm = 0.0;
  double z_cm = 0.0;
  double time_ns = 0.0;
};

std::vector<LevelPopulation> populations_from_sites(
    const std::vector<LevelInfo>& levels,
    const std::vector<PhotonSourceSite>& sites);

std::vector<const PhotonSourceSite*> matching_source_sites(
    const KineticComponent& component, const std::vector<LevelInfo>& levels,
    const std::vector<PhotonSourceSite>& sites);

const PhotonSourceSite* sample_source_site(
    const std::vector<const PhotonSourceSite*>& matching_sites,
    TRandom3& random);

// Compatibility overload. Prefer building matching_source_sites once per
// kinetic component and sampling from that pool for every emitted photon.
const PhotonSourceSite* sample_source_site(
    const KineticComponent& component, const std::vector<LevelInfo>& levels,
    const std::vector<PhotonSourceSite>& sites, TRandom3& random);

double sample_wavelength_nm(const KineticComponent& component,
                            TRandom3& random, double minimum_nm,
                            double maximum_nm);

EmissionDelaySample sample_emission_delay(
    const KineticComponent& component, TRandom3& random);

double sample_emission_delay_ns(const KineticComponent& component,
                                TRandom3& random);

long long number_of_mc_photons(double expected_photons,
                               long long mc_samples, TRandom3& random);

constexpr double kHcEvNm = 1239.8419843320026;
inline double photon_energy_ev(const double wavelength_nm) {
  return std::isfinite(wavelength_nm) && wavelength_nm > 0.0
             ? kHcEvNm / wavelength_nm
             : 0.0;
}

}  // namespace photonemission

#endif  // SECONDARY_AVALANCHES_PHOTON_MODEL_HH
