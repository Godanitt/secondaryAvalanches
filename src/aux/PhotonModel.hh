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

struct KineticComponent {
  std::string name;
  std::string region;
  SourceSelector source;
  std::vector<SpectralPeak> peaks;
  double source_population = 0.0;
  double expected_photons = 0.0;
  double probability = 0.0;
  double mean_lifetime_ns = 0.0;
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

  const std::size_t i_name = optional_col("name") != std::numeric_limits<std::size_t>::max()
                               ? optional_col("name")
                               : optional_col("parameter");
  const std::size_t i_value = optional_col("value") != std::numeric_limits<std::size_t>::max()
                                ? optional_col("value")
                                : optional_col("secondary_optimal");
  if (i_name == std::numeric_limits<std::size_t>::max()) {
    throw std::runtime_error("Columna ausente en parameter CSV: name/parameter");
  }
  if (i_value == std::numeric_limits<std::size_t>::max()) {
    throw std::runtime_error("Columna ausente en parameter CSV: value/secondary_optimal");
  }

  std::map<std::string, double> out;
  std::string line;
  while (std::getline(in, line)) {
    if (trim(line).empty()) continue;
    const auto fields = split_csv_line(line);
    if (fields.size() <= std::max(i_name, i_value)) continue;
    const std::string name = trim(fields[i_name]);
    const std::string value_text = trim(fields[i_value]);
    if (name.empty() || value_text.empty()) continue;
    try {
      out[name] = std::stod(value_text);
    } catch (...) {
      // Ignore NaN/empty/non-numeric fields; only the central value column matters here.
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

inline void push_component(std::vector<KineticComponent>& out,
                           std::string name,
                           std::string region,
                           SourceSelector source,
                           std::vector<SpectralPeak> peaks,
                           const double source_population,
                           const double expected_photons_one_sample,
                           const double probability,
                           const double lifetime_ns,
                           const int /*mc_samples*/) {
  if (!std::isfinite(expected_photons_one_sample) || expected_photons_one_sample <= 0.0) return;
  KineticComponent c;
  c.name = std::move(name);
  c.region = std::move(region);
  c.source = std::move(source);
  c.peaks = std::move(peaks);
  c.source_population = source_population;
  c.expected_photons = expected_photons_one_sample;
  c.probability = probability;
  c.mean_lifetime_ns = std::max(0.0, lifetime_ns);
  out.push_back(std::move(c));
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

  // Disjoint source populations for the second-continuum cascade.  The
  // intervals must never overlap, otherwise the 11.7--12.0 eV 4s population
  // would also be counted as Ar**.
  //   Ar(1s4,1s5): 11.5--11.7 eV
  //   Ar(1s2,1s3): 11.7--12.0 eV
  //   Ar**:         12.0--100 eV
  const SourceSelector sel_ar_precursor{"Ar", {"EXC"}, 11.50, 11.70};
  const SourceSelector sel_ar_high_4s{"Ar", {"EXC"}, 11.70, 12.00};
  const SourceSelector sel_ar_upper{"Ar", {"EXC"}, 12.00, 100.0};

  const double p_precursor_population = population_for_selector(sel_ar_precursor, populations);
  const double p_high_4s_population = population_for_selector(sel_ar_high_4s, populations);
  const double p_upper_population = population_for_selector(sel_ar_upper, populations);
  if (p_precursor_population <= 0.0 && p_high_4s_population <= 0.0 &&
      p_upper_population <= 0.0) {
    return;
  }

  const double tau_upper = std::max(
      pget(p2nd, "tau_Ar_dbleStar_ns", 30.0), 1.0e-30);
  const double gamma_upper = 1.0 / tau_upper;
  const double upper_productive = gamma_upper + n * f_ar * rates.upper_ar;
  const double p_upper = safe_ratio(
      upper_productive,
      upper_productive + n * rates.upper_additives);

  const double high_4s_productive = n * f_ar * rates.four_s_ar;
  const double p_high_4s = safe_ratio(
      high_4s_productive,
      high_4s_productive + n * rates.four_s_additives);

  const double formation_productive =
      n * n * f_ar * f_ar * rates.excimer_formation;
  const double p_form = safe_ratio(
      formation_productive,
      formation_productive + n * rates.four_s_additives);
  if (p_form <= 0.0) return;

  const double tau_s = std::max(pget(p2nd, "tau_S_ns", 11.3), 1.0e-30);
  const double tau_t = std::max(pget(p2nd, "tau_T_ns", 3140.0), 1.0e-30);
  const double inv_tau_s = 1.0 / tau_s;
  const double inv_tau_t = 1.0 / tau_t;
  const double p_rad_s = safe_ratio(
      inv_tau_s,
      inv_tau_s + n * rates.excimer_additives);
  const double p_rad_t = safe_ratio(
      inv_tau_t,
      inv_tau_t + n * rates.excimer_additives);

  double f_s = 0.1;
  double f_t = 0.9;
  ar2nd_formation_fractions(p2nd, f_s, f_t);

  const double scale = pget(p2nd, "scale_Ar2nd", 1.0);
  const double w_upper_to_1s = std::max(
      0.0, pget(p2nd, "W_Ar_dbleStar_to_1s", 1.0));
  const double lambda = pget(p2nd, "lambda_Ar2nd_nm", 128.0);
  const double fwhm = pget(p2nd, "fwhm_Ar2nd_nm", 10.0);
  const auto peaks = peaks_ar_second_continuum(lambda, fwhm);

  const double path_precursor = p_form;
  const double path_high_4s = p_high_4s * p_form;
  const double path_upper = w_upper_to_1s * p_upper * p_form;

  const double singlet_factor = scale * f_s * p_rad_s;
  // Both excimer branches are always generated.  The slow component is not
  // disabled by hand: it disappears naturally in mixtures through p_rad_t,
  // whose long triplet lifetime makes additive quenching dominant.
  const double triplet_weight = clamp_value(pget(p2nd, "triplet_weight", 1.0), 0.0, 1.0);
  const double triplet_factor = scale * triplet_weight * f_t * p_rad_t;
  const std::string prefix = mixture_name(gas) + "_VUV_Ar2nd";

  push_component(out, prefix + "_1s4_1s5_singlet", "VUV", sel_ar_precursor, peaks,
                 p_precursor_population,
                 p_precursor_population * path_precursor * singlet_factor,
                 path_precursor * singlet_factor, tau_s, mc_samples);
  push_component(out, prefix + "_1s2_1s3_singlet", "VUV", sel_ar_high_4s, peaks,
                 p_high_4s_population,
                 p_high_4s_population * path_high_4s * singlet_factor,
                 path_high_4s * singlet_factor, tau_s, mc_samples);
  push_component(out, prefix + "_upper_singlet", "VUV", sel_ar_upper, peaks,
                 p_upper_population,
                 p_upper_population * path_upper * singlet_factor,
                 path_upper * singlet_factor, tau_s, mc_samples);

  push_component(out, prefix + "_1s4_1s5_triplet", "VUV", sel_ar_precursor, peaks,
                 p_precursor_population,
                 p_precursor_population * path_precursor * triplet_factor,
                 path_precursor * triplet_factor, tau_t, mc_samples);
  push_component(out, prefix + "_1s2_1s3_triplet", "VUV", sel_ar_high_4s, peaks,
                 p_high_4s_population,
                 p_high_4s_population * path_high_4s * triplet_factor,
                 path_high_4s * triplet_factor, tau_t, mc_samples);
  push_component(out, prefix + "_upper_triplet", "VUV", sel_ar_upper, peaks,
                 p_upper_population,
                 p_upper_population * path_upper * triplet_factor,
                 path_upper * triplet_factor, tau_t, mc_samples);
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

  const auto p2nd = read_optional_parameter_csv_values(ar2nd_parameter_csv(parameter_dir));
  const double br_cf4_d_to_x = clamp_value(pget(p2nd, "Br_CF4_D_to_X", pget(p, "Br_CF4_D_to_X", 0.1)), 0.0, 1.0);
  const double lambda_cf4_d_to_x_nm = pget(p2nd, "lambda_CF4_D_to_X_nm", 155.0);
  const double fwhm_cf4_d_to_x_nm = pget(p2nd, "fwhm_CF4_D_to_X_nm", 10.0);

  const double denom_dbl = n * f * K2 + n * (1.0 - f) * K + 1.0 / 30.0;
  const double frac_dbl_to_cf4 = safe_ratio(K2 * n * f, denom_dbl);

  const double frac1 = safe_ratio(f * n, f * n + K1);
  const double frac2 = safe_ratio(1.0, 1.0 + K3 * n * f);
  const double denom3 = (1.0 / tau_3rd) + f * n * K4;
  const double frac3 = safe_ratio(f * n * K4, denom3);
  const double frac4 = safe_ratio(1.0 / tau_3rd, denom3);

  std::vector<KineticComponent> out;

  const double vis_cf3 = N * p_CF3_vis * P_CF3;
  const double vis_ar = N * frac_dbl_to_cf4 * p_DbleStar * P_Ar_dbl;
  push_component(out, "ArCF4_VIS_CF3_direct", "VIS", sel_cf3, peaks_cf3_visible(), P_CF3,
                 vis_cf3, p_CF3_vis, 0.0, mc_samples);
  push_component(out, "ArCF4_VIS_ArDbleStar_transfer", "VIS", sel_ar_dbl, peaks_cf3_visible(), P_Ar_dbl,
                 vis_ar, frac_dbl_to_cf4 * p_DbleStar, 30.0, mc_samples);

  // CF4+(D)->CF4+(X) branch from the current scintillation scheme:
  // Y_150 = Br_D_to_X * Y_CF4,UV.  It is an added VUV branch relative to the
  // fitted CF4 ionic UV channel; do not reduce the ordinary CF4 UV component.
  const double uv_cf3 = p_CF3_uv * N * p_CF3_vis * P_CF3;
  const double uv_ar_dbl = p_CF3_uv * N * frac_dbl_to_cf4 * p_DbleStar * P_Ar_dbl;
  const double uv_cf4_total = N * (frac1 * frac2) * (p_CF4_dir * P_CF4);
  const double uv_cf4_150 = br_cf4_d_to_x * uv_cf4_total;
  const double uv_ar3_transfer = N * (frac1 * frac2) * (frac3 * P_Ar3rd * PAr_3rd);
  const double uv_ar3_direct = N * (tercer_continuo * frac4 * P_Ar3rd);

  push_component(out, "ArCF4_UV_CF3_direct", "UV", sel_cf3, peaks_cf3_uv(), P_CF3,
                 uv_cf3, p_CF3_uv * p_CF3_vis, 0.0, mc_samples);
  push_component(out, "ArCF4_UV_ArDbleStar_transfer", "UV", sel_ar_dbl, peaks_cf3_uv(), P_Ar_dbl,
                 uv_ar_dbl, p_CF3_uv * frac_dbl_to_cf4 * p_DbleStar, 30.0, mc_samples);
  push_component(out, "ArCF4_UV_CF4_direct", "UV", sel_cf4, peaks_cf4_uv(), P_CF4,
                 uv_cf4_total, frac1 * frac2 * p_CF4_dir, 0.0, mc_samples);
  push_component(out, "ArCF4_VUV155_CF4_direct", "VUV", sel_cf4,
                 peaks_cf4_vuv_branch(lambda_cf4_d_to_x_nm, fwhm_cf4_d_to_x_nm), P_CF4,
                 uv_cf4_150, br_cf4_d_to_x * frac1 * frac2 * p_CF4_dir, 0.0, mc_samples);
  push_component(out, "ArCF4_UV_Ar3rd_transfer", "UV", sel_ar3rd, peaks_ar_3rd_uv(), P_Ar3rd,
                 uv_ar3_transfer, frac1 * frac2 * frac3 * PAr_3rd, tau_3rd, mc_samples);
  push_component(out, "ArCF4_UV_Ar3rd_third_continuum", "UV", sel_ar3rd, peaks_ar_3rd_uv(), P_Ar3rd,
                 uv_ar3_direct, tercer_continuo * frac4, tau_3rd, mc_samples);

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
    push_component(out, "ArCF4_IR_" + s, "IR", sel, peaks_single(line.lambda, 2.5), P,
                   yield, prob, tau, mc_samples);
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
  const double factor_N2 = safe_ratio(radiative_n2,
                                      radiative_n2 + n * f * K_N2_Q_N2 + n * (1.0 - f) * K_N2_Q_Ar);
  const double factor_Ar_meta = safe_ratio(K_ArMeta_Q_N2c * f * n,
      (K_ArMeta_Q_N2b + K_ArMeta_Q_N2c) * f * n + K_ArMeta_Q_2Ar * std::pow(1.0 - f, 2) * n * n);
  const double factor_Ar_res = safe_ratio(K_ArRes_Q_N2c * f * n,
      (K_ArRes_Q_N2b + K_ArRes_Q_N2c) * f * n + K_ArRes_Q_2Ar * std::pow(1.0 - f, 2) * n * n);

  std::vector<KineticComponent> out;
  const auto n2_peaks = peaks_n2_second_positive();

  push_component(out, "ArN2_UV_N2_C_direct", "UV", sel_n2, n2_peaks, P_N2_pop,
                 Nnorm * factor_N2 * (P_N2_pop * P_N2),
                 factor_N2 * P_N2, tau_N2, mc_samples);
  push_component(out, "ArN2_UV_ArMeta_to_N2C", "UV", sel_ar_meta, n2_peaks, P_Ar_meta,
                 Nnorm * factor_N2 * (P_Ar_meta * factor_Ar_meta),
                 factor_N2 * factor_Ar_meta, tau_N2, mc_samples);
  push_component(out, "ArN2_UV_ArRes_to_N2C", "UV", sel_ar_res, n2_peaks, P_Ar_res,
                 Nnorm * factor_N2 * (P_Ar_res * factor_Ar_res),
                 factor_N2 * factor_Ar_res, tau_N2, mc_samples);
  const double ar_dbl_mix = frac_Ar_dbleStar * factor_Ar_meta + (1.0 - frac_Ar_dbleStar) * factor_Ar_res;
  push_component(out, "ArN2_UV_ArDbleStar_to_N2C", "UV", sel_ar_dbl, n2_peaks, P_Ar_dbl,
                 Nnorm * factor_N2 * P_Ar_dbl * P_Ar_dbleStar * ar_dbl_mix,
                 factor_N2 * P_Ar_dbleStar * ar_dbl_mix, tau_N2, mc_samples);

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
    push_component(out, "ArN2_IR_" + s, "IR", sel, peaks_single(line.lambda, 2.8), P,
                   yield, prob, tau, mc_samples);
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
  const double frac1 = safe_ratio(f * n, f * n + K1);
  const double frac2 = safe_ratio(1.0, 1.0 + K3 * n * f);
  const auto p2nd = read_optional_parameter_csv_values(ar2nd_parameter_csv(parameter_dir));
  const double br_cf4_d_to_x = clamp_value(pget(p2nd, "Br_CF4_D_to_X", pget(p, "Br_CF4_D_to_X", 0.1)), 0.0, 1.0);
  const double lambda_cf4_d_to_x_nm = pget(p2nd, "lambda_CF4_D_to_X_nm", 155.0);
  const double fwhm_cf4_d_to_x_nm = pget(p2nd, "fwhm_CF4_D_to_X_nm", 10.0);

  std::vector<KineticComponent> out;
  push_component(out, "HeCF4_VIS_CF3_direct_ArCF4Fit", "VIS", sel_cf3, peaks_cf3_visible(), P_CF3,
                 N * p_CF3_vis * P_CF3, p_CF3_vis, 0.0, mc_samples);
  push_component(out, "HeCF4_UV_CF3_direct_ArCF4Fit", "UV", sel_cf3, peaks_cf3_uv(), P_CF3,
                 p_CF3_uv * N * p_CF3_vis * P_CF3, p_CF3_uv * p_CF3_vis, 0.0, mc_samples);
  const double uv_cf4_total = N * frac1 * frac2 * p_CF4_dir * P_CF4;
  push_component(out, "HeCF4_UV_CF4_direct_ArCF4Fit", "UV", sel_cf4, peaks_cf4_uv(), P_CF4,
                 uv_cf4_total, frac1 * frac2 * p_CF4_dir, 0.0, mc_samples);
  push_component(out, "HeCF4_VUV155_CF4_direct_ArCF4Fit", "VUV", sel_cf4,
                 peaks_cf4_vuv_branch(lambda_cf4_d_to_x_nm, fwhm_cf4_d_to_x_nm), P_CF4,
                 br_cf4_d_to_x * uv_cf4_total,
                 br_cf4_d_to_x * frac1 * frac2 * p_CF4_dir, 0.0, mc_samples);
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

const PhotonSourceSite* sample_source_site(
    const KineticComponent& component, const std::vector<LevelInfo>& levels,
    const std::vector<PhotonSourceSite>& sites, TRandom3& random);

double sample_wavelength_nm(const KineticComponent& component,
                            TRandom3& random, double minimum_nm,
                            double maximum_nm);

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
