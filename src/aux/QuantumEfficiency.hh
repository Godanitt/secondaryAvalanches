#ifndef SECONDARY_AVALANCHES_QUANTUM_EFFICIENCY_HH
#define SECONDARY_AVALANCHES_QUANTUM_EFFICIENCY_HH

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <TRandom3.h>

namespace photonfeedback {

constexpr double kHcEvNm = 1239.8419843320026;

struct QeConfig {
  std::string material = "Ti";
  std::string csv_path = "data/qe/qe_materials.csv";
  std::string model_mode = "measured_extended";
  double fallback_qe = 1.0e-4;
  double fallback_threshold_ev = 3.0;
  double qe0 = std::numeric_limits<double>::quiet_NaN();
  double lambda0_nm = std::numeric_limits<double>::quiet_NaN();
  double work_function_override_ev = std::numeric_limits<double>::quiet_NaN();
  double extrapolate_min_nm = 100.0;
  double extrapolate_max_nm = 400.0;
};

struct QePoint {
  double lambda_nm = 0.0;
  double qe = 0.0;
};

struct QeMaterial {
  bool loaded = false;
  std::string material;
  std::string group;
  double work_function_ev = std::numeric_limits<double>::quiet_NaN();
  double threshold_ev = std::numeric_limits<double>::quiet_NaN();
  std::vector<QePoint> points;
};

inline std::string trim_qe(const std::string& input) {
  const auto first = input.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  const auto last = input.find_last_not_of(" \t\r\n");
  return input.substr(first, last - first + 1);
}

inline std::string lower_qe(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

inline std::vector<std::string> split_qe_csv(const std::string& line) {
  std::vector<std::string> fields;
  std::string field;
  bool in_quotes = false;
  for (std::size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (c == '"') {
      if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
        field.push_back('"');
        ++i;
      } else {
        in_quotes = !in_quotes;
      }
    } else if (c == ',' && !in_quotes) {
      fields.push_back(trim_qe(field));
      field.clear();
    } else {
      field.push_back(c);
    }
  }
  fields.push_back(trim_qe(field));
  return fields;
}

inline double parse_qe_double(const std::string& text) {
  if (trim_qe(text).empty()) return std::numeric_limits<double>::quiet_NaN();
  try {
    return std::stod(text);
  } catch (...) {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

inline std::string infer_qe_group(const std::string& material) {
  const std::string name = lower_qe(material);
  const std::vector<std::string> photocathodes = {
      "csi", "cs2te", "cste", "k2cssb", "na2ksb", "bialkali",
      "multialkali", "s20", "gaasp", "gaas", "ingaasp"};
  for (const auto& token : photocathodes) {
    if (name.find(token) != std::string::npos) return "photocathode";
  }
  return "metal";
}

inline bool is_photocathode(const std::string& group) {
  return lower_qe(trim_qe(group)) == "photocathode";
}

inline QeMaterial load_qe_material(const QeConfig& config) {
  if (config.material.empty() || config.model_mode == "constant_threshold") return {};
  std::ifstream input(config.csv_path);
  if (!input) throw std::runtime_error("No puedo abrir QE CSV: " + config.csv_path);

  std::string line;
  std::map<std::string, std::size_t> columns;
  while (std::getline(input, line)) {
    line = trim_qe(line);
    if (line.empty() || line[0] == '#') continue;
    const auto header = split_qe_csv(line);
    for (std::size_t i = 0; i < header.size(); ++i) columns[lower_qe(header[i])] = i;
    break;
  }
  if (columns.empty()) throw std::runtime_error("QE CSV sin cabecera: " + config.csv_path);

  auto field = [&](const std::vector<std::string>& row, const std::string& name) -> std::string {
    const auto it = columns.find(lower_qe(name));
    if (it == columns.end() || it->second >= row.size()) return "";
    return row[it->second];
  };

  QeMaterial model;
  const std::string wanted = lower_qe(trim_qe(config.material));
  while (std::getline(input, line)) {
    const std::string stripped = trim_qe(line);
    if (stripped.empty() || stripped[0] == '#') continue;
    const auto row = split_qe_csv(line);
    const std::string material = field(row, "material");
    if (lower_qe(material) != wanted) continue;

    std::string wavelength_text = field(row, "wavelength_nm");
    if (wavelength_text.empty()) wavelength_text = field(row, "lambda_nm");
    const double wavelength = parse_qe_double(wavelength_text);
    const double qe = parse_qe_double(field(row, "qe"));
    if (!std::isfinite(wavelength) || wavelength <= 0.0 || !std::isfinite(qe) || qe < 0.0) continue;

    if (!model.loaded) {
      model.loaded = true;
      model.material = material;
      model.group = field(row, "group");
      if (model.group.empty()) model.group = infer_qe_group(material);
      model.work_function_ev = parse_qe_double(field(row, "work_function_ev"));
    }
    model.points.push_back({wavelength, std::clamp(qe, 0.0, 1.0)});
  }

  if (!model.loaded || model.points.empty()) {
    throw std::runtime_error("No encuentro QE para el material '" + config.material + "' en " + config.csv_path);
  }
  std::sort(model.points.begin(), model.points.end(),
            [](const QePoint& a, const QePoint& b) { return a.lambda_nm < b.lambda_nm; });

  std::vector<QePoint> unique;
  for (std::size_t i = 0; i < model.points.size();) {
    const double wavelength = model.points[i].lambda_nm;
    double sum = 0.0;
    std::size_t count = 0;
    std::size_t j = i;
    while (j < model.points.size() && std::abs(model.points[j].lambda_nm - wavelength) < 1.0e-9) {
      sum += model.points[j].qe;
      ++count;
      ++j;
    }
    unique.push_back({wavelength, count > 0 ? sum / static_cast<double>(count) : 0.0});
    i = j;
  }
  model.points = std::move(unique);
  model.threshold_ev = (std::isfinite(model.work_function_ev) && model.work_function_ev > 0.0)
      ? model.work_function_ev
      : kHcEvNm / model.points.back().lambda_nm;
  return model;
}

inline double photon_energy_ev(const double lambda_nm) {
  return std::isfinite(lambda_nm) && lambda_nm > 0.0 ? kHcEvNm / lambda_nm : 0.0;
}

inline double interpolate_linear_qe(const std::vector<QePoint>& points,
                                    const double lambda_nm) {
  if (points.empty() || lambda_nm < points.front().lambda_nm || lambda_nm > points.back().lambda_nm) return 0.0;
  auto upper = std::lower_bound(points.begin(), points.end(), lambda_nm,
                                [](const QePoint& point, double value) { return point.lambda_nm < value; });
  if (upper == points.begin()) return upper->qe;
  if (upper == points.end()) return points.back().qe;
  const auto& p1 = *(upper - 1);
  const auto& p2 = *upper;
  const double fraction = (lambda_nm - p1.lambda_nm) / (p2.lambda_nm - p1.lambda_nm);
  return std::clamp(p1.qe + fraction * (p2.qe - p1.qe), 0.0, 1.0);
}

inline double interpolate_log_qe(const std::vector<QePoint>& points,
                                 const double lambda_nm) {
  if (points.empty() || lambda_nm < points.front().lambda_nm || lambda_nm > points.back().lambda_nm) return 0.0;
  auto upper = std::lower_bound(points.begin(), points.end(), lambda_nm,
                                [](const QePoint& point, double value) { return point.lambda_nm < value; });
  if (upper == points.begin()) return upper->qe;
  if (upper == points.end()) return points.back().qe;
  const auto& p1 = *(upper - 1);
  const auto& p2 = *upper;
  if (p1.qe <= 0.0 || p2.qe <= 0.0) return interpolate_linear_qe(points, lambda_nm);
  const double fraction = (lambda_nm - p1.lambda_nm) / (p2.lambda_nm - p1.lambda_nm);
  return std::clamp(std::exp(std::log(p1.qe) + fraction * (std::log(p2.qe) - std::log(p1.qe))), 0.0, 1.0);
}

inline double effective_work_function(const QeConfig& config, const QeMaterial& material,
                                      const double anchor_lambda_nm) {
  double phi = config.work_function_override_ev;
  if (!std::isfinite(phi) || phi <= 0.0) phi = material.work_function_ev;
  if (!std::isfinite(phi) || phi <= 0.0) phi = material.threshold_ev;
  if (!std::isfinite(phi) || phi <= 0.0) phi = config.fallback_threshold_ev;
  const double anchor_energy = photon_energy_ev(anchor_lambda_nm);
  if (anchor_energy > 0.0 && phi >= anchor_energy) phi = 0.95 * anchor_energy;
  return phi;
}

inline double metal_shape(const double lambda_nm, const double phi_ev) {
  const double energy = photon_energy_ev(lambda_nm);
  if (energy <= phi_ev || phi_ev <= 0.0) return 0.0;
  const double excess = (energy - phi_ev) / energy;
  return excess * excess;
}

inline double metal_tail(const double lambda_nm, const QePoint& anchor,
                         const double phi_ev) {
  const double anchor_shape = metal_shape(anchor.lambda_nm, phi_ev);
  const double shape = metal_shape(lambda_nm, phi_ev);
  if (anchor_shape <= 0.0 || shape <= 0.0) return 0.0;
  return std::clamp(anchor.qe * shape / anchor_shape, 0.0, 1.0);
}

inline double qe_for_lambda(const double lambda_nm, const QeConfig& config,
                            const QeMaterial& material) {
  if (config.model_mode == "constant_threshold" || !material.loaded) {
    return photon_energy_ev(lambda_nm) > config.fallback_threshold_ev ? config.fallback_qe : 0.0;
  }
  if (config.model_mode == "measured_table") {
    return interpolate_linear_qe(material.points, lambda_nm);
  }
  if (lambda_nm < config.extrapolate_min_nm || lambda_nm > config.extrapolate_max_nm) return 0.0;
  if (lambda_nm >= material.points.front().lambda_nm && lambda_nm <= material.points.back().lambda_nm) {
    return interpolate_log_qe(material.points, lambda_nm);
  }
  if (is_photocathode(material.group)) return 0.0;
  const QePoint& anchor = lambda_nm < material.points.front().lambda_nm
      ? material.points.front() : material.points.back();
  return metal_tail(lambda_nm, anchor, effective_work_function(config, material, anchor.lambda_nm));
}

inline double photoelectron_threshold_ev(const QeConfig& config,
                                         const QeMaterial& material) {
  if (std::isfinite(config.work_function_override_ev) && config.work_function_override_ev > 0.0) {
    return config.work_function_override_ev;
  }
  if (material.loaded && std::isfinite(material.threshold_ev) && material.threshold_ev > 0.0) {
    return material.threshold_ev;
  }
  return config.fallback_threshold_ev;
}


inline double photoelectron_energy_ev(const double wavelength_nm,
                                      const QeConfig& config,
                                      const QeMaterial& material) {
  return std::max(0.0, photon_energy_ev(wavelength_nm) -
                        photoelectron_threshold_ev(config, material));
}

struct PhotoElectronResult {
  bool emitted = false;
  double energy_ev = 0.0;
};

QeMaterial load_photoelectron_material(const QeConfig& config);
double quantum_efficiency(double wavelength_nm, const QeConfig& config,
                          const QeMaterial& material);
PhotoElectronResult emit_photoelectron(double wavelength_nm,
                                       double extraction_efficiency,
                                       const QeConfig& config,
                                       const QeMaterial& material,
                                       TRandom3& random);

}  // namespace photonfeedback

#endif  // SECONDARY_AVALANCHES_QUANTUM_EFFICIENCY_HH
