#include "QuantumEfficiency.hh"

#include <algorithm>

namespace photonfeedback {

QeMaterial load_photoelectron_material(const QeConfig& config) {
  return load_qe_material(config);
}

double quantum_efficiency(const double wavelength_nm, const QeConfig& config,
                          const QeMaterial& material) {
  return qe_for_lambda(wavelength_nm, config, material);
}

PhotoElectronResult emit_photoelectron(
    const double wavelength_nm, const double extraction_efficiency,
    const QeConfig& config, const QeMaterial& material, TRandom3& random) {
  PhotoElectronResult result;
  const double qe = quantum_efficiency(wavelength_nm, config, material);
  if (qe <= 0.0 || random.Uniform() >= qe) return result;
  if (extraction_efficiency <= 0.0 || random.Uniform() >= extraction_efficiency) {
    return result;
  }
  result.emitted = true;
  result.energy_ev = photoelectron_energy_ev(wavelength_nm, config, material);
  return result;
}

}  // namespace photonfeedback
