#include "PhotonModel.hh"

#include <algorithm>
#include <cmath>
#include <limits>

namespace photonemission {

std::vector<KineticComponent> create_photon_emission(
    const GasData& gas, const std::vector<LevelPopulation>& populations,
    const std::string& parameters_dir, const int mc_samples) {
  return build_kinetic_components(gas, populations, parameters_dir, mc_samples);
}

double slowest_emission_lifetime_ns(
    const std::vector<KineticComponent>& components) {
  double lifetime_ns = 0.0;
  for (const auto& component : components) {
    lifetime_ns =
        std::max(lifetime_ns, std::max(0.0, component.mean_total_delay_ns));
  }
  return lifetime_ns;
}

}  // namespace photonemission

namespace photonemission {
namespace {

const SpectralPeak& sample_peak(const KineticComponent& component,
                                TRandom3& random) {
  if (component.peaks.empty()) {
    throw std::runtime_error("Photon component without spectral peaks.");
  }
  double total = 0.0;
  for (const auto& peak : component.peaks) total += std::max(0.0, peak.weight);
  if (total <= 0.0) return component.peaks.front();

  const double target = random.Uniform(0.0, total);
  double accumulated = 0.0;
  for (const auto& peak : component.peaks) {
    accumulated += std::max(0.0, peak.weight);
    if (target <= accumulated) return peak;
  }
  return component.peaks.back();
}

}  // namespace

std::vector<LevelPopulation> populations_from_sites(
    const std::vector<LevelInfo>& levels,
    const std::vector<PhotonSourceSite>& sites) {
  std::vector<LevelPopulation> populations(levels.size());
  for (std::size_t i = 0; i < levels.size(); ++i) {
    populations[i].level = static_cast<int>(i);
    populations[i].known_level = true;
    populations[i].level_info = levels[i];
  }
  for (const auto& site : sites) {
    if (site.level < 0 || site.level >= static_cast<int>(populations.size())) {
      continue;
    }
    populations[static_cast<std::size_t>(site.level)].n_events += 1.0;
  }
  return populations;
}

std::vector<const PhotonSourceSite*> matching_source_sites(
    const KineticComponent& component, const std::vector<LevelInfo>& levels,
    const std::vector<PhotonSourceSite>& sites) {
  std::vector<const PhotonSourceSite*> matching;
  matching.reserve(sites.size());
  for (const auto& site : sites) {
    if (site.level < 0 || site.level >= static_cast<int>(levels.size())) continue;
    if (selector_matches_level(
            component.source, levels[static_cast<std::size_t>(site.level)])) {
      matching.push_back(&site);
    }
  }
  return matching;
}

const PhotonSourceSite* sample_source_site(
    const std::vector<const PhotonSourceSite*>& matching_sites,
    TRandom3& random) {
  if (matching_sites.empty()) return nullptr;
  return matching_sites[static_cast<std::size_t>(
      random.Integer(matching_sites.size()))];
}

const PhotonSourceSite* sample_source_site(
    const KineticComponent& component, const std::vector<LevelInfo>& levels,
    const std::vector<PhotonSourceSite>& sites, TRandom3& random) {
  const auto matching = matching_source_sites(component, levels, sites);
  return sample_source_site(matching, random);
}

double sample_wavelength_nm(const KineticComponent& component,
                            TRandom3& random, const double minimum_nm,
                            const double maximum_nm) {
  const auto& peak = sample_peak(component, random);
  if (!std::isfinite(peak.sigma_nm) || peak.sigma_nm <= 0.0) {
    return std::clamp(peak.center_nm, minimum_nm, maximum_nm);
  }
  for (int attempt = 0; attempt < 100; ++attempt) {
    const double wavelength = random.Gaus(peak.center_nm, peak.sigma_nm);
    if (wavelength >= minimum_nm && wavelength <= maximum_nm) {
      return wavelength;
    }
  }
  return std::clamp(peak.center_nm, minimum_nm, maximum_nm);
}

EmissionDelaySample sample_emission_delay(
    const KineticComponent& component, TRandom3& random) {
  EmissionDelaySample sample;
  for (const auto& stage : component.kinetic_stages) {
    if (!std::isfinite(stage.tau_total_ns) || stage.tau_total_ns <= 0.0) {
      continue;
    }
    const double waiting_time_ns = random.Exp(stage.tau_total_ns);
    if (stage.kind == KineticStageKind::Emission) {
      sample.emission_delay_ns += waiting_time_ns;
    } else {
      sample.excitation_delay_ns += waiting_time_ns;
    }
  }
  return sample;
}

double sample_emission_delay_ns(const KineticComponent& component,
                                TRandom3& random) {
  return sample_emission_delay(component, random).total_delay_ns();
}

long long number_of_mc_photons(const double expected_photons,
                               const long long mc_samples,
                               TRandom3& random) {
  if (!std::isfinite(expected_photons) || expected_photons <= 0.0 ||
      mc_samples <= 0) {
    return 0;
  }
  const long double exact = static_cast<long double>(expected_photons) *
                            static_cast<long double>(mc_samples);
  if (exact >= static_cast<long double>(
                   std::numeric_limits<long long>::max())) {
    throw std::runtime_error(
        "The number of Monte Carlo photons exceeds Long64_t; reduce mcSamples.");
  }
  const long long integer = static_cast<long long>(std::floor(exact));
  const double fraction = static_cast<double>(exact - integer);
  return integer + (random.Uniform() < fraction ? 1LL : 0LL);
}

}  // namespace photonemission
