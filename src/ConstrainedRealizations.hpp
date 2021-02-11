#ifndef CORAS_CONSTRAINED_REALIZATIONS_H
#define CORAS_CONSTRAINED_REALIZATIONS_H

#include <cstddef>
#include <functional>
#include <vector>

#include <gsl/gsl_rng.h>

#include "LinearRSDCorrection.hpp"
#include "LogNormalPoissonRealization.hpp"
#include "transformations.hpp"
#include "SphericalFourierBesselDecomposition.hpp"
#include "SphericalGridFunction.hpp"
#include "WienerFilter.hpp"

/**
 * \defgroup RECONSTRUCTION Reconstruction
 *
 * \brief Main tools for the reconstruction of density and velocity fields.
 * @{
 */

/**
 * \brief Class implementing a generator of constrained realizations of the density contrast and peculiar velocity
 * fields, based on a Wiener filter and random log-normal Poisson realizations.
 */
class ConstrainedRealizations
{
public:
	ConstrainedRealizations(const std::vector<double> &radialCoordinates, const std::vector<double> &thetaCoordinates, const std::vector<double> &phiCoordinates, ReferenceFrame redshiftReferenceFrame,
							double maxRadius, double radialResolution, std::size_t maxMultipole, std::size_t radialBinNumber, std::size_t thetaBinNumber, std::size_t phiBinNumber,
							double fftPeriodicBoundaryDistance, std::size_t fftBinNumber, double logTransformSmoothingScale,
							double meanGalaxyDensity, const std::function<double(double)> &selectionFunction, const std::function<double(double)> &selectionFunctionLogDerivative, const std::function<double(double)> &sigmaGalaxy, const std::function<double(double)> &normalizedPowerSpectrum,
							const std::string &directory, const std::string &identifier,
							bool usePrecomputedCouplingAndNoiseMatrices = true,
							const gsl_rng_type *randomNumberGeneratorType = gsl_rng_ranlxd2);

	void generate(std::size_t realizationNumber);

	std::vector<double> &get_random_galaxy_radial_coordinates();
	std::vector<double> &get_random_galaxy_theta_coordinates();
	std::vector<double> &get_random_galaxy_phi_coordinates();

	SphericalGridFunction get_constrained_density_contrast(double sigmaMatter, double normalizedGrowthRate, double smoothingScale);
	SphericalGridFunction get_constrained_radial_velocity(double normalizedGrowthRate, double smoothingScale);
	SphericalGridFunction get_constrained_theta_velocity(double normalizedGrowthRate, double smoothingScale);
	SphericalGridFunction get_constrained_phi_velocity(double normalizedGrowthRate, double smoothingScale);

	SphericalGridFunction get_noise_density_contrast(double sigmaMatter, double smoothingScale);
	SphericalGridFunction get_noise_radial_velocity(double normalizedGrowthRate, double smoothingScale);
	SphericalGridFunction get_noise_theta_velocity(double normalizedGrowthRate, double smoothingScale);
	SphericalGridFunction get_noise_phi_velocity(double normalizedGrowthRate, double smoothingScale);

	SphericalGridFunction get_survey_estimate_density_contrast(double sigmaMatter, double normalizedGrowthRate, double smoothingScale);
	SphericalGridFunction get_survey_estimate_radial_velocity(double normalizedGrowthRate, double smoothingScale);
	SphericalGridFunction get_survey_estimate_theta_velocity(double normalizedGrowthRate, double smoothingScale);
	SphericalGridFunction get_survey_estimate_phi_velocity(double normalizedGrowthRate, double smoothingScale);

	SphericalGridFunction get_random_estimate_density_contrast(double sigmaMatter, double smoothingScale);
	SphericalGridFunction get_random_estimate_radial_velocity(double normalizedGrowthRate, double smoothingScale);
	SphericalGridFunction get_random_estimate_theta_velocity(double normalizedGrowthRate, double smoothingScale);
	SphericalGridFunction get_random_estimate_phi_velocity(double normalizedGrowthRate, double smoothingScale);

	SphericalGridFunction get_random_signal_density_contrast(double sigmaMatter, double smoothingScale);
	SphericalGridFunction get_random_signal_radial_velocity(double normalizedGrowthRate, double smoothingScale);
	SphericalGridFunction get_random_signal_theta_velocity(double normalizedGrowthRate, double smoothingScale);
	SphericalGridFunction get_random_signal_phi_velocity(double normalizedGrowthRate, double smoothingScale);

private:
	double MaxRadius;
	double RadialResolution;
	std::size_t MaxMultipole;
	std::size_t RadialBinNumber;
	std::size_t ThetaBinNumber;
	std::size_t PhiBinNumber;

	double FFTPeriodicBoundaryDistance;
	std::size_t FFTBinNumber;
	double LogTransformSmoothingScale;

	double MeanGalaxyDensity;

	std::function<double(double)> SelectionFunction;
	std::function<double(double)> SigmaGalaxy;
	std::function<double(double)> NormalizedPowerSpectrum;

	LogNormalPoissonRealization RandomRealizations;

	LinearRSDCorrection SurveyRSDCorrection;

	WienerFilter SurveyWienerFilter;
	WienerFilter RandomWienerFilter;

	SphericalFourierBesselDecomposition RawSurveySFBD;
	SphericalFourierBesselDecomposition SurveySFBD;
	SphericalFourierBesselDecomposition RandomSFBD;

	const gsl_rng_type *RandomNumberGeneratorType;

	std::size_t RealizationNumber;

	std::vector<double> RandomGalaxyRadialCoordinates;
	std::vector<double> RandomGalaxyThetaCoordinates;
	std::vector<double> RandomGalaxyPhiCoordinates;

	SphericalGridFunction NormalizedSurveyEstimateDensityContrast;
	SphericalGridFunction NormalizedSurveyEstimateRadialVelocity;
	SphericalGridFunction NormalizedSurveyEstimateThetaVelocity;
	SphericalGridFunction NormalizedSurveyEstimatePhiVelocity;

	SphericalGridFunction NormalizedRandomEstimateDensityContrast;
	SphericalGridFunction NormalizedRandomEstimateRadialVelocity;
	SphericalGridFunction NormalizedRandomEstimateThetaVelocity;
	SphericalGridFunction NormalizedRandomEstimatePhiVelocity;

	SphericalGridFunction NormalizedRandomSignalDensityContrast;
	SphericalGridFunction NormalizedRandomSignalRadialVelocity;
	SphericalGridFunction NormalizedRandomSignalThetaVelocity;
	SphericalGridFunction NormalizedRandomSignalPhiVelocity;

	double SurveyEstimateGrowthRate;
	double SurveyEstimateDensityContrastSmoothingScale;
	double SurveyEstimateRadialVelocitySmoothingScale;
	double SurveyEstimateThetaVelocitySmoothingScale;
	double SurveyEstimatePhiVelocitySmoothingScale;

	double RandomEstimateDensityContrastSmoothingScale;
	double RandomEstimateRadialVelocitySmoothingScale;
	double RandomEstimateThetaVelocitySmoothingScale;
	double RandomEstimatePhiVelocitySmoothingScale;

	double RandomSignalSmoothingScale;

	void reset_survey_fields();
	void reset_random_fields();

	void update_estimate_density_contrast(double smoothingScale,
										  SphericalGridFunction &normalizedDensityContrast, double &normalizedDensityContrastSmoothingScale,
										  const SphericalFourierBesselDecomposition &sfbd);
	void update_estimate_radial_velocity(double smoothingScale,
										 SphericalGridFunction &normalizedRadialVelocity, double &normalizedRadialVelocitySmoothingScale,
										 const SphericalFourierBesselDecomposition &sfbd);
	void update_estimate_theta_velocity(double smoothingScale,
										SphericalGridFunction &normalizedThetaVelocity, double &normalizedThetaVelocitySmoothingScale,
										const SphericalFourierBesselDecomposition &sfbd);
	void update_estimate_phi_velocity(double smoothingScale,
									  SphericalGridFunction &normalizedPhiVelocity, double &normalizedPhiVelocitySmoothingScale,
									  const SphericalFourierBesselDecomposition &sfbd);

	void update_random_signal_fields(double smoothingScale);

	void update_survey_estimate_decomposition(double normalizedGrowthRate);
};

/** @} */

#endif