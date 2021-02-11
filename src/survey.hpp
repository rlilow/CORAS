#ifndef CORAS_SURVEY_H
#define CORAS_SURVEY_H

#include <cstddef>
#include <functional>
#include <vector>

#include <gsl/gsl_rng.h>

#include "Cartesian1DGridFunction.hpp"
#include "cosmology.hpp"
#include "transformations.hpp"

/**
 * \addtogroup RECONSTRUCTION
 * @{
 */

/**
 * \name Survey
 *
 * Preparation and analysis of survey data.
 */
///@{

void apply_partial_volume_limit(const std::vector<double> &inputRedshiftVelocities, const std::vector<double> &inputThetaCoordinates, const std::vector<double> &inputPhiCoordinates, const std::vector<double> &inputApparentMagnitudes,
                                ReferenceFrameChange referenceFrameChange, double maxRadius, double volumeLimitRadius, double maxApparentMagnitude, double omegaMatter,
                                std::vector<double> &outputRedshiftVelocities, std::vector<double> &outputRadialCoordinates, std::vector<double> &outputThetaCoordinates, std::vector<double> &outputPhiCoordinates, std::vector<double> &outputApparentMagnitudes,
                                const std::function<double(double redshiftVelocity)> &luminosityEvolutionCorrection = luminosity_evolution_correction_2MRS, const std::function<double(double redshiftVelocity)> &kCorrection = k_correction_2MRS);

/**
 * Prepare the 2MRS data for the reconstruction process. The input data consist of the redshift velocities \a
 * inputRedshiftVelocities, galactic latitudes \a inputLatitudes, galactic longitudes \a inputLongitudes and apparent
 * magnitudes \a inputApparentMagnitudes of all galaxies. The reference frame used for the redshifts is transformed
 * according to \a referenceFrameChange, and all those galaxies are removed which are blueshifted or lie beyond the
 * comoving radius \a maxRadius. Afterwards a partial volume limit is applied by removing all galaxies within the
 * comoving radius \a volumeLimitRadius which would be too faint to be detectable beyond that radius, given the maximal
 * apparent magnitude \a maxApparentMagnitude. Both of these steps assume a flat LambdaCDM cosmology with dimensionless
 * matter density parameter \a omegaMatter, and the functions of redshift velocity \a luminosityEvolutionCorrection
 * default: luminosity_evolution_correction_2MRS) and \a kCorrection (default: k_correction_2MRS) accounting for galaxy
 * luminosity evolution and k-correction, respectively. By default these are set to the functions found to be suitable
 * for 2MRS in [Branchini+ MNRAS 424 (2012) 472]. After that, if \a fillZOA (default: \c true) is \c true, the zone of
 * avoidance in the galactic plane is filled up with mock galaxies sampled from the adjacent regions, using a GSL random
 * number generator of type \a randomNumberGeneratorType (default: gsl_rng_ranlxd2) with seed \a
 * randomNumberGeneratorSeed (default: 1) (see GSL documentation for details). The output consists of the redshift
 * velocities \a outputRedshiftVelocities, comoving radial redshift coordinates \a outputRadialCoordinates, theta
 * coordinates \a outputThetaCoordinates, phi coordinates \a outputPhiCoordinates and apparent magnitudes \a
 * outputApparentMagnitudes of the resulting set of galaxies.
 */
void prepare_2MRS_data(const std::vector<double> &inputRedshiftVelocities, const std::vector<double> &inputLatitudes, const std::vector<double> &inputLongitudes, const std::vector<double> &inputApparentMagnitudes,
                       ReferenceFrameChange referenceFrameChange, double maxRadius, double volumeLimitRadius, double maxApparentMagnitude, double omegaMatter,
                       std::vector<double> &outputRedshiftVelocities, std::vector<double> &outputRadialCoordinates, std::vector<double> &outputThetaCoordinates, std::vector<double> &outputPhiCoordinates, std::vector<double> &outputApparentMagnitudes,
                       const std::function<double(double redshiftVelocity)> &luminosityEvolutionCorrection = luminosity_evolution_correction_2MRS, const std::function<double(double redshiftVelocity)> &kCorrection = k_correction_2MRS,
                       bool fillZOA = true, std::size_t randomNumberGeneratorSeed = 1, const gsl_rng_type *randomNumberGeneratorType = gsl_rng_ranlxd2);

void prepare_distance_catalog_data(const std::vector<double> &inputRedshiftVelocities, const std::vector<double> &inputLatitudes, const std::vector<double> &inputLongitudes, const std::vector<double> &inputDistanceModuli, const std::vector<double> &inputDistanceModulusErrors,
                                   ReferenceFrameChange inputToReferenceFrame, double minRedshiftVelocity, double maxRedshiftVelocity, double omegaMatter,
                                   const SphericalGridFunction &fiducialRadialVelocity, ReferenceFrameChange cmbToReferenceFrame, double fiducialHubble, double maxSigmaDeviation,
                                   std::vector<double> &outputRedshiftVelocities, std::vector<double> &outputRadialCoordinates, std::vector<double> &outputThetaCoordinates, std::vector<double> &outputPhiCoordinates, std::vector<double> &outputDistanceModuli, std::vector<double> &outputDistanceModulusErrors);

/**
 * Estimate the mean galaxy number density within the radius \a radius given the distances \a galaxyDistances to the
 * observed galaxies and their radial selection function \c selectionFunction.
 */
double estimate_mean_density(double radius, const std::vector<double> &observedGalaxyDistances, const std::function<double(double)> &selectionFunction);

/**
 * Use the F/T estimator [Davis+ ApJ 254 (1982) 437, Branchini+ MNRAS 424 (2012) 472] to infer the selection function \a
 * selectionFunction and its logarithmic derivative \a selectionFunctionLogDerivative of a flux-limited redshift survey,
 * specified by the observed redshift velocities and apparent magnitudes of all galaxies, \a redshiftVelocities and \a
 * apparentMagnitudes, as well as the maximal detectable apparent magnitude \a maxApparentMagnitude. In addition, the
 * radial distribution function \a radialDistributionFunction of the survey is computed, too. All three functions are
 * estimated up to a maximal comoving radius \a maxRadius, divided into \a radialBinNumber bins, assuming a flat
 * LambdaCDM cosmology with dimensionless matter density parameter \a omegaMatter. The functions of redshift velocity \a
 * luminosityEvolutionCorrection (default: luminosity_evolution_correction_2MRS) and \a kCorrection (default:
 * k_correction_2MRS) account for galaxy luminosity evolution and k-correction, respectively. TBy default these are set
 * to the functions found to be suitable for 2MRS in [Branchini+ MNRAS 424 (2012) 472]
 */
void estimate_selection_function(const std::vector<double> &redshiftVelocities, const std::vector<double> &apparentMagnitudes, double maxApparentMagnitude,
                                 double maxRadius, std::size_t radialBinNumber, double omegaMatter,
                                 Cartesian1DGridFunction &selectionFunction, Cartesian1DGridFunction &selectionFunctionLogDerivative, Cartesian1DGridFunction &radialDistributionFunction,
                                 const std::function<double(double redshiftVelocity)> &luminosityEvolutionCorrection = luminosity_evolution_correction_2MRS, const std::function<double(double redshiftVelocity)> &kCorrection = k_correction_2MRS);

/**
 * Estimate the RMS of the galaxy number density contrast, smoothed with a top-hat filter of width \a smoothingScale,
 * for a set of galaxies at redshift velocities \a redshiftVelocities, and angular positions \a thetaCoordinates, \a
 * phiCoordinates, with apparent magnitudes \a apparentMagnitudes. It is computed from volume-limited galaxy subsamples,
 * given the maximal observable magnitude \a maxApparentMagnitude, in \a radialBinNumber bins between the distances \a
 * minRadius and \a maxRadius. The result is written to \a sigmaGalaxy. In addition, the number of galaxies contained in
 * each subsample is written to \a galaxyNumbers. The dimensionless matter density parameter \a omegaMatter is used to
 * translate between redshifts and distances. The functions of redshift velocity \a luminosityEvolutionCorrection
 * (default: luminosity_evolution_correction_2MRS) and \a kCorrection (default: k_correction_2MRS) account for galaxy
 * luminosity evolution and k-correction, respectively. By default these are set to the functions found to be suitable
 * for 2MRS in [Branchini+ MNRAS 424 (2012) 472]
 */
void estimate_sigma_galaxy(const std::vector<double> &redshiftVelocities, const std::vector<double> &thetaCoordinates, const std::vector<double> &phiCoordinates, const std::vector<double> &apparentMagnitudes,
                           const double minRadius, const double maxRadius, const std::size_t radialBinNumber, const double maxApparentMagnitude, const double omegaMatter, const double smoothingScale,
                           Cartesian1DGridFunction &sigmaGalaxy, std::vector<std::size_t> &galaxyNumbers,
                           const std::function<double(double redshift)> &luminosityEvolutionCorrection = luminosity_evolution_correction_2MRS, const std::function<double(double redshift)> &kCorrection = k_correction_2MRS);

///@}

/** @} */

#endif