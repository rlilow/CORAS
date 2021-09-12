#ifndef CORAS_COSMOLOGY_H
#define CORAS_COSMOLOGY_H

#include <functional>

/**
 * \defgroup COSMOLOGY Cosmology
 *
 * \brief Distances and magnitudes in a Lambda-CDM cosmology.
 * @{
 */

/**
 * Speed of light in km/s.
 */
constexpr double SPEED_OF_LIGHT = 2.99792458e5;

/**
 * Hubble constant in (km/s)/(Mpc/h).
 */
constexpr double HUBBLE_NORMALIZATION = 100.0;

constexpr double HUBBLE_RADIUS = SPEED_OF_LIGHT / HUBBLE_NORMALIZATION;

double expansion_function(double redshiftVelocity, double omegaMatter);

double proper_distance(double redshiftVelocity, double omegaMatter);
double comoving_distance(double redshiftVelocity, double omegaMatter);
double luminosity_distance(double redshiftVelocity, double omegaMatter);

double proper_distance_redshift_velocity_derivative(double redshiftVelocity, double omegaMatter);
double comoving_distance_redshift_velocity_derivative(double redshiftVelocity, double omegaMatter);
double luminosity_distance_redshift_velocity_derivative(double redshiftVelocity, double omegaMatter);

double redshift_velocity_from_proper_distance(double comovingDistance, double omegaMatter);
double redshift_velocity_from_comoving_distance(double comovingDistance, double omegaMatter);
double redshift_velocity_from_luminosity_distance(double luminosityDistance, double omegaMatter);

double distance_modulus(double redshiftVelocity, double omegaMatter, double hubble);
double distance_modulus_redshift_velocity_derivative(double redshiftVelocity, double omegaMatter);

double k_correction_2MRS(double redshiftVelocity);
double luminosity_evolution_correction_2MRS(double redshiftVelocity);

/**
 * Return the absolute magnitude of an object with cosmological and observed redshift velocities \a
 * observedRedshiftVelocity and \a cosmologicalRedshiftVelocity, respectively, and apparent magnitude \a
 * apparentMagnitude, using the dimensionless matter density parameter \a omegaMatter, the dimensionless Hubble constant
 * \a hubble, and the functions \a luminosityEvolutionCorrection and \a kCorrection accounting for galaxy luminosity
 * evolution and k-correction, depending on the cosmological and observed redshift velocities, respectively.
 */
double absolute_magnitude(double cosmologicalRedshiftVelocity, double observedRedshiftVelocity, double apparentMagnitude, double omegaMatter, double hubble,
                          const std::function<double(double cosmologicalRedshiftVelocity)> &luminosityEvolutionCorrection, const std::function<double(double observedRedshiftVelocity)> &kCorrection);

/**
 * Return the apparent magnitude of an object with cosmological and observed redshift velocities \a
 * observedRedshiftVelocity and \a cosmologicalRedshiftVelocity, respectively, and apparent magnitude \a
 * apparentMagnitude, using the dimensionless matter density parameter \a omegaMatter, the dimensionless Hubble constant
 * \a hubble, and the functions \a luminosityEvolutionCorrection and \a kCorrection accounting for galaxy luminosity
 * evolution and k-correction, depending on the cosmological and observed redshift velocities, respectively.
 */
double apparent_magnitude(double cosmologicalRedshiftVelocity, double observedRedshiftVelocity, double absoluteMagnitude, double omegaMatter, double hubble,
                          const std::function<double(double cosmologicalRedshiftVelocity)> &luminosityEvolutionCorrection, const std::function<double(double observedRedshiftVelocity)> &kCorrection);

/** @} */

#endif