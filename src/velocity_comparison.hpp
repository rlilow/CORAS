#ifndef CORAS_VELOCITY_COMPARISON_H
#define CORAS_VELOCITY_COMPARISON_H

#include <cstddef>
#include <functional>
#include <vector>

#include <gsl/gsl_matrix.h>

#include "../src/Cartesian1DGridFunction.hpp"
#include "../src/SphericalGridFunction.hpp"
#include "../src/transformations.hpp"

/**
 * \addtogroup RECONSTRUCTION
 * @{
 */

/**
 * \name Velocity Comparison
 *
 * Methods for the comparison of reconstructed and observed velocities.
 */
///@{

void estimate_parameters_via_radial_velocity_comparison(const std::function<double(double normalizedGrowthRate, double radius, double theta, double phi)> &normalizedWienerRadialVelocity, const SphericalGridFunction &normalizedNoiseRadialVelocity, const ReferenceFrameChange cmbToReferenceFrame,
                                                        const std::vector<double> &redshiftVelocities, const std::vector<double> &thetaCoordinates, const std::vector<double> &phiCoordinates, const std::vector<double> &distanceModuli, const std::vector<double> &distanceModulusErrors,
                                                        double omegaMatter, double normalizedGrowthRateGuess, double normalizedGrowthRateMin, double normalizedGrowthRateMax, double hubbleGuess,
                                                        double &normalizedGrowthRateEstimate, double &externalBulkXVelocityEstimate, double &externalBulkYVelocityEstimate, double &externalBulkZVelocityEstimate, double &hubbleEstimate,
                                                        gsl_matrix *covarianceMatrix, double &reducedChiSquared);

void compute_radial_velocity_correlation_functions(const std::vector<double> &observedRedshiftVelocities, const std::vector<double> &observedThetaCoordinates, const std::vector<double> &observedPhiCoordinates, const std::vector<double> &observedDistanceModuli,
                                                   const SphericalGridFunction &reconstructedRadialVelocity,
                                                   ReferenceFrameChange refToCMBFrame, double omegaMatter, double hubble,
                                                   double maxDistance, std::size_t distanceBinNumber,
                                                   Cartesian1DGridFunction &observedRadialVelocityCorrFunc, Cartesian1DGridFunction &reconstructedRadialVelocityCorrFunc, Cartesian1DGridFunction &residualRadialVelocityCorrFunc);

void compute_tensor_smoothed_radial_velocity_points(const std::vector<double> &observedRedshiftVelocities, const std::vector<double> &observedThetaCoordinates, const std::vector<double> &observedPhiCoordinates, const std::vector<double> &observedDistanceModuli, const std::vector<double> &observedDistanceModulusErrors,
                                                    const SphericalGridFunction &reconstructedRadialVelocity,
                                                    ReferenceFrameChange refToCMBFrame, double omegaMatter, double hubble, double minSmoothingScale,
                                                    std::vector<double> &smoothObservedRadialVelocities, std::vector<double> &smoothReconstructedRadialVelocities, std::vector<double> &adaptiveSmoothingScales,
                                                    std::size_t minNeighbourNumber = 8);

///@}

/** @} */

#endif