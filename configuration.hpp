#ifndef CORAS_CONFIGURATION_H
#define CORAS_CONFIGURATION_H

#include <cmath>
#include <string>

#include "src/transformations.hpp"

/**
 * \defgroup CONFIGURATION Configuration
 * 
 * \brief Specification of the different parameters of CORAS.
 * @{
 */

const std::string CONFIGURATION_COMMENT = "";

const std::string DATA_DIRECTORY = "data/";

const std::string FIDUCIAL_POWER_SPECTRUM_FILE_NAME = "data/CosmicEmu_spectrum_Planck18.dat";
const std::string REDSHIFT_CATALOG_FILE_NAME = "data/2MRS_group_member_catalog.dat";
const std::string DISTANCE_CATALOG_FILE_NAME = "data/CF3_group_catalog.dat";

const ReferenceFrame REDSHIFT_CATALOG_REDSHIFT_REFERENCE_FRAME = HELIOCENTRIC_FRAME;
constexpr bool REDSHIFT_CATALOG_COLLAPSE_FINGERS_OF_GOD = true;
constexpr double REDSHIFT_CATALOG_MAX_APPARENT_MAGNITUDE = 11.754; // 2MRS
constexpr double REDSHIFT_CATALOG_VOLUME_LIMIT_RADIUS = 30.0;      // [Mpc/h]

const ReferenceFrame DISTANCE_CATALOG_REDSHIFT_REFERENCE_FRAME = HELIOCENTRIC_FRAME;
constexpr double DISTANCE_CATALOG_MAX_REDSHIFT_VELOCITY = 16000.0; // [km/s]
constexpr double DISTANCE_CATALOG_OUTLIER_MAX_SIGMA = 5.0;

constexpr double FIDUCIAL_HUBBLE = 0.6736; // Planck 2018 values [Aghanim et al., A&A 641 (2020) A1]
constexpr double FIDUCIAL_OMEGA_MATTER = 0.3153;
const double FIDUCIAL_GROWTH_RATE = std::pow(FIDUCIAL_OMEGA_MATTER, 0.55);
constexpr double FIDUCIAL_FIELD_SMOOTHING_SCALE = 5.0; // [Mpc/h]

const std::vector<ReferenceFrame> RECONSTRUCTION_REDSHIFT_REFERENCE_FRAMES = {CMB_FRAME, LOCAL_GROUP_FRAME};
constexpr double RECONSTRUCTION_MAX_RADIUS = 200.0;        // [Mpc/h]
constexpr double RECONSTRUCTION_RADIAL_RESOLUTION = 120.0; // max. radial Fourier mode times max. radius
constexpr std::size_t RECONSTRUCTION_MAX_MULTIPOLE = 60;
constexpr std::size_t RECONSTRUCTION_RADIAL_BIN_NUMBER = 101;
constexpr std::size_t RECONSTRUCTION_THETA_BIN_NUMBER = 51;
constexpr std::size_t RECONSTRUCTION_PHI_BIN_NUMBER = 100;
constexpr std::size_t RECONSTRUCTION_FFT_PERIODIC_BOUNDARY_DISTANCE = 2.0 * RECONSTRUCTION_MAX_RADIUS;
constexpr std::size_t RECONSTRUCTION_FFT_BIN_NUMBER = 150;
constexpr double RECONSTRUCTION_LOG_TRAFO_SMOOTHING_SCALE = 1.0; // [Mpc/h]
constexpr bool RECONSTRUCTION_USE_PRECOMPUTED_RSD_CORRECTION_AND_WIENER_FILTER = true;

constexpr std::size_t INITIAL_CONSTRAINED_REALIZATION = 1;
constexpr std::size_t FINAL_CONSTRAINED_REALIZATION = 50;

constexpr double SIGMA_SCALE = 8.0; // [Mpc/h]
constexpr std::size_t SIGMA_GALAXY_RADIAL_BIN_NUMBER = 18;

constexpr std::size_t SELECTION_FUNCTION_BIN_NUMBER = 101;

const std::vector<double> PARAMETER_ESTIMATE_SMOOTHING_SCALES = {5.0, 10.0, 15.0, 20.0, 25.0, 30.0};                  // [Mpc/h]
const std::vector<double> PARAMETER_ESTIMATE_MIN_REDSHIFT_VELOCITIES = {0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0}; // [km/s]
constexpr double PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_MIN = 0.0;
constexpr double PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_MAX = 2.0;
constexpr std::size_t PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_BIN_NUMBER = 21;

constexpr double ESTIMATED_NORMALIZED_GROWTH_RATE = 0.402;    // these parameters have been inferred by comparing the reconstructed velocities from 2MRS with the observed velocities from Cosmicflows-3
constexpr double ESTIMATED_EXTERNAL_BULK_X_VELOCITY = 107.0;  // [km/s]
constexpr double ESTIMATED_EXTERNAL_BULK_Y_VELOCITY = -124.0; // [km/s]
constexpr double ESTIMATED_EXTERNAL_BULK_Z_VELOCITY = -1.0;   // [km/s]
constexpr double ESTIMATED_DISTANCE_CATALOG_HUBBLE = 0.7272;

constexpr double RECONSTRUCTED_FIELD_SMOOTHING_SCALE = 5.0;     // [Mpc/h]
constexpr double RECONSTRUCTED_FIELD_MIN_SMOOTHING_SCALE = 1.0; // [Mpc/h]

constexpr std::size_t CARTESIAN_GRID_BIN_NUMBER = 201;

constexpr std::size_t SGP_BIN_NUMBER = 201;
const std::vector<std::size_t> SGP_NOISE_REALIZATIONS = {1, 2};

const std::size_t GAUSSIAN_VOLUME_AVERAGE_RADIAL_BIN_NUMBER = 401;

const std::vector<double> TENSOR_SMOOTHING_SCALES = {15.0, 30.0}; // [Mpc/h]
constexpr std::size_t TENSOR_SMOOTHING_RADIAL_BIN_NUMBER = 21;
constexpr std::size_t TENSOR_SMOOTHING_THETA_BIN_NUMBER = 61;
constexpr std::size_t TENSOR_SMOOTHING_PHI_BIN_NUMBER = 121;

constexpr double CORRELATION_FUNCTION_MAX_DISTANCE = 100.0; // [Mpc/h]
constexpr std::size_t CORRELATION_FUNCTION_DISTANCE_BIN_NUMBER = 21;
constexpr double CORRELATION_FUNCTION_MAX_REDSHIFT_VELOCITY = 10000.0; // [km/s]

/** @} */

#endif