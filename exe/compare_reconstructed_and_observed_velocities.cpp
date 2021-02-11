#include <chrono>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

#include <fftw3.h>
#include <omp.h>

#include "../src/Cartesian1DGridFunction.hpp"
#include "../src/ConstrainedRealizations.hpp"
#include "../src/FileTable.hpp"
#include "../src/NormalizedPowerSpectrum.hpp"
#include "../src/SphericalGridFunction.hpp"
#include "../src/survey.hpp"
#include "../src/transformations.hpp"
#include "../src/velocity_comparison.hpp"

#include "../configuration.hpp"

int main(int argc, char **argv)
{
    fftw_init_threads(); // initialize parallel execution of the FFTs used to generate random realizations
    fftw_plan_with_nthreads(omp_get_max_threads());

    const FileTable inputPowerSpectrumFileTable(FIDUCIAL_POWER_SPECTRUM_FILE_NAME, 2); // read the power spectrum from file

    const auto inputWavenumbers = inputPowerSpectrumFileTable.column<double>(0);
    const auto inputPowerSpectrumValues = inputPowerSpectrumFileTable.column<double>(1);

    NormalizedPowerSpectrum normalizedPowerSpectrum(inputWavenumbers, inputPowerSpectrumValues,
                                                    FIDUCIAL_HUBBLE, SIGMA_SCALE);

    const double fiducialNormalizedGrowthRate = FIDUCIAL_GROWTH_RATE * std::sqrt(normalizedPowerSpectrum.variance());

    const FileTable inputRedshiftCatalogFile(REDSHIFT_CATALOG_FILE_NAME, 10, '\t'); // read the redshift catalog properties from file

    const auto inputRedshiftCatalogLatitudes = inputRedshiftCatalogFile.column<double>(1);
    const auto inputRedshiftCatalogLongitudes = inputRedshiftCatalogFile.column<double>(2);
    const auto inputRedshiftCatalogRedshiftVelocities = REDSHIFT_CATALOG_COLLAPSE_FINGERS_OF_GOD ? inputRedshiftCatalogFile.column<double>(9)
                                                                                                 : inputRedshiftCatalogFile.column<double>(3);
    const auto inputRedshiftCatalogApparentMagnitudes = inputRedshiftCatalogFile.column<double>(4);

    const FileTable inputDistanceCatalogFile(DISTANCE_CATALOG_FILE_NAME, 7); // read the distance catalog properties from file

    const auto inputDistanceCatalogLatitudes = inputDistanceCatalogFile.column<double>(2);
    const auto inputDistanceCatalogLongitudes = inputDistanceCatalogFile.column<double>(3);
    const auto inputDistanceCatalogRedshiftVelocities = inputDistanceCatalogFile.column<double>(4);
    const auto inputDistanceCatalogDistanceModuli = inputDistanceCatalogFile.column<double>(5);
    const auto inputDistanceCatalogDistanceModulusErrors = inputDistanceCatalogFile.column<double>(6);

    for (std::size_t i_ref = 0; i_ref < RECONSTRUCTION_REDSHIFT_REFERENCE_FRAMES.size(); ++i_ref)
    {
        const ReferenceFrame redshiftReferenceFrame = RECONSTRUCTION_REDSHIFT_REFERENCE_FRAMES[i_ref];

        const ReferenceFrameChange redshiftCatalogInputToReferenceFrame = get_reference_frame_change(REDSHIFT_CATALOG_REDSHIFT_REFERENCE_FRAME, redshiftReferenceFrame);
        const ReferenceFrameChange distanceCatalogInputToReferenceFrame = get_reference_frame_change(DISTANCE_CATALOG_REDSHIFT_REFERENCE_FRAME, redshiftReferenceFrame);
        const ReferenceFrameChange cmbToReferenceFrame = get_reference_frame_change(CMB_FRAME, redshiftReferenceFrame);
        const ReferenceFrameChange referenceToCMBFrame = get_reference_frame_change(redshiftReferenceFrame, CMB_FRAME);

        std::string redshiftReferenceFrameName;

        if (redshiftReferenceFrame == CMB_FRAME)
        {
            redshiftReferenceFrameName = "CMB";
        }
        else if (redshiftReferenceFrame == LOCAL_GROUP_FRAME)
        {
            redshiftReferenceFrameName = "LG";
        }
        else
        {
            std::cout << std::endl
                      << " Error: Invalid choice of reference frame" << std::endl
                      << std::endl;

            exit(EXIT_FAILURE);
        }

        const auto time1 = std::chrono::high_resolution_clock::now();
        std::cout << "Preparing " << redshiftReferenceFrameName << " reference frame..." << std::flush;

        std::vector<double> reconstructionRedshiftVelocities, reconstructionRadialCoordinates, reconstructionThetaCoordinates, reconstructionPhiCoordinates, reconstructionApparentKsMagnitudes;

        prepare_2MRS_data(inputRedshiftCatalogRedshiftVelocities, inputRedshiftCatalogLatitudes, inputRedshiftCatalogLongitudes, inputRedshiftCatalogApparentMagnitudes,
                          redshiftCatalogInputToReferenceFrame, RECONSTRUCTION_MAX_RADIUS, REDSHIFT_CATALOG_VOLUME_LIMIT_RADIUS, REDSHIFT_CATALOG_MAX_APPARENT_MAGNITUDE, FIDUCIAL_OMEGA_MATTER,
                          reconstructionRedshiftVelocities, reconstructionRadialCoordinates, reconstructionThetaCoordinates, reconstructionPhiCoordinates, reconstructionApparentKsMagnitudes);

        Cartesian1DGridFunction sigmaGalaxyAboveVolumeLimit;
        std::vector<std::size_t> volumeLimitGalaxyNumbers;

        estimate_sigma_galaxy(reconstructionRedshiftVelocities, reconstructionThetaCoordinates, reconstructionPhiCoordinates, reconstructionApparentKsMagnitudes,
                              REDSHIFT_CATALOG_VOLUME_LIMIT_RADIUS, RECONSTRUCTION_MAX_RADIUS, SIGMA_GALAXY_RADIAL_BIN_NUMBER, REDSHIFT_CATALOG_MAX_APPARENT_MAGNITUDE, FIDUCIAL_OMEGA_MATTER, SIGMA_SCALE,
                              sigmaGalaxyAboveVolumeLimit, volumeLimitGalaxyNumbers);

        auto sigmaGalaxy = [&](double radius) {
            if (radius < sigmaGalaxyAboveVolumeLimit.coordinate(0))
            {
                return sigmaGalaxyAboveVolumeLimit.value(0);
            }
            else
            {
                return sigmaGalaxyAboveVolumeLimit(radius);
            }
        };

        Cartesian1DGridFunction selectionFunction, selectionFunctionLogDerivative, radialDistributionFunction;

        estimate_selection_function(reconstructionRedshiftVelocities, reconstructionApparentKsMagnitudes, REDSHIFT_CATALOG_MAX_APPARENT_MAGNITUDE,
                                    RECONSTRUCTION_MAX_RADIUS, SELECTION_FUNCTION_BIN_NUMBER, FIDUCIAL_OMEGA_MATTER,
                                    selectionFunction, selectionFunctionLogDerivative, radialDistributionFunction);

        const double meanGalaxyDensity = estimate_mean_density(RECONSTRUCTION_MAX_RADIUS, reconstructionRadialCoordinates, selectionFunction);

        const std::string redshiftReferenceFrameComment = "_z" + redshiftReferenceFrameName;
        const std::string precomputedRSDCorrectionAndWienerFilterComment = redshiftReferenceFrameComment + CONFIGURATION_COMMENT;

        ConstrainedRealizations constrainedRealizations(reconstructionRadialCoordinates, reconstructionThetaCoordinates, reconstructionPhiCoordinates, redshiftReferenceFrame,
                                                        RECONSTRUCTION_MAX_RADIUS, RECONSTRUCTION_RADIAL_RESOLUTION, RECONSTRUCTION_MAX_MULTIPOLE, RECONSTRUCTION_RADIAL_BIN_NUMBER, RECONSTRUCTION_THETA_BIN_NUMBER, RECONSTRUCTION_PHI_BIN_NUMBER,
                                                        RECONSTRUCTION_FFT_PERIODIC_BOUNDARY_DISTANCE, RECONSTRUCTION_FFT_BIN_NUMBER, RECONSTRUCTION_LOG_TRAFO_SMOOTHING_SCALE, // set FFT bin number to 1, as no random realizations will be generated
                                                        meanGalaxyDensity, selectionFunction, selectionFunctionLogDerivative, sigmaGalaxy, normalizedPowerSpectrum,
                                                        DATA_DIRECTORY, precomputedRSDCorrectionAndWienerFilterComment,
                                                        RECONSTRUCTION_USE_PRECOMPUTED_RSD_CORRECTION_AND_WIENER_FILTER);

        SphericalGridFunction fiducialRadialVelocity = constrainedRealizations.get_survey_estimate_radial_velocity(fiducialNormalizedGrowthRate, FIDUCIAL_FIELD_SMOOTHING_SCALE);

        const double fiducialOriginXVelocity = fiducialRadialVelocity(0.0, M_PI_2, 0.0);
        const double fiducialOriginYVelocity = fiducialRadialVelocity(0.0, M_PI_2, M_PI_2);
        const double fiducialOriginZVelocity = fiducialRadialVelocity(0.0, 0.0, 0.0);

        double localGroupXVelocity, localGroupYVelocity, localGroupZVelocity;

        transform_spherical_to_cartesian_coordinates(LOCAL_GROUP_TO_CMB.relativeVelocityAmplitude, LOCAL_GROUP_TO_CMB.relativeVelocityTheta, LOCAL_GROUP_TO_CMB.relativeVelocityPhi,
                                                     localGroupXVelocity, localGroupYVelocity, localGroupZVelocity);

        const double fiducialExternalBulkXVelocity = localGroupXVelocity - fiducialOriginXVelocity;
        const double fiducialExternalBulkYVelocity = localGroupYVelocity - fiducialOriginYVelocity;
        const double fiducialExternalBulkZVelocity = localGroupZVelocity - fiducialOriginZVelocity;

        double fiducialExternalBulkVelocityAmplitude, fiducialExternalBulkVelocityTheta, fiducialExternalBulkVelocityPhi;

        transform_cartesian_to_spherical_coordinates(fiducialExternalBulkXVelocity, fiducialExternalBulkYVelocity, fiducialExternalBulkZVelocity,
                                                     fiducialExternalBulkVelocityAmplitude, fiducialExternalBulkVelocityTheta, fiducialExternalBulkVelocityPhi);

        const ReferenceFrameChange fiducialExternalBulkVelocityShift = {fiducialExternalBulkVelocityAmplitude, fiducialExternalBulkVelocityTheta, fiducialExternalBulkVelocityPhi};

        fiducialRadialVelocity.apply([&](double velocity, double radius, double theta, double phi) {
            return change_reference_frame(velocity, theta, phi, fiducialExternalBulkVelocityShift);
        });

        SphericalGridFunction reconstructedRadialVelocityMinSmooth = constrainedRealizations.get_survey_estimate_radial_velocity(ESTIMATED_NORMALIZED_GROWTH_RATE, RECONSTRUCTED_FIELD_MIN_SMOOTHING_SCALE);

        double estimatedExternalBulkVelocityAmplitude, estimatedExternalBulkVelocityTheta, estimatedExternalBulkVelocityPhi;

        transform_cartesian_to_spherical_coordinates(ESTIMATED_EXTERNAL_BULK_X_VELOCITY, ESTIMATED_EXTERNAL_BULK_Y_VELOCITY, ESTIMATED_EXTERNAL_BULK_Z_VELOCITY,
                                                     estimatedExternalBulkVelocityAmplitude, estimatedExternalBulkVelocityTheta, estimatedExternalBulkVelocityPhi);

        const ReferenceFrameChange estimatedExternalBulkVelocityShift = {estimatedExternalBulkVelocityAmplitude, estimatedExternalBulkVelocityTheta, estimatedExternalBulkVelocityPhi};

        reconstructedRadialVelocityMinSmooth.apply([&](double velocity, double radius, double theta, double phi) {
            return change_reference_frame(velocity, theta, phi, estimatedExternalBulkVelocityShift);
        });

        std::vector<double> tensorSmoothingComparisonRedshiftVelocities, tensorSmoothingComparisonRadialCoordinates, tensorSmoothingComparisonThetaCoordinates, tensorSmoothingComparisonPhiCoordinates, tensorSmoothingComparisonDistanceModuli, tensorSmoothingComparisonDistanceModulusErrors;

        prepare_distance_catalog_data(inputDistanceCatalogRedshiftVelocities, inputDistanceCatalogLatitudes, inputDistanceCatalogLongitudes, inputDistanceCatalogDistanceModuli, inputDistanceCatalogDistanceModulusErrors,
                                      distanceCatalogInputToReferenceFrame, 0.0, DISTANCE_CATALOG_MAX_REDSHIFT_VELOCITY, FIDUCIAL_OMEGA_MATTER,
                                      fiducialRadialVelocity, cmbToReferenceFrame, FIDUCIAL_HUBBLE, DISTANCE_CATALOG_OUTLIER_MAX_SIGMA,
                                      tensorSmoothingComparisonRedshiftVelocities, tensorSmoothingComparisonRadialCoordinates, tensorSmoothingComparisonThetaCoordinates, tensorSmoothingComparisonPhiCoordinates, tensorSmoothingComparisonDistanceModuli, tensorSmoothingComparisonDistanceModulusErrors);

        const auto time2 = std::chrono::high_resolution_clock::now();
        std::cout << " done (" << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "s)" << std::endl;
        std::cout << "Comparing reconstructed and observed velocities... " << std::flush;

        for (std::size_t i_s = 0; i_s < TENSOR_SMOOTHING_SCALES.size(); ++i_s)
        {
            const double minTensorSmoothingScale = TENSOR_SMOOTHING_SCALES[i_s];

            const std::string smoothComment = "_rSmoothMin" + std::to_string(minTensorSmoothingScale);
            const std::string smoothComparisonComment = redshiftReferenceFrameComment + smoothComment + CONFIGURATION_COMMENT;

            const std::string velocityComparisonTensorSmoothedPointsFileName = DATA_DIRECTORY + "velocity_comparison_tensor_smoothed_points" + smoothComparisonComment + ".dat";
            const std::string velocityComparisonTensorSmoothedFieldsFileName = DATA_DIRECTORY + "velocity_comparison_tensor_smoothed_fields" + smoothComparisonComment + ".dat";

            std::vector<double> smoothedObservedRadialVelocities, smoothedReconstructedRadialVelocities, adaptiveSmoothingScales;

            compute_tensor_smoothed_radial_velocity_points(tensorSmoothingComparisonRedshiftVelocities, tensorSmoothingComparisonThetaCoordinates, tensorSmoothingComparisonPhiCoordinates, tensorSmoothingComparisonDistanceModuli, tensorSmoothingComparisonDistanceModulusErrors,
                                                           reconstructedRadialVelocityMinSmooth,
                                                           referenceToCMBFrame, FIDUCIAL_OMEGA_MATTER, ESTIMATED_DISTANCE_CATALOG_HUBBLE, minTensorSmoothingScale,
                                                           smoothedObservedRadialVelocities, smoothedReconstructedRadialVelocities, adaptiveSmoothingScales);

            std::ofstream velocityComparisonTensorSmoothedPointsFile(velocityComparisonTensorSmoothedPointsFileName);

            velocityComparisonTensorSmoothedPointsFile.setf(std::ios::fixed);

            velocityComparisonTensorSmoothedPointsFile << "# "
                                                       << "radius[Mpc/h]\t"
                                                       << "theta\t"
                                                       << "phi\t"
                                                       << "vObs[km/s]\t"
                                                       << "vRec[km/s]\t"
                                                       << "rSmooth[Mpc/h]"
                                                       << std::endl;

            for (std::size_t i_g = 0; i_g < tensorSmoothingComparisonRadialCoordinates.size(); ++i_g)
            {
                velocityComparisonTensorSmoothedPointsFile << std::setprecision(2)
                                                           << tensorSmoothingComparisonRadialCoordinates[i_g] << "\t"
                                                           << std::setprecision(6)
                                                           << tensorSmoothingComparisonThetaCoordinates[i_g] << "\t"
                                                           << tensorSmoothingComparisonPhiCoordinates[i_g] << "\t"
                                                           << std::setprecision(2)
                                                           << smoothedObservedRadialVelocities[i_g] << "\t"
                                                           << smoothedReconstructedRadialVelocities[i_g] << "\t"
                                                           << adaptiveSmoothingScales[i_g]
                                                           << std::endl;
            }

            velocityComparisonTensorSmoothedPointsFile.close();

            SphericalGridFunction smoothedObservedRadialVelocityField, smoothedReconstructedRadialVelocityField, adaptiveSmoothingScaleField;

            compute_tensor_smoothed_radial_velocity_fields(tensorSmoothingComparisonRedshiftVelocities, tensorSmoothingComparisonThetaCoordinates, tensorSmoothingComparisonPhiCoordinates, tensorSmoothingComparisonDistanceModuli, tensorSmoothingComparisonDistanceModulusErrors,
                                                           reconstructedRadialVelocityMinSmooth,
                                                           referenceToCMBFrame, FIDUCIAL_OMEGA_MATTER, ESTIMATED_DISTANCE_CATALOG_HUBBLE, minTensorSmoothingScale,
                                                           RECONSTRUCTION_MAX_RADIUS, TENSOR_SMOOTHING_RADIAL_BIN_NUMBER, TENSOR_SMOOTHING_THETA_BIN_NUMBER, TENSOR_SMOOTHING_PHI_BIN_NUMBER,
                                                           smoothedObservedRadialVelocityField, smoothedReconstructedRadialVelocityField, adaptiveSmoothingScaleField);

            std::ofstream velocityComparisonTensorSmoothedFields(velocityComparisonTensorSmoothedFieldsFileName);

            velocityComparisonTensorSmoothedFields.setf(std::ios::fixed);

            velocityComparisonTensorSmoothedFields << "# "
                                                   << "radius[Mpc/h]\t"
                                                   << "theta\t"
                                                   << "phi\t"
                                                   << "vObs[km/s]\t"
                                                   << "vRec[km/s]\t"
                                                   << "rSmooth[Mpc/h]"
                                                   << std::endl;

            smoothedObservedRadialVelocityField.evaluate_for_all_grid_points([&](std::size_t i_r, std::size_t i_t, std::size_t i_p) {
                const double radius = smoothedObservedRadialVelocityField.radial_coordinate(i_r);
                const double theta = smoothedObservedRadialVelocityField.theta_coordinate(i_t);
                const double phi = smoothedObservedRadialVelocityField.phi_coordinate(i_p);

                const double observedVelocity = smoothedObservedRadialVelocityField.value(i_r, i_t, i_p);
                const double reconstructedVelocity = smoothedReconstructedRadialVelocityField.value(i_r, i_t, i_p);
                const double adaptiveSmoothingScale = adaptiveSmoothingScaleField.value(i_r, i_t, i_p);

                velocityComparisonTensorSmoothedFields << std::setprecision(2)
                                                       << radius << "\t"
                                                       << std::setprecision(6)
                                                       << theta << "\t"
                                                       << phi << "\t"
                                                       << std::setprecision(2)
                                                       << observedVelocity << "\t"
                                                       << reconstructedVelocity << "\t"
                                                       << adaptiveSmoothingScale
                                                       << std::endl;
            });

            velocityComparisonTensorSmoothedFields.close();
        }

        std::vector<double> correlationFunctionComparisonRedshiftVelocities, correlationFunctionComparisonRadialCoordinates, correlationFunctionComparisonThetaCoordinates, correlationFunctionComparisonPhiCoordinates, correlationFunctionComparisonDistanceModuli, correlationFunctionComparisonDistanceModulusErrors;

        prepare_distance_catalog_data(inputDistanceCatalogRedshiftVelocities, inputDistanceCatalogLatitudes, inputDistanceCatalogLongitudes, inputDistanceCatalogDistanceModuli, inputDistanceCatalogDistanceModulusErrors,
                                      distanceCatalogInputToReferenceFrame, 0.0, CORRELATION_FUNCTION_MAX_REDSHIFT_VELOCITY, FIDUCIAL_OMEGA_MATTER,
                                      fiducialRadialVelocity, cmbToReferenceFrame, FIDUCIAL_HUBBLE, DISTANCE_CATALOG_OUTLIER_MAX_SIGMA,
                                      correlationFunctionComparisonRedshiftVelocities, correlationFunctionComparisonRadialCoordinates, correlationFunctionComparisonThetaCoordinates, correlationFunctionComparisonPhiCoordinates, correlationFunctionComparisonDistanceModuli, correlationFunctionComparisonDistanceModulusErrors);

        const std::string correlationComparisonComment = redshiftReferenceFrameComment + CONFIGURATION_COMMENT;

        const std::string velocityComparisonCorrelationFunctionsFileName = DATA_DIRECTORY + "velocity_comparison_correlation_functions" + correlationComparisonComment + ".dat";

        Cartesian1DGridFunction observedRadialVelocityCorrelationFunction, reconstructedRadialVelocityCorrelationFunction, residualRadialVelocityCorrelationFunction;

        compute_radial_velocity_correlation_functions(correlationFunctionComparisonRedshiftVelocities, correlationFunctionComparisonThetaCoordinates, correlationFunctionComparisonPhiCoordinates, correlationFunctionComparisonDistanceModuli,
                                                      reconstructedRadialVelocityMinSmooth,
                                                      referenceToCMBFrame, FIDUCIAL_OMEGA_MATTER, ESTIMATED_DISTANCE_CATALOG_HUBBLE,
                                                      CORRELATION_FUNCTION_MAX_DISTANCE, CORRELATION_FUNCTION_DISTANCE_BIN_NUMBER,
                                                      observedRadialVelocityCorrelationFunction, reconstructedRadialVelocityCorrelationFunction, residualRadialVelocityCorrelationFunction);

        std::ofstream velocityComparisonCorrelationFunctionsFile(velocityComparisonCorrelationFunctionsFileName);

        velocityComparisonCorrelationFunctionsFile.setf(std::ios::fixed);

        velocityComparisonCorrelationFunctionsFile << "# "
                                                   << "dist[Mpc/h]\t"
                                                   << "psiObs[km^2/s^2]\t"
                                                   << "psiRec[km^2/s^2]\t"
                                                   << "psiRes[km^2/s^2]"
                                                   << std::endl;

        for (std::size_t i_g = 0; i_g < observedRadialVelocityCorrelationFunction.bin_number(); ++i_g)
        {
            velocityComparisonCorrelationFunctionsFile << std::setprecision(2)
                                                       << observedRadialVelocityCorrelationFunction.coordinate(i_g) << "\t"
                                                       << observedRadialVelocityCorrelationFunction.value(i_g) << "\t"
                                                       << reconstructedRadialVelocityCorrelationFunction.value(i_g) << "\t"
                                                       << residualRadialVelocityCorrelationFunction.value(i_g)
                                                       << std::endl;
        }

        velocityComparisonCorrelationFunctionsFile.close();

        const auto time3 = std::chrono::high_resolution_clock::now();
        std::cout << " done (" << std::chrono::duration_cast<std::chrono::seconds>(time3 - time2).count() << "s)" << std::endl;
    }

    fftw_cleanup_threads();

    return 0;
}