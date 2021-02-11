#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <omp.h>

#include "../src/ConstrainedRealizations.hpp"
#include "../src/FileTable.hpp"
#include "../src/NormalizedPowerSpectrum.hpp"
#include "../src/SphericalGridFunction.hpp"
#include "../src/SphericalGridFunctionInterpolator.hpp"
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
    const auto inputDistanceCatalogRedshiftVelocitities = inputDistanceCatalogFile.column<double>(4);
    const auto inputDistanceCatalogDistanceModuli = inputDistanceCatalogFile.column<double>(5);
    const auto inputDistanceCatalogDistanceModulusErrors = inputDistanceCatalogFile.column<double>(6);

    for (std::size_t i_ref = 0; i_ref < RECONSTRUCTION_REDSHIFT_REFERENCE_FRAMES.size(); ++i_ref)
    {
        const ReferenceFrame redshiftReferenceFrame = RECONSTRUCTION_REDSHIFT_REFERENCE_FRAMES[i_ref];

        const ReferenceFrameChange redshiftCatalogInputToReferenceFrame = get_reference_frame_change(REDSHIFT_CATALOG_REDSHIFT_REFERENCE_FRAME, redshiftReferenceFrame);
        const ReferenceFrameChange distanceCatalogInputToReferenceFrame = get_reference_frame_change(DISTANCE_CATALOG_REDSHIFT_REFERENCE_FRAME, redshiftReferenceFrame);
        const ReferenceFrameChange cmbToReferenceFrame = get_reference_frame_change(CMB_FRAME, redshiftReferenceFrame);

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

        std::vector<double> reconstructionRedshiftVelocities, reconstructionRadialCoordinates, reconstructionThetaCoordinates, reconstructionPhiCoordinates, reconstructionApparentMagnitudes;

        prepare_2MRS_data(inputRedshiftCatalogRedshiftVelocities, inputRedshiftCatalogLatitudes, inputRedshiftCatalogLongitudes, inputRedshiftCatalogApparentMagnitudes,
                          redshiftCatalogInputToReferenceFrame, RECONSTRUCTION_MAX_RADIUS, REDSHIFT_CATALOG_VOLUME_LIMIT_RADIUS, REDSHIFT_CATALOG_MAX_APPARENT_MAGNITUDE, FIDUCIAL_OMEGA_MATTER,
                          reconstructionRedshiftVelocities, reconstructionRadialCoordinates, reconstructionThetaCoordinates, reconstructionPhiCoordinates, reconstructionApparentMagnitudes);

        Cartesian1DGridFunction sigmaGalaxyAboveVolumeLimit;
        std::vector<std::size_t> volumeLimitGalaxyNumbers;

        estimate_sigma_galaxy(reconstructionRedshiftVelocities, reconstructionThetaCoordinates, reconstructionPhiCoordinates, reconstructionApparentMagnitudes,
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

        estimate_selection_function(reconstructionRedshiftVelocities, reconstructionApparentMagnitudes, REDSHIFT_CATALOG_MAX_APPARENT_MAGNITUDE,
                                    RECONSTRUCTION_MAX_RADIUS, SELECTION_FUNCTION_BIN_NUMBER, FIDUCIAL_OMEGA_MATTER,
                                    selectionFunction, selectionFunctionLogDerivative, radialDistributionFunction);

        const double meanGalaxyDensity = estimate_mean_density(RECONSTRUCTION_MAX_RADIUS, reconstructionRadialCoordinates, selectionFunction);

        const std::string redshiftReferenceFrameComment = "_z" + redshiftReferenceFrameName;
        const std::string precomputedRSDCorrectionAndWienerFilterComment = redshiftReferenceFrameComment + CONFIGURATION_COMMENT;

        ConstrainedRealizations constrainedRealizations(reconstructionRadialCoordinates, reconstructionThetaCoordinates, reconstructionPhiCoordinates, redshiftReferenceFrame,
                                                        RECONSTRUCTION_MAX_RADIUS, RECONSTRUCTION_RADIAL_RESOLUTION, RECONSTRUCTION_MAX_MULTIPOLE, RECONSTRUCTION_RADIAL_BIN_NUMBER, RECONSTRUCTION_THETA_BIN_NUMBER, RECONSTRUCTION_PHI_BIN_NUMBER,
                                                        RECONSTRUCTION_FFT_PERIODIC_BOUNDARY_DISTANCE, RECONSTRUCTION_FFT_BIN_NUMBER, RECONSTRUCTION_LOG_TRAFO_SMOOTHING_SCALE,
                                                        meanGalaxyDensity, selectionFunction, selectionFunctionLogDerivative, sigmaGalaxy, normalizedPowerSpectrum,
                                                        DATA_DIRECTORY, precomputedRSDCorrectionAndWienerFilterComment,
                                                        RECONSTRUCTION_USE_PRECOMPUTED_RSD_CORRECTION_AND_WIENER_FILTER);

        const double normalizedGrowthRateBinWidth = (PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_MAX - PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_MIN) / static_cast<double>(PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_BIN_NUMBER - 1);

        std::vector<std::vector<SphericalGridFunction>> normalizedSurveyEstimateRadialVelocities(PARAMETER_ESTIMATE_SMOOTHING_SCALES.size(), std::vector<SphericalGridFunction>(PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_BIN_NUMBER));

        for (std::size_t i_growth = 0; i_growth < PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_BIN_NUMBER; ++i_growth)
        {
            const double normalizedGrowthRate = (i_growth == PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_BIN_NUMBER - 1) ? PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_MAX
                                                                                                                       : PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_MIN + static_cast<double>(i_growth) * normalizedGrowthRateBinWidth;

            for (std::size_t i_rSmooth = 0; i_rSmooth < PARAMETER_ESTIMATE_SMOOTHING_SCALES.size(); ++i_rSmooth)
            {
                const double smoothingScale = PARAMETER_ESTIMATE_SMOOTHING_SCALES[i_rSmooth];

                normalizedSurveyEstimateRadialVelocities[i_rSmooth][i_growth] = constrainedRealizations.get_survey_estimate_radial_velocity(normalizedGrowthRate, smoothingScale) / normalizedGrowthRate;
            }
        }

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

        const auto time2 = std::chrono::high_resolution_clock::now();
        std::cout << " done (" << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "s)" << std::endl;

        const std::size_t cutNumber = PARAMETER_ESTIMATE_SMOOTHING_SCALES.size() * PARAMETER_ESTIMATE_MIN_REDSHIFT_VELOCITIES.size();

        std::vector<double> averageNormalizedGrowthRateEstimates(cutNumber, 0.0);
        std::vector<double> averageNormalizedGrowthRateCRVariances(cutNumber, 0.0);
        std::vector<double> averageNormalizedGrowthRateMLVariances(cutNumber, 0.0);
        std::vector<double> averageExternalBulkXVelocityEstimates(cutNumber, 0.0);
        std::vector<double> averageExternalBulkXVelocityCRVariances(cutNumber, 0.0);
        std::vector<double> averageExternalBulkXVelocityMLVariances(cutNumber, 0.0);
        std::vector<double> averageExternalBulkYVelocityEstimates(cutNumber, 0.0);
        std::vector<double> averageExternalBulkYVelocityCRVariances(cutNumber, 0.0);
        std::vector<double> averageExternalBulkYVelocityMLVariances(cutNumber, 0.0);
        std::vector<double> averageExternalBulkZVelocityEstimates(cutNumber, 0.0);
        std::vector<double> averageExternalBulkZVelocityCRVariances(cutNumber, 0.0);
        std::vector<double> averageExternalBulkZVelocityMLVariances(cutNumber, 0.0);
        std::vector<double> averageHubbleEstimates(cutNumber, 0.0);
        std::vector<double> averageHubbleCRVariances(cutNumber, 0.0);
        std::vector<double> averageHubbleMLVariances(cutNumber, 0.0);
        std::vector<double> averageReducedChiSquares(cutNumber, 0.0);

        for (std::size_t i_CR = INITIAL_CONSTRAINED_REALIZATION; i_CR <= FINAL_CONSTRAINED_REALIZATION; ++i_CR)
        {
            const auto timeCR1 = std::chrono::high_resolution_clock::now();
            std::cout << "Analyzing CR " << i_CR << "..." << std::flush;

            constrainedRealizations.generate(i_CR);

            const std::string realizationComment = "_CR" + std::to_string(i_CR);
            const std::string individualComment = redshiftReferenceFrameComment + realizationComment + CONFIGURATION_COMMENT;

            const std::string parameterFitIndividualFileName = DATA_DIRECTORY + "parameters" + individualComment + ".dat";

            std::ofstream parameterFitIndividualFile = std::ofstream(parameterFitIndividualFileName);

            parameterFitIndividualFile.setf(std::ios::fixed);

            parameterFitIndividualFile << "# "
                                       << "rSmooth[Mpc/h]\t"
                                       << "czMin[km/s]\t"
                                       << "fsig8\t"
                                       << "fsig8_e\t"
                                       << "vExtX[km/s]\t"
                                       << "vExtX_e[km/s]\t"
                                       << "vExtY[km/s]\t"
                                       << "vExtY_e[km/s]\t"
                                       << "vExtZ[km/s]\t"
                                       << "vExtZ_e[km/s]\t"
                                       << "h\t"
                                       << "h_e\t"
                                       << "chiSq/dof"
                                       << std::endl;

            for (std::size_t i_rSmooth = 0; i_rSmooth < PARAMETER_ESTIMATE_SMOOTHING_SCALES.size(); ++i_rSmooth)
            {
                const double smoothingScale = PARAMETER_ESTIMATE_SMOOTHING_SCALES[i_rSmooth];

                SphericalGridFunction normalizedRandomNoiseRadialVelocity = constrainedRealizations.get_noise_radial_velocity(1.0, smoothingScale);

                SphericalGridFunctionInterpolator normalizedSurveyEstimateRadialVelocityInterpolator(PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_MIN, PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_MAX, normalizedSurveyEstimateRadialVelocities[i_rSmooth]);

                for (std::size_t i_czMin = 0; i_czMin < PARAMETER_ESTIMATE_MIN_REDSHIFT_VELOCITIES.size(); ++i_czMin)
                {
                    const double minRedshiftVelocity = PARAMETER_ESTIMATE_MIN_REDSHIFT_VELOCITIES[i_czMin];

                    const std::size_t i_cut = i_rSmooth * PARAMETER_ESTIMATE_MIN_REDSHIFT_VELOCITIES.size() + i_czMin;

                    std::vector<double> comparisonRedshiftVelocities, comparisonRadialCoordinates, comparisonThetaCoordinates, comparisonPhiCoordinates, comparisonDistanceModuli, comparisonDistanceModulusErrors;

                    prepare_distance_catalog_data(inputDistanceCatalogRedshiftVelocitities, inputDistanceCatalogLatitudes, inputDistanceCatalogLongitudes, inputDistanceCatalogDistanceModuli, inputDistanceCatalogDistanceModulusErrors,
                                                  distanceCatalogInputToReferenceFrame, minRedshiftVelocity, DISTANCE_CATALOG_MAX_REDSHIFT_VELOCITY, FIDUCIAL_OMEGA_MATTER,
                                                  fiducialRadialVelocity, cmbToReferenceFrame, FIDUCIAL_HUBBLE, DISTANCE_CATALOG_OUTLIER_MAX_SIGMA,
                                                  comparisonRedshiftVelocities, comparisonRadialCoordinates, comparisonThetaCoordinates, comparisonPhiCoordinates, comparisonDistanceModuli, comparisonDistanceModulusErrors);

                    double normalizedGrowthRateEstimate, externalBulkXVelocityEstimate, externalBulkYVelocityEstimate, externalBulkZVelocityEstimate, reducedChiSquare, hubbleEstimate;

                    gsl_matrix *covarianceMatrix = gsl_matrix_alloc(5, 5);

                    estimate_parameters_via_radial_velocity_comparison(normalizedSurveyEstimateRadialVelocityInterpolator, normalizedRandomNoiseRadialVelocity, cmbToReferenceFrame,
                                                                       comparisonRedshiftVelocities, comparisonThetaCoordinates, comparisonPhiCoordinates, comparisonDistanceModuli, comparisonDistanceModulusErrors,
                                                                       FIDUCIAL_OMEGA_MATTER, fiducialNormalizedGrowthRate, PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_MIN, PARAMETER_ESTIMATE_NORMALIZED_GROWTH_RATE_MAX, FIDUCIAL_HUBBLE,
                                                                       normalizedGrowthRateEstimate, externalBulkXVelocityEstimate, externalBulkYVelocityEstimate, externalBulkZVelocityEstimate, hubbleEstimate,
                                                                       covarianceMatrix, reducedChiSquare);

                    const double normalizedGrowthRateError = std::sqrt(gsl_matrix_get(covarianceMatrix, 0, 0));
                    const double externalBulkXVelocityError = std::sqrt(gsl_matrix_get(covarianceMatrix, 1, 1));
                    const double externalBulkYVelocityError = std::sqrt(gsl_matrix_get(covarianceMatrix, 2, 2));
                    const double externalBulkZVelocityError = std::sqrt(gsl_matrix_get(covarianceMatrix, 3, 3));
                    const double hubbleError = std::sqrt(gsl_matrix_get(covarianceMatrix, 4, 4));

                    gsl_matrix_free(covarianceMatrix);

                    averageNormalizedGrowthRateEstimates[i_cut] += normalizedGrowthRateEstimate;
                    averageNormalizedGrowthRateCRVariances[i_cut] += gsl_pow_2(normalizedGrowthRateEstimate);
                    averageNormalizedGrowthRateMLVariances[i_cut] += gsl_pow_2(normalizedGrowthRateError);
                    averageExternalBulkXVelocityEstimates[i_cut] += externalBulkXVelocityEstimate;
                    averageExternalBulkXVelocityCRVariances[i_cut] += gsl_pow_2(externalBulkXVelocityEstimate);
                    averageExternalBulkXVelocityMLVariances[i_cut] += gsl_pow_2(externalBulkXVelocityError);
                    averageExternalBulkYVelocityEstimates[i_cut] += externalBulkYVelocityEstimate;
                    averageExternalBulkYVelocityCRVariances[i_cut] += gsl_pow_2(externalBulkYVelocityEstimate);
                    averageExternalBulkYVelocityMLVariances[i_cut] += gsl_pow_2(externalBulkYVelocityError);
                    averageExternalBulkZVelocityEstimates[i_cut] += externalBulkZVelocityEstimate;
                    averageExternalBulkZVelocityCRVariances[i_cut] += gsl_pow_2(externalBulkZVelocityEstimate);
                    averageExternalBulkZVelocityMLVariances[i_cut] += gsl_pow_2(externalBulkZVelocityError);
                    averageHubbleEstimates[i_cut] += hubbleEstimate;
                    averageHubbleCRVariances[i_cut] += gsl_pow_2(hubbleEstimate);
                    averageHubbleMLVariances[i_cut] += gsl_pow_2(hubbleError);
                    averageReducedChiSquares[i_cut] += reducedChiSquare;

                    parameterFitIndividualFile << std::setprecision(1)
                                               << smoothingScale << "\t"
                                               << minRedshiftVelocity << "\t"
                                               << std::setprecision(5)
                                               << normalizedGrowthRateEstimate << "\t"
                                               << normalizedGrowthRateError << "\t"
                                               << std::setprecision(2)
                                               << externalBulkXVelocityEstimate << "\t"
                                               << externalBulkXVelocityError << "\t"
                                               << externalBulkYVelocityEstimate << "\t"
                                               << externalBulkYVelocityError << "\t"
                                               << externalBulkZVelocityEstimate << "\t"
                                               << externalBulkZVelocityError << "\t"
                                               << std::setprecision(5)
                                               << hubbleEstimate << "\t"
                                               << hubbleError << "\t"
                                               << std::setprecision(3)
                                               << reducedChiSquare
                                               << std::endl;
                }
            }

            parameterFitIndividualFile.close();

            const auto timeCR2 = std::chrono::high_resolution_clock::now();
            std::cout << " done (" << std::chrono::duration_cast<std::chrono::seconds>(timeCR2 - timeCR1).count() << "s)" << std::endl;
        }

        const std::string realizationRangeComment = "_CR" + std::to_string(INITIAL_CONSTRAINED_REALIZATION) + "-" + std::to_string(FINAL_CONSTRAINED_REALIZATION);
        const std::string averageComment = redshiftReferenceFrameComment + realizationRangeComment + CONFIGURATION_COMMENT;

        const std::string parameterFitAverageFileName = DATA_DIRECTORY + "parameters" + averageComment + ".dat";

        std::ofstream parameterFitAverageFile = std::ofstream(parameterFitAverageFileName);

        parameterFitAverageFile << "# "
                                << "rSmooth[Mpc/h]\t"
                                << "czMin[km/s]\t"
                                << "fsig8\t"
                                << "fsig8_eShot\t"
                                << "fsig8_eDist\t"
                                << "vExtX[km/s]\t"
                                << "vExtX_eShot[km/s]\t"
                                << "vExtX_eDist[km/s]\t"
                                << "vExtY[km/s]\t"
                                << "vExtY_eShot[km/s]\t"
                                << "vExtY_eDist[km/s]\t"
                                << "vExtZ[km/s]\t"
                                << "vExtZ_eShot[km/s]\t"
                                << "vExtZ_eDist[km/s]\t"
                                << "h\t"
                                << "h_eShot\t"
                                << "h_eDist\t"
                                << "chiSq/dof"
                                << std::endl;

        parameterFitAverageFile.setf(std::ios::fixed);

        for (std::size_t i_rSmooth = 0; i_rSmooth < PARAMETER_ESTIMATE_SMOOTHING_SCALES.size(); ++i_rSmooth)
        {
            const double smoothingScale = PARAMETER_ESTIMATE_SMOOTHING_SCALES[i_rSmooth];

            for (std::size_t i_czMin = 0; i_czMin < PARAMETER_ESTIMATE_MIN_REDSHIFT_VELOCITIES.size(); ++i_czMin)
            {
                const double minRedshiftVelocity = PARAMETER_ESTIMATE_MIN_REDSHIFT_VELOCITIES[i_czMin];

                const std::size_t i_cut = i_rSmooth * PARAMETER_ESTIMATE_MIN_REDSHIFT_VELOCITIES.size() + i_czMin;

                const double constrainedRealizationNumber = static_cast<double>(FINAL_CONSTRAINED_REALIZATION - INITIAL_CONSTRAINED_REALIZATION + 1);

                const double averageNormalizedGrowthRateEstimate = averageNormalizedGrowthRateEstimates[i_cut] / constrainedRealizationNumber;
                const double averageNormalizedGrowthRateCRError = std::sqrt(averageNormalizedGrowthRateCRVariances[i_cut] / constrainedRealizationNumber - gsl_pow_2(averageNormalizedGrowthRateEstimate));
                const double averageNormalizedGrowthRateMLError = std::sqrt(averageNormalizedGrowthRateMLVariances[i_cut] / constrainedRealizationNumber);
                const double averageExternalBulkXVelocityEstimate = averageExternalBulkXVelocityEstimates[i_cut] / constrainedRealizationNumber;
                const double averageExternalBulkXVelocityCRError = std::sqrt(averageExternalBulkXVelocityCRVariances[i_cut] / constrainedRealizationNumber - gsl_pow_2(averageExternalBulkXVelocityEstimate));
                const double averageExternalBulkXVelocityMLError = std::sqrt(averageExternalBulkXVelocityMLVariances[i_cut] / constrainedRealizationNumber);
                const double averageExternalBulkYVelocityEstimate = averageExternalBulkYVelocityEstimates[i_cut] / constrainedRealizationNumber;
                const double averageExternalBulkYVelocityCRError = std::sqrt(averageExternalBulkYVelocityCRVariances[i_cut] / constrainedRealizationNumber - gsl_pow_2(averageExternalBulkYVelocityEstimate));
                const double averageExternalBulkYVelocityMLError = std::sqrt(averageExternalBulkYVelocityMLVariances[i_cut] / constrainedRealizationNumber);
                const double averageExternalBulkZVelocityEstimate = averageExternalBulkZVelocityEstimates[i_cut] / constrainedRealizationNumber;
                const double averageExternalBulkZVelocityCRError = std::sqrt(averageExternalBulkZVelocityCRVariances[i_cut] / constrainedRealizationNumber - gsl_pow_2(averageExternalBulkZVelocityEstimate));
                const double averageExternalBulkZVelocityMLError = std::sqrt(averageExternalBulkZVelocityMLVariances[i_cut] / constrainedRealizationNumber);
                const double averageHubbleEstimate = averageHubbleEstimates[i_cut] / constrainedRealizationNumber;
                const double averageHubbleCRError = std::sqrt(averageHubbleCRVariances[i_cut] / constrainedRealizationNumber - gsl_pow_2(averageHubbleEstimate));
                const double averageHubbleMLError = std::sqrt(averageHubbleMLVariances[i_cut] / constrainedRealizationNumber);
                const double averageReducedChiSquare = averageReducedChiSquares[i_cut] / constrainedRealizationNumber;

                parameterFitAverageFile << std::setprecision(1)
                                        << smoothingScale << "\t"
                                        << minRedshiftVelocity << "\t"
                                        << std::setprecision(5)
                                        << averageNormalizedGrowthRateEstimate << "\t"
                                        << averageNormalizedGrowthRateCRError << "\t"
                                        << averageNormalizedGrowthRateMLError << "\t"
                                        << std::setprecision(2)
                                        << averageExternalBulkXVelocityEstimate << "\t"
                                        << averageExternalBulkXVelocityCRError << "\t"
                                        << averageExternalBulkXVelocityMLError << "\t"
                                        << averageExternalBulkYVelocityEstimate << "\t"
                                        << averageExternalBulkYVelocityCRError << "\t"
                                        << averageExternalBulkYVelocityMLError << "\t"
                                        << averageExternalBulkZVelocityEstimate << "\t"
                                        << averageExternalBulkZVelocityCRError << "\t"
                                        << averageExternalBulkZVelocityMLError << "\t"
                                        << std::setprecision(5)
                                        << averageHubbleEstimate << "\t"
                                        << averageHubbleCRError << "\t"
                                        << averageHubbleMLError << "\t"
                                        << std::setprecision(3)
                                        << averageReducedChiSquare
                                        << std::endl;
            }
        }

        parameterFitAverageFile.close();
    }

    fftw_cleanup_threads();

    return 0;
}