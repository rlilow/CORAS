#ifndef CORAS_MDPL2_H
#define CORAS_MDPL2_H

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "Cartesian3DGridFunction.hpp"
#include "transformations.hpp"

constexpr double BoxLength = 1000.0;
constexpr std::size_t BoxBinNumber = 512;
constexpr std::size_t GridBinNumber = BoxBinNumber / 2;
constexpr double BoxBinWidth = BoxLength / static_cast<double>(BoxBinNumber);
constexpr double GridBinWidth = BoxLength / static_cast<double>(GridBinNumber);
constexpr double HubbleSim = 0.6777;
constexpr double OmegaMatterSim = 0.307115;
constexpr double OmegaBaryonSim = 0.048206;
constexpr double OmegaLambdaSim = 0.692885;
constexpr double LinearSigma8Sim = 0.8228;
constexpr double SpectralIndexSim = 0.96;
const double GrowthRateSim = std::pow(OmegaMatterSim, 0.55);
const double LinearNormalizedGrowthRateSim = GrowthRateSim * LinearSigma8Sim;

double luminosity_evolution_correction_MDPL2_SAGE(double redshift);

std::size_t read_MDPL2_galaxy_data(const std::string &fileName,
                                   std::vector<std::size_t> &galaxy_rowId,
                                   std::vector<std::size_t> &galaxy_rockstarId,
                                   std::vector<std::size_t> &galaxy_hostHaloId,
                                   std::vector<std::size_t> &galaxy_mainHaloId,
                                   std::vector<std::size_t> &galaxy_galaxyType,
                                   std::vector<double> &galaxy_haloMass,
                                   std::vector<double> &galaxy_maxVelocity,
                                   std::vector<double> &galaxy_peakVelocity,
                                   std::vector<double> &galaxy_xCoordinate,
                                   std::vector<double> &galaxy_yCoordinate,
                                   std::vector<double> &galaxy_zCoordinate,
                                   std::vector<double> &galaxy_xVelocity,
                                   std::vector<double> &galaxy_yVelocity,
                                   std::vector<double> &galaxy_zVelocity,
                                   std::vector<double> &galaxy_stellarMass,
                                   std::vector<double> &galaxy_meanStarAge);

std::size_t read_MDPL2_galaxy_sage_data(const std::string &fileName,
                                        std::vector<std::size_t> &galaxy_rowId,
                                        std::vector<std::size_t> &galaxy_rockstarId,
                                        std::vector<std::size_t> &galaxy_hostHaloId,
                                        std::vector<std::size_t> &galaxy_mainHaloId,
                                        std::vector<std::size_t> &galaxy_galaxyType,
                                        std::vector<double> &galaxy_haloMass,
                                        std::vector<double> &galaxy_maxVelocity,
                                        std::vector<double> &galaxy_xCoordinate,
                                        std::vector<double> &galaxy_yCoordinate,
                                        std::vector<double> &galaxy_zCoordinate,
                                        std::vector<double> &galaxy_xVelocity,
                                        std::vector<double> &galaxy_yVelocity,
                                        std::vector<double> &galaxy_zVelocity,
                                        std::vector<double> &galaxy_stellarMass,
                                        std::vector<double> &galaxy_meanStarAge);

std::size_t read_MDPL2_galaxy_sage_data_new(const std::string &fileName,
                                            std::vector<std::size_t> &galaxy_galaxyId,
                                            std::vector<std::size_t> &galaxy_rockstarId,
                                            std::vector<std::size_t> &galaxy_mainHaloId,
                                            std::vector<double> &galaxy_haloMass,
                                            std::vector<double> &galaxy_xCoordinate,
                                            std::vector<double> &galaxy_yCoordinate,
                                            std::vector<double> &galaxy_zCoordinate,
                                            std::vector<double> &galaxy_xVelocity,
                                            std::vector<double> &galaxy_yVelocity,
                                            std::vector<double> &galaxy_zVelocity,
                                            std::vector<double> &galaxy_stellarMass);

std::size_t read_MDPL2_halo_data(const std::string &fileName,
                                 std::vector<std::size_t> &halo_rowId,
                                 std::vector<std::size_t> &halo_rockstarId,
                                 std::vector<std::intmax_t> &halo_pId,
                                 std::vector<std::intmax_t> &halo_upId,
                                 std::vector<double> &halo_virialMass,
                                 std::vector<double> &halo_virialRadius,
                                 std::vector<double> &halo_scaleRadius,
                                 std::vector<double> &halo_rmsVelocity,
                                 std::vector<double> &halo_maxVelocity,
                                 std::vector<double> &halo_xCoordinate,
                                 std::vector<double> &halo_yCoordinate,
                                 std::vector<double> &halo_zCoordinate,
                                 std::vector<double> &halo_xVelocity,
                                 std::vector<double> &halo_yVelocity,
                                 std::vector<double> &halo_zVelocity,
                                 std::vector<double> &halo_200cMass,
                                 std::vector<double> &halo_500cMass,
                                 std::vector<double> &halo_accretionMass,
                                 std::vector<double> &halo_peakMass,
                                 std::vector<double> &halo_accretionVelocity,
                                 std::vector<double> &halo_peakVelocity,
                                 std::vector<double> &halo_peakMassMaxVelocity);

void ignore_fortran_unformatted_boundary_field(std::ifstream &file);

template <typename T>
void read_fortran_unformatted_value(std::ifstream &file, T &variable);

Cartesian3DGridFunction read_MDPL2_field(const std::string &fileName);
Cartesian3DGridFunction read_MDPL2_field(const std::string &fileName, double smoothingScale, bool useTophatKernel = false);
Cartesian3DGridFunction read_MDPL2_velocity_divergence_field(const std::string &xVelocityFileName, const std::string &yVelocityFileName, const std::string &zVelocityFileName,
                                                             double smoothingScale, bool useTophatKernel = false);

void read_MDPL2_density_field_and_get_linear_velocity(const std::string &fileName, double smoothingScale, double beta,
                                                      Cartesian3DGridFunction &densityContrast, Cartesian3DGridFunction &xVelocity, Cartesian3DGridFunction &yVelocity, Cartesian3DGridFunction &zVelocity,
                                                      bool useTophatKernel = false);

CoordinateSystemRotation MDPL2_mock_galactic_coordinate_system_change(double thetaMockVirgo, double phiMockVirgo, double thetaMockLGVelocity, double phiMockLGVelocity);

double evaluate_within_sphere_in_periodic_box(double xCenter, double yCenter, double zCenter, double boxLength,
                                              double radius, double theta, double phi,
                                              const Cartesian3DGridFunction &func);
//   const std::function<double(double x, double y, double z)> &func);

#endif