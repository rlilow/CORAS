#ifndef CORAS_SPHERICAL_FOURIER_BESSEL_DECOMPOSITION_H
#define CORAS_SPHERICAL_FOURIER_BESSEL_DECOMPOSITION_H

#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "SphericalGridFunction.hpp"

/**
 * \addtogroup RECONSTRUCTION
 * @{
 */

/**
 * \brief Class implementing the spherical Fourier Bessel (SFB) decomposition of some scalar field.
 *
 * It provides methods to compute the radial wavenumbers and normalization constants, to access the SFB coefficients,
 * and to evaluate the decomposed field, its inverse divergence as well as its inverse laplacian. Furthermore, the
 * explicit computation of the SFB coefficients of the density of a given set of point sources is implemented.
 */
class SphericalFourierBesselDecomposition
{
public:
	/**
	 * Constructor instantiating the SFB decomposition of a zero-valued field for a maximal radius \a maxRadius, radial
	 * resolution \a radialResolution (product of maximal radius and maximal considered radial wavenumber) and maximal
	 * multipole order \a maxMultipole. This means that all SFB coefficients are set to zero.
	 */
	SphericalFourierBesselDecomposition(double maxRadius, double radialResolution, std::size_t maxMultipole);

	/**
	 * Constructor instantiating the SFB decomposition of the density of a set of point sources at the radial and
	 * angular coordinates \a  radialCoordinates, \a thetaCoordinates and \a phiCoordinates, for a maximal radius \a
	 * maxRadius, \a radialResolution (product of maximal radius and maximal considered radial wavenumber) and maximal
	 * multipole order \a maxMultipole. The individual particle-contributions are weighted with the radially dependent
	 * function \a weightingFunction (default: 1), and the radially dependent function \a isotropicContribution
	 * (default: 0) is subtracted after the weighting, e.g. to set the mean field to zero.
	 */
	SphericalFourierBesselDecomposition(
		const std::vector<double> &radialCoordinates, const std::vector<double> &thetaCoordinates, const std::vector<double> &phiCoordinates,
		double maxRadius, double radialResolution, std::size_t maxMultipole,
		const std::function<double(double radius)> &weightingFunction = [](double) { return 1.0; }, const std::function<double(double radius)> &isotropicContribution = [](double) { return 0.0; });

	/**
	 * Constructor instantiating a previously saved SphericalFourierBesselDecomposition from the binary file \a
	 * fileName.
	 *
	 * The file must have been created using SphericalFourierBesselDecomposition::save_object_to_file.
	 */
	SphericalFourierBesselDecomposition(const std::string &fileName);

	/**
	 * Copy-constructor instantiating the same SFB decomposition as \a otherSphericalFourierBesselDecomposition.
	 */
	SphericalFourierBesselDecomposition(const SphericalFourierBesselDecomposition &otherSphericalFourierBesselDecomposition) = default;

	/**
	 * Constructor instantiating an empty SphericalFourierBesselDecomposition object, intended to be properly
	 * initialized at a later point. For example, this can be used to choose a constructor based on some conditional
	 * statement.
	 */
	SphericalFourierBesselDecomposition();

	/**
	 * Return the value of the radial wavenumber k_ln corresponding to the angular and radial modes \a  l and \a n.
	 */
	double radial_wavenumber(std::size_t l, std::size_t n) const;

	/**
	 * Return the value of the radial normalization constant C_ln corresponding to the angular and radial modes \a  l
	 * and \a n.
	 */
	double radial_normalization_constant(std::size_t l, std::size_t n) const;

	/**
	 * Return a reference to the complex value of the SFB coefficient of the decomposed field at the angular and radial
	 * mode numbers \a  l, \a m and \a n.
	 *
	 * Note that m needs to be non-negative since only those SFB coefficients are stored.
	 */
	std::complex<double> &SFB_coefficient(std::size_t l, std::size_t m, std::size_t n);

	/**
	 * Return the complex value of the SFB coefficient of the decomposed field at the angular and radial mode numbers \a
	 * l, \a m and \a n.
	 */
	std::complex<double> SFB_coefficient(std::size_t l, std::intmax_t m, std::size_t n) const;

	/**
	 * Return the value of the decomposed field at the coordinates \a  radius, \a theta and \a phi, smoothed with a
	 * Gaussian window of width \a smoothingScale (default: 0).
	 */
	double field(double radius, double theta, double phi, double smoothingScale = 0.0) const;

	/**
	 * Return the value of the radial component of the inverse divergence of the decomposed field at the coordinates \a
	 * radius, \a theta and \a phi, smoothed with a Gaussian window of width \a smoothingScale (default: 0).
	 */
	double radial_inverse_divergence(double radius, double theta, double phi, double smoothingScale = 0.0) const;

	/**
	 * Return the value of the theta component of the inverse divergence of the decomposed field at the coordinates \a
	 * radius, \a theta and \a phi, smoothed with a Gaussian window of width \a smoothingScale (default: 0).
	 */
	double theta_inverse_divergence(double radius, double theta, double phi, double smoothingScale = 0.0) const;

	/**
	 * Return the value of the phi component of the inverse divergence of the decomposed field at the coordinates \a
	 * radius, \a theta and \a phi, smoothed with a Gaussian window of width \a smoothingScale (default: 0).
	 */
	double phi_inverse_divergence(double radius, double theta, double phi, double smoothingScale = 0.0) const;

	/**
	 * Return the value of the inverse laplacian of the decomposed field at the coordinates \a  radius, \a theta and \a
	 * phi, smoothed with a Gaussian window of width \a smoothingScale (default: 0).
	 */
	double inverse_laplacian(double radius, double theta, double phi, double smoothingScale = 0.0) const;

	/**
	 * Return the \c SphericalGridFunction obtained by evaluating the decomposed field on a regular spherical grid with
	 * \a radialBinNumber radial bins between zero and the maximal radius of the SFB decomposition, \a thetaBinNumber
	 * theta bins and \a phiBinNumber phi bins, smoothed with a Gaussian window of width \a smoothingScale (default: 0).
	 */
	SphericalGridFunction field_grid(std::size_t radialBinNumber, std::size_t thetaBinNumber, std::size_t phiBinNumber, double smoothingScale = 0.0) const;

	/**
	 * Return the \c SphericalGridFunction obtained by evaluating the radial component of the inverse divergence of the
	 * decomposed field on a regular spherical grid with \a radialBinNumber radial bins between zero and the maximal
	 * radius of the SFB decomposition, \a thetaBinNumber theta bins and \a phiBinNumber phi bins, smoothed with a
	 * Gaussian window of width \a smoothingScale (default: 0).
	 */
	SphericalGridFunction radial_inverse_divergence_grid(std::size_t radialBinNumber, std::size_t thetaBinNumber, std::size_t phiBinNumber, double smoothingScale = 0.0) const;

	/**
	 * Return the \c SphericalGridFunction obtained by evaluating the theta component of the inverse divergence of the
	 * decomposed field on a regular spherical grid with \a radialBinNumber radial bins between zero and the maximal
	 * radius of the SFB decomposition, \a thetaBinNumber theta bins and \a phiBinNumber phi bins, smoothed with a
	 * Gaussian window of width \a smoothingScale (default: 0).
	 */
	SphericalGridFunction theta_inverse_divergence_grid(std::size_t radialBinNumber, std::size_t thetaBinNumber, std::size_t phiBinNumber, double smoothingScale = 0.0) const;

	/**
	 * Return the \c SphericalGridFunction obtained by evaluating the phi component of the inverse divergence of the
	 * decomposed field on a regular spherical grid with \a radialBinNumber radial bins between zero and the maximal
	 * radius of the SFB decomposition, \a thetaBinNumber theta bins and \a phiBinNumber phi bins, smoothed with a
	 * Gaussian window of width \a smoothingScale (default: 0).
	 */
	SphericalGridFunction phi_inverse_divergence_grid(std::size_t radialBinNumber, std::size_t thetaBinNumber, std::size_t phiBinNumber, double smoothingScale = 0.0) const;

	/**
	 * Return the \c SphericalGridFunction obtained by evaluating the inverse laplacian of the decomposed field on a
	 * regular spherical grid with \a radialBinNumber radial bins between zero and the maximal radius of the SFB
	 * decomposition, \a thetaBinNumber theta bins and \a phiBinNumber phi bins, smoothed with a Gaussian window of
	 * width \a smoothingScale (default: 0).
	 */
	SphericalGridFunction inverse_laplacian_grid(std::size_t radialBinNumber, std::size_t thetaBinNumber, std::size_t phiBinNumber, double smoothingScale = 0.0) const;

	/**
	 * Return the maximal radius used in the decomposition.
	 */
	double maximal_radius() const;

	/**
	 * Return the radial resolution used in the decomposition.
	 */
	double radial_resolution() const;

	/**
	 * Return the maximal multipole order used in the decomposition.
	 */
	std::size_t maximal_multipole() const;

	/**
	 * Return the maximal radial wavenumber used in the decomposition.
	 */
	double maximal_radial_wavenumber() const;

	/**
	 * Return the number of radial modes used in the decomposition.
	 */
	std::size_t radial_mode_number() const;

	/**
	 * Save this SphericalGridFunction to the binary file \a fileName.
	 *
	 * This allows to load this object at a later point using SphericalGridFunction::SphericalGridFunction(const
	 * std::string &) or SphericalGridFunction::load_object_from_file.
	 */
	void save_object_to_file(const std::string &fileName) const;

private:
	/**
	 * Maximal radius used in the decomposition.
	 */
	double MaxRadius;

	/**
	 * Radial resolution used in the decomposition (product of maximal radius and maximal considered wavenumber).
	 */
	double RadialResolution;

	/**
	 * Maximal multipole order used in the decomposition.
	 */
	std::size_t MaxMultipole;

	/**
	 * Maximal radial wavenumber used in the decomposition.
	 */
	double MaxRadialWavenumber;

	/**
	 * Number of radial modes used in the decomposition.
	 */
	std::size_t RadialModeNumber;

	/**
	 * Vector containing the complex SFB coefficients.
	 */
	std::vector<std::complex<double>> SFBCoefficients;

	/**
	 * Vector containing the radial wavenumbers k_ln.
	 */
	std::vector<double> RadialWavenumbers;

	/**
	 * Vector containing the radial normalization constants C_ln.
	 */
	std::vector<double> RadialNormalizationConstants;

	/**
	 * Initialize SphericalFourierBesselDecomposition::RadialModeNumber,
	 * SphericalFourierBesselDecomposition::SFBCoefficients, SphericalFourierBesselDecomposition::RadialWavenumbers and
	 * SphericalFourierBesselDecomposition::RadialNormalizationConstants.
	 */
	void initialize_SFB_decomposition();

	/**
	 * Compute the SFB coefficients in SphericalFourierBesselDecomposition::SFBCoefficients corresponding to the set
	 * of point sources at the radial and angular coordinates \a radialCoordinates, \a thetaCoordinates and \a
	 * phiCoordinates, using radial weighting function \a weightingFunction, and subtracting the radial function \a
	 * isotropicContribution.
	 */
	void compute_spherical_fourier_bessel_decomposition_of_discrete_field(const std::vector<double> &radialCoordinates, const std::vector<double> &thetaCoordinates, const std::vector<double> &phiCoordinates,
																		  const std::function<double(double radius)> &weightingFunction, const std::function<double(double radius)> &isotropicContribution);

	/**
	 * Loads a previously saved SphericalFourierBesselDecomposition from the binary file \a fileName.
	 * 
	 * The file must have been created using SphericalFourierBesselDecomposition::save_object_to_file.
	 */
	void load_object_from_file(const std::string &fileName);

	/**
	 * Auxiliary method that performs the actual computation of the value obtained by evaluating some linear functional
	 * of the decomposed field at the coordinates \a  radius, \a theta and \a phi, smoothed by a Gaussian on the scale
	 * \a smoothingScale. The linear functional is fixed by specifying the radial and angular kernels \a
	 * radialKernelFunction, \a thetaKernelFunction and \a phiKernelFunction, appearing in the SFB expansion of the
	 * linear functional of the decomposed field.
	 */
	double compute_value_for_given_kernels(double radius, double theta, double phi, double smoothingScale,
										   const std::function<double(std::size_t l, double k_ln, double radius)> &radialKernelFunction,
										   const std::function<double(std::size_t l, std::size_t m, double theta)> &thetaKernelFunction,
										   const std::function<std::complex<double>(std::size_t m, double phi)> &phiKernelFunction) const;

	/**
	 * Auxiliary method that performs the actual computation of the \c SphericalGridFunction obtained by evaluating some
	 * linear functional of the decomposed field on a regular spherical grid with \a radialBinNumber radial bins between
	 * zero and the maximal radius of the SFB decomposition, \a thetaBinNumber theta bins and \a phiBinNumber phi bins,
	 * smoothed by a Gaussian on the scale \a smoothingScale. The linear functional is fixed by specifying the radial
	 * and angular kernels \a radialKernelFunction, \a thetaKernelFunction and \a phiKernelFunction, appearing in the
	 * SFB expansion of the linear functional of the decomposed field.
	 */
	SphericalGridFunction compute_grid_function_for_given_kernels(std::size_t radialBinNumber, std::size_t thetaBinNumber, std::size_t phiBinNumber, double smoothingScale,
																  const std::function<double(std::size_t l, double k_ln, double radius)> &radialKernelFunction,
																  const std::function<double(std::size_t l, std::size_t m, double theta)> &thetaKernelFunction,
																  const std::function<std::complex<double>(std::size_t m, double phi)> &phiKernelFunction) const;

	/**
	 * Radial kernel appearing in the SFB expansion of the decomposed field.
	 */
	static double radial_kernel_field(std::size_t l, double k_ln, double radius);

	/**
	 * Radial kernel appearing in the SFB expansion of the radial component of the inverse divergence of the decomposed
	 * field.
	 */
	static double radial_kernel_radial_inverse_divergence(std::size_t l, double k_ln, double radius);

	/**
	 * Radial kernel appearing in the SFB expansion of the theta and phi components of the inverse divergence of the
	 * decomposed field.
	 */
	static double radial_kernel_angular_inverse_divergence(std::size_t l, double k_ln, double radius);

	/**
	 * Radial kernel appearing in the SFB expansion of the inverse laplacian of the decomposed field.
	 */
	static double radial_kernel_inverse_laplacian(std::size_t l, double k_ln, double radius);

	/**
	 * Theta kernel appearing in the SFB expansion of the decomposed field, the radial component of its inverse
	 * divergence and its inverse laplacian.
	 */
	static double theta_kernel_scalar_and_radial_component(std::size_t l, std::size_t m, double theta);

	/**
	 * Theta kernel appearing in the SFB expansion of the theta component of the inverse divergence of the decomposed
	 * field.
	 */
	static double theta_kernel_theta_inverse_divergence(std::size_t l, std::size_t m, double theta);

	/**
	 * Theta kernel appearing in the SFB expansion of the phi component of the inverse divergence of the decomposed
	 * field.
	 */
	static double theta_kernel_phi_inverse_divergence(std::size_t l, std::size_t m, double theta);

	/**
	 * Phi kernel appearing in the SFB expansion of the decomposed field, the radial and theta components of its inverse
	 * divergence, and its inverse laplacian.
	 */
	static std::complex<double> phi_kernel_scalar_radial_and_theta_component(std::size_t m, double phi);

	/**
	 * Phi kernel appearing in the SFB expansion of the phi components of the inverse divergence of the decomposed
	 * field.
	 */
	static std::complex<double> phi_kernel_phi_component(std::size_t m, double phi);
};

/** @} */

#endif