#ifndef CORAS_WIENER_FILTER_H
#define CORAS_WIENER_FILTER_H

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "SphericalFourierBesselDecomposition.hpp"

/**
 * \addtogroup RECONSTRUCTION
 * @{
 */

/**
 * \brief Class implementing the Wiener filter for the spherical Fourier Bessel (SFB) decomposition of the density
 * contrast of a set of points.
 *
 * It provides methods to apply the Wiener filter, and to compute the expected variance of the residual between Wiener
 * filtered and original density contrast as well as radial velocity field.
 */
class WienerFilter
{
public:
	/**
	 * Constructor instantiating the Wiener filter for the SFB decomposed density contrast of a set of points, using a
	 * maximal radius \a maxRadius, radial resolution \a radialResolution (product of maximal radius and maximal
	 * considered radial wavenumber) and maximal multipole order \a maxMultipole. It also depends on the selection
	 * function \a selectionFunction. 
	 *
	 * Note that this constructor precomputes the entries of the noise matrix, which is computationally the most
	 * expensive step of the Wiener filter application. Therefore, the constructor should only be called once for the
	 * same arguments.
	 */
	WienerFilter(double maxRadius, double radialResolution, std::size_t maxMultipole,
				 const std::function<double(double)> &selectionFunction);

	/**
	 * Constructor instantiating a previously saved WienerFilter from the binary file \a fileName.
	 *
	 * The file must have been created using WienerFilter::save_object_to_file.
	 */
	WienerFilter(const std::string &fileName);

	/**
	 * Constructor instantiating an empty WienerFilter object, intended to be properly initialized at a later point. For
	 * example, this can be used to choose a constructor based on some conditional statement.
	 */
	WienerFilter();

	/**
	 * Apply the Wiener filter to the SFB decomposition of the density contrast \a densityContrastSFBDecomposition,
	 * using the mean number density \a meanDensity, as well as the power spectra of the density contrast data-data and
	 * signal-data correlations, \a dataDataPowerSpectrum and \a signalDataPowerSpectrum, respectively.
	 */
	void apply_to(SphericalFourierBesselDecomposition &densityContrastSFBDecomposition, double meanDensity, const std::function<double(double)> &dataDataPowerSpectrum, const std::function<double(double)> &signalDataPowerSpectrum) const;

	/**
	 * Save this WienerFilter to the binary file \a fileName.
	 *
	 * This allows to load this object at a later point using WienerFilter::WienerFilter(const std::string &) or
	 * WienerFilter::load_object_from_file.
	 */
	void save_object_to_file(const std::string &fileName) const;

private:
	/**
	 * Maximal radius used in the SFB decomposition.
	 */
	double MaxRadius;

	/**
	 * Radial resolution used in the SFB decomposition (product of maximal radius and maximal considered wavenumber).
	 */
	double RadialResolution;

	/**
	 * Maximal multipole used in the SFB decomposition.
	 */
	std::size_t MaxMultipole;

	/**
	 * Maximal radial wavenumber used in the SFB decomposition.
	 */
	double MaxRadialWavenumber;

	/**
	 * Maximal number of radial modes used in the SFB decomposition.
	 */
	std::size_t RadialModeNumber;

	/**
	 * Entries of the reduced noise matrix, i.e. the noise matrix multiplied by the mean number density.
	 */
	std::vector<double> ReducedNoiseMatrix;

	/**
	 * Compute the entries of the reduced noise matrix, using the selection function \a selectionFunction, and write
	 * them into WienerFilter::ReducedNoiseMatrix.
	 */
	void compute_noise_matrix(const std::function<double(double)> &selectionFunction);

	/**
	 * Loads a previously saved WienerFilter from the binary file \a fileName.
	 * 
	 * The file must have been created using WienerFilter::save_object_to_file.
	 */
	void load_object_from_file(const std::string &fileName);

	/**
	 * Return the value of the integral kernel used in the computation of the noise matrix. The kernel depends on the
	 * multipole order \a l, the two radial wavenumbers \a k1 and \a k2, the radius \a r, and the selection function \a
	 * selectionFunction.
	 */
	double noise_matrix_kernel(std::size_t l, double k1, double k2, double r, const std::function<double(double)> &selectionFunction) const;
};

/** @} */

#endif