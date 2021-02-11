#ifndef CORAS_LINEAR_RSD_CORRECTION_H
#define CORAS_LINEAR_RSD_CORRECTION_H

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "SphericalFourierBesselDecomposition.hpp"
#include "transformations.hpp"

/**
 * \addtogroup RECONSTRUCTION
 * @{
 */

/**
 * \brief Class implementing the correction of linear redshift-space distortions (RSD) for the spherical Fourier Bessel
 * (SFB) decomposition of the density contrast.
 *
 * It provides methods to transform from redshift- to configuration space and vice versa.
 */
class LinearRSDCorrection
{
public:
	/**
	 * Constructor instantiating the linear RSD correction for the SFB decomposed density contrast, using a maximal
	 * radius \a maxRadius, radial resolution \a radialResolution (product of maximal radius and maximal considered
	 * radial wavenumber) and maximal multipole order \a maxMultipole. It also depends on the logarithmic derivative of
	 * the selection function \a selectionFunctionLogDerivative, the product of selection and weighting function \a
	 * selectionTimesWeightingFunction, and on the reference frame \a referenceFrame (default: LOCAL_GROUP_FRAME) in
	 * which the redshifts are specified. For the latter, only the choices LOCAL_GROUP_FRAME and CMB_FRAME allowed.
	 *
	 * Note that this constructor precomputes the entries of the radial coupling matrix, which is computationally the
	 * most expensive step of the RSD correction. Therefore, the constructor should only be called once for the same
	 * arguments.
	 */
	LinearRSDCorrection(double maxRadius, double radialResolution, std::size_t maxMultipole,
						const std::function<double(double)> &selectionFunctionLogDerivative, const std::function<double(double)> &selectionTimesWeightingFunction, ReferenceFrame referenceFrame = LOCAL_GROUP_FRAME);

	/**
	 * Constructor instantiating a previously saved LinearRSDCorrection from the binary file \a fileName.
	 * 
	 * The file must have been created using LinearRSDCorrection::save_object_to_file.
	 */
	LinearRSDCorrection(const std::string &fileName);

	/**
	 * Constructor instantiating an empty LinearRSDCorrection object, intended to be properly initialized at a later
	 * point. For example, this can be used to choose a constructor based on some conditional statement.
	 */
	LinearRSDCorrection();

	/**
	 * Transform the SFB decomposition of the density contrast \a densityContrastSFBDecomposition from redshift to
	 * configuration space, using the following linear relation between velocity field and density contrast, div v = -
	 * \a beta delta.
	 */
	void transform_redshift_to_configuration_space(SphericalFourierBesselDecomposition &densityContrastSFBDecomposition, double beta) const;

	/**
	 * Transform the SFB decomposition of the density contrast \a densityContrastSFBDecomposition from configuration to
	 * redshift space, using the following linear relation between velocity field and density contrast, div v = - \a
	 * beta delta
	 */
	void transform_configuration_to_redshift_space(SphericalFourierBesselDecomposition &densityContrastSFBDecomposition, double beta) const;

	/**
	 * Save this LinearRSDCorrection to the binary file \a fileName.
	 *
	 * This allows to load this object at a later point using LinearRSDCorrection::LinearRSDCorrection(const std::string
	 * &) or LinearRSDCorrection::load_object_from_file.
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
	 * Entries of the reduced radial coupling matrix, i.e. the radial coupling matrix for beta = 1 and with the identity
	 * matrix subtracted.
	 */
	std::vector<double> ReducedCouplingMatrix;

	/**
	 * Compute the entries of the reduced radial coupling matrix, using the logarithm of the selection function \a
	 * selectionFunctionLogDerivative, the product of selection and weighting function \a
	 * selectionTimesWeightingFunction, and the redshift reference frame \a referenceFrame, and write them into
	 * LinearRSDCorrection::ReducedCouplingMatrix.
	 */
	void compute_coupling_matrix(const std::function<double(double)> &selectionFunctionLogDerivative, const std::function<double(double)> &selectionTimesWeightingFunction, ReferenceFrame referenceFrame);

	/**
	 * Loads a previously saved LinearRSDCorrection from the binary file \a fileName.
	 * 
	 * The file must have been created using LinearRSDCorrection::save_object_to_file.
	 */
	void load_object_from_file(const std::string &fileName);

	/**
	 * Return the value of the integral kernel used in the computation of the radial coupling matrix. The kernel depends
	 * on the multipole order \a l, the two radial wavenumbers \a k1 and \a k2, the radius \a r, the logarithm of the
	 * selection function \a selectionFunctionLogDerivative, and the product of selection and weighting function \a
	 * selectionTimesWeightingFunction. If Local Group frame redshifts are used, \a subtractObserverVelocity needs to be
	 * set to \c true, to subtract the observer velocity. If CMB frame redshifts are used, it needs to be set to \c
	 * false.
	 */
	double coupling_matrix_kernel(std::size_t l, double k1, double k2, double r, const std::function<double(double)> &selectionFunctionLogDerivative, const std::function<double(double)> &selectionTimesWeightingFunction, bool subtractObserverVelocity) const;
};

/** @} */

#endif