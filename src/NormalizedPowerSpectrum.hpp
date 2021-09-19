#ifndef CORAS_NORMALISED_POWER_SPECTRUM_H
#define CORAS_NORMALISED_POWER_SPECTRUM_H

#include <functional>
#include <vector>

#include <gsl/gsl_interp.h>

#include "SplineInterpolator.hpp"

/**
 * \addtogroup RECONSTRUCTION
 * @{
 */

/**
 * \brief Class implementing a density contrast power spectrum normalised to the density variance on some smoothing
 * scale.
 */
class NormalizedPowerSpectrum
{
public:
	NormalizedPowerSpectrum(const std::vector<double> &wavenumbers, const std::vector<double> &powerSpectrumValues,
							double lengthRescaling, double smoothingScale,
							bool useTophatFilter = true,
							double integrationAbsoluteError = 0.0, double integrationRelativeError = 1e-5,
							const gsl_interp_type *gslInterpolationType = gsl_interp_cspline);

	double evaluate(double k) const;
	double operator()(double k) const;
	double smoothing_scale() const;
	double variance() const;
	double variance(const std::function<double(double k)> &kernel,
					bool useTophatFilter = true,
					double integrationAbsoluteError = 0.0, double integrationRelativeError = 1e-5) const;

private:
	double SmoothingScale;
	SplineInterpolator Interpolator;
	double Variance;
	double MaxWavenumber;

	double tophat_filter_integrand(double k) const;
	double gaussian_filter_integrand(double k) const;
};

/** @} */

#endif