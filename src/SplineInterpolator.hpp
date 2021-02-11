#ifndef CORAS_SPLINE_INTERPOLATOR_H
#define CORAS_SPLINE_INTERPOLATOR_H

#include <vector>

#include <gsl/gsl_interp.h>

/**
 * \addtogroup AUXILIARY
 * @{
 */

/**
 * \brief Class implementing a wrapper for the one-dimensional spline interpolation of the GSL.
 */
class SplineInterpolator
{
public:
	/**
	 * Constructor instantiating a spline interpolator, computing a continouos interpolating function passing through
	 * the data points specified by the arguments \a arguments and their respective function values \a functionValues,
	 * using the GSL interpolation type \a gslInterpolationType (see GSL documentation for details).
	 */
	SplineInterpolator(const std::vector<double> &arguments, const std::vector<double> &functionValues, const gsl_interp_type *gslInterpolationType = gsl_interp_cspline);
	SplineInterpolator(const SplineInterpolator &otherSplineInterpolator);
	SplineInterpolator();
	SplineInterpolator &operator=(const SplineInterpolator &otherSplineInterpolator);

	/**
	 * Return the value of the interpolating function at the argument \a argument.
	 */
	double evaluate(double argument) const;

	/**
	 * Overloads the parentheses operator to return the value of the interpolating function at the argument \a argument.
     *
     * This provides the same functionality as SplineInterpolator::evaluate.
	 */
	double operator()(double argument) const;

	/**
	 * Return the derivative of the interpolating function at the argument \a argument.
	 */
	double derivative(double argument) const;

	/**
	 * Return the value of the second derivative of the interpolating function at the argument \a argument.
	 */
	double second_derivative(double argument) const;

	/**
	 * Return the value of the integral of the interpolating function between the arguments \a
	 * lowerlowerIntegrationLimit and \a upperIntegrationLimit.
	 */
	double integral(double lowerIntegrationLimit, double upperIntegrationLimit) const;

	/**
	 * Deconstructor taking care of freeing the GSL interpolation objects SplineInterpolator::Accelerator and
	 * SplineInterpolator::Spline.
	 */
	~SplineInterpolator();

private:
	std::vector<double> Arguments;
	std::vector<double> FunctionValues;

	const gsl_interp_type *InterpolationType;

	/**
	 * GSL interpolation accelerator.
	 */
	gsl_interp_accel *Accelerator;
	gsl_interp *Interpolator;
};

/** @} */

#endif