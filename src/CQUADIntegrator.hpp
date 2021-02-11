#ifndef CORAS_CQUAD_INTEGRATOR_H
#define CORAS_CQUAD_INTEGRATOR_H

#include <cstddef>
#include <functional>

#include <gsl/gsl_integration.h>

/**
 * \addtogroup AUXILIARY
 * @{
 */

/**
 * \brief Class implementing a wrapper for the one-dimensional CQUAD integration algorithm of the GSL.
 */
class CQUADIntegrator
{
public:
	/**
	 * Constructor instantiating a CQUAD integrator which is integrating the function \a integrand, using a CQUAD
	 * workspace size \a workSpaceSize (see GSL documentation for details).
	 */
	CQUADIntegrator(const std::function<double(double)> &integrand, std::size_t workSpaceSize = 100);

	/**
	 * Perform the integral from \a lowerIntegrationLimit to \a upperIntegrationLimit, using absolute and relative error
	 * bounds \a absoluteErrorBound and \a relativeErrorBound, respectively, and write the resulting value into \a
	 * result, the estimated absolute error into \a error and the number of integrand evaluations needed into \a
	 * evaluationNumber.
	 */
	void integrate(double lowerIntegrationLimit, double upperIntegrationLimit,
				   double absoluteErrorBound, double relativeErrorBound,
				   double &result, double &error, std::size_t &evaluationNumber) const;

	/**
	 * Perform the integral from \a lowerIntegrationLimit to \a upperIntegrationLimit, using absolute and relative error
	 * bounds \a absoluteErrorBound and \a relativeErrorBound, respectively, and write the resulting value into \a
	 * result and the estimated absolute error into \a error.
	 */
	void integrate(double lowerIntegrationLimit, double upperIntegrationLimit,
				   double absoluteErrorBound, double relativeErrorBound,
				   double &result, double &error) const;

	/**
	 * Perform the integral from \a lowerIntegrationLimit to \a upperIntegrationLimit, using absolute and relative error
	 * bounds \a absoluteErrorBound and \a relativeErrorBound, respectively, and write the resulting value into \a
	 * result.
	 */
	void integrate(double lowerIntegrationLimit, double upperIntegrationLimit,
				   double absoluteErrorBound, double relativeErrorBound,
				   double &result) const;

	/**
	 * Deconstructor taking care of freeing the GSL integration workspace CQUADIntegrator::Workspace.
	 */
	~CQUADIntegrator();

private:
	/**
	 * Function to be integrated.
	 */
	std::function<double(double)> Integrand;

	/**
	 * Size of the GSL integration workspace used.
	 */
	std::size_t WorkSpaceSize;

	/**
	 * GSL integration workspace.
	 */
	gsl_integration_cquad_workspace *Workspace;

	/**
	 * GSL function \c struct used by the GSL integrator.
	 */
	gsl_function GSLFunction;

	/**
	 * Auxiliary function that provides the argument structure expected by the GSL.
	 */
	static double gsl_integrand(double variable, void *parameters);
};

/** @} */

#endif