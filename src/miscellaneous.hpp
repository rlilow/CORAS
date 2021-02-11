#ifndef CORAS_MISCELLANEOUS_H
#define CORAS_MISCELLANEOUS_H

#include <cstddef>
#include <functional>
#include <string>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

/**
 * \addtogroup AUXILIARY
 * @{
 */

/**
 * \name Miscellaneous
 */
///@{

/**
 * \brief Class implementing a wrapper for the function used in the multi-dimensional minimizer of the GSL.
 *
 * It is derived from GSL's gsl_multimin_function, and allows the use of more general function types for the
 * minimization, e.g. lambda expressions or member functions.
 */
class GSLMultiminFunctionWrapper : public gsl_multimin_function
{
public:
	/**
	 * Constructor instantiating a wrapper for GSL's gsl_multimin_function, to minimize the function \a func, depending on
	 * a gsl_vector of dimension \a dimension (see GSL documentation for details).
	 */
	GSLMultiminFunctionWrapper(const std::size_t dimension, const std::function<double(const gsl_vector *minimizerArguments)> &func);

private:
	/**
	 * Function to be minimized.
	 */
	std::function<double(const gsl_vector *minimizerArguments)> Function;

	/**
	 * Wrapper for calling GSLMultiminFunctionWrapper::Function using the function signature expected by the GSL.
	 */
	static double invoke(const gsl_vector *minimizerArguments, void *parameters);
};

/**
 * Invert the general invertible square gsl_matrix \a matrix using singular value decomposition, and write the result to
 * \a inverse.
 *
 * \a inverse needs to be allocated beforehand and share the same dimensions with \a matrix. If both ar the same, the
 * inversion is performed in-place. Otherwise, \a matrix is not modified.
 */
void invert_matrix_svd(gsl_matrix *matrix, gsl_matrix *inverse);

/**
 * Return \c true if the file \a fileName already exists and \c false otherwise.
 */
bool file_exists(const std::string &fileName);

///@}

/** @} */

#endif