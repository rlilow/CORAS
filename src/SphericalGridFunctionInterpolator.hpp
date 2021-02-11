#ifndef CORAS_SPHERICAL_GRID_FUNCTION_INTERPOLATOR_H
#define CORAS_SPHERICAL_GRID_FUNCTION_INTERPOLATOR_H

#include <cstddef>
#include <vector>

#include "SphericalGridFunction.hpp"

/**
 * \addtogroup GRIDFUNCTIONS
 * @{
 */

/**
 * \brief Class implementing a linear interpolator between a 1-parameter family of SphericalGridFunction objects.
 */
class SphericalGridFunctionInterpolator
{
public:
	/**
	 * Constructor instantiating a linear interpolator between the 1-parameter family of SphericalGridFunction objects
	 * contained in the vector \a gridFunctions. The parameter is discretized into equidistant bins between \a
	 * minParameter for the first vector entry to \a maxParameter for the last entry.
	 */
	SphericalGridFunctionInterpolator(double minParameter, double maxParameter, const std::vector<SphericalGridFunction> &gridFunctions);

	/**
	 * Return the linearly interpolated SphericalGridFunction at the parameter \a parameter.
	 */
	SphericalGridFunction interpolate(double parameter) const;

	/**
	 * Overloads the parentheses operator to return the linearly interpolated SphericalGridFunction at the parameter \a
	 * parameter.
	 *
	 * This provides the same functionality as SphericalGridFunctionInterpolator::interpolate(double).
	 */
	SphericalGridFunction operator()(double parameter) const;

	/**
	 * Return the linearly interpolated function value at the parameter \a parameter, and the coordinates \a radius, \a
	 * theta and \a phi.
	 */
	double interpolate(double parameter, double radius, double theta, double phi) const;

	/**
	 * Overloads the parentheses operator to return the linearly interpolated function value at the parameter \a
	 * parameter, and the coordinates \a radius, \a theta and \a phi.
	 *
	 * This provides the same functionality as SphericalGridFunctionInterpolator::interpolate(double, double, double,
	 * double).
	 */
	double operator()(double parameter, double radius, double theta, double phi) const;

	/**
	 * Return the SphericalGridFunction at the parameter bin \a bin.
	 */
	SphericalGridFunction at(std::size_t bin) const;

	/**
	 * Return the linearly interpolated function value of the SphericalGridFunction at the parameter bin \a bin.
	 */
	double at(std::size_t bin, double radius, double theta, double phi) const;

	/**
	 * Return the parameter value at the parameter bin \a bin.
	 */
	double parameter(std::size_t bin) const;

	/**
	 * Return the minimal parameter value.
	 */
	double min_parameter() const;

	/**
	 * Return the minimal parameter value.
	 */
	double max_parameter() const;

	/**
	 * Return the parameter bin width.
	 */
	double bin_width() const;

	/**
	 * Return the number of parameter bins.
	 */
	std::size_t bin_number() const;

private:
	/**
	 * Minimal parameter value.
	 */
	double MinParameter;

	/**
	 * Maximal parameter value.
	 */
	double MaxParameter;

	/**
	 * Reference to the vector containing the SphericalGridFunction objects.
	 */
	const std::vector<SphericalGridFunction> &GridFunctions;

	/**
	 * Vector containing the parameter values.
	 */
	std::vector<double> ParameterValues;

	/**
	 * Number of parameter bins.
	 */
	std::size_t BinNumber;

	/**
	 * Width of parameter bins.
	 */
	double BinWidth;
};

/** @} */

#endif