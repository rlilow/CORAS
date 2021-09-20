#ifndef CORAS_SPHERICAL_GRID_FUNCTION_H
#define CORAS_SPHERICAL_GRID_FUNCTION_H

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "Cartesian1DGridFunction.hpp"

/**
 * \defgroup GRIDFUNCTIONS Grid Functions
 *
 * \brief Operations on functions discretized on a regular grid.
 * @{
 */

/**
 * \brief Class implementing a function discretized on a regular spherical grid.
 */
class SphericalGridFunction
{
public:
	/**
	 * Constructor instantiating a constant function with value \a value discretized on a regular spherical grid with
	 * maximal radius \a maxRadius, and \a radialBinNumber, \a thetaBinNumber and \a phiBinNumber bins on the radial,
	 * theta and phi axes, respectively.
	 */
	SphericalGridFunction(double maxRadius, std::size_t radialBinNumber, std::size_t thetaBinNumber, std::size_t phiBinNumber,
						  double value = 0.0);

	/**
	 * Constructor instantiating the function \a func discretized on a regular spherical grid with maximal radius \a
	 * maxRadius, and \a radialBinNumber, \a thetaBinNumber and \a phiBinNumber bins on the radial, theta and phi axes,
	 * respectively.
	 */
	SphericalGridFunction(double maxRadius, std::size_t radialBinNumber, std::size_t thetaBinNumber, std::size_t phiBinNumber,
						  const std::function<double(double radius, double theta, double phi)> &func);

	/**
	 * Constructor instantiating a previously saved SphericalGridFunction from the binary file \a fileName.
	 * 
	 * The file must have been created using SphericalGridFunction::save_object_to_file.
	 */
	SphericalGridFunction(const std::string &fileName);

	/**
	 * Copy-constructor instantiating the same grid function as \a otherSphericalGridFunction.
	 */
	SphericalGridFunction(const SphericalGridFunction &otherSphericalGridFunction) = default;

	/**
	 * Copy-constructor instantiating the same grid function as \a otherSphericalGridFunction and applying the function
	 * \a func, depending on the grid function value, afterwards.
	 */
	SphericalGridFunction(const SphericalGridFunction &otherSphericalGridFunction, const std::function<double(double value)> &func);

	/**
	 * Copy-constructor instantiating the same grid function as \a otherSphericalGridFunction and applying the function
	 * \a func, depending on the grid function value and the coordinates, afterwards.
	 */
	SphericalGridFunction(const SphericalGridFunction &otherSphericalGridFunction, const std::function<double(double value, double radius, double theta, double phi)> &func);

	/**
	 * Constructor instantiating an empty SphericalGridFunction object, intended to be properly initialized at a later
	 * point. For example, this can be used to choose a constructor based on some conditional statement.
	 */
	SphericalGridFunction();

	/**
	 * Return a reference to the function value at the grid point specified by \a radialBin, \a thetaBin and \a phiBin.
	 */
	double &value(std::size_t radialBin, std::size_t thetaBin, std::size_t phiBin);

	/**
	 * Return the function value at the grid point specified by \a radialBin, \a thetaBin and \a phiBin.
	 */
	double value(std::size_t radialBin, std::size_t thetaBin, std::size_t phiBin) const;

	/**
	 * Return a reference to the vector of all function values.
	 */
	std::vector<double> &value_vector();

	/**
	 * Return a reference to the \c const vector of all function values.
	 */
	const std::vector<double> &value_vector() const;

	/**
	 * Return the linearly interpolated function value at the coordinates \a  radius, \a theta and \a phi.
	 */
	double interpolate(double radius, double theta, double phi) const;

	/**
	 * Overloads the parentheses operator to return the linearly interpolated function value at the coordinates \a
	 * radius, \a theta and \a phi.
	 *
	 * This provides the same functionality as SphericalGridFunction::interpolate.
	 */
	double operator()(double radius, double theta, double phi) const;

	/**
	 * Return the volume element at the grid point specified by the radial and theta bins \a radialBin and \a thetaBin,
	 * respectively.
	 */
	double volume_element(std::size_t radialBin, std::size_t thetaBin) const;

	/**
	 * Return the solid angle element at the grid point specified by the theta bin \a thetaBin.
	 */
	double solid_angle_element(std::size_t thetaBin) const;

	/**
	 * Return the radial coordinate value at the radial-axis point specified by \a radialBin.
	 */
	double radial_coordinate(std::size_t radialBin) const;

	/**
	 * Return the theta coordinate value at the theta-axis point specified by \a thetaBin.
	 */
	double theta_coordinate(std::size_t thetaBin) const;

	/**
	 * Return the phi coordinate value at the phi-axis point specified by \a phiBin.
	 */
	double phi_coordinate(std::size_t phiBin) const;

	/**
	 * Return the vector of all radial-coordinate values.
	 */
	std::vector<double> radial_coordinate_vector() const;

	/**
	 * Return the vector of all theta-coordinate values.
	 */
	std::vector<double> theta_coordinate_vector() const;

	/**
	 * Return the vector of all phi-coordinate values.
	 */
	std::vector<double> phi_coordinate_vector() const;

	/**
	 * Return the maximal radial-coordinate value.
	 */
	double maximal_radius() const;

	/**
	 * Return the radial-axis bin width.
	 */
	double radial_bin_width() const;

	/**
	 * Return the theta-axis bin width.
	 */
	double theta_bin_width() const;

	/**
	 * Return the theta-axis bin width.
	 */
	double phi_bin_width() const;

	/**
	 * Return the number of radial-axis bins.
	 */
	std::size_t radial_bin_number() const;

	/**
	 * Return the number of theta-axis bins.
	 */
	std::size_t theta_bin_number() const;

	/**
	 * Return the number of phi-axis bins.
	 */
	std::size_t phi_bin_number() const;

	/**
	 * Apply the function \a func, depending on the grid function value, to all grid points. If \a parallelized is
	 * set to \c true, this is performed for different grid points in parallel.
	 */
	void apply(const std::function<double(double value)> &func, const bool parallelized = false);

	/**
	 * Apply the function \a func, depending on the grid function value and coordinates, to all grid points. If \a
	 * parallelized is set to \c true, this is performed for different grid points in parallel.
	 */
	void apply(const std::function<double(double value, double radius, double theta, double phi)> &func, const bool parallelized = false);

	/**
	 * Evaluate the function \a func, depending on the grid bins, for all grid points. If \a parallelized is set to
	 * \c true, this is performed for different grid points in parallel.
	 */
	void evaluate_for_all_grid_points(const std::function<void(std::size_t radialBin, std::size_t thetaBin, std::size_t phiBin)> &func, const bool parallelized = false) const;

	/**
	 * Return the average of the function \a func, depending on the grid function value, over the whole grid volume.
	 *
	 * If no argument is given, the volume average of the grid function itself is returned.
	 */
	double average(const std::function<double(double value)> &func = [](double value)
				   { return value; }) const;

	/**
	 * Return a Cartesian1DGridFunction describing the averages of the function \a func, depending on the grid function
	 * value, over spheres of different radii. The averages are computed for the radii specified by the radial grid
	 * coordinate axis.
	 */
	Cartesian1DGridFunction averages(const std::function<double(double value)> &func = [](double value)
									 { return value; }) const;

	/**
	 * Return the average of the function \a func, depending on the grid function value and the coordinates, over the
	 * whole grid volume.
	 */
	double average(const std::function<double(double value, double radius, double theta, double phi)> &func) const;

	/**
	 * Return a Cartesian1DGridFunction describing the averages of the function \a func, depending on the grid function
	 * value and the coordinates, over spheres of different radii. The averages are computed for the radii specified by
	 * the radial grid coordinate axis.
	 */
	Cartesian1DGridFunction averages(const std::function<double(double value, double radius, double theta, double phi)> &func) const;

	/**
	 * Return a Cartesian1DGridFunction describing the radially dependent average of the function \a func, depending on
	 * the grid function value, over the whole solid angle.
	 */
	Cartesian1DGridFunction angular_average(const std::function<double(double value)> &func = [](double value)
											{ return value; }) const;

	/**
	 * Return a Cartesian1DGridFunction describing the radially dependent average of the function \a func, depending on
	 * the grid function value and coordinates, over the whole solid angle.
	 */
	Cartesian1DGridFunction angular_average(const std::function<double(double value, double radius, double theta, double phi)> &func) const;

	/**
	 * Write the grid function values into the text file \a fileName.
	 */
	void write_to_file(const std::string &fileName) const;

	/**
	 * Save this SphericalGridFunction to the binary file \a fileName.
	 *
	 * This allows to load this object at a later point using SphericalGridFunction::SphericalGridFunction(const
	 * std::string &) or SphericalGridFunction::load_object_from_file.
	 *
	 * For a human-readable text file output use SphericalGridFunction::write_to_file instead.
	 */
	void save_object_to_file(const std::string &fileName) const;

	/**
	 * Assign this SphericalGridFunction to the SphericalGridFunction \a otherGridFunction.
	 */
	SphericalGridFunction &operator=(const SphericalGridFunction &otherGridFunction) = default;

	/**
	 * Assign this SphericalGridFunction to the the elementwise sum of itself and the SphericalGridFunction \a
	 * otherGridFunction, which needs to share the same numbers of bins as well as maximal radius.
	 */
	SphericalGridFunction &operator+=(const SphericalGridFunction &otherGridFunction);

	/**
	 * Assign this SphericalGridFunction to the elementwise difference of itself and the SphericalGridFunction \a
	 * otherGridFunction, which needs to share the same numbers of bins as well as maximal radius.
	 */
	SphericalGridFunction &operator-=(const SphericalGridFunction &otherGridFunction);

	/**
	 * Assign this SphericalGridFunction to the elementwise product of itself and the SphericalGridFunction \a
	 * otherGridFunction, which needs to share the same numbers of bins as well as maximal radius.
	 */
	SphericalGridFunction &operator*=(const SphericalGridFunction &otherGridFunction);

	/**
	 * Assign this SphericalGridFunction to the elementwise quotient of itself and the SphericalGridFunction \a
	 * otherGridFunction, which needs to share the same numbers of bins as well as maximal radius.
	 */
	SphericalGridFunction &operator/=(const SphericalGridFunction &otherGridFunction);

	/**
	 * Assign this SphericalGridFunction to the elementwise sum of itself and the scalar \a otherValue.
	 */
	SphericalGridFunction &operator+=(double otherValue);

	/**
	 * Assign this SphericalGridFunction to the elementwise difference of itself and the scalar \a otherValue.
	 */
	SphericalGridFunction &operator-=(double otherValue);

	/**
	 * Assign this SphericalGridFunction to the elementwise product of itself and the scalar \a otherValue.
	 */
	SphericalGridFunction &operator*=(double otherValue);

	/**
	 * Assign this SphericalGridFunction to the elementwise quotient of itself and the scalar \a otherValue.
	 */
	SphericalGridFunction &operator/=(double otherValue);

	/**
	 * Return \c true if this SphericalGridFunction is defined on the same spherical grid as \a otherGridFunction,
	 * otherwise returns \c false.
	 */
	bool has_same_grid_as(const SphericalGridFunction &otherGridFunction) const;

private:
	/**
	 * Maximal radial-coordinate value.
	 */
	double MaxRadius;

	/**
	 * Number of radial-axis bins.
	 */
	std::size_t RadialBinNumber;

	/**
	 * Number of theta-axis bins.
	 */
	std::size_t ThetaBinNumber;

	/**
	 * Number of phi-axis bins.
	 */
	std::size_t PhiBinNumber;

	/**
	 * Width of radial-axis bins.
	 */
	double RadialBinWidth;

	/**
	 * Width of theta-axis bins.
	 */
	double ThetaBinWidth;

	/**
	 * Width of phi-axis bins.
	 */
	double PhiBinWidth;

	/**
	 * Vector containing all the radial-coordinate values.
	 */
	std::vector<double> RadialCoordinates;

	/**
	 * Vector containing all the theta-coordinate values.
	 */
	std::vector<double> ThetaCoordinates;

	/**
	 * Vector containing all the phi-coordinate values.
	 */
	std::vector<double> PhiCoordinates;

	/**
	 * Vector containing all the grid function values.
	 */
	std::vector<double> FunctionValues;

	/**
	 * Initializes the vectors SphericalGridFunction::RadialCoordinates, SphericalGridFunction::ThetaCoordinates and
	 * SphericalGridFunction::PhiCoordinates.
	 */
	void initialize_coordinates();

	/**
	 * Loads a previously saved SphericalGridFunction from the binary file \a fileName.
	 * 
	 * The file must have been created using SphericalGridFunction::save_object_to_file.
	 */
	void load_object_from_file(const std::string &fileName);
};

/**
 * \name SphericalGridFunction Binary Operators
 *
 * Binary operators between SphericalGridFunction objects and scalars.
 */
///@{

/**
 * Return the elementwise sum of the SphericalGridFunction objects \a leftGridFunction and \a rightGridFunction, which
 * need to share the same numbers of bins as well as maximal radius.
 */
SphericalGridFunction operator+(const SphericalGridFunction &leftGridFunction, const SphericalGridFunction &rightGridFunction);

/**
 * Return the elementwise difference of the SphericalGridFunction objects \a leftGridFunction and \a rightGridFunction,
 * which need to share the same coordinate grid.
 */
SphericalGridFunction operator-(const SphericalGridFunction &leftGridFunction, const SphericalGridFunction &rightGridFunction);

/**
 * Return the elementwise product of the SphericalGridFunction objects \a leftGridFunction and \a rightGridFunction,
 * which need to share the same coordinate grid.
 */
SphericalGridFunction operator*(const SphericalGridFunction &leftGridFunction, const SphericalGridFunction &rightGridFunction);

/**
 * Return the elementwise quotient of the SphericalGridFunction objects \a leftGridFunction and \a rightGridFunction,
 * which need to share the same coordinate grid.
 */
SphericalGridFunction operator/(const SphericalGridFunction &leftGridFunction, const SphericalGridFunction &rightGridFunction);

/**
 * Return the elementwise sum of the SphericalGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
SphericalGridFunction operator+(const SphericalGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise difference of the SphericalGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
SphericalGridFunction operator-(const SphericalGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise product of the SphericalGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
SphericalGridFunction operator*(const SphericalGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise quotient of the SphericalGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
SphericalGridFunction operator/(const SphericalGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise sum of the scalar \a leftValue and the SphericalGridFunction \a rightGridFunction.
 */
SphericalGridFunction operator+(double leftValue, const SphericalGridFunction &rightGridFunction);

/**
 * Return the elementwise difference of the scalar \a leftValue and the SphericalGridFunction \a rightGridFunction.
 */
SphericalGridFunction operator-(double leftValue, const SphericalGridFunction &rightGridFunction);

/**
 * Return the elementwise product of the scalar \a leftValue and the SphericalGridFunction \a rightGridFunction.
 */
SphericalGridFunction operator*(double leftValue, const SphericalGridFunction &rightGridFunction);

/**
 * Return the elementwise quotient of the scalar \a leftValue and the SphericalGridFunction \a rightGridFunction.
 */
SphericalGridFunction operator/(double leftValue, const SphericalGridFunction &rightGridFunction);

///@}

/**
 * \name SphericalGridFunction Transformations
 *
 * General transformations between up to three SphericalGridFunction objects.
 */
///@{

/**
 * Compute the transformation \a trafo of the SphericalGridFunction \a input, and write the result into \a output. \a
 * trafo is a function of the grid bins.
 *
 * The transformation can be performed in-place. Otherwise, \a input is not modified.
 */
void transform_spherical_grid_functions(const SphericalGridFunction &input,
										SphericalGridFunction &output,
										const std::function<void(std::size_t radialBin, std::size_t thetaBin, std::size_t phiBin,
																 double &output)> &trafo);

/**
 * Compute the transformation \a trafo of the SphericalGridFunction objects \a input1 and \a input2, which need to share
 * the same coordinate grid, and write the result into \a output1 and \a output2. \a trafo is a function of the grid
 * bins.
 *
 * The transformation can be performed (partially) in-place. Otherwise, \a input1 and \a input2 are not modified.
 */
void transform_spherical_grid_functions(const SphericalGridFunction &input1, const SphericalGridFunction &input2,
										SphericalGridFunction &output1, SphericalGridFunction &output2,
										const std::function<void(std::size_t radialBin, std::size_t thetaBin, std::size_t phiBin,
																 double &output1, double &output2)> &trafo);

/**
 * Compute the transformation \a trafo of the SphericalGridFunction objects \a input1, \a input2 and \a input3, which
 * need to share the same coordinate grid, and write the result into \a output1, \a output2 and \a output3. \a trafo is
 * a function of the grid bins.
 *
 * The transformation can be performed (partially) in-place. Otherwise, \a input1, \a input2 and \a input3 are not
 * modified.
 */
void transform_spherical_grid_functions(const SphericalGridFunction &input1, const SphericalGridFunction &input2, const SphericalGridFunction &input3,
										SphericalGridFunction &output1, SphericalGridFunction &output2, SphericalGridFunction &output3,
										const std::function<void(std::size_t radialBin, std::size_t thetaBin, std::size_t phiBin,
																 double &output1, double &output2, double &output3)> &trafo);

///@}

/** @} */

#endif