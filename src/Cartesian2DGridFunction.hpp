#ifndef CORAS_CARTESIAN_2D_GRID_H
#define CORAS_CARTESIAN_2D_GRID_H

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

/**
 * \addtogroup GRIDFUNCTIONS
 * @{
 */

/**
 * \brief Class implementing a function discretized on a regular two-dimensional Cartesian grid.
 */
class Cartesian2DGridFunction
{
public:
	/**
	 * Constructor instantiating a constant function with value \a value discretized on a regular two-dimensional
	 * Cartesian grid with \a xBinNumber bins from \a minXCoordinate to \a maxXCoordinate on the x-axis and \a
	 * yBinNumber bins from \a minYCoordinate to \a maxYCoordinate on the y-axis.
	 */
	Cartesian2DGridFunction(double minXCoordinate, double maxXCoordinate, std::size_t xBinNumber,
							double minYCoordinate, double maxYCoordinate, std::size_t yBinNumber,
							double value = 0.0);

	/**
	 * Constructor instantiating the function \a func discretized on a regular two-dimensional Cartesian grid with \a
	 * xBinNumber bins from \a minXCoordinate to \a maxXCoordinate on the x-axis and \a yBinNumber bins from \a
	 * minYCoordinate to \a maxYCoordinate on the y-axis.
	 */
	Cartesian2DGridFunction(double minXCoordinate, double maxXCoordinate, std::size_t xBinNumber,
							double minYCoordinate, double maxYCoordinate, std::size_t yBinNumber,
							const std::function<double(double x, double y)> &func);

	/**
	 * Constructor instantiating a previously saved Cartesian2DGridFunction from the binary file \a fileName.
	 * 
	 * The file must have been created using Cartesian2DGridFunction::save_object_to_file.
	 */
	Cartesian2DGridFunction(const std::string &fileName);

	/**
	 * Copy-constructor instantiating the same grid function as \a otherGridFunction.
	 */
	Cartesian2DGridFunction(const Cartesian2DGridFunction &otherGridFunction) = default;

	/**
	 * Copy-constructor instantiating the same grid function as \a otherGridFunction and applying the function \a func,
	 * depending on the grid function value, afterwards.
	 */
	Cartesian2DGridFunction(const Cartesian2DGridFunction &otherGridFunction, const std::function<double(double value)> &func);

	/**
	 * Copy-constructor instantiating the same grid function as \a otherGridFunction and applying the function \a func,
	 * depending on the grid function value and coordinates, afterwards.
	 */
	Cartesian2DGridFunction(const Cartesian2DGridFunction &otherGridFunction, const std::function<double(double value, double x, double y)> &func);

	/**
	 * Constructor instantiating an empty Cartesian2DGridFunction object, intended to be properly initialized at a later
	 * point. For example, this can be used to choose a constructor based on some conditional statement.
	 */
	Cartesian2DGridFunction();

	/**
	 * Return a reference to the function value at the grid point specified by the bins \a xBin and \a yBin.
	 */
	double &value(std::size_t xBin, std::size_t yBin);

	/**
	 * Return the function value at the grid point specified by the bins \a xBin and \a yBin.
	 */
	double value(std::size_t xBin, std::size_t yBin) const;

	/**
	 * Return a reference to the vector of all function values.
	 */
	std::vector<double> &value_vector();

	/**
	 * Return a reference to the \c const vector of all function values.
	 */
	const std::vector<double> &value_vector() const;

	/**
	 * Return the linearly interpolated function value at the coordinates \a  x and \a y.
	 */
	double interpolate(double x, double y) const;

	/**
	 * Overloads the parentheses operator to return the linearly interpolated function value at the coordinates \a  x
	 * and \a y.
	 *
	 * This provides the same functionality as Cartesian2DGridFunction::interpolate.
	 */
	double operator()(double x, double y) const;

	/**
	 * Return the volume element at the grid point specified by the bins \a xBin and \a yBin.
	 */
	double volume_element(std::size_t xBin, std::size_t yBin) const;

	/**
	 * Return the x-coordinate value of the x-axis bin \a xBin.
	 */
	double x_coordinate(std::size_t xBin) const;

	/**
	 * Return the y-coordinate value of the y-axis bin \a yBin.
	 */
	double y_coordinate(std::size_t yBin) const;

	/**
	 * Return the vector of all x-coordinate values.
	 */
	std::vector<double> x_coordinate_vector() const;

	/**
	 * Return the vector of all y-coordinate values.
	 */
	std::vector<double> y_coordinate_vector() const;

	/**
	 * Return the minimal x-coordinate value.
	 */
	double minimal_x_coordinate() const;

	/**
	 * Return the maximal x-coordinate value.
	 */
	double maximal_x_coordinate() const;

	/**
	 * Return the minimal y-coordinate value.
	 */
	double minimal_y_coordinate() const;

	/**
	 * Return the maximal y-coordinate value.
	 */
	double maximal_y_coordinate() const;

	/**
	 * Return the x-axis bin width.
	 */
	double x_bin_width() const;

	/**
	 * Return the y-axis bin width.
	 */
	double y_bin_width() const;

	/**
	 * Return the number of x-axis bins.
	 */
	std::size_t x_bin_number() const;

	/**
	 * Return the number of y-axis bins.
	 */
	std::size_t y_bin_number() const;

	/**
	 * Apply the function \a func, depending on the grid function value, to all grid points. If \a parallelized is
	 * set to \c true, this is performed for different grid points in parallel.
	 */
	void apply(const std::function<double(double value)> &func, const bool parallelized = false);

	/**
	 * Apply the function \a func, depending on the grid function value and coordinates, to all grid points. If \a
	 * parallelized is set to \c true, this is performed for different grid points in parallel.
	 */
	void apply(const std::function<double(double value, double x, double y)> &func, const bool parallelized = false);

	/**
	 * Evaluate the function \a func, depending on the grid bins, for all grid points. If \a parallelized is set to
	 * \c true, this is performed for different grid points in parallel.
	 */
	void evaluate_for_all_grid_points(const std::function<void(std::size_t xBin, std::size_t yBin)> &func, const bool parallelized = false) const;

	/**
	 * Return the average of the function \a func, depending on the grid function value, over the whole grid area.
	 *
	 * If no argument is given, the volume average of the grid function itself is returned.
	 */
	double average(const std::function<double(double value)> &func = [](double value) { return value; }) const;

	/**
	 * Return the average of the function \a func, depending on the grid function value and the coordinates, over
	 * the whole grid area.
	 */
	double average(const std::function<double(double value, double x, double y)> &func) const;

	/**
	 * Write the grid function values into the file \a fileName.
	 */
	void write_to_file(const std::string &fileName) const;

	/**
	 * Save this Cartesian2DGridFunction to the binary file \a fileName.
	 *
	 * This allows to load this object at a later point using Cartesian2DGridFunction::Cartesian2DGridFunction(const
	 * std::string &) or Cartesian2DGridFunction::load_object_from_file.
	 *
	 * For a human-readable text file output use Cartesian2DGridFunction::write_to_file instead.
	 */
	void save_object_to_file(const std::string &fileName) const;

	/**
	 * Assign this Cartesian2DGridFunction to the Cartesian2DGridFunction \a otherGridFunction.
	 */
	Cartesian2DGridFunction &operator=(const Cartesian2DGridFunction &otherGridFunction) = default;

	/**
	 * Assign this Cartesian2DGridFunction to the the elementwise sum of itself and the Cartesian2DGridFunction \a
	 * otherGridFunction, which needs to share the same numbers of bins as well as minimal and maximal coordinate
	 * values.
	 */
	Cartesian2DGridFunction &operator+=(const Cartesian2DGridFunction &otherGridFunction);

	/**
	 * Assign this Cartesian2DGridFunction to the elementwise difference of itself and the Cartesian2DGridFunction \a
	 * otherGridFunction, which needs to share the same numbers of bins as well as minimal and maximal coordinate
	 * values.
	 */
	Cartesian2DGridFunction &operator-=(const Cartesian2DGridFunction &otherGridFunction);

	/**
	 * Assign this Cartesian2DGridFunction to the elementwise product of itself and the Cartesian2DGridFunction \a
	 * otherGridFunction, which needs to share the same numbers of bins as well as minimal and maximal coordinate
	 * values.
	 */
	Cartesian2DGridFunction &operator*=(const Cartesian2DGridFunction &otherGridFunction);

	/**
	 * Assign this Cartesian2DGridFunction to the elementwise quotient of itself and the Cartesian2DGridFunction \a
	 * otherGridFunction, which needs to share the same numbers of bins as well as minimal and maximal coordinate
	 * values.
	 */
	Cartesian2DGridFunction &operator/=(const Cartesian2DGridFunction &otherGridFunction);

	/**
	 * Assign this Cartesian2DGridFunction to the elementwise sum of itself and the scalar \a otherValue.
	 */
	Cartesian2DGridFunction &operator+=(double otherValue);

	/**
	 * Assign this Cartesian2DGridFunction to the elementwise difference of itself and the scalar \a otherValue.
	 */
	Cartesian2DGridFunction &operator-=(double otherValue);

	/**
	 * Assign this Cartesian2DGridFunction to the elementwise product of itself and the scalar \a otherValue.
	 */
	Cartesian2DGridFunction &operator*=(double otherValue);

	/**
	 * Assign this Cartesian2DGridFunction to the elementwise quotient of itself and the scalar \a otherValue.
	 */
	Cartesian2DGridFunction &operator/=(double otherValue);

	/**
	 * Return \c true if this Cartesian2DGridFunction is defined on the same spherical grid as \a otherGridFunction,
	 * otherwise returns \c false.
	 */
	bool has_same_grid_as(const Cartesian2DGridFunction &otherGridFunction) const;

private:
	/**
	 * Minimal x-coordinate value.
	 */
	double MinXCoordinate;

	/**
	 * Maximal x-coordinate value.
	 */
	double MaxXCoordinate;

	/**
	 * Minimal y-coordinate value.
	 */
	double MinYCoordinate;

	/**
	 * Maximal y-coordinate value.
	 */
	double MaxYCoordinate;

	/**
	 * Number of x-axis bins.
	 */
	std::size_t XBinNumber;

	/**
	 * Number of y-axis bins.
	 */
	std::size_t YBinNumber;

	/**
	 * Width of x-axis bins.
	 */
	double XBinWidth;

	/**
	 * Width of y-axis bins.
	 */
	double YBinWidth;

	/**
	 * Vector containing all the x-coordinate values.
	 */
	std::vector<double> XCoordinates;

	/**
	 * Vector containing all the y-coordinate values.
	 */
	std::vector<double> YCoordinates;

	/**
	 * Vector containing all the grid function values.
	 */
	std::vector<double> FunctionValues;

	/**
	 * Initializes the vectors Cartesian2DGridFunction::XCoordinates and Cartesian2DGridFunction::YCoordinates.
	 */
	void initialize_coordinates();

	/**
	 * Loads a previously saved Cartesian2DGridFunction from the binary file \a fileName.
	 * 
	 * The file must have been created using Cartesian2DGridFunction::save_object_to_file.
	 */
	void load_object_from_file(const std::string &fileName);
};

/**
 * \name Cartesian2DGridFunction Binary Operators
 *
 * Binary operators between Cartesian2DGridFunction objects and scalars.
 */
///@{

/**
 * Return the elementwise sum of the Cartesian2DGridFunction objects \a leftGridFunction and \a rightGridFunction, which
 * need to share the same numbers of bins as well as minimal and maximal coordinate values.
 */
Cartesian2DGridFunction operator+(const Cartesian2DGridFunction &leftGridFunction, const Cartesian2DGridFunction &rightGridFunction);

/**
 * Return the elementwise difference of the Cartesian2DGridFunction objects \a leftGridFunction and \a
 * rightGridFunction, which need to share the same coordinate grid.
 */
Cartesian2DGridFunction operator-(const Cartesian2DGridFunction &leftGridFunction, const Cartesian2DGridFunction &rightGridFunction);

/**
 * Return the elementwise product of the Cartesian2DGridFunction objects \a leftGridFunction and \a rightGridFunction,
 * which need to share the same coordinate grid.
 */
Cartesian2DGridFunction operator*(const Cartesian2DGridFunction &leftGridFunction, const Cartesian2DGridFunction &rightGridFunction);

/**
 * Return the elementwise quotient of the Cartesian2DGridFunction objects \a leftGridFunction and \a rightGridFunction,
 * which need to share the same coordinate grid.
 */
Cartesian2DGridFunction operator/(const Cartesian2DGridFunction &leftGridFunction, const Cartesian2DGridFunction &rightGridFunction);

/**
 * Return the elementwise sum of the Cartesian2DGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
Cartesian2DGridFunction operator+(const Cartesian2DGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise difference of the Cartesian2DGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
Cartesian2DGridFunction operator-(const Cartesian2DGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise product of the Cartesian2DGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
Cartesian2DGridFunction operator*(const Cartesian2DGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise quotient of the Cartesian2DGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
Cartesian2DGridFunction operator/(const Cartesian2DGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise sum of the scalar \a leftValue and the Cartesian2DGridFunction \a rightGridFunction.
 */
Cartesian2DGridFunction operator+(double leftValue, const Cartesian2DGridFunction &rightGridFunction);

/**
 * Return the elementwise difference of the scalar \a leftValue and the Cartesian2DGridFunction \a rightGridFunction.
 */
Cartesian2DGridFunction operator-(double leftValue, const Cartesian2DGridFunction &rightGridFunction);

/**
 * Return the elementwise product of the scalar \a leftValue and the Cartesian2DGridFunction \a rightGridFunction.
 */
Cartesian2DGridFunction operator*(double leftValue, const Cartesian2DGridFunction &rightGridFunction);

/**
 * Return the elementwise quotient of the scalar \a leftValue and the Cartesian2DGridFunction \a rightGridFunction.
 */
Cartesian2DGridFunction operator/(double leftValue, const Cartesian2DGridFunction &rightGridFunction);

///@}

/**
 * \name Cartesian1DGridFunction Transformations
 *
 * General transformations between up to three Cartesian1DGridFunction objects.
 */
///@{

/**
 * Compute the transformation \a trafo of the Cartesian2DGridFunction \a input, and write the result into \a output. \a
 * trafo is a function of the grid bins.
 *
 * The transformation can be performed in-place. Otherwise, \a input is not modified.
 */
void transform_cartesian_2D_grid_functions(const Cartesian2DGridFunction &input,
										   Cartesian2DGridFunction &output,
										   const std::function<void(std::size_t xBin, std::size_t yBin,
																	double &output)> &trafo);

/**
 * Compute the transformation \a trafo of the Cartesian2DGridFunction objects \a input1 and \a input2, which need to
 * share the same coordinate grid, and write the result into \a output1 and \a output2. \a trafo is a function of the
 * grid bins.
 *
 * The transformation can be performed (partially) in-place. Otherwise, \a input1 and \a input2 are not modified.
 */
void transform_cartesian_2D_grid_functions(const Cartesian2DGridFunction &input1, const Cartesian2DGridFunction &input2,
										   Cartesian2DGridFunction &output1, Cartesian2DGridFunction &output2,
										   const std::function<void(std::size_t xBin, std::size_t yBin,
																	double &output1, double &output2)> &trafo);

/**
 * Compute the transformation \a trafo of the Cartesian2DGridFunction objects \a input1, \a input2 and \a input3, which
 * need to share the same coordinate grid, and write the result into \a output1, \a output2 and \a output3. \a trafo is
 * a function of the grid bins.
 *
 * The transformation can be performed (partially) in-place. Otherwise, \a input1, \a input2 and \a input3 are not
 * modified.
 */
void transform_cartesian_2D_grid_functions(const Cartesian2DGridFunction &input1, const Cartesian2DGridFunction &input2, const Cartesian2DGridFunction &input3,
										   Cartesian2DGridFunction &output1, Cartesian2DGridFunction &output2, Cartesian2DGridFunction &output3,
										   const std::function<void(std::size_t xBin, std::size_t yBin,
																	double &output1, double &output2, double &output3)> &trafo);

///@}

/** @} */

#endif