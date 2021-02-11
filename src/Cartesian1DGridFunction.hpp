#ifndef CORAS_CARTESIAN_1D_GRID_H
#define CORAS_CARTESIAN_1D_GRID_H

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

/**
 * \addtogroup GRIDFUNCTIONS
 * @{
 */

/**
 * \brief Class implementing a function discretized on a regular one-dimensional Cartesian grid.
 */
class Cartesian1DGridFunction
{
public:
	/**
	 * Constructor instantiating a constant function with value \a value discretized on a regular one-dimensional
	 * Cartesian grid with \a binNumber bins from \a minCoordinate to \a maxCoordinate.
	 */
	Cartesian1DGridFunction(double minCoordinate, double maxCoordinate, std::size_t binNumber,
							double value = 0.0);

	/**
	 * Constructor instantiating the function \a func discretized on a regular one-dimensional Cartesian grid with \a
	 * binNumber bins from \a minCoordinate to \a maxCoordinate.
	 */
	Cartesian1DGridFunction(double minCoordinate, double maxCoordinate, std::size_t binNumber,
							const std::function<double(double x)> &func);

	/**
	 * Constructor instantiating a previously saved Cartesian1DGridFunction from the binary file \a fileName.
	 * 
	 * The file must have been created using Cartesian1DGridFunction::save_object_to_file.
	 */
	Cartesian1DGridFunction(const std::string &fileName);

	/**
	 * Copy-constructor instantiating the same grid function as \a otherGridFunction.
	 */
	Cartesian1DGridFunction(const Cartesian1DGridFunction &otherGridFunction) = default;

	/**
	 * Copy-constructor instantiating the same grid function as \a otherGridFunction and applying the function \a func,
	 * depending on the grid function value, afterwards.
	 */
	Cartesian1DGridFunction(const Cartesian1DGridFunction &otherGridFunction, const std::function<double(double value)> &func);

	/**
	 * Copy-constructor instantiating the same grid function as \a otherGridFunction and applying the function \a func,
	 * depending on the grid function value and the coordinate value, afterwards.
	 */
	Cartesian1DGridFunction(const Cartesian1DGridFunction &otherGridFunction, const std::function<double(double value, double x)> &func);

	/**
	 * Constructor instantiating an empty Cartesian1DGridFunction object, intended to be properly initialized at a later
	 * point. For example, this can be used to choose a constructor based on some conditional statement.
	 */
	Cartesian1DGridFunction();

	/**
	 * Return a reference to the function value at the grid point specified by the bin \a bin.
	 */
	double &value(std::size_t bin);

	/**
	 * Return the function value at the grid point specified by the bin \a bin.
	 */
	double value(std::size_t bin) const;

	/**
	 * Return a reference to the vector of all function values.
	 */
	std::vector<double> &value_vector();

	/**
	 * Return a reference to the \c const vector of all function values.
	 */
	const std::vector<double> &value_vector() const;

	/**
	 * Return the linearly interpolated function value at the coordinate \a x.
	 */
	double interpolate(double x) const;

	/**
	 * Overloads the parentheses operator to return the linearly interpolated function value at the coordinate \a x.
	 *
	 * This provides the same functionality as Cartesian1DGridFunction::interpolate.
	 */
	double operator()(double x) const;

	/**
	 * Return the line element at the grid point specified by the bin \a bin.
	 */
	double line_element(std::size_t bin) const;

	/**
	 * Return the coordinate value of the bin \a bin.
	 */
	double coordinate(std::size_t bin) const;

	/**
	 * Return the vector of all coordinate values.
	 */
	std::vector<double> coordinate_vector() const;

	/**
	 * Return the minimal coordinate value.
	 */
	double minimal_coordinate() const;

	/**
	 * Return the maximal coordinate value.
	 */
	double maximal_coordinate() const;

	/**
	 * Return the bin width.
	 */
	double bin_width() const;

	/**
	 * Return the number of bins.
	 */
	std::size_t bin_number() const;

	/**
	 * Apply the function \a func, depending on the grid function value, to all grid points. If \a parallelized is
	 * set to \c true, this is performed for different grid points in parallel.
	 */
	void apply(const std::function<double(double value)> &func, const bool parallelized = false);

	/**
	 * Apply the function \a func, depending on the grid function value and coordinate, to all grid points. If \a
	 * parallelized is set to \c true, this is performed for different grid points in parallel.
	 */
	void apply(const std::function<double(double value, double x)> &func, const bool parallelized = false);

	/**
	 * Evaluate the function \a func, depending on the grid bin, for all grid points. If \a parallelized is set to
	 * \c true, this is performed for different grid points in parallel.
	 */
	void evaluate_for_all_grid_points(const std::function<void(std::size_t bin)> &func, const bool parallelized = false) const;

	/**
	 * Return the average of the function \a func, depending on the grid function value, over the whole grid length.
	 *
	 * If no argument is given, the average of the grid function itself is returned.
	 */
	double average(const std::function<double(double value)> &func = [](double value) { return value; }) const;

	/**
	 * Return the average of the function \a func, depending on the grid function value and the coordinate value,
	 * over the whole grid length.
	 */
	double average(const std::function<double(double value, double x)> &func) const;

	/**
	 * Write the grid function values into the file \a fileName.
	 */
	void write_to_file(const std::string &fileName) const;

	/**
	 * Save this Cartesian1DGridFunction to the binary file \a fileName.
	 *
	 * This allows to load this object at a later point using Cartesian1DGridFunction::Cartesian1DGridFunction(const
	 * std::string &) or Cartesian1DGridFunction::load_object_from_file.
	 *
	 * For a human-readable text file output use Cartesian1DGridFunction::write_to_file instead.
	 */
	void save_object_to_file(const std::string &fileName) const;

	/**
	 * Assign this Cartesian1DGridFunction to the Cartesian1DGridFunction \a otherGridFunction.
	 */
	Cartesian1DGridFunction &operator=(const Cartesian1DGridFunction &otherGridFunction) = default;

	/**
	 * Assign this Cartesian1DGridFunction to the the elementwise sum of itself and the Cartesian1DGridFunction \a
	 * otherGridFunction, which needs to share the same number of bins as well as minimal and maximal coordinate values.
	 */
	Cartesian1DGridFunction &operator+=(const Cartesian1DGridFunction &otherGridFunction);

	/**
	 * Assign this Cartesian1DGridFunction to the elementwise difference of itself and the Cartesian1DGridFunction \a
	 * otherGridFunction, which needs to share the same number of bins as well as minimal and maximal coordinate values.
	 */
	Cartesian1DGridFunction &operator-=(const Cartesian1DGridFunction &otherGridFunction);

	/**
	 * Assign this Cartesian1DGridFunction to the elementwise product of itself and the Cartesian1DGridFunction \a
	 * otherGridFunction, which needs to share the same number of bins as well as minimal and maximal coordinate values.
	 */
	Cartesian1DGridFunction &operator*=(const Cartesian1DGridFunction &otherGridFunction);

	/**
	 * Assign this Cartesian1DGridFunction to the elementwise quotient of itself and the Cartesian1DGridFunction \a
	 * otherGridFunction, which needs to share the same number of bins as well as minimal and maximal coordinate values.
	 */
	Cartesian1DGridFunction &operator/=(const Cartesian1DGridFunction &otherGridFunction);

	/**
	 * Assign this Cartesian1DGridFunction to the elementwise sum of itself and the scalar \a otherValue.
	 */
	Cartesian1DGridFunction &operator+=(double otherValue);

	/**
	 * Assign this Cartesian1DGridFunction to the elementwise difference of itself and the scalar \a otherValue.
	 */
	Cartesian1DGridFunction &operator-=(double otherValue);

	/**
	 * Assign this Cartesian1DGridFunction to the elementwise product of itself and the scalar \a otherValue.
	 */
	Cartesian1DGridFunction &operator*=(double otherValue);

	/**
	 * Assign this Cartesian1DGridFunction to the elementwise quotient of itself and the scalar \a otherValue.
	 */
	Cartesian1DGridFunction &operator/=(double otherValue);

	/**
	 * Return \c true if this Cartesian1DGridFunction is defined on the same spherical grid as \a otherGridFunction,
	 * otherwise returns \c false.
	 */
	bool has_same_grid_as(const Cartesian1DGridFunction &otherGridFunction) const;

private:
	/**
	 * Minimal coordinate value.
	 */
	double MinCoordinate;

	/**
	 * Maximal coordinate value.
	 */
	double MaxCoordinate;

	/**
	 * Number of bins.
	 */
	std::size_t BinNumber;

	/**
	 * Width of bins.
	 */
	double BinWidth;

	/**
	 * Vector containing all the coordinate values.
	 */
	std::vector<double> Coordinates;

	/**
	 * Vector containing all the grid function values.
	 */
	std::vector<double> FunctionValues;

	/**
	 * Initializes the vector Cartesian1DGridFunction::Coordinates.
	 */
	void initialize_coordinates();

	/**
	 * Loads a previously saved Cartesian1DGridFunction from the binary file \a fileName.
	 * 
	 * The file must have been created using Cartesian1DGridFunction::save_object_to_file.
	 */
	void load_object_from_file(const std::string &fileName);
};

/**
 * \name Cartesian1DGridFunction Binary Operators
 *
 * Binary operators between Cartesian1DGridFunction objects and scalars.
 */
///@{

/**
 * Return the elementwise sum of the Cartesian1DGridFunction objects \a leftGridFunction and \a rightGridFunction, which
 * need to share the same number of bins as well as minimal and maximal coordinate values.
 */
Cartesian1DGridFunction operator+(const Cartesian1DGridFunction &leftGridFunction, const Cartesian1DGridFunction &rightGridFunction);

/**
 * Return the elementwise difference of the Cartesian1DGridFunction objects \a leftGridFunction and \a
 * rightGridFunction, which need to share the same coordinate grid.
 */
Cartesian1DGridFunction operator-(const Cartesian1DGridFunction &leftGridFunction, const Cartesian1DGridFunction &rightGridFunction);

/**
 * Return the elementwise product of the Cartesian1DGridFunction objects \a leftGridFunction and \a rightGridFunction,
 * which need to share the same coordinate grid.
 */
Cartesian1DGridFunction operator*(const Cartesian1DGridFunction &leftGridFunction, const Cartesian1DGridFunction &rightGridFunction);

/**
 * Return the elementwise quotient of the Cartesian1DGridFunction objects \a leftGridFunction and \a rightGridFunction,
 * which need to share the same coordinate grid.
 */
Cartesian1DGridFunction operator/(const Cartesian1DGridFunction &leftGridFunction, const Cartesian1DGridFunction &rightGridFunction);

/**
 * Return the elementwise sum of the Cartesian1DGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
Cartesian1DGridFunction operator+(const Cartesian1DGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise difference of the Cartesian1DGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
Cartesian1DGridFunction operator-(const Cartesian1DGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise product of the Cartesian1DGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
Cartesian1DGridFunction operator*(const Cartesian1DGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise quotient of the Cartesian1DGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
Cartesian1DGridFunction operator/(const Cartesian1DGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise sum of the scalar \a leftValue and the Cartesian1DGridFunction \a rightGridFunction.
 */
Cartesian1DGridFunction operator+(double leftValue, const Cartesian1DGridFunction &rightGridFunction);

/**
 * Return the elementwise difference of the scalar \a leftValue and the Cartesian1DGridFunction \a rightGridFunction.
 */
Cartesian1DGridFunction operator-(double leftValue, const Cartesian1DGridFunction &rightGridFunction);

/**
 * Return the elementwise product of the scalar \a leftValue and the Cartesian1DGridFunction \a rightGridFunction.
 */
Cartesian1DGridFunction operator*(double leftValue, const Cartesian1DGridFunction &rightGridFunction);

/**
 * Return the elementwise quotient of the scalar \a leftValue and the Cartesian1DGridFunction \a rightGridFunction.
 */
Cartesian1DGridFunction operator/(double leftValue, const Cartesian1DGridFunction &rightGridFunction);

///@}

/**
 * \name Cartesian1DGridFunction Transformations
 *
 * General transformations between up to three Cartesian1DGridFunction objects.
 */
///@{

/**
 * Compute the transformation \a trafo of the Cartesian1DGridFunction \a input, and write the result into \a output. \a
 * trafo is a function of the grid bins.
 *
 * The transformation can be performed in-place. Otherwise, \a input is not modified.
 */
void transform_cartesian_1D_grid_functions(const Cartesian1DGridFunction &input,
										   Cartesian1DGridFunction &output,
										   const std::function<void(std::size_t bin,
																	double &output)> &trafo);

/**
 * Compute the transformation \a trafo of the Cartesian1DGridFunction objects \a input1 and \a input2, which need to
 * share the same coordinate grid, and write the result into \a output1 and \a output2. \a trafo is a function of the
 * grid bins.
 *
 * The transformation can be performed (partially) in-place. Otherwise, \a input1 and \a input2 are not modified.
 */
void transform_cartesian_1D_grid_functions(const Cartesian1DGridFunction &input1, const Cartesian1DGridFunction &input2,
										   Cartesian1DGridFunction &output1, Cartesian1DGridFunction &output2,
										   const std::function<void(std::size_t Bin,
																	double &output1, double &output2)> &trafo);

/**
 * Compute the transformation \a trafo of the Cartesian1DGridFunction objects \a input1, \a input2 and \a input3, which
 * need to share the same coordinate grid, and write the result into \a output1, \a output2 and \a output3. \a trafo is
 * a function of the grid bins.
 *
 * The transformation can be performed (partially) in-place. Otherwise, \a input1, \a input2 and \a input3 are not
 * modified.
 */
void transform_cartesian_1D_grid_functions(const Cartesian1DGridFunction &input1, const Cartesian1DGridFunction &input2, const Cartesian1DGridFunction &input3,
										   Cartesian1DGridFunction &output1, Cartesian1DGridFunction &output2, Cartesian1DGridFunction &output3,
										   const std::function<void(std::size_t Bin,
																	double &output1, double &output2, double &output3)> &trafo);

///@}

/** @} */

#endif