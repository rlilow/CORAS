#ifndef CORAS_ANGULAR_GRID_FUNCTION_H
#define CORAS_ANGULAR_GRID_FUNCTION_H

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

/**
 * \addtogroup GRIDFUNCTIONS
 * @{
 */

/**
 * \brief Class implementing a function discretized on a regular two-dimensional angular grid.
 */
class AngularGridFunction
{
public:
	/**
	 * Constructor instantiating a constant function with value \a value discretized on a regular two-dimensional
	 * angular grid with \a thetaBinNumber and \a phiBinNumber bins on the theta and phi axes, respectively.
	 */
	AngularGridFunction(std::size_t thetaBinNumber, std::size_t phiBinNumber,
						double value = 0.0);

	/**
	 * Constructor instantiating the function \a func discretized on a regular two-dimensional angular grid with \a
	 * thetaBinNumber and \a phiBinNumber bins on the theta and phi axes, respectively.
	 */
	AngularGridFunction(std::size_t thetaBinNumber, std::size_t phiBinNumber,
						const std::function<double(double theta, double phi)> &func);

	/**
	 * Constructor instantiating a previously saved AngularGridFunction from the binary file \a fileName.
	 * 
	 * The file must have been created using AngularGridFunction::save_object_to_file.
	 */
	AngularGridFunction(const std::string &fileName);

	/**
	 * Copy-constructor instantiating the same grid function as \a otherAngularGridFunction.
	 */
	AngularGridFunction(const AngularGridFunction &otherAngularGridFunction) = default;

	/**
	 * Copy-constructor instantiating the same grid function as \a otherAngularGridFunction and applying the function \a
	 * func, depending on the grid function value, afterwards.
	 */
	AngularGridFunction(const AngularGridFunction &otherAngularGridFunction, const std::function<double(double value)> &func);

	/**
	 * Copy-constructor instantiating the same grid function as \a otherAngularGridFunction and applying the function \a
	 * func, depending on the grid function value and coordinates, afterwards.
	 */
	AngularGridFunction(const AngularGridFunction &otherAngularGridFunction, const std::function<double(double value, double theta, double phi)> &func);

	/**
	 * Constructor instantiating an empty AngularGridFunction object, intended to be properly initialized at a later
	 * point. For example, this can be used to choose a constructor based on some conditional statement.
	 */
	AngularGridFunction();

	/**
	 * Return a reference to the function value at the grid point specified by \a thetaBin and \a phiBin.
	 */
	double &value(std::size_t thetaBin, std::size_t phiBin);

	/**
	 * Return the function value at the grid point specified by \a thetaBin and \a phiBin.
	 */
	double value(std::size_t thetaBin, std::size_t phiBin) const;

	/**
	 * Return a reference to the vector of all function values.
	 */
	std::vector<double> &value_vector();

	/**
	 * Return a reference to the \c const vector of all function values.
	 */
	const std::vector<double> &value_vector() const;

	/**
	 * Return the linearly interpolated function value at the coordinates \a theta and \a phi.
	 */
	double interpolate(double theta, double phi) const;

	/**
	 * Overloads the parentheses operator to return the linearly interpolated function value at the coordinates \a theta
	 * and \a phi.
	 *
	 * This provides the same functionality as AngularGridFunction::interpolate.
	 */
	double operator()(double theta, double phi) const;

	/**
	 * Return the solid angle element at the grid point specified by the theta bin \a thetaBin.
	 */
	double solid_angle_element(std::size_t thetaBin) const;

	/**
	 * Return the theta coordinate value at the theta-axis point specified by \a thetaBin.
	 */
	double theta_coordinate(std::size_t thetaBin) const;

	/**
	 * Return the phi coordinate value at the phi-axis point specified by \a phiBin.
	 */
	double phi_coordinate(std::size_t phiBin) const;

	/**
	 * Return the vector of all theta-coordinate values.
	 */
	std::vector<double> theta_coordinate_vector() const;

	/**
	 * Return the vector of all phi-coordinate values.
	 */
	std::vector<double> phi_coordinate_vector() const;

	/**
	 * Return the theta-axis bin width.
	 */
	double theta_bin_width() const;

	/**
	 * Return the theta-axis bin width.
	 */
	double phi_bin_width() const;

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
	void apply(const std::function<double(double value, double theta, double phi)> &func, const bool parallelized = false);

	/**
	 * Evaluate the function \a func, depending on the grid bins, for all grid points. If \a parallelized is set to
	 * \c true, this is performed for different grid points in parallel.
	 */
	void evaluate_for_all_grid_points(const std::function<void(std::size_t thetaBin, std::size_t phiBin)> &func, const bool parallelized = false) const;

	/**
	 * Return the average of the function \a func, depending on the grid function value, over the whole solid angle.
	 *
	 * If no argument is given, the average of the grid function itself is returned.
	 */
	double average(const std::function<double(double value)> &func = [](double value) { return value; }) const;

	/**
	 * Return the average of the function \a func, depending on the grid function value and the coordinates, over
	 * the whole solid angle.
	 */
	double average(const std::function<double(double value, double theta, double phi)> &func) const;

	/**
	 * Write the grid function values into the text file \a fileName.
	 */
	void write_to_file(const std::string &fileName) const;

	/**
	 * Save this AngularGridFunction to the binary file \a fileName.
	 *
	 * This allows to load this object at a later point using AngularGridFunction::AngularGridFunction(const std::string
	 * &) or AngularGridFunction::load_object_from_file.
	 *
	 * For a human-readable text file output use AngularGridFunction::write_to_file instead.
	 */
	void save_object_to_file(const std::string &fileName) const;

	/**
	 * Assign this AngularGridFunction to the AngularGridFunction \a otherGridFunction.
	 */
	AngularGridFunction &operator=(const AngularGridFunction &otherGridFunction) = default;

	/**
	 * Assign this AngularGridFunction to the the elementwise sum of itself and the AngularGridFunction \a
	 * otherGridFunction, which needs to share the same numbers of bins as well as maximal radius.
	 */
	AngularGridFunction &operator+=(const AngularGridFunction &otherGridFunction);

	/**
	 * Assign this AngularGridFunction to the elementwise difference of itself and the AngularGridFunction \a
	 * otherGridFunction, which needs to share the same numbers of bins as well as maximal radius.
	 */
	AngularGridFunction &operator-=(const AngularGridFunction &otherGridFunction);

	/**
	 * Assign this AngularGridFunction to the elementwise product of itself and the AngularGridFunction \a
	 * otherGridFunction, which needs to share the same numbers of bins as well as maximal radius.
	 */
	AngularGridFunction &operator*=(const AngularGridFunction &otherGridFunction);

	/**
	 * Assign this AngularGridFunction to the elementwise quotient of itself and the AngularGridFunction \a
	 * otherGridFunction, which needs to share the same numbers of bins as well as maximal radius.
	 */
	AngularGridFunction &operator/=(const AngularGridFunction &otherGridFunction);

	/**
	 * Assign this AngularGridFunction to the elementwise sum of itself and the scalar \a otherValue.
	 */
	AngularGridFunction &operator+=(double otherValue);

	/**
	 * Assign this AngularGridFunction to the elementwise difference of itself and the scalar \a otherValue.
	 */
	AngularGridFunction &operator-=(double otherValue);

	/**
	 * Assign this AngularGridFunction to the elementwise product of itself and the scalar \a otherValue.
	 */
	AngularGridFunction &operator*=(double otherValue);

	/**
	 * Assign this AngularGridFunction to the elementwise quotient of itself and the scalar \a otherValue.
	 */
	AngularGridFunction &operator/=(double otherValue);

	/**
	 * Return \c true if this AngularGridFunction is defined on the same angular grid as \a otherGridFunction,
	 * otherwise returns \c false.
	 */
	bool has_same_grid_as(const AngularGridFunction &otherGridFunction) const;

private:
	/**
	 * Number of theta-axis bins.
	 */
	std::size_t ThetaBinNumber;

	/**
	 * Number of phi-axis bins.
	 */
	std::size_t PhiBinNumber;

	/**
	 * Width of theta-axis bins.
	 */
	double ThetaBinWidth;

	/**
	 * Width of phi-axis bins.
	 */
	double PhiBinWidth;

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
	 * Initializes the vectors AngularGridFunction::ThetaCoordinates and AngularGridFunction::PhiCoordinates.
	 */
	void initialize_coordinates();

	/**
	 * Loads a previously saved AngularGridFunction from the binary file \a fileName.
	 * 
	 * The file must have been created using AngularGridFunction::save_object_to_file.
	 */
	void load_object_from_file(const std::string &fileName);
};

/**
 * \name AngularGridFunction Binary Operators
 *
 * Binary operators between AngularGridFunction objects and scalars.
 */
///@{

/**
 * Return the elementwise sum of the AngularGridFunction objects \a leftGridFunction and \a rightGridFunction, which
 * need to share the same numbers of bins as well as maximal radius.
 */
AngularGridFunction operator+(const AngularGridFunction &leftGridFunction, const AngularGridFunction &rightGridFunction);

/**
 * Return the elementwise difference of the AngularGridFunction objects \a leftGridFunction and \a rightGridFunction,
 * which need to share the same coordinate grid.
 */
AngularGridFunction operator-(const AngularGridFunction &leftGridFunction, const AngularGridFunction &rightGridFunction);

/**
 * Return the elementwise product of the AngularGridFunction objects \a leftGridFunction and \a rightGridFunction, which
 * need to share the same numbers of bins as well as maximal radius.
 */
AngularGridFunction operator*(const AngularGridFunction &leftGridFunction, const AngularGridFunction &rightGridFunction);

/**
 * Return the elementwise quotient of the AngularGridFunction objects \a leftGridFunction and \a rightGridFunction,
 * which need to share the same coordinate grid.
 */
AngularGridFunction operator/(const AngularGridFunction &leftGridFunction, const AngularGridFunction &rightGridFunction);

/**
 * Return the elementwise sum of the AngularGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
AngularGridFunction operator+(const AngularGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise difference of the AngularGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
AngularGridFunction operator-(const AngularGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise product of the AngularGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
AngularGridFunction operator*(const AngularGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise quotient of the AngularGridFunction \a leftGridFunction and the scalar \a rightValue.
 */
AngularGridFunction operator/(const AngularGridFunction &leftGridFunction, double rightValue);

/**
 * Return the elementwise sum of the scalar \a leftValue and the AngularGridFunction \a rightGridFunction.
 */
AngularGridFunction operator+(double leftValue, const AngularGridFunction &rightGridFunction);

/**
 * Return the elementwise difference of the scalar \a leftValue and the AngularGridFunction \a rightGridFunction.
 */
AngularGridFunction operator-(double leftValue, const AngularGridFunction &rightGridFunction);

/**
 * Return the elementwise product of the scalar \a leftValue and the AngularGridFunction \a rightGridFunction.
 */
AngularGridFunction operator*(double leftValue, const AngularGridFunction &rightGridFunction);

/**
 * Return the elementwise quotient of the scalar \a leftValue and the AngularGridFunction \a rightGridFunction.
 */
AngularGridFunction operator/(double leftValue, const AngularGridFunction &rightGridFunction);

///@}

/**
 * \name AngularGridFunction Transformations
 *
 * General transformations between up to three AngularGridFunction objects.
 */
///@{

/**
 * Compute the transformation \a trafo of the AngularGridFunction \a input, and write the result into \a output. \a
 * trafo is a function of the grid bins.
 *
 * The transformation can be performed in-place. Otherwise, \a input is not modified.
 */
void transform_angular_grid_functions(const AngularGridFunction &input,
									  AngularGridFunction &output,
									  const std::function<void(std::size_t thetaBin, std::size_t phiBin,
															   double &output)> &trafo);

/**
 * Compute the transformation \a trafo of the AngularGridFunction objects \a input1 and \a input2, which need to share
 * the same coordinate grid, and write the result into \a output1 and \a output2. \a trafo is a function of the grid
 * bins.
 *
 * The transformation can be performed (partially) in-place. Otherwise, \a input1 and \a input2 are not modified.
 */
void transform_angular_grid_functions(const AngularGridFunction &input1, const AngularGridFunction &input2,
									  AngularGridFunction &output1, AngularGridFunction &output2,
									  const std::function<void(std::size_t thetaBin, std::size_t phiBin,
															   double &output1, double &output2)> &trafo);

/**
 * Compute the transformation \a trafo of the AngularGridFunction objects \a input1, \a input2 and \a input3, which need
 * to share the same coordinate grid, and write the result into \a output1, \a output2 and \a output3. \a trafo is a
 * function of the grid bins.
 *
 * The transformation can be performed (partially) in-place. Otherwise, \a input1, \a input2 and \a input3 are not
 * modified.
 */
void transform_angular_grid_functions(const AngularGridFunction &input1, const AngularGridFunction &input2, const AngularGridFunction &input3,
									  AngularGridFunction &output1, AngularGridFunction &output2, AngularGridFunction &output3,
									  const std::function<void(std::size_t thetaBin, std::size_t phiBin,
															   double &output1, double &output2, double &output3)> &trafo);

///@}

/** @} */

#endif