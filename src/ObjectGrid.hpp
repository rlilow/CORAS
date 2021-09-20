#ifndef CORAS_OBJECT_GRID_H
#define CORAS_OBJECT_GRID_H

#include <functional>
#include <vector>

#include "Cartesian3DGridFunction.hpp"

/**
 * \addtogroup AUXILIARY
 * @{
 */

/**
 * \brief Class implementing the association of a set of objects in 3D space with cells in a regular Cartesian grid.
 *
 * It allows the faster operation on objects within a given distance by only iterating through the objects in the
 * relevant grid cells. The grid can have either periodic or non-periodic boundary conditions.
 */
class ObjectGrid
{
public:
    ObjectGrid(
        double gridLength, std::size_t gridBinNumber,
        double centralXCoordinate, double centralYCoordinate, double centralZCoordinate,
        const std::vector<double> &xCoordinates, const std::vector<double> &yCoordinates, const std::vector<double> &zCoordinates,
        bool usePeriodicBoundaryConditions = true,
        const std::function<bool(std::size_t index)> &condition = [](std::size_t)
        { return true; });

    std::size_t object_number() const;
    bool uses_periodic_boundary_conditions() const;

    double x_coordinate(std::size_t index) const;
    double y_coordinate(std::size_t index) const;
    double z_coordinate(std::size_t index) const;

    double x_distance(std::size_t index1, std::size_t index2) const;
    double y_distance(std::size_t index1, std::size_t index2) const;
    double z_distance(std::size_t index1, std::size_t index2) const;

    double coordinate_distance(double coordinateDifference) const;

    double distance(double xDifference, double yDifference, double zDifference) const;
    double distance(std::size_t index1, std::size_t index2) const;

    void evaluate_for_all_objects(const std::function<void(std::size_t index)> &function) const;

    void evaluate_for_all_objects_within_distance(double xCenter, double yCenter, double zCenter, double maxDistance,
                                                  const std::function<void(double relativeX, double relativeY, double relativeZ, double distance, std::size_t index)> &function) const;

    void evaluate_for_all_objects_within_distance(std::size_t centralIndex, double maxDistance,
                                                  const std::function<void(double relativeX, double relativeY, double relativeZ, double distance, std::size_t index)> &function) const;

private:
    double GridLength;
    std::size_t GridBinNumber;
    double GridBinWidth;
    double XCoordinateShift;
    double YCoordinateShift;
    double ZCoordinateShift;
    const std::vector<double> *XCoordinates;
    const std::vector<double> *YCoordinates;
    const std::vector<double> *ZCoordinates;
    bool UsePeriodicBoundaryConditions;
    std::vector<std::vector<std::vector<std::vector<std::size_t>>>> IndexGrid;
    std::size_t ObjectNumber;

    void initialize_grid(const std::function<bool(std::size_t index)> &condition);

    std::size_t bin(double shiftedCoordinate) const;

    double shift_x_coordinate(double x) const;
    double shift_y_coordinate(double y) const;
    double shift_z_coordinate(double z) const;

    void evaluate_for_all_objects_in_sub_grid(double xCenter, double yCenter, double zCenter, double minDistance,
                                              const std::function<void(std::size_t index)> &function) const;
};

/** @} */

#endif