#ifndef CORAS_TRANSFORMATIONS_H
#define CORAS_TRANSFORMATIONS_H

#include <cmath>

#include "Cartesian3DGridFunction.hpp"
#include "SphericalGridFunction.hpp"

/**
 * \defgroup AUXILIARY Auxiliary
 * 
 * \brief Auxiliary tools.
 * @{
 */

/**
 * \name Coordinate Transformations
 *
 * Transformations of coordinates, scalar fields and vector fields between different coordinate systems and reference
 * frames.
 */
///@{

/**
 * Radians per degree.
 */
constexpr double RadiansPerDegree = M_PI / 180.0;

/**
 * Compute the cross product between two vectors with Cartesian components \a x1Input, \a y1Input, \a z1Input and \a
 * x2Input, \a y2Input, \a z2Input, and write the result into \a xOutput, \a yOutput, \a zOutput.
 */
void cartesian_cross_product(double x1Input, double y1Input, double z1Input,
                             double x2Input, double y2Input, double z2Input,
                             double &xOutput, double &yOutput, double &zOutput);

/**
 * Return the scalar product between two vectors with Cartesian components \a x1, \a y1, \a z1 and \a x2, \a y2, \a z2.
 */
double cartesian_scalar_product(double x1, double y1, double z1,
                                double x2, double y2, double z2);

/**
 * Return the norm of a vector with Cartesian components \a x, \a y, \a z.
 */
double cartesian_norm(double x, double y, double z);

/**
 * Compute the cross product between two vectors with Spherical components \a radius1Input, \a theta1Input, \a phi1Input
 * and \a radius2Input, \a theta2Input, \a phi2Input, and write the result into \a radiusOutput, \a thetaOutput, \a
 * phiOutput.
 */
void spherical_cross_product(double radius1Input, double theta1Input, double phi1Input,
                             double radius2Input, double theta2Input, double phi2Input,
                             double &radiusOutput, double &thetaOutput, double &phiOutput);

/**
 * Return the scalar product between two vectors with Spherical components \a radius1, \a theta1, \a phi1 and \a
 * radius2, \a theta2, \a phi2.
 */
double spherical_scalar_product(double radius1, double theta1, double phi1,
                                double radius2, double theta2, double phi2);

/**
 * Return the norm of a vector with Spherical components \a radius, \a theta, \a phi.
 */
double spherical_norm(double radius, double theta, double phi);

/**
 * Compute the cross product between two vectors with celestial components \a radius1Input, \a latitude1Input, \a
 * longitude1Input and \a radius2Input, \a latitude2Input, \a longitude2Input, and write the result into \a
 * radiusOutput, \a latitudeOutput, \a longitudeOutput.
 */
void celestial_cross_product(double radius1Input, double latitude1Input, double longitude1Input,
                             double radius2Input, double latitude2Input, double longitude2Input,
                             double &radiusOutput, double &latitudeOutput, double &longitudeOutput);

/**
 * Return the scalar product between two vectors with celestial components \a radius1, \a latitude1, \a longitude1 and
 * \a radius2, \a latitude2, \a longitude2.
 */
double celestial_scalar_product(double radius1, double latitude1, double longitude1,
                                double radius2, double latitude2, double longitude2);

/**
 * Return the norm of a vector with celestial components \a radius, \a latitude, \a longitude.
 */
double celestial_norm(double radius, double latitude, double longitude);

/**
 * Transform the spherical coordinates \a radius, \a theta, \a phi to the Cartesian coordinates \a x, \a y, \a z.
 */
void transform_spherical_to_cartesian_coordinates(double radius, double theta, double phi,
                                                  double &x, double &y, double &z);

/**
 * Transform the Cartesian coordinates \a x, \a y, \a z to the spherical coordinates \a radius, \a theta, \a phi.
 */
void transform_cartesian_to_spherical_coordinates(double x, double y, double z,
                                                  double &radius, double &theta, double &phi);

/**
 * Transform the celestial angular coordinates \a latitude, \a longitude to the spherical angular coordinates \a theta,
 * \a phi.
 */
void transform_celestial_to_spherical_coordinates(double latitude, double longitude,
                                                  double &theta, double &phi);

/**
 * Transform the spherical angular coordinates \a theta, \a phi to the celestial angular coordinates \a latitude, \a
 * longitude.
 */
void transform_spherical_to_celestial_coordinates(double theta, double phi,
                                                  double &latitude, double &longitude);

/**
 * Transform the celestial coordinates \a radius, \a latitude, \a longitude to the Cartesian coordinates \a x, \a y, \a
 * z.
 */
void transform_celestial_to_cartesian_coordinates(double radius, double latitude, double longitude,
                                                  double &x, double &y, double &z);

/**
 * Transform the Cartesian coordinates \a x, \a y, \a z to the celestial coordinates \a radius, \a latitude, \a
 * longitude.
 */
void transform_cartesian_to_celestial_coordinates(double x, double y, double z,
                                                  double &radius, double &latitude, double &longitude);

/**
 * Transform the values of the spherical vector field components \a radialComponent, \a thetaComponent, \a phiComponent,
 * evaluated at the angular coordinates \a theta, \a phi, to the Cartesian components \a xComponent, \a yComponent, \a
 * zComponent.
 */
void transform_spherical_to_cartesian_vector_field_value(double theta, double phi,
                                                         double radialComponent, double thetaComponent, double phiComponent,
                                                         double &xComponent, double &yComponent, double &zComponent);

/**
 * Transform the values of the Cartesian vector field components \a xComponent, \a yComponent, \a zComponent, evaluated
 * at the angular coordinates \a theta, \a phi, to the Cartesian components \a radialComponent, \a thetaComponent, \a
 * phiComponent.
 */
void transform_cartesian_to_spherical_vector_field_value(double theta, double phi,
                                                         double xComponent, double yComponent, double zComponent,
                                                         double &radialComponent, double &thetaComponent, double &phiComponent);

/**
 * Transform the spherical vector field components \a radialComponent, \a thetaComponent, \a phiComponent to the
 * Cartesian components \a xComponent, \a yComponent, \a zComponent. All components are stored as SphericalGridFunction
 * objects.
 */
void transform_spherical_to_cartesian_vector_field(const SphericalGridFunction &radialComponent, const SphericalGridFunction &thetaComponent, const SphericalGridFunction &phiComponent,
                                                   SphericalGridFunction &xComponent, SphericalGridFunction &yComponent, SphericalGridFunction &zComponent);

/**
 * Transform the Cartesian vector field components \a xComponent, \a yComponent, \a zComponent to the Cartesian
 * components \a radialComponent, \a thetaComponent, \a phiComponent. All components are stored as SphericalGridFunction
 * objects.
 */
void transform_cartesian_to_spherical_vector_field(const SphericalGridFunction &xComponent, const SphericalGridFunction &yComponent, const SphericalGridFunction &zComponent,
                                                   SphericalGridFunction &radialComponent, SphericalGridFunction &thetaComponent, SphericalGridFunction &phiComponent);

/**
 * Transform the spherical vector field components \a radialComponent, \a thetaComponent, \a phiComponent to the
 * Cartesian components \a xComponent, \a yComponent, \a zComponent. All components are stored as
 * Cartesian3DGridFunction objects.
 */
void transform_spherical_to_cartesian_vector_field(const Cartesian3DGridFunction &radialComponent, const Cartesian3DGridFunction &thetaComponent, const Cartesian3DGridFunction &phiComponent,
                                                   Cartesian3DGridFunction &xComponent, Cartesian3DGridFunction &yComponent, Cartesian3DGridFunction &zComponent);

/**
 * Transform the Cartesian vector field components \a xComponent, \a yComponent, \a zComponent to the Cartesian
 * components \a radialComponent, \a thetaComponent, \a phiComponent. All components are stored as
 * Cartesian3DGridFunction objects.
 */
void transform_cartesian_to_spherical_vector_field(const Cartesian3DGridFunction &xComponent, const Cartesian3DGridFunction &yComponent, const Cartesian3DGridFunction &zComponent,
                                                   Cartesian3DGridFunction &radialComponent, Cartesian3DGridFunction &thetaComponent, Cartesian3DGridFunction &phiComponent);

/**
 * \brief Structure defining a rotation of the coordinate system by specifying the angular coordinates of the target z-
 * and x-axes in the current coordinate system.
 */
struct CoordinateSystemRotation
{
  /**
   * Theta coordinate of the target z-axis.
   */
  double ThetaZ;

  /**
   * Phi coordinate of the target z-axis.
   */
  double PhiZ;

  /**
   * Theta coordinate of the target x-axis.
   */
  double ThetaX;

  /**
   * Phi coordinate of the target x-axis.
   */
  double PhiX;
};

extern const CoordinateSystemRotation
    /**
     * Do not rotate the coordinate system.
     */
    NoCoordinateSystemRotation,

    /**
     * Rotate from the equatorial to the galactic coordinate system.
     */
    EquatorialToGalactic,

    /**
     * Rotate from the galactic to the equatorial coordinate system.
     */
    GalacticToEquatorial,

    /**
     * Rotate from the equatorial to the supergalactic coordinate system.
     */
    EquatorialToSupergalactic,

    /**
     * Rotate from the supergalactic to the equatorial coordinate system.
     */
    SupergalacticToEquatorial,

    /**
     * Rotate from the galactic to the supergalactic coordinate system.
     */
    GalacticToSupergalactic,

    /**
     * Rotate from the supergalactic to the galactic coordinate system.
     */
    SupergalacticToGalactic;

/**
 * Apply the rotation \a coordinateSystemRotation to the Cartesian coordinates \a xInput, \a yInput, \a zInput, and
 * write the result to \a xOutput, \a yOutput, \a zOutput. If \a invert (default: \c false) is set to \c true,
 * the inverse rotation is performed.
 */
void rotate_cartesian_coordinates(double xInput, double yInput, double zInput,
                                  double &xOutput, double &yOutput, double &zOutput,
                                  CoordinateSystemRotation coordinateSystemRotation, bool invert = false);

/**
 * Apply the rotation \a coordinateSystemRotation to the spherical angular coordinates \a thetaInput, \a phiInput, and
 * write the result to \a thetaOutput, \a phiOutput. If \a invert (default: \c false) is set to \c true, the
 * inverse rotation is performed.
 */
void rotate_spherical_coordinates(double thetaInput, double phiInput,
                                  double &thetaOutput, double &phiOutput,
                                  CoordinateSystemRotation coordinateSystemRotation, bool invert = false);

/**
 * Apply the rotation \a coordinateSystemRotation to the celestial angular coordinates \a latitudeInput, \a
 * longitudeInput, and write the result to \a latitudeOutput, \a longitudeOutput. If \a invert (default: \c
 * false) is set to \c true, the inverse rotation is performed.
 */
void rotate_celestial_coordinates(double latitudeInput, double longitudeInput,
                                  double &latitudeOutput, double &longitudeOutput,
                                  CoordinateSystemRotation coordinateSystemRotation, bool invert = false);

/**
 * Apply the rotation \a coordinateSystemRotation to the scalar field \a fieldInput, and write the result to \a
 * fieldOutput. If \a invert (default: \c false) is set to \c true, the inverse rotation is performed.
 */
void rotate_scalar_field(const SphericalGridFunction &fieldInput, SphericalGridFunction &fieldOutput,
                         CoordinateSystemRotation coordinateSystemRotation, bool invert = false);

/**
 * Apply the rotation \a coordinateSystemRotation to the vector field with Cartesian components \a xFieldInput, \a
 * yFieldInput, \a zFieldInput, and write the result to \a xFieldOutput, \a yFieldOutput, \a zFieldOutput. If \a
 * invert (default: \c false) is set to \c true, the inverse rotation is performed.
 */
void rotate_cartesian_vector_field(const SphericalGridFunction &xFieldInput, const SphericalGridFunction &yFieldInput, const SphericalGridFunction &zFieldInput,
                                   SphericalGridFunction &xFieldOutput, SphericalGridFunction &yFieldOutput, SphericalGridFunction &zFieldOutput,
                                   CoordinateSystemRotation coordinateSystemRotation, bool invert = false);

/**
 * Apply the rotation \a coordinateSystemRotation to the vector field with spherical components \a radialFieldInput, \a
 * thetaFieldInput, \a phiFieldInput, and write the result to \a radialFieldOutput, \a thetaFieldOutput, \a
 * phiFieldOutput. If \a invert (default: \c false) is set to \c true, the inverse rotation is performed.
 */
void rotate_spherical_vector_field(const SphericalGridFunction &radialFieldInput, const SphericalGridFunction &thetaFieldInput, const SphericalGridFunction &phiFieldInput,
                                   SphericalGridFunction &radialFieldOutput, SphericalGridFunction &thetaFieldOutput, SphericalGridFunction &phiFieldOutput,
                                   CoordinateSystemRotation coordinateSystemRotation, bool invert = false);

/**
 * Return the CoordinateSystemRotation that satisfies two conditions: 1) Rotate the initial reference vector towards \a
 * thetaReferenceInitial, \a phiReferenceInitial onto the final reference vector towards \a thetaReferenceFinal, \a
 * phiReferenceFinal. 2) Rotate the plane spanned by the initial reference vector and the vector towards \a
 * thetaPlaneInitial, \a phiPlaneInitial onto the plane spanned by the final reference vector and the vector towards \a
 * thetaPlaneFinal, \a phiPlaneFinal.
 */
CoordinateSystemRotation reference_vector_and_plane_coordinate_system_rotation(double thetaReferenceInitial, double phiReferenceInitial, double thetaPlaneInitial, double phiPlaneInitial,
                                                                               double thetaReferenceFinal, double phiReferenceFinal, double thetaPlaneFinal, double phiPlaneFinal);

/**
 * \brief Enumeration defining various astronomical velocity and redshift reference frames.
 */
enum ReferenceFrame
{
  /**
   * Heliocentric reference frame.
   */
  HELIOCENTRIC_FRAME,

  /**
   * Local Group reference frame.
   */
  LOCAL_GROUP_FRAME,

  /**
   * CMB reference frame.
   */
  CMB_FRAME
};

/**
 * \brief Structure defining a change of reference frame by specifying the relative velocity between the final and
 * initial frames.
 */
struct ReferenceFrameChange
{
  /**
   * Amplitude of the relative velocity between the final and initial reference frames.
   */
  double relativeVelocityAmplitude;

  /**
   * Theta direction of the relative velocity between the final and initial reference frames.
   */
  double relativeVelocityTheta;

  /**
   * Phi direction of the relative velocity between the final and initial reference frames.
   */
  double relativeVelocityPhi;
};

extern const ReferenceFrameChange
    /**
     * Do not change the reference frame.
     */
    NO_REFERENCE_FRAME_CHANGE,

    /**
     * Transform from the heliocentric to the Local Group reference frame.
     */
    HELIOCENTRIC_TO_LOCAL_GROUP,

    /**
     * Transform from the Local Group to the heliocentric reference frame.
     */
    LOCAL_GROUP_TO_HELIOCENTRIC,

    /**
     * Transform from the heliocentric to the CMB reference frame.
     */
    HELIOCENTRIC_TO_CMB,

    /**
     * Transform from the CMB to the heliocentric reference frame.
     */
    CMB_TO_HELIOCENTRIC,

    /**
     * Transform from the Local Group to the CMB reference frame.
     */
    LOCAL_GROUP_TO_CMB,

    /**
     * Transform from the CMB to the Local Group reference frame.
     */
    CMB_TO_LOCAL_GROUP;

/**
 * Return the abbreviated name of the ReferenceFrame \a referenceFrameChange as a string.
 */
std::string get_reference_frame_name(ReferenceFrame referenceFrame);

ReferenceFrameChange get_reference_frame_change(ReferenceFrame inputFrame, ReferenceFrame outputFrame);

/**
 * Apply the change of reference frame \a referenceFrameChange to the radial velocity \a velocity pointing towards the
 * spherical angular coordinates \a theta, \a phi, and return the resulting radial velocity. If \a invert (default: \c
 * false) is set to \c true, the inverse reference frame change is performed.
 */
double change_reference_frame(double velocity, double theta, double phi, ReferenceFrameChange referenceFrameChange, bool invert = false);

///@}

/** @} */

#endif