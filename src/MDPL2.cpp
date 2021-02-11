#include "MDPL2.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <omp.h>

#include "Cartesian3DGridFunction.hpp"
#include "cosmology.hpp"
#include "transformations.hpp"

double luminosity_evolution_correction_MDPL2_SAGE(const double redshift)
{
    return 0.0;
}

void ignore_fortran_unformatted_boundary_field(std::ifstream &file)
{
    char ignore[4];

    file.read(ignore, 4);
}

template <typename T>
void read_fortran_unformatted_value(std::ifstream &file, T &variable)
{
    file.read(reinterpret_cast<char *>(&variable), sizeof(T));
}

Cartesian3DGridFunction read_MDPL2_field(const std::string &fileName)
{
    Cartesian3DGridFunction gridFunction(0.0, BoxLength, BoxBinNumber,
                                         0.0, BoxLength, BoxBinNumber,
                                         0.0, BoxLength, BoxBinNumber);

    std::ifstream inputFile(fileName, std::ios::binary);

    if (inputFile.fail())
    {
        std::cout << std::endl
                  << " read_MDPL2_field Error: File can not be opened" << std::endl
                  << std::endl;

        exit(EXIT_FAILURE);
    }

    int xNum, yNum, zNum;

    ignore_fortran_unformatted_boundary_field(inputFile);

    read_fortran_unformatted_value(inputFile, xNum);
    read_fortran_unformatted_value(inputFile, yNum);
    read_fortran_unformatted_value(inputFile, zNum);

    ignore_fortran_unformatted_boundary_field(inputFile);

    float gridValue;

    for (int i_z = 0; i_z < zNum; ++i_z)
    {
        ignore_fortran_unformatted_boundary_field(inputFile);

        for (int i_y = 0; i_y < yNum; ++i_y)
        {
            for (int i_x = 0; i_x < xNum; ++i_x)
            {
                read_fortran_unformatted_value(inputFile, gridValue);

                gridFunction.value(i_x, i_y, i_z) = static_cast<double>(gridValue);
            }
        }

        ignore_fortran_unformatted_boundary_field(inputFile);
    }

    inputFile.close();

    return gridFunction;
}

Cartesian3DGridFunction read_MDPL2_field(const std::string &fileName, double smoothingScale, bool useTophatKernel)
{
    // fftw_init_threads();
    // fftw_plan_with_nthreads(omp_get_max_threads());

    const int BinNumber1D = static_cast<int>(BoxBinNumber);
    const int BoxPointNumber = BinNumber1D * BinNumber1D * BinNumber1D;
    const double RealBinWidth3D = BoxLength / static_cast<double>(BinNumber1D);
    const double FourierBinWidth3D = 2.0 * M_PI / RealBinWidth3D / static_cast<double>(BinNumber1D); // 3D bin width in Fourier space chosen to match the definition of the discrete Fourier transform (DFT) in FFTW3 (see FFTW3 documentation for details)
    const double BackwardTrafoNormalization3D = gsl_pow_3(FourierBinWidth3D / (2.0 * M_PI));         // standard normalizations of 3D Fourier transform
    const double ForwardTrafoNormalization3D = gsl_pow_3(RealBinWidth3D);
    fftw_complex *FourierField = fftw_alloc_complex(BoxPointNumber);
    fftw_complex *RealField = fftw_alloc_complex(BoxPointNumber);
    fftw_plan BackwardTrafoPlan3D = fftw_plan_dft_3d(BinNumber1D, BinNumber1D, BinNumber1D, // FFT plan for transforming the field from Fourier to real space (see FFTW3 documentation for details)
                                                     FourierField, RealField,
                                                     FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan ForwardTrafoPlan3D = fftw_plan_dft_3d(BinNumber1D, BinNumber1D, BinNumber1D, // FFT plan for transforming the field from Fourier to real space (see FFTW3 documentation for details)
                                                    RealField, FourierField,
                                                    FFTW_FORWARD, FFTW_ESTIMATE);

    const auto fourier_coordinate_3D = [&](const int bin) {
        if (bin <= BinNumber1D / 2) // discretization chosen to match the definition of the DFT in FFTW3; the first half of the bins corresponds to positive, the second half to negative coordinates (see FFTW3 documentation for details)
        {
            return static_cast<double>(bin) * FourierBinWidth3D;
        }
        else
        {
            return static_cast<double>(bin - BinNumber1D) * FourierBinWidth3D;
        }
    };

    std::function<double(double)> window_function;

    if (useTophatKernel)
    {
        window_function = [&](double k) {
            return (k > 0.0) ? 3.0 / (k * smoothingScale) * gsl_sf_bessel_j1(k * smoothingScale)
                             : 1.0;
        };
    }
    else
    {
        window_function = [&](double k) {
            return std::exp(-gsl_pow_2(k * smoothingScale) / 2.0);
        };
    }

    Cartesian3DGridFunction gridFunction(0.0, BoxLength, BoxBinNumber,
                                         0.0, BoxLength, BoxBinNumber,
                                         0.0, BoxLength, BoxBinNumber);

    std::ifstream inputFile(fileName, std::ios::binary);

    if (inputFile.fail())
    {
        std::cout << std::endl
                  << " read_MDPL2_field Error: File can not be opened" << std::endl
                  << std::endl;

        exit(EXIT_FAILURE);
    }

    int xNum, yNum, zNum;

    ignore_fortran_unformatted_boundary_field(inputFile);

    read_fortran_unformatted_value(inputFile, xNum);
    read_fortran_unformatted_value(inputFile, yNum);
    read_fortran_unformatted_value(inputFile, zNum);

    ignore_fortran_unformatted_boundary_field(inputFile);

    float gridValue;

    for (int i_z = 0; i_z < zNum; ++i_z)
    {
        ignore_fortran_unformatted_boundary_field(inputFile);

        for (int i_y = 0; i_y < yNum; ++i_y)
        {
            for (int i_x = 0; i_x < xNum; ++i_x)
            {
                read_fortran_unformatted_value(inputFile, gridValue);

                const std::size_t i_xyz = (i_x * BinNumber1D + i_y) * BinNumber1D + i_z;

                RealField[i_xyz][0] = static_cast<double>(gridValue) * ForwardTrafoNormalization3D;
                RealField[i_xyz][1] = 0.0;
            }
        }

        ignore_fortran_unformatted_boundary_field(inputFile);
    }

    inputFile.close();

    fftw_execute(ForwardTrafoPlan3D);

    for (int i_x = 0; i_x < BinNumber1D; ++i_x) // apply the smoothing
    {
        const double kx = fourier_coordinate_3D(i_x);

        for (int i_y = 0; i_y < BinNumber1D; ++i_y)
        {
            const double ky = fourier_coordinate_3D(i_y);

            for (int i_z = 0; i_z < BinNumber1D; ++i_z)
            {
                const double kz = fourier_coordinate_3D(i_z);

                const double k = std::sqrt(kx * kx + ky * ky + kz * kz);

                const std::size_t i_xyz = (i_x * BinNumber1D + i_y) * BinNumber1D + i_z;

                const double windowFunction = window_function(k);

                FourierField[i_xyz][0] *= windowFunction * BackwardTrafoNormalization3D;
                FourierField[i_xyz][1] *= windowFunction * BackwardTrafoNormalization3D;
            }
        }
    }

    fftw_execute(BackwardTrafoPlan3D);

    for (int i_x = 0; i_x < BinNumber1D; ++i_x)
    {
        for (int i_y = 0; i_y < BinNumber1D; ++i_y)
        {
            for (int i_z = 0; i_z < BinNumber1D; ++i_z)
            {
                const std::size_t i_xyz = (i_x * BinNumber1D + i_y) * BinNumber1D + i_z;

                gridFunction.value(i_x, i_y, i_z) = RealField[i_xyz][0];
            }
        }
    }

    fftw_free(FourierField);
    fftw_free(RealField);

    fftw_destroy_plan(BackwardTrafoPlan3D);
    fftw_destroy_plan(ForwardTrafoPlan3D);

    // fftw_cleanup_threads();

    return gridFunction;
}

Cartesian3DGridFunction read_MDPL2_velocity_divergence_field(const std::string &xVelocityFileName, const std::string &yVelocityFileName, const std::string &zVelocityFileName,
                                                             double smoothingScale, bool useTophatKernel)
{
    // fftw_init_threads();
    // fftw_plan_with_nthreads(omp_get_max_threads());

    const int BinNumber3D = static_cast<int>(BoxBinNumber);
    const int BoxPointNumber = BinNumber3D * BinNumber3D * BinNumber3D;
    const int BinNumbers3D[3] = {BinNumber3D, BinNumber3D, BinNumber3D};
    const double RealBinWidth3D = BoxLength / static_cast<double>(BinNumber3D);
    const double FourierBinWidth3D = 2.0 * M_PI / RealBinWidth3D / static_cast<double>(BinNumber3D); // 3D bin width in Fourier space chosen to match the definition of the discrete Fourier transform (DFT) in FFTW3 (see FFTW3 documentation for details)
    const double BackwardTrafoNormalization3D = gsl_pow_3(FourierBinWidth3D / (2.0 * M_PI));         // standard normalizations of 3D Fourier transform
    const double ForwardTrafoNormalization3D = gsl_pow_3(RealBinWidth3D);
    fftw_complex *FourierVelocity = fftw_alloc_complex(3 * BoxPointNumber);
    fftw_complex *RealVelocity = fftw_alloc_complex(3 * BoxPointNumber);
    fftw_complex *FourierVelocityDivergence = fftw_alloc_complex(BoxPointNumber);
    fftw_complex *RealVelocityDivergence = fftw_alloc_complex(BoxPointNumber);
    fftw_plan VelocityForwardTrafoPlan3D = fftw_plan_many_dft(3, BinNumbers3D, 3, // FFT plan for simultaneously transforming all three components of the velocity field from real to Fourier space (see FFTW3 documentation for details)
                                                              RealVelocity, NULL,
                                                              1, BoxPointNumber,
                                                              FourierVelocity, NULL,
                                                              1, BoxPointNumber,
                                                              FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan VelocityDivergenceBackwardTrafoPlan3D = fftw_plan_dft_3d(BinNumber3D, BinNumber3D, BinNumber3D, // FFT plan for transforming the velocity divergence from Fourier to real space (see FFTW3 documentation for details)
                                                                       FourierVelocityDivergence, RealVelocityDivergence,
                                                                       FFTW_BACKWARD, FFTW_ESTIMATE);

    const auto fourier_coordinate_3D = [&](const int bin) {
        if (bin <= BinNumber3D / 2) // discretization chosen to match the definition of the DFT in FFTW3; the first half of the bins corresponds to positive, the second half to negative coordinates (see FFTW3 documentation for details)
        {
            return static_cast<double>(bin) * FourierBinWidth3D;
        }
        else
        {
            return static_cast<double>(bin - BinNumber3D) * FourierBinWidth3D;
        }
    };

    std::function<double(double)> window_function;

    if (useTophatKernel)
    {
        window_function = [&](double k) {
            return (k > 0.0) ? 3.0 / (k * smoothingScale) * gsl_sf_bessel_j1(k * smoothingScale)
                             : 1.0;
        };
    }
    else
    {
        window_function = [&](double k) {
            return std::exp(-gsl_pow_2(k * smoothingScale) / 2.0);
        };
    }

    Cartesian3DGridFunction gridFunction(0.0, BoxLength, BoxBinNumber,
                                         0.0, BoxLength, BoxBinNumber,
                                         0.0, BoxLength, BoxBinNumber);

    std::ifstream xVelocityFile(xVelocityFileName, std::ios::binary);
    std::ifstream yVelocityFile(yVelocityFileName, std::ios::binary);
    std::ifstream zVelocityFile(zVelocityFileName, std::ios::binary);

    if (xVelocityFile.fail() or
        yVelocityFile.fail() or
        zVelocityFile.fail())
    {
        std::cout << std::endl
                  << " read_MDPL2_velocity_divergence_field Error: File can not be opened" << std::endl
                  << std::endl;

        exit(EXIT_FAILURE);
    }

    int xVelocityXNum, xVelocityYNum, xVelocityZNum,
        yVelocityXNum, yVelocityYNum, yVelocityZNum,
        zVelocityXNum, zVelocityYNum, zVelocityZNum;

    ignore_fortran_unformatted_boundary_field(xVelocityFile);
    ignore_fortran_unformatted_boundary_field(yVelocityFile);
    ignore_fortran_unformatted_boundary_field(zVelocityFile);

    read_fortran_unformatted_value(xVelocityFile, xVelocityXNum);
    read_fortran_unformatted_value(xVelocityFile, xVelocityYNum);
    read_fortran_unformatted_value(xVelocityFile, xVelocityZNum);

    read_fortran_unformatted_value(yVelocityFile, yVelocityXNum);
    read_fortran_unformatted_value(yVelocityFile, yVelocityYNum);
    read_fortran_unformatted_value(yVelocityFile, yVelocityZNum);

    read_fortran_unformatted_value(zVelocityFile, zVelocityXNum);
    read_fortran_unformatted_value(zVelocityFile, zVelocityYNum);
    read_fortran_unformatted_value(zVelocityFile, zVelocityZNum);

    ignore_fortran_unformatted_boundary_field(xVelocityFile);
    ignore_fortran_unformatted_boundary_field(yVelocityFile);
    ignore_fortran_unformatted_boundary_field(zVelocityFile);

    if ((xVelocityXNum != yVelocityXNum) or (xVelocityXNum != zVelocityXNum) or
        (xVelocityYNum != yVelocityYNum) or (xVelocityYNum != zVelocityYNum) or
        (xVelocityZNum != yVelocityZNum) or (xVelocityZNum != zVelocityZNum))
    {
        std::cout << std::endl
                  << " read_MDPL2_velocity_divergence_field Error: Bin numbers differ between files" << std::endl
                  << std::endl;

        exit(EXIT_FAILURE);
    }

    float xVelocityGridValue, yVelocityGridValue, zVelocityGridValue;

    for (int i_z = 0; i_z < xVelocityZNum; ++i_z)
    {
        ignore_fortran_unformatted_boundary_field(xVelocityFile);
        ignore_fortran_unformatted_boundary_field(yVelocityFile);
        ignore_fortran_unformatted_boundary_field(zVelocityFile);

        for (int i_y = 0; i_y < xVelocityYNum; ++i_y)
        {
            for (int i_x = 0; i_x < xVelocityXNum; ++i_x)
            {
                read_fortran_unformatted_value(xVelocityFile, xVelocityGridValue);
                read_fortran_unformatted_value(yVelocityFile, yVelocityGridValue);
                read_fortran_unformatted_value(zVelocityFile, zVelocityGridValue);

                const std::size_t i_xyz = (i_x * BinNumber3D + i_y) * BinNumber3D + i_z;

                RealVelocity[i_xyz][0] = static_cast<double>(xVelocityGridValue) * ForwardTrafoNormalization3D;
                RealVelocity[i_xyz][1] = 0.0;
                RealVelocity[i_xyz + BoxPointNumber][0] = static_cast<double>(yVelocityGridValue) * ForwardTrafoNormalization3D;
                RealVelocity[i_xyz + BoxPointNumber][1] = 0.0;
                RealVelocity[i_xyz + 2 * BoxPointNumber][0] = static_cast<double>(zVelocityGridValue) * ForwardTrafoNormalization3D;
                RealVelocity[i_xyz + 2 * BoxPointNumber][1] = 0.0;
            }
        }

        ignore_fortran_unformatted_boundary_field(xVelocityFile);
        ignore_fortran_unformatted_boundary_field(yVelocityFile);
        ignore_fortran_unformatted_boundary_field(zVelocityFile);
    }

    xVelocityFile.close();
    yVelocityFile.close();
    zVelocityFile.close();

    fftw_execute(VelocityForwardTrafoPlan3D);

    for (int i_x = 0; i_x < BinNumber3D; ++i_x) // apply the smoothing and compute the divergence
    {
        const double kx = fourier_coordinate_3D(i_x);

        for (int i_y = 0; i_y < BinNumber3D; ++i_y)
        {
            const double ky = fourier_coordinate_3D(i_y);

            for (int i_z = 0; i_z < BinNumber3D; ++i_z)
            {
                const double kz = fourier_coordinate_3D(i_z);

                const double k = std::sqrt(kx * kx + ky * ky + kz * kz);

                const std::size_t i_xyz = (i_x * BinNumber3D + i_y) * BinNumber3D + i_z;

                const double windowFunction = window_function(k);

                const double fourierXVelocityRe = FourierVelocity[i_xyz][0] * windowFunction;
                const double fourierXVelocityIm = FourierVelocity[i_xyz][1] * windowFunction;
                const double fourierYVelocityRe = FourierVelocity[i_xyz + BoxPointNumber][0] * windowFunction;
                const double fourierYVelocityIm = FourierVelocity[i_xyz + BoxPointNumber][1] * windowFunction;
                const double fourierZVelocityRe = FourierVelocity[i_xyz + 2 * BoxPointNumber][0] * windowFunction;
                const double fourierZVelocityIm = FourierVelocity[i_xyz + 2 * BoxPointNumber][1] * windowFunction;

                FourierVelocityDivergence[i_xyz][0] = (kx * fourierXVelocityIm + ky * fourierYVelocityIm + kz * fourierZVelocityIm) * BackwardTrafoNormalization3D; // div(v)(k) = - i k * v(k)
                FourierVelocityDivergence[i_xyz][1] = -(kx * fourierXVelocityRe + ky * fourierYVelocityRe + kz * fourierZVelocityRe) * BackwardTrafoNormalization3D;
            }
        }
    }

    fftw_execute(VelocityDivergenceBackwardTrafoPlan3D);

    for (int i_x = 0; i_x < BinNumber3D; ++i_x)
    {
        for (int i_y = 0; i_y < BinNumber3D; ++i_y)
        {
            for (int i_z = 0; i_z < BinNumber3D; ++i_z)
            {
                const std::size_t i_xyz = (i_x * BinNumber3D + i_y) * BinNumber3D + i_z;

                gridFunction.value(i_x, i_y, i_z) = RealVelocityDivergence[i_xyz][0];
            }
        }
    }

    fftw_free(FourierVelocity);
    fftw_free(RealVelocity);
    fftw_free(FourierVelocityDivergence);
    fftw_free(RealVelocityDivergence);

    fftw_destroy_plan(VelocityForwardTrafoPlan3D);
    fftw_destroy_plan(VelocityDivergenceBackwardTrafoPlan3D);

    // fftw_cleanup_threads();

    return gridFunction;
}

void read_MDPL2_density_field_and_get_linear_velocity(const std::string &fileName, double smoothingScale, double beta,
                                                      Cartesian3DGridFunction &densityContrast, Cartesian3DGridFunction &xVelocity, Cartesian3DGridFunction &yVelocity, Cartesian3DGridFunction &zVelocity,
                                                      bool useTophatKernel)
{
    // fftw_init_threads();
    // fftw_plan_with_nthreads(omp_get_max_threads());

    const int BinNumber3D = static_cast<int>(BoxBinNumber);
    const int BoxPointNumber = BinNumber3D * BinNumber3D * BinNumber3D;
    const int BinNumbers3D[3] = {BinNumber3D, BinNumber3D, BinNumber3D};
    const double RealBinWidth3D = BoxLength / static_cast<double>(BinNumber3D);
    const double FourierBinWidth3D = 2.0 * M_PI / RealBinWidth3D / static_cast<double>(BinNumber3D); // 3D bin width in Fourier space chosen to match the definition of the discrete Fourier transform (DFT) in FFTW3 (see FFTW3 documentation for details)
    const double BackwardTrafoNormalization3D = gsl_pow_3(FourierBinWidth3D / (2.0 * M_PI));         // standard normalizations of 3D Fourier transform
    const double ForwardTrafoNormalization3D = gsl_pow_3(RealBinWidth3D);
    fftw_complex *FourierDensity = fftw_alloc_complex(BoxPointNumber);
    fftw_complex *RealDensity = fftw_alloc_complex(BoxPointNumber);
    fftw_complex *FourierVelocity(fftw_alloc_complex(3 * BoxPointNumber));
    fftw_complex *RealVelocity(fftw_alloc_complex(3 * BoxPointNumber));
    fftw_plan DensityBackwardTrafoPlan3D = fftw_plan_dft_3d(BinNumber3D, BinNumber3D, BinNumber3D, // FFT plan for transforming the field from Fourier to real space (see FFTW3 documentation for details)
                                                            FourierDensity, RealDensity,
                                                            FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan DensityForwardTrafoPlan3D = fftw_plan_dft_3d(BinNumber3D, BinNumber3D, BinNumber3D, // FFT plan for transforming the field from Fourier to real space (see FFTW3 documentation for details)
                                                           RealDensity, FourierDensity,
                                                           FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan VelocityBackwardTrafoPlan3D(fftw_plan_many_dft(3, BinNumbers3D, 3, // FFT plan for simultaneously transforming all three components of the velocity field from Fourier to real space (see FFTW3 documentation for details)
                                                             FourierVelocity, NULL,
                                                             1, BoxPointNumber,
                                                             RealVelocity, NULL,
                                                             1, BoxPointNumber,
                                                             FFTW_BACKWARD, FFTW_ESTIMATE));

    const auto fourier_coordinate_3D = [&](const int bin) {
        if (bin <= BinNumber3D / 2) // discretization chosen to match the definition of the DFT in FFTW3; the first half of the bins corresponds to positive, the second half to negative coordinates (see FFTW3 documentation for details)
        {
            return static_cast<double>(bin) * FourierBinWidth3D;
        }
        else
        {
            return static_cast<double>(bin - BinNumber3D) * FourierBinWidth3D;
        }
    };

    std::function<double(double)> window_function;

    if (useTophatKernel)
    {
        window_function = [&](double k) {
            return (k > 0.0) ? 3.0 / (k * smoothingScale) * gsl_sf_bessel_j1(k * smoothingScale)
                             : 1.0;
        };
    }
    else
    {
        window_function = [&](double k) {
            return std::exp(-gsl_pow_2(k * smoothingScale) / 2.0);
        };
    }

    densityContrast = Cartesian3DGridFunction(0.0, BoxLength, BoxBinNumber,
                                              0.0, BoxLength, BoxBinNumber,
                                              0.0, BoxLength, BoxBinNumber);

    xVelocity = densityContrast;
    yVelocity = densityContrast;
    zVelocity = densityContrast;

    std::ifstream inputFile(fileName, std::ios::binary);

    if (inputFile.fail())
    {
        std::cout << std::endl
                  << " read_MDPL2_density_field_and_get_linear_velocity Error: File can not be opened" << std::endl
                  << std::endl;

        exit(EXIT_FAILURE);
    }

    int xNum, yNum, zNum;

    ignore_fortran_unformatted_boundary_field(inputFile);

    read_fortran_unformatted_value(inputFile, xNum);
    read_fortran_unformatted_value(inputFile, yNum);
    read_fortran_unformatted_value(inputFile, zNum);

    ignore_fortran_unformatted_boundary_field(inputFile);

    float gridValue;

    for (int i_z = 0; i_z < zNum; ++i_z)
    {
        ignore_fortran_unformatted_boundary_field(inputFile);

        for (int i_y = 0; i_y < yNum; ++i_y)
        {
            for (int i_x = 0; i_x < xNum; ++i_x)
            {
                read_fortran_unformatted_value(inputFile, gridValue);

                const std::size_t i_xyz = (i_x * BinNumber3D + i_y) * BinNumber3D + i_z;

                RealDensity[i_xyz][0] = (static_cast<double>(gridValue) - 1.0) * ForwardTrafoNormalization3D;
                RealDensity[i_xyz][1] = 0.0;
            }
        }

        ignore_fortran_unformatted_boundary_field(inputFile);
    }

    inputFile.close();

    fftw_execute(DensityForwardTrafoPlan3D);

    for (int i_x = 0; i_x < BinNumber3D; ++i_x) // apply the smoothing and compute the linear velocity
    {
        const double kx = fourier_coordinate_3D(i_x);

        for (int i_y = 0; i_y < BinNumber3D; ++i_y)
        {
            const double ky = fourier_coordinate_3D(i_y);

            for (int i_z = 0; i_z < BinNumber3D; ++i_z)
            {
                if ((i_x == 0) and (i_y == 0) and (i_z == 0)) // corresponds to the zero-mode, which vanishes for zero-mean density contrast and velocity fields
                {
                    FourierDensity[0][0] = 0.0;
                    FourierDensity[0][1] = 0.0;

                    FourierVelocity[0][0] = 0.0;
                    FourierVelocity[0][1] = 0.0;
                    FourierVelocity[BoxPointNumber][0] = 0.0;
                    FourierVelocity[BoxPointNumber][1] = 0.0;
                    FourierVelocity[2 * BoxPointNumber][0] = 0.0;
                    FourierVelocity[2 * BoxPointNumber][1] = 0.0;
                }
                else
                {
                    const double kz = fourier_coordinate_3D(i_z);

                    const double kSquared = kx * kx + ky * ky + kz * kz;
                    const double k = std::sqrt(kSquared);

                    const std::size_t i_xyz = (i_x * BinNumber3D + i_y) * BinNumber3D + i_z;

                    const double windowFunction = window_function(k);

                    const double fourierDensityRe = FourierDensity[i_xyz][0] * windowFunction;
                    const double fourierDensityIm = FourierDensity[i_xyz][1] * windowFunction;

                    FourierDensity[i_xyz][0] = fourierDensityRe * BackwardTrafoNormalization3D;
                    FourierDensity[i_xyz][1] = fourierDensityIm * BackwardTrafoNormalization3D;

                    FourierVelocity[i_xyz][0] = -beta * HUBBLE_NORMALIZATION * kx / kSquared * fourierDensityIm * BackwardTrafoNormalization3D; // v(k) = beta H_0 delta_g(k) i k / |k|^2
                    FourierVelocity[i_xyz][1] = beta * HUBBLE_NORMALIZATION * kx / kSquared * fourierDensityRe * BackwardTrafoNormalization3D;
                    FourierVelocity[i_xyz + BoxPointNumber][0] = -beta * HUBBLE_NORMALIZATION * ky / kSquared * fourierDensityIm * BackwardTrafoNormalization3D;
                    FourierVelocity[i_xyz + BoxPointNumber][1] = beta * HUBBLE_NORMALIZATION * ky / kSquared * fourierDensityRe * BackwardTrafoNormalization3D;
                    FourierVelocity[i_xyz + 2 * BoxPointNumber][0] = -beta * HUBBLE_NORMALIZATION * kz / kSquared * fourierDensityIm * BackwardTrafoNormalization3D;
                    FourierVelocity[i_xyz + 2 * BoxPointNumber][1] = beta * HUBBLE_NORMALIZATION * kz / kSquared * fourierDensityRe * BackwardTrafoNormalization3D;
                }
            }
        }
    }

    fftw_execute(DensityBackwardTrafoPlan3D);
    fftw_execute(VelocityBackwardTrafoPlan3D);

    for (int i_x = 0; i_x < BinNumber3D; ++i_x)
    {
        for (int i_y = 0; i_y < BinNumber3D; ++i_y)
        {
            for (int i_z = 0; i_z < BinNumber3D; ++i_z)
            {
                const std::size_t i_xyz = (i_x * BinNumber3D + i_y) * BinNumber3D + i_z;

                densityContrast.value(i_x, i_y, i_z) = RealDensity[i_xyz][0];
                xVelocity.value(i_x, i_y, i_z) = RealVelocity[i_xyz][0];
                yVelocity.value(i_x, i_y, i_z) = RealVelocity[i_xyz + BoxPointNumber][0];
                zVelocity.value(i_x, i_y, i_z) = RealVelocity[i_xyz + 2 * BoxPointNumber][0];
            }
        }
    }

    fftw_free(FourierDensity);
    fftw_free(RealDensity);
    fftw_free(FourierVelocity);
    fftw_free(RealVelocity);

    fftw_destroy_plan(DensityBackwardTrafoPlan3D);
    fftw_destroy_plan(DensityForwardTrafoPlan3D);
    fftw_destroy_plan(VelocityBackwardTrafoPlan3D);

    // fftw_cleanup_threads();
}

std::size_t read_MDPL2_galaxy_data(const std::string &fileName,
                                   std::vector<std::size_t> &galaxy_rowId,
                                   std::vector<std::size_t> &galaxy_rockstarId,
                                   std::vector<std::size_t> &galaxy_hostHaloId,
                                   std::vector<std::size_t> &galaxy_mainHaloId,
                                   std::vector<std::size_t> &galaxy_galaxyType,
                                   std::vector<double> &galaxy_haloMass,
                                   std::vector<double> &galaxy_maxVelocity,
                                   std::vector<double> &galaxy_peakVelocity,
                                   std::vector<double> &galaxy_xCoordinate,
                                   std::vector<double> &galaxy_yCoordinate,
                                   std::vector<double> &galaxy_zCoordinate,
                                   std::vector<double> &galaxy_xVelocity,
                                   std::vector<double> &galaxy_yVelocity,
                                   std::vector<double> &galaxy_zVelocity,
                                   std::vector<double> &galaxy_stellarMass,
                                   std::vector<double> &galaxy_meanStarAge)
{
    galaxy_rowId.clear();
    galaxy_rockstarId.clear();
    galaxy_hostHaloId.clear();
    galaxy_mainHaloId.clear();
    galaxy_galaxyType.clear();
    galaxy_haloMass.clear();
    galaxy_maxVelocity.clear();
    galaxy_peakVelocity.clear();
    galaxy_xCoordinate.clear();
    galaxy_yCoordinate.clear();
    galaxy_zCoordinate.clear();
    galaxy_xVelocity.clear();
    galaxy_yVelocity.clear();
    galaxy_zVelocity.clear();
    galaxy_stellarMass.clear();
    galaxy_meanStarAge.clear();

    std::size_t galaxyNumber, rockstarId, rowId, hostHaloId, mainHaloId, galaxyType;
    double haloMass, maxVelocity, peakVelocity, xCoordinate, yCoordinate, zCoordinate, xVelocity, yVelocity, zVelocity, stellarMass, meanStarAge;

    const std::string binaryFileName = std::string(fileName).replace(fileName.find_last_of('.') + 1, std::string::npos, "bin");

    std::ifstream binaryGalaxyInputFile(binaryFileName, std::ios::binary);

    if (binaryGalaxyInputFile.fail())
    {
        std::ifstream formattedGalaxyInputFile(fileName, std::ios::in);

        if (formattedGalaxyInputFile.fail())
        {
            std::cout << std::endl
                      << " Error: Formatted galaxy input file can not be opened" << std::endl
                      << std::endl;

            exit(EXIT_FAILURE);
        }

        std::string line;

        std::getline(formattedGalaxyInputFile, line); // header

        char delimiter;

        while (std::getline(formattedGalaxyInputFile, line))
        {
            std::istringstream lineStream(line);

            lineStream >>
                rowId >> delimiter >>
                rockstarId >> delimiter >>
                hostHaloId >> delimiter >>
                mainHaloId >> delimiter >>
                galaxyType >> delimiter >>
                haloMass >> delimiter >>
                maxVelocity >> delimiter >>
                peakVelocity >> delimiter >>
                xCoordinate >> delimiter >>
                yCoordinate >> delimiter >>
                zCoordinate >> delimiter >>
                xVelocity >> delimiter >>
                yVelocity >> delimiter >>
                zVelocity >> delimiter >>
                stellarMass >> delimiter >>
                meanStarAge;

            // galaxy_rowId.push_back(rowId);
            galaxy_rockstarId.push_back(rockstarId);
            galaxy_hostHaloId.push_back(hostHaloId);
            galaxy_mainHaloId.push_back(mainHaloId);
            // galaxy_galaxyType.push_back(galaxyType);
            galaxy_haloMass.push_back(haloMass / HubbleSim);
            // galaxy_maxVelocity.push_back(maxVelocity);
            // galaxy_peakVelocity.push_back(peakVelocity);
            galaxy_xCoordinate.push_back(xCoordinate);
            galaxy_yCoordinate.push_back(yCoordinate);
            galaxy_zCoordinate.push_back(zCoordinate);
            galaxy_xVelocity.push_back(xVelocity);
            galaxy_yVelocity.push_back(yVelocity);
            galaxy_zVelocity.push_back(zVelocity);
            galaxy_stellarMass.push_back(stellarMass / HubbleSim);
            // galaxy_meanStarAge.push_back(meanStarAge);
        }

        formattedGalaxyInputFile.close();

        galaxyNumber = galaxy_xCoordinate.size();

        std::ofstream binaryGalaxyOutputFile(binaryFileName, std::ios::binary);

        binaryGalaxyOutputFile.write((char *)&galaxyNumber, sizeof(std::size_t));

        for (std::size_t i_g = 0; i_g < galaxyNumber; ++i_g)
        {
            // binaryGalaxyOutputFile.write((char *)&galaxy_rowId[i_g], sizeof(std::size_t));
            binaryGalaxyOutputFile.write((char *)&galaxy_rockstarId[i_g], sizeof(std::size_t));
            binaryGalaxyOutputFile.write((char *)&galaxy_hostHaloId[i_g], sizeof(std::size_t));
            binaryGalaxyOutputFile.write((char *)&galaxy_mainHaloId[i_g], sizeof(std::size_t));
            // binaryGalaxyOutputFile.write((char *)&galaxy_galaxyType[i_g], sizeof(std::size_t));
            binaryGalaxyOutputFile.write((char *)&galaxy_haloMass[i_g], sizeof(double));
            // binaryGalaxyOutputFile.write((char *)&galaxy_maxVelocity[i_g], sizeof(double));
            // binaryGalaxyOutputFile.write((char *)&galaxy_peakVelocity[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_xCoordinate[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_yCoordinate[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_zCoordinate[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_xVelocity[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_yVelocity[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_zVelocity[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_stellarMass[i_g], sizeof(double));
            // binaryGalaxyOutputFile.write((char *)&galaxy_meanStarAge[i_g], sizeof(double));
        }

        binaryGalaxyOutputFile.close();
    }
    else
    {
        binaryGalaxyInputFile.read((char *)&galaxyNumber, sizeof(std::size_t));

        // galaxy_rowId.resize(galaxyNumber);
        galaxy_rockstarId.resize(galaxyNumber);
        galaxy_hostHaloId.resize(galaxyNumber);
        galaxy_mainHaloId.resize(galaxyNumber);
        // galaxy_galaxyType.resize(galaxyNumber);
        galaxy_haloMass.resize(galaxyNumber);
        // galaxy_maxVelocity.resize(galaxyNumber);
        // galaxy_peakVelocity.resize(galaxyNumber);
        galaxy_xCoordinate.resize(galaxyNumber);
        galaxy_yCoordinate.resize(galaxyNumber);
        galaxy_zCoordinate.resize(galaxyNumber);
        galaxy_xVelocity.resize(galaxyNumber);
        galaxy_yVelocity.resize(galaxyNumber);
        galaxy_zVelocity.resize(galaxyNumber);
        galaxy_stellarMass.resize(galaxyNumber);
        // galaxy_meanStarAge.resize(galaxyNumber);

        for (std::size_t i_g = 0; i_g < galaxyNumber; ++i_g)
        {
            // binaryGalaxyInputFile.read((char *)&galaxy_rowId[i_g], sizeof(std::size_t));
            binaryGalaxyInputFile.read((char *)&galaxy_rockstarId[i_g], sizeof(std::size_t));
            binaryGalaxyInputFile.read((char *)&galaxy_hostHaloId[i_g], sizeof(std::size_t));
            binaryGalaxyInputFile.read((char *)&galaxy_mainHaloId[i_g], sizeof(std::size_t));
            // binaryGalaxyInputFile.read((char *)&galaxy_galaxyType[i_g], sizeof(std::size_t));
            binaryGalaxyInputFile.read((char *)&galaxy_haloMass[i_g], sizeof(double));
            // binaryGalaxyInputFile.read((char *)&galaxy_maxVelocity[i_g], sizeof(double));
            // binaryGalaxyInputFile.read((char *)&galaxy_peakVelocity[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_xCoordinate[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_yCoordinate[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_zCoordinate[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_xVelocity[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_yVelocity[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_zVelocity[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_stellarMass[i_g], sizeof(double));
            // binaryGalaxyInputFile.read((char *)&galaxy_meanStarAge[i_g], sizeof(double));
        }
    }

    binaryGalaxyInputFile.close();

    return galaxy_xCoordinate.size();
}

std::size_t read_MDPL2_galaxy_sage_data(const std::string &fileName,
                                        std::vector<std::size_t> &galaxy_rowId,
                                        std::vector<std::size_t> &galaxy_rockstarId,
                                        std::vector<std::size_t> &galaxy_hostHaloId,
                                        std::vector<std::size_t> &galaxy_mainHaloId,
                                        std::vector<std::size_t> &galaxy_galaxyType,
                                        std::vector<double> &galaxy_haloMass,
                                        std::vector<double> &galaxy_maxVelocity,
                                        std::vector<double> &galaxy_xCoordinate,
                                        std::vector<double> &galaxy_yCoordinate,
                                        std::vector<double> &galaxy_zCoordinate,
                                        std::vector<double> &galaxy_xVelocity,
                                        std::vector<double> &galaxy_yVelocity,
                                        std::vector<double> &galaxy_zVelocity,
                                        std::vector<double> &galaxy_stellarMass,
                                        std::vector<double> &galaxy_meanStarAge)
{
    galaxy_rowId.clear();
    galaxy_rockstarId.clear();
    galaxy_hostHaloId.clear();
    galaxy_mainHaloId.clear();
    galaxy_galaxyType.clear();
    galaxy_haloMass.clear();
    galaxy_maxVelocity.clear();
    galaxy_xCoordinate.clear();
    galaxy_yCoordinate.clear();
    galaxy_zCoordinate.clear();
    galaxy_xVelocity.clear();
    galaxy_yVelocity.clear();
    galaxy_zVelocity.clear();
    galaxy_stellarMass.clear();
    galaxy_meanStarAge.clear();

    std::size_t galaxyNumber, rockstarId, rowId, hostHaloId, mainHaloId, galaxyType;
    double haloMass, maxVelocity, xCoordinate, yCoordinate, zCoordinate, xVelocity, yVelocity, zVelocity, stellarMass, meanStarAge;

    const std::string binaryFileName = std::string(fileName).replace(fileName.find_last_of('.') + 1, std::string::npos, "bin");

    std::ifstream binaryGalaxyInputFile(binaryFileName, std::ios::binary);

    if (binaryGalaxyInputFile.fail())
    {
        std::ifstream formattedGalaxyInputFile(fileName, std::ios::in);

        if (formattedGalaxyInputFile.fail())
        {
            std::cout << std::endl
                      << " Error: Formatted galaxy input file can not be opened" << std::endl
                      << std::endl;

            exit(EXIT_FAILURE);
        }

        std::string line;

        std::getline(formattedGalaxyInputFile, line); // header

        char delimiter;

        while (std::getline(formattedGalaxyInputFile, line))
        {
            std::istringstream lineStream(line);

            lineStream >>
                rowId >> delimiter >>
                rockstarId >> delimiter >>
                hostHaloId >> delimiter >>
                mainHaloId >> delimiter >>
                galaxyType >> delimiter >>
                haloMass >> delimiter >>
                maxVelocity >> delimiter >>
                xCoordinate >> delimiter >>
                yCoordinate >> delimiter >>
                zCoordinate >> delimiter >>
                xVelocity >> delimiter >>
                yVelocity >> delimiter >>
                zVelocity >> delimiter >>
                stellarMass >> delimiter >>
                meanStarAge;

            // galaxy_rowId.push_back(rowId);
            galaxy_rockstarId.push_back(rockstarId);
            galaxy_hostHaloId.push_back(hostHaloId);
            galaxy_mainHaloId.push_back(mainHaloId);
            // galaxy_galaxyType.push_back(galaxyType);
            galaxy_haloMass.push_back(haloMass / HubbleSim);
            // galaxy_maxVelocity.push_back(maxVelocity);
            galaxy_xCoordinate.push_back(xCoordinate);
            galaxy_yCoordinate.push_back(yCoordinate);
            galaxy_zCoordinate.push_back(zCoordinate);
            galaxy_xVelocity.push_back(xVelocity);
            galaxy_yVelocity.push_back(yVelocity);
            galaxy_zVelocity.push_back(zVelocity);
            galaxy_stellarMass.push_back(stellarMass / HubbleSim);
            // galaxy_meanStarAge.push_back(meanStarAge);
        }

        formattedGalaxyInputFile.close();

        galaxyNumber = galaxy_xCoordinate.size();

        std::ofstream binaryGalaxyOutputFile(binaryFileName, std::ios::binary);

        binaryGalaxyOutputFile.write((char *)&galaxyNumber, sizeof(std::size_t));

        for (std::size_t i_g = 0; i_g < galaxyNumber; ++i_g)
        {
            // binaryGalaxyOutputFile.write((char *)&galaxy_rowId[i_g], sizeof(std::size_t));
            binaryGalaxyOutputFile.write((char *)&galaxy_rockstarId[i_g], sizeof(std::size_t));
            binaryGalaxyOutputFile.write((char *)&galaxy_hostHaloId[i_g], sizeof(std::size_t));
            binaryGalaxyOutputFile.write((char *)&galaxy_mainHaloId[i_g], sizeof(std::size_t));
            // binaryGalaxyOutputFile.write((char *)&galaxy_galaxyType[i_g], sizeof(std::size_t));
            binaryGalaxyOutputFile.write((char *)&galaxy_haloMass[i_g], sizeof(double));
            // binaryGalaxyOutputFile.write((char *)&galaxy_maxVelocity[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_xCoordinate[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_yCoordinate[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_zCoordinate[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_xVelocity[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_yVelocity[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_zVelocity[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_stellarMass[i_g], sizeof(double));
            // binaryGalaxyOutputFile.write((char *)&galaxy_meanStarAge[i_g], sizeof(double));
        }

        binaryGalaxyOutputFile.close();
    }
    else
    {
        binaryGalaxyInputFile.read((char *)&galaxyNumber, sizeof(std::size_t));

        // galaxy_rowId.resize(galaxyNumber);
        galaxy_rockstarId.resize(galaxyNumber);
        galaxy_hostHaloId.resize(galaxyNumber);
        galaxy_mainHaloId.resize(galaxyNumber);
        // galaxy_galaxyType.resize(galaxyNumber);
        galaxy_haloMass.resize(galaxyNumber);
        // galaxy_maxVelocity.resize(galaxyNumber);
        galaxy_xCoordinate.resize(galaxyNumber);
        galaxy_yCoordinate.resize(galaxyNumber);
        galaxy_zCoordinate.resize(galaxyNumber);
        galaxy_xVelocity.resize(galaxyNumber);
        galaxy_yVelocity.resize(galaxyNumber);
        galaxy_zVelocity.resize(galaxyNumber);
        galaxy_stellarMass.resize(galaxyNumber);
        // galaxy_meanStarAge.resize(galaxyNumber);

        for (std::size_t i_g = 0; i_g < galaxyNumber; ++i_g)
        {
            // binaryGalaxyInputFile.read((char *)&galaxy_rowId[i_g], sizeof(std::size_t));
            binaryGalaxyInputFile.read((char *)&galaxy_rockstarId[i_g], sizeof(std::size_t));
            binaryGalaxyInputFile.read((char *)&galaxy_hostHaloId[i_g], sizeof(std::size_t));
            binaryGalaxyInputFile.read((char *)&galaxy_mainHaloId[i_g], sizeof(std::size_t));
            // binaryGalaxyInputFile.read((char *)&galaxy_galaxyType[i_g], sizeof(std::size_t));
            binaryGalaxyInputFile.read((char *)&galaxy_haloMass[i_g], sizeof(double));
            // binaryGalaxyInputFile.read((char *)&galaxy_maxVelocity[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_xCoordinate[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_yCoordinate[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_zCoordinate[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_xVelocity[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_yVelocity[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_zVelocity[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_stellarMass[i_g], sizeof(double));
            // binaryGalaxyInputFile.read((char *)&galaxy_meanStarAge[i_g], sizeof(double));
        }
    }

    binaryGalaxyInputFile.close();

    return galaxy_xCoordinate.size();
}

std::size_t read_MDPL2_galaxy_sage_data_new(const std::string &fileName,
                                            std::vector<std::size_t> &galaxy_galaxyId,
                                            std::vector<std::size_t> &galaxy_rockstarId,
                                            std::vector<std::size_t> &galaxy_mainHaloId,
                                            std::vector<double> &galaxy_haloMass,
                                            std::vector<double> &galaxy_xCoordinate,
                                            std::vector<double> &galaxy_yCoordinate,
                                            std::vector<double> &galaxy_zCoordinate,
                                            std::vector<double> &galaxy_xVelocity,
                                            std::vector<double> &galaxy_yVelocity,
                                            std::vector<double> &galaxy_zVelocity,
                                            std::vector<double> &galaxy_stellarMass)
{
    galaxy_galaxyId.clear();
    galaxy_rockstarId.clear();
    galaxy_mainHaloId.clear();
    galaxy_haloMass.clear();
    galaxy_xCoordinate.clear();
    galaxy_yCoordinate.clear();
    galaxy_zCoordinate.clear();
    galaxy_xVelocity.clear();
    galaxy_yVelocity.clear();
    galaxy_zVelocity.clear();
    galaxy_stellarMass.clear();

    std::size_t galaxyNumber, rowId, galaxyId, rockstarId, hostHaloId, mainHaloId, galaxyType;
    double haloMass, maxVelocity, xCoordinate, yCoordinate, zCoordinate, xVelocity, yVelocity, zVelocity, stellarMass, meanStarAge;

    const std::string binaryFileName = std::string(fileName).replace(fileName.find_last_of('.') + 1, std::string::npos, "bin");

    std::ifstream binaryGalaxyInputFile(binaryFileName, std::ios::binary);

    if (binaryGalaxyInputFile.fail())
    {
        std::ifstream formattedGalaxyInputFile(fileName, std::ios::in);

        if (formattedGalaxyInputFile.fail())
        {
            std::cout << std::endl
                      << " Error: Formatted galaxy input file can not be opened" << std::endl
                      << std::endl;

            exit(EXIT_FAILURE);
        }

        std::string line;

        std::getline(formattedGalaxyInputFile, line); // header

        char delimiter;

        while (std::getline(formattedGalaxyInputFile, line))
        {
            std::istringstream lineStream(line);

            std::string ignore;

            lineStream >>
                rowId >> delimiter >>
                galaxyId >> delimiter >>
                rockstarId >> delimiter >>
                hostHaloId >> delimiter >>
                mainHaloId >> delimiter >>
                galaxyType >> delimiter >>
                haloMass >> delimiter >>
                maxVelocity >> delimiter >>
                xCoordinate >> delimiter >>
                yCoordinate >> delimiter >>
                zCoordinate >> delimiter >>
                xVelocity >> delimiter >>
                yVelocity >> delimiter >>
                zVelocity >> delimiter >>
                stellarMass >> delimiter >>
                meanStarAge;

            galaxy_galaxyId.push_back(galaxyId);
            galaxy_rockstarId.push_back(rockstarId);
            galaxy_mainHaloId.push_back(mainHaloId);
            galaxy_haloMass.push_back(haloMass / HubbleSim);
            galaxy_xCoordinate.push_back(xCoordinate);
            galaxy_yCoordinate.push_back(yCoordinate);
            galaxy_zCoordinate.push_back(zCoordinate);
            galaxy_xVelocity.push_back(xVelocity);
            galaxy_yVelocity.push_back(yVelocity);
            galaxy_zVelocity.push_back(zVelocity);
            galaxy_stellarMass.push_back(stellarMass / HubbleSim);
        }

        formattedGalaxyInputFile.close();

        galaxyNumber = galaxy_xCoordinate.size();

        std::ofstream binaryGalaxyOutputFile(binaryFileName, std::ios::binary);

        binaryGalaxyOutputFile.write((char *)&galaxyNumber, sizeof(std::size_t));

        for (std::size_t i_g = 0; i_g < galaxyNumber; ++i_g)
        {
            binaryGalaxyOutputFile.write((char *)&galaxy_galaxyId[i_g], sizeof(std::size_t));
            binaryGalaxyOutputFile.write((char *)&galaxy_rockstarId[i_g], sizeof(std::size_t));
            binaryGalaxyOutputFile.write((char *)&galaxy_mainHaloId[i_g], sizeof(std::size_t));
            binaryGalaxyOutputFile.write((char *)&galaxy_haloMass[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_xCoordinate[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_yCoordinate[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_zCoordinate[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_xVelocity[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_yVelocity[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_zVelocity[i_g], sizeof(double));
            binaryGalaxyOutputFile.write((char *)&galaxy_stellarMass[i_g], sizeof(double));
        }

        binaryGalaxyOutputFile.close();
    }
    else
    {
        binaryGalaxyInputFile.read((char *)&galaxyNumber, sizeof(std::size_t));

        galaxy_galaxyId.resize(galaxyNumber);
        galaxy_rockstarId.resize(galaxyNumber);
        galaxy_mainHaloId.resize(galaxyNumber);
        galaxy_haloMass.resize(galaxyNumber);
        galaxy_xCoordinate.resize(galaxyNumber);
        galaxy_yCoordinate.resize(galaxyNumber);
        galaxy_zCoordinate.resize(galaxyNumber);
        galaxy_xVelocity.resize(galaxyNumber);
        galaxy_yVelocity.resize(galaxyNumber);
        galaxy_zVelocity.resize(galaxyNumber);
        galaxy_stellarMass.resize(galaxyNumber);

        for (std::size_t i_g = 0; i_g < galaxyNumber; ++i_g)
        {
            binaryGalaxyInputFile.read((char *)&galaxy_galaxyId[i_g], sizeof(std::size_t));
            binaryGalaxyInputFile.read((char *)&galaxy_rockstarId[i_g], sizeof(std::size_t));
            binaryGalaxyInputFile.read((char *)&galaxy_mainHaloId[i_g], sizeof(std::size_t));
            binaryGalaxyInputFile.read((char *)&galaxy_haloMass[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_xCoordinate[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_yCoordinate[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_zCoordinate[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_xVelocity[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_yVelocity[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_zVelocity[i_g], sizeof(double));
            binaryGalaxyInputFile.read((char *)&galaxy_stellarMass[i_g], sizeof(double));
        }
    }

    binaryGalaxyInputFile.close();

    return galaxy_xCoordinate.size();
}

std::size_t read_MDPL2_halo_data(const std::string &fileName,
                                 std::vector<std::size_t> &halo_rowId,
                                 std::vector<std::size_t> &halo_rockstarId,
                                 std::vector<std::intmax_t> &halo_pId,
                                 std::vector<std::intmax_t> &halo_upId,
                                 std::vector<double> &halo_virialMass,
                                 std::vector<double> &halo_virialRadius,
                                 std::vector<double> &halo_scaleRadius,
                                 std::vector<double> &halo_rmsVelocity,
                                 std::vector<double> &halo_maxVelocity,
                                 std::vector<double> &halo_xCoordinate,
                                 std::vector<double> &halo_yCoordinate,
                                 std::vector<double> &halo_zCoordinate,
                                 std::vector<double> &halo_xVelocity,
                                 std::vector<double> &halo_yVelocity,
                                 std::vector<double> &halo_zVelocity,
                                 std::vector<double> &halo_200cMass,
                                 std::vector<double> &halo_500cMass,
                                 std::vector<double> &halo_accretionMass,
                                 std::vector<double> &halo_peakMass,
                                 std::vector<double> &halo_accretionVelocity,
                                 std::vector<double> &halo_peakVelocity,
                                 std::vector<double> &halo_peakMassMaxVelocity)
{
    halo_rowId.clear();
    halo_rockstarId.clear();
    halo_pId.clear();
    halo_upId.clear();
    halo_virialMass.clear();
    halo_virialRadius.clear();
    halo_scaleRadius.clear();
    halo_rmsVelocity.clear();
    halo_maxVelocity.clear();
    halo_xCoordinate.clear();
    halo_yCoordinate.clear();
    halo_zCoordinate.clear();
    halo_xVelocity.clear();
    halo_yVelocity.clear();
    halo_zVelocity.clear();
    halo_200cMass.clear();
    halo_500cMass.clear();
    halo_accretionMass.clear();
    halo_peakMass.clear();
    halo_accretionVelocity.clear();
    halo_peakVelocity.clear();
    halo_peakMassMaxVelocity.clear();

    std::size_t haloNumber, halo_rowId_value, halo_rockstarId_value;
    std::intmax_t halo_pId_value, halo_upId_value;
    double halo_virialMass_value, halo_virialRadius_value, halo_scaleRadius_value, halo_rmsVelocity_value, halo_maxVelocity_value, halo_xCoordinate_value, halo_yCoordinate_value, halo_zCoordinate_value, halo_xVelocity_value, halo_yVelocity_value, halo_zVelocity_value, halo_200cMass_value, halo_500cMass_value, halo_accretionMass_value, halo_peakMass_value, halo_accretionVelocity_value, halo_peakVelocity_value, halo_peakMassMaxVelocity_value;

    const std::string binaryFileName = std::string(fileName).replace(fileName.find_last_of('.') + 1, std::string::npos, "bin");

    std::ifstream binaryHaloInputFile(binaryFileName, std::ios::binary);

    if (binaryHaloInputFile.fail())
    {
        std::ifstream formattedHaloInputFile(fileName, std::ios::in);

        if (formattedHaloInputFile.fail())
        {
            std::cout << std::endl
                      << " Error: Formatted halo input file can not be opened" << std::endl
                      << std::endl;

            exit(EXIT_FAILURE);
        }

        std::string line;

        std::getline(formattedHaloInputFile, line); // header

        char delimiter;

        while (std::getline(formattedHaloInputFile, line))
        {
            std::istringstream lineStream(line);

            lineStream >>
                halo_rowId_value >> delimiter >>
                halo_rockstarId_value >> delimiter >>
                halo_pId_value >> delimiter >>
                halo_upId_value >> delimiter >>
                halo_virialMass_value >> delimiter >>
                halo_virialRadius_value >> delimiter >>
                halo_scaleRadius_value >> delimiter >>
                halo_rmsVelocity_value >> delimiter >>
                halo_maxVelocity_value >> delimiter >>
                halo_xCoordinate_value >> delimiter >>
                halo_yCoordinate_value >> delimiter >>
                halo_zCoordinate_value >> delimiter >>
                halo_xVelocity_value >> delimiter >>
                halo_yVelocity_value >> delimiter >>
                halo_zVelocity_value >> delimiter >>
                halo_200cMass_value >> delimiter >>
                halo_500cMass_value >> delimiter >>
                halo_accretionMass_value >> delimiter >>
                halo_peakMass_value >> delimiter >>
                halo_accretionVelocity_value >> delimiter >>
                halo_peakVelocity_value >> delimiter >>
                halo_peakMassMaxVelocity_value;

            // halo_rowId.push_back(halo_rowId_value);
            halo_rockstarId.push_back(halo_rockstarId_value);
            halo_pId.push_back(halo_pId_value);
            halo_upId.push_back(halo_upId_value);
            halo_virialMass.push_back(halo_virialMass_value / HubbleSim);
            halo_virialRadius.push_back(halo_virialRadius_value);
            // halo_scaleRadius.push_back(halo_scaleRadius_value);
            halo_rmsVelocity.push_back(halo_rmsVelocity_value);
            // halo_maxVelocity.push_back(halo_maxVelocity_value);
            halo_xCoordinate.push_back(halo_xCoordinate_value);
            halo_yCoordinate.push_back(halo_yCoordinate_value);
            halo_zCoordinate.push_back(halo_zCoordinate_value);
            halo_xVelocity.push_back(halo_xVelocity_value);
            halo_yVelocity.push_back(halo_yVelocity_value);
            halo_zVelocity.push_back(halo_zVelocity_value);
            halo_200cMass.push_back(halo_200cMass_value / HubbleSim);
            // halo_500cMass.push_back(halo_500cMass_value / HubbleSim);
            // halo_accretionMass.push_back(halo_accretionMass_value / HubbleSim);
            // halo_peakMass.push_back(halo_peakMass_value / HubbleSim);
            // halo_accretionVelocity.push_back(halo_accretionVelocity_value);
            // halo_peakVelocity.push_back(halo_peakVelocity_value);
            // halo_peakMassMaxVelocity.push_back(halo_peakMassMaxVelocity_value);

            formattedHaloInputFile.peek();
        }

        formattedHaloInputFile.close();

        haloNumber = halo_xCoordinate.size();

        std::ofstream binaryHaloOutputFile(binaryFileName, std::ios::binary);

        binaryHaloOutputFile.write((char *)&haloNumber, sizeof(std::size_t));

        for (std::size_t i_h = 0; i_h < haloNumber; ++i_h)
        {
            // binaryHaloOutputFile.write((char *)&halo_rowId[i_h], sizeof(std::size_t));
            binaryHaloOutputFile.write((char *)&halo_rockstarId[i_h], sizeof(std::size_t));
            binaryHaloOutputFile.write((char *)&halo_pId[i_h], sizeof(std::intmax_t));
            binaryHaloOutputFile.write((char *)&halo_upId[i_h], sizeof(std::intmax_t));
            binaryHaloOutputFile.write((char *)&halo_virialMass[i_h], sizeof(double));
            binaryHaloOutputFile.write((char *)&halo_virialRadius[i_h], sizeof(double));
            // binaryHaloOutputFile.write((char *)&halo_scaleRadius[i_h], sizeof(double));
            binaryHaloOutputFile.write((char *)&halo_rmsVelocity[i_h], sizeof(double));
            // binaryHaloOutputFile.write((char *)&halo_maxVelocity[i_h], sizeof(double));
            binaryHaloOutputFile.write((char *)&halo_xCoordinate[i_h], sizeof(double));
            binaryHaloOutputFile.write((char *)&halo_yCoordinate[i_h], sizeof(double));
            binaryHaloOutputFile.write((char *)&halo_zCoordinate[i_h], sizeof(double));
            binaryHaloOutputFile.write((char *)&halo_xVelocity[i_h], sizeof(double));
            binaryHaloOutputFile.write((char *)&halo_yVelocity[i_h], sizeof(double));
            binaryHaloOutputFile.write((char *)&halo_zVelocity[i_h], sizeof(double));
            binaryHaloOutputFile.write((char *)&halo_200cMass[i_h], sizeof(double));
            // binaryHaloOutputFile.write((char *)&halo_500cMass[i_h], sizeof(double));
            // binaryHaloOutputFile.write((char *)&halo_accretionMass[i_h], sizeof(double));
            // binaryHaloOutputFile.write((char *)&halo_peakMass[i_h], sizeof(double));
            // binaryHaloOutputFile.write((char *)&halo_accretionVelocity[i_h], sizeof(double));
            // binaryHaloOutputFile.write((char *)&halo_peakVelocity[i_h], sizeof(double));
            // binaryHaloOutputFile.write((char *)&halo_peakMassMaxVelocity[i_h], sizeof(double));
        }

        binaryHaloOutputFile.close();
    }
    else
    {
        binaryHaloInputFile.read((char *)&haloNumber, sizeof(std::size_t));

        // halo_rowId.resize(haloNumber);
        halo_rockstarId.resize(haloNumber);
        halo_pId.resize(haloNumber);
        halo_upId.resize(haloNumber);
        halo_virialMass.resize(haloNumber);
        halo_virialRadius.resize(haloNumber);
        // halo_scaleRadius.resize(haloNumber);
        halo_rmsVelocity.resize(haloNumber);
        // halo_maxVelocity.resize(haloNumber);
        halo_xCoordinate.resize(haloNumber);
        halo_yCoordinate.resize(haloNumber);
        halo_zCoordinate.resize(haloNumber);
        halo_xVelocity.resize(haloNumber);
        halo_yVelocity.resize(haloNumber);
        halo_zVelocity.resize(haloNumber);
        halo_200cMass.resize(haloNumber);
        // halo_500cMass.resize(haloNumber);
        // halo_accretionMass.resize(haloNumber);
        // halo_peakMass.resize(haloNumber);
        // halo_accretionVelocity.resize(haloNumber);
        // halo_peakVelocity.resize(haloNumber);
        // halo_peakMassMaxVelocity.resize(haloNumber);

        for (std::size_t i_h = 0; i_h < haloNumber; ++i_h)
        {
            // binaryHaloInputFile.read((char *)&halo_rowId[i_h], sizeof(std::size_t));
            binaryHaloInputFile.read((char *)&halo_rockstarId[i_h], sizeof(std::size_t));
            binaryHaloInputFile.read((char *)&halo_pId[i_h], sizeof(std::intmax_t));
            binaryHaloInputFile.read((char *)&halo_upId[i_h], sizeof(std::intmax_t));
            binaryHaloInputFile.read((char *)&halo_virialMass[i_h], sizeof(double));
            binaryHaloInputFile.read((char *)&halo_virialRadius[i_h], sizeof(double));
            // binaryHaloInputFile.read((char *)&halo_scaleRadius[i_h], sizeof(double));
            binaryHaloInputFile.read((char *)&halo_rmsVelocity[i_h], sizeof(double));
            // binaryHaloInputFile.read((char *)&halo_maxVelocity[i_h], sizeof(double));
            binaryHaloInputFile.read((char *)&halo_xCoordinate[i_h], sizeof(double));
            binaryHaloInputFile.read((char *)&halo_yCoordinate[i_h], sizeof(double));
            binaryHaloInputFile.read((char *)&halo_zCoordinate[i_h], sizeof(double));
            binaryHaloInputFile.read((char *)&halo_xVelocity[i_h], sizeof(double));
            binaryHaloInputFile.read((char *)&halo_yVelocity[i_h], sizeof(double));
            binaryHaloInputFile.read((char *)&halo_zVelocity[i_h], sizeof(double));
            binaryHaloInputFile.read((char *)&halo_200cMass[i_h], sizeof(double));
            // binaryHaloInputFile.read((char *)&halo_500cMass[i_h], sizeof(double));
            // binaryHaloInputFile.read((char *)&halo_accretionMass[i_h], sizeof(double));
            // binaryHaloInputFile.read((char *)&halo_peakMass[i_h], sizeof(double));
            // binaryHaloInputFile.read((char *)&halo_accretionVelocity[i_h], sizeof(double));
            // binaryHaloInputFile.read((char *)&halo_peakVelocity[i_h], sizeof(double));
            // binaryHaloInputFile.read((char *)&halo_peakMassMaxVelocity[i_h], sizeof(double));
        }
    }

    binaryHaloInputFile.close();

    return halo_xCoordinate.size();
}

CoordinateSystemRotation MDPL2_mock_galactic_coordinate_system_change(const double thetaMockVirgo, const double phiMockVirgo, const double thetaMockLGVelocity, const double phiMockLGVelocity)
{
    constexpr double latitudeRealVirgo = 74.436845; // NASA/IPAC Extragalactic Database
    constexpr double longitudeRealVirgo = 283.807431;

    double thetaRealVirgo, phiRealVirgo;

    transform_celestial_to_spherical_coordinates(latitudeRealVirgo, longitudeRealVirgo,
                                                 thetaRealVirgo, phiRealVirgo);

    const double thetaRealLGVelocity = CMB_TO_LOCAL_GROUP.relativeVelocityTheta;
    const double phiRealLGVelocity = CMB_TO_LOCAL_GROUP.relativeVelocityPhi;

    return reference_vector_and_plane_coordinate_system_rotation(thetaMockLGVelocity, phiMockLGVelocity, thetaMockVirgo, phiMockVirgo,
                                                                 thetaRealLGVelocity, phiRealLGVelocity, thetaRealVirgo, phiRealVirgo);
}

// double evaluate_field_in_sphere_around_point(const Cartesian3DGridFunction &field, const double xCenter, const double yCenter, const double zCenter,
//                                              const double radius, const double theta, const double phi)
// {
//     double x, y, z;

//     transform_spherical_to_cartesian_coordinates(radius, theta, phi,
//                                                  x, y, z);

//     const double xShifted = x + xCenter - BoxLength * std::floor((x + xCenter) / BoxLength);
//     const double yShifted = y + yCenter - BoxLength * std::floor((y + yCenter) / BoxLength);
//     const double zShifted = z + zCenter - BoxLength * std::floor((z + zCenter) / BoxLength);

//     return field(xShifted, yShifted, zShifted);
// }

double evaluate_within_sphere_in_periodic_box(const double xCenter, const double yCenter, const double zCenter, const double boxLength,
                                              const double radius, const double theta, const double phi,
                                              const Cartesian3DGridFunction &func)
//   const std::function<double(double x, double y, double z)> &func)
{
    double x, y, z;

    transform_spherical_to_cartesian_coordinates(radius, theta, phi,
                                                 x, y, z);

    const double xShifted = x + xCenter - boxLength * std::floor((x + xCenter) / boxLength);
    const double yShifted = y + yCenter - boxLength * std::floor((y + yCenter) / boxLength);
    const double zShifted = z + zCenter - boxLength * std::floor((z + zCenter) / boxLength);

    return func(xShifted, yShifted, zShifted);
}