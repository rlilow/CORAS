# CORAS

A framework for the fast generation of **CO**nstrained **R**ealizations from **A**ll-sky **S**urveys

## Description

CORAS reconstructs the density and velocity fields from a given all-sky redshift survey, using a Wiener filter estimator.
By combining this estimator with random log-normal Poisson realizations of the density field and galaxy distribution, constrained realizations can be generated, which sample the space of field realizations compatible with the observed data.
By comparing the reconstructed velocities with those obtained from a galaxy distance catalog, CORAS can furthermore infer the normalized growth rate of structures and the external bulk flow contribution from sources beyond the survey volume. As a nuisance parameter, it also fixes the Hubble constant used to translate observed distance moduli to velocities.

## Data

### Input data

The `data` directory contains all input data needed to apply CORAS to 2MRS (Huchra et al., [ApJS](https://iopscience.iop.org/article/10.1088/0067-0049/199/2/26) 199 (2012) 26; Macri et al., [ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab465a) 245 (2019) 6), and to compare the reconstructed velocities to those obtained from Cosmicflows-3 (Tully et al., [AJ](https://iopscience.iop.org/article/10.3847/0004-6256/152/2/50) 152 (2016) 50). The fiducial power spectrum used in the Wiener filter was computed with the [Cosmic Emu](https://github.com/lanl/CosmicEmu) (Heitmann et al., [AJ](https://iopscience.iop.org/article/10.3847/0004-637X/820/2/108) 820 (2016) 108).

### Reconstructed fields

For easy access, the normalized density contrast (density contrast divided by sigma_8) and peculiar velocity fields (with respect to the CMB frame) reconstructed with CORAS using the parameter values inferred from the velocity-velocity comparison are directly available in this [Dropbox folder](https://www.dropbox.com/sh/3nebvt1lskxshtu/AAByegavgA_-l1x118tZkaSAa?dl=0). These fields have been smoothed with a 5 Mpc/h Gaussian, and were reconstructed using observed redshifts in the LG frame. They are discretized on a regular 201^3 Cartesian grid in galactic coordinates, with each coordinate running from -200 Mpc/h to +200 Mpc/h in steps of 2 Mpc/h. The field values at the coordinates (x_i, y_j, z_k) with

> x/y/z_i = 2 \* (i - 100), 0 <= i <= 200

are written to the line

> l = (i \* 201 + j) \* 201 + k

of the file (excluding the header line and starting at l = 0).

These files can also be generated using the `compute_reconstructed_fields_on_cartesian_grid.x` executable, as described in more detail below in the section **Running the code**.

## Installation

### Install dependencies
Before compiling the code, the following dependencies need to be installed:

- compiler supporting at least C++11 (e.g. gcc 7 or clang 5)
- [GNU Scientific Library](https://www.gnu.org/software/gsl/) (GSL)
- [FFTW3](http://www.fftw.org/)

### Download and compilation

Clone the repository into the desired location and change into the root directory by running

```bash
git clone https://github.com/rlilow/CORAS.git
cd CORAS
```

If the GSL or FFTW3 are not in the standard paths, specify their include and library paths in the beginning of the `Makefile` by changing the corresponding variables `GSL_INCLUDE_PATH`, `GSL_LIB_PATH`, `FFTW3_INCLUDE_PATH` and `FFTW3_LIB_PATH`.
Afterwards run

```bash
make
```

This will create a number of executables in the `exe` directory.

## Documentation 

If you have Doxygen (https://www.doxygen.nl/index.html) installed, you can build a detailed documentation of the different classes and functions in CORAS by running

```bash
make doc
```

from within the root directory. This will create the directory `doc` containing the documentation.
To view it, open the file `doc/documentation.html` in the browser.

## Running the code

To run CORAS in its default settings, simply execute the desired `*.x` file in the `exe` directory, e.g.

```
exe/compute_reconstructed_fields_on_cartesian_grid.x
```

By default, all output is written to the `data` directory.

To change the settings, modify the parameters in `configuration.hpp` and re-run `make` before executing a `*.x` file.

### Executables 

The different executables in the `exe` directory can be run independently. They do the following:

- `compute_reconstructed_fields_on_cartesian_grid.x`  
  Compute the reconstructed density and velocity fields for given normalized growth rate and external bulk flow contribution on a Cartesian grid. The reconstructed fields described above in the **Data** section have been computed with this executable and the default settings.
    
- `analyze_reconstructed_fields.x`  
  Analyze the reconstructed density and velocity fields for given normalized growth rate and external bulk flow contribution. Among other things, this computes both fields on a slice through the supergalactic plane, and their standard deviations as a function of radius.
    
- `compare_reconstructed_and_observed_velocities.x`  
  Compare the reconstructed and observed radial velocities for given normalized growth rate, external bulk flow contribution and distance catalog Hubble constant. This computes the tensor-smoothed values of both velocity fields as well as their correlation functions.
    
- `estimate_parameters.x`  
  Compute the normalized growth rate, the external bulk flow contribution and the distance catalog Hubble constant via a maximum-likelihood comparison of the reconstructed and observed radial velocities.

## Development and contributing
The Doxygen documentation of the code is still incomplete and will be extended in the near future.

Contributions are very welcome! If you have any suggestions for new features or improvements, or if you encounter any errors, please create a GitHub issue.