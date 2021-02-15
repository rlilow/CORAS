# CORAS

A framework for the fast generation of **CO**nstrained **R**ealizations from **A**ll-sky **S**urveys

## Description

CORAS reconstructs the density and velocity fields from a given all-sky redshift survey, using a Wiener filter estimator.
By combining this estimator with random log-normal Poisson realizations of the density field and galaxy distribution, constrained realizations can be generated, which sample the space of field realizations compatible with the observed data.
By comparing the reconstructed velocities with those obtained from a galaxy distance catalog, CORAS can furthermore infer the normalized growth rate of structures and the external bulk flow contribution from sources beyond the survey volume.
As a nuisance parameter, it also fixes the Hubble constant used to translate observed distance moduli to velocities.

## Data

### Input data

The `data` directory contains all input data needed to apply CORAS to 2MRS and perform the velocity comparison with Cosmicflows-3.

2MRS is described in

> John P. Huchra et al., [ApJS](https://iopscience.iop.org/article/10.1088/0067-0049/199/2/26) 199 (2012) 26  
> Lucas M. Macri et al., [ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab465a) 245 (2019) 6

the 2MRS group catalog in

> R. Brent Tully, [AJ](https://iopscience.iop.org/article/10.1088/0004-6256/149/5/171) 149 (2015) 171

Cosmicflows-3 in

> R. Brent Tully et al., [AJ](https://iopscience.iop.org/article/10.3847/0004-6256/152/2/50) 152 (2016) 50

and the power spectrum emulator Cosmic Emu in

> Katrin Heitmann et al., [AJ](https://iopscience.iop.org/article/10.3847/0004-637X/820/2/108) 820 (2016) 108

The original 2MRS data are available alongside their publications, the 2MRS group catalog and Cosmicflows-3 are available in the [Extragalactic Distance Database](http://edd.ifa.hawaii.edu), and Cosmic Emu is available in [this repository](https://github.com/lanl/CosmicEmu).

### Reconstructed fields on a grid

The reconstructed normalized density contrast and peculiar velocity fields on a grid are available in this [Dropbox folder](https://www.dropbox.com/sh/3nebvt1lskxshtu/AAByegavgA_-l1x118tZkaSAa?dl=0).
The normalized density contrast is the density contrast divided by sigma_8.
The peculiar velocity field is with respect to the CMB frame, and its components point along the Galactic coordinate axes.
Both fields have been smoothed with a 5 Mpc/h Gaussian, and were reconstructed using observed redshifts in the Local Group (LG) frame.

They are discretized on a regular 201^3 Cartesian grid in comoving Galactic coordinates, with each coordinate running from -200 Mpc/h to +200 Mpc/h in steps of 2 Mpc/h.
The field values at the coordinates (x_i, y_j, z_k) with

> x/y/z_i = 2 \* (i - 100), 0 <= i <= 200

are written to the line

> l = (i \* 201 + j) \* 201 + k

of the file (excluding the header line and starting at l = 0).

These files can also be generated using the `compute_reconstructed_fields_on_cartesian_grid.x` executable, as described in more detail below in the section **Running the code**.

### Reconstructed velocities of Cosmicflows-3 galaxies and groups

The reconstructed peculiar velocities at the positions of the Cosmicflows-3 galaxies and groups are available in the [Extragalactic Distance Database](http://edd.ifa.hawaii.edu) in the table _Lilow-Nusser CF3 Peculiar Velocities_.
Both the galaxy and group velocities have been obtained by evaluating the velocity field at the comoving LG-frame redshift distance of the group.
This minimizes the Malmquist bias and reduces the contamination by incoherent small-scale motions (e.g. fingers-of-god).
Thus, the only difference between galaxy and group velocities is the angular position at which they have been evaluated.
The table also contains the uncertainty in the peculiar velocity components, estimated from the scatter between 50 constrained realizations.
This uncertainty is the same for each velocity component and only depends on distance.

A small number of galaxies has negative LG-frame group redshifts.
Since all of their actual distances are small compared to the 5 Mpc/h smoothing scale, their peculiar velocities have simply been evaluated at zero distance.
They additionally carry the flag -1.
Another small number of galaxies have LG-frame group redshift distances beyond the reconstruction boundary at 200 Mpc/h (corresponding to cz = 20,300 km/s).
They carry the flag 1, and their peculiar velocities are set to zero.
All remaining galaxies carry the flag 0.

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

from within the root directory.
This will create the directory `doc` containing the documentation.
To view it, open the file `doc/documentation.html` in the browser.

## Running the code

To run CORAS in its default settings, simply execute the desired `*.x` file in the `exe` directory, e.g.

```
exe/compute_reconstructed_fields_on_cartesian_grid.x
```

By default, all output is written to the `data` directory.

To change the settings, modify the parameters in `configuration.hpp` and re-run `make` before executing a `*.x` file.

### Executables 

The different executables in the `exe` directory can be run independently.
They do the following:

- `compute_reconstructed_fields_on_cartesian_grid.x`  
  Compute the reconstructed density and velocity fields for given normalized growth rate and external bulk flow contribution on a Cartesian grid.
  The reconstructed fields described above in the **Data** section have been computed with this executable and the default settings.
    
- `analyze_reconstructed_fields.x`  
  Analyze the reconstructed density and velocity fields for given normalized growth rate and external bulk flow contribution.
  Among other things, this computes both fields on a slice through the supergalactic plane, and their standard deviations as a function of radius.
    
- `compare_reconstructed_and_observed_velocities.x`  
  Compare the reconstructed and observed radial velocities for given normalized growth rate, external bulk flow contribution and distance catalog Hubble constant.
  This computes the tensor-smoothed values of both velocity fields as well as their correlation functions.
    
- `estimate_parameters.x`  
  Compute the normalized growth rate, the external bulk flow contribution and the distance catalog Hubble constant via a maximum-likelihood comparison of the reconstructed and observed radial velocities.

## Development and contributing
The Doxygen documentation of the code is still incomplete and will be extended in the near future.

Contributions are very welcome! If you have any suggestions for new features or improvements, or if you encounter any errors, please create a GitHub issue.