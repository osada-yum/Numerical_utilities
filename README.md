# Numerical_utilities
Items for numerical simulation implemented by Fortran.

## Purpose
Collect the algorithms for numerical computation.

## Usage
If you use fpm(Fortran package manager)[https://github.com/fortran-lang/fpm],
add
```
    [dependencies]
    numeric_simulation = { git = "https://github.com/osada-yum/Numerical_utilities.git" }
```
to your `fpm.toml`.

See `demo/demo_variance_covariance_kahan.f90` to understand usage of fortran modules.

# TODO

## Algorithm for correction to numerical error
[x] Kahan summation algorithm
[x] One-pass algorithm for variance and covariance computed by the Kahan summation algorithm
  [x] MPI subroutine for user-defined-type
[x] Welford's online algorithm
