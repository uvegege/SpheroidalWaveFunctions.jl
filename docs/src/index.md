# SpheroidalWaveFunctions.jl Documentation

A Julia implementation for computing angular and radial spheroidal wave functions and their associated eigenvalues.

## Background

Spheroidal wave functions are solutions to the Helmholtz equation in spheroidal coordinates. They appear in electromagnetic and acoustic scattering problems, among other applications.

This package was developed to support [AnalyticEMModes.jl](https://github.com/uvegege/AnalyticEMModes.jl), where these functions were needed for electromagnetic mode calculations. Since no native Julia implementation was available at the time (or at least none could be found), this package provides the necessary functionality based on standard mathematical references and established numerical methods.

The goal is not to provide cutting-edge or highly optimized algorithms, but simply to offer a working set of functionalities for practical use in Julia.


## Overview

This package provides functions to compute:
- **Angular spheroidal wave functions** $S_{mn}(c, x)$ and their derivatives
- **Radial spheroidal wave functions** $R_{mn}^{(i)}(c, \xi)$ and their derivatives
- **Characteristic values (eigenvalues)** $\lambda_{mn}(c)$

Both **prolate** and **oblate** cases are supported.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/uvegege/SpheroidalWaveFunctions.jl")
```

## Quick Example

```julia
using SpheroidalWaveFunctions

# Parameters
m, n, c = 1, 2, 0.5
x = 0.3
ξ = 1.5

# Compute angular function (prolate)
S, dS = prolate_angular_leg(m, n, c, x)

# Compute characteristic value
λ = prolate_cv(m, n, c)

# Compute radial function (first kind)
R, dR = prolate_radial1(m, n, c, ξ)
```


## API Reference

### Exported Functions

#### Angular Functions
- `prolate_angular_leg(m, n, c, x)` - Prolate angular function (Legendre expansion)
- `prolate_angular_ps(m, n, c, x)` - Prolate angular function (power series)
- `oblate_angular_leg(m, n, c, x)` - Oblate angular function (Legendre expansion)
- `oblate_angular_ps(m, n, c, x)` - Oblate angular function (power series)

#### Radial Functions
- `prolate_radial1(m, n, c, ξ)` - Prolate radial function, first kind
- `prolate_radial2(m, n, c, ξ)` - Prolate radial function, second kind
- `oblate_radial1(m, n, c, ξ)` - Oblate radial function, first kind
- `oblate_radial2(m, n, c, ξ)` - Oblate radial function, second kind

#### Characteristic Values
- `prolate_cv(m, n, c)` - Single prolate characteristic value
- `oblate_cv(m, n, c)` - Single oblate characteristic value
- `prolate_cv_seq(m, n, c)` - Sequence of prolate characteristic values
- `oblate_cv_seq(m, n, c)` - Sequence of oblate characteristic values

All angular and radial functions return `(value, derivative)` tuples.

See the [API Reference](@ref) for detailed function documentation.

## Limitations and Known Issues

- Limited testing coverage compared to mature packages like SciPy's `special` module
- May not handle all edge cases or extreme parameter values
- Numerical stability for very large parameters has not been extensively tested

## Contributing

Contributions are welcome and appreciated! Whether you find a bug, have an idea for improvement, or want to add new functionality, please feel free to:

- Open an issue to report bugs or discuss enhancements
- Submit a pull request with improvements or fixes
- Share your use cases and feedback

Any help in expanding test coverage, improving documentation, or optimizing performance is especially appreciated.

## Related Packages

- **Python**: [scipy.special.pro_ang1](https://docs.scipy.org/doc/scipy/reference/special.html) and related functions
- **Mathematica**: Built-in `SpheroidalPS`, `SpheroidalQS` functions
- **Fortran**: Various implementations in SLATEC and elsewhere


## Mathematical References

The implementation follows standard definitions and algorithms from:

- **NIST Digital Library of Mathematical Functions (DLMF)**: Chapter 30, "Spheroidal Wave Functions"
- **Abramowitz & Stegun**: *Handbook of Mathematical Functions*, Chapter 21
- **Zhang, S., & Jin, J.** (1996): *Computation of Special Functions*, Wiley (Fortran implementation used as reference by SciPy)
- **Falloon, P. E., Abbott, P. C., & Wang, J. B.** (2003): "Theory and computation of spheroidal wavefunctions", *Journal of Physics A: Mathematical and General*, 36(20), 5477
- **Li, L. W., Kang, X. K., & Leong, M. S.** (2002): *Spheroidal Wave Functions in Electromagnetic Theory*, Wiley
- **Adelman, R., Gumerov, N. A., & Duraiswami, R.** (2014): "Software for Computing the Spheroidal Wave Functions Using Arbitrary Precision Arithmetic", arXiv:1408.0074

The numerical methods implemented are standard approaches including:
- Eigenvalue computation via symmetric tridiagonal matrix methods
- Newton-Raphson refinement with continued fractions
- Mixed forward-backward recursion for expansion coefficients

Some implementation details (such as matrix sizes and initial values) were informed by examining SciPy's routines and Zhang & Jin's Fortran code for guidance on numerical stability and validation.


