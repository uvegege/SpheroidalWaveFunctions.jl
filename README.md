# SpheroidalWaveFunctions.jl

A Julia implementation for computing angular and radial spheroidal wave functions ($S_{mn}$, $R_{mn}^{(i)}$) and their associated eigenvalues ($\lambda_{mn}$).

## Background

This package was created to support the development of [AnalyticEMModes.jl](https://github.com/your-repo/AnalyticEMModes.jl), where spheroidal wave functions are needed for electromagnetic mode calculations. Since no native Julia implementation was available at the time, this package provides the necessary functionality based on standard references and established numerical methods.

## Features

- **Angular spheroidal wave functions**: Both Legendre expansion and power series methods
- **Radial spheroidal wave functions**: First and second kind
- **Characteristic values (eigenvalues)**: Individual and sequential computation
- **Prolate and oblate cases**: Support for both coordinate systems
- **Pure Julia implementation**: No external dependencies for core computations

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/your-username/SpheroidalWaveFunctions.jl")
```

## Quick Start

```julia
using SpheroidalWaveFunctions

# Compute angular function (prolate case)
m, n, c = 1, 2, 0.5
x = 0.3
S, dS = prolate_angular_leg(m, n, c, x)

# Compute characteristic value
λ = prolate_cv(m, n, c)

# Compute radial function (first kind)
ξ = 1.5
R, dR = prolate_radial1(m, n, c, ξ)

# Oblate case
S_obl, dS_obl = oblate_angular_ps(m, n, c, x)
λ_obl = oblate_cv(m, n, c)
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

## Mathematical References

The implementation follows standard definitions and algorithms from:

- **NIST Digital Library of Mathematical Functions (DLMF)**: Chapter 30, "Spheroidal Wave Functions"
- **Abramowitz & Stegun**: *Handbook of Mathematical Functions*, Chapter 21
- **Falloon, P. E., Abbott, P. C., & Wang, J. B.** (2003): "Theory and computation of spheroidal wavefunctions", *Journal of Physics A: Mathematical and General*, 36(20), 5477
- **Li, L. W., Kang, X. K., & Leong, M. S.** (2002): *Spheroidal Wave Functions in Electromagnetic Theory*, Wiley
- **Adelman, R., Gumerov, N. A., & Duraiswami, R.** (2014): "Software for Computing the Spheroidal Wave Functions Using Arbitrary Precision Arithmetic", arXiv:1408.0074

The numerical methods used include:
- Eigenvalue computation via tridiagonal matrix methods
- Newton-Raphson refinement with continued fractions
- Mixed forward-backward recursion for expansion coefficients

## Limitations and Known Issues

- This is a straightforward implementation prioritizing correctness over performance
- Not optimized for very large values of parameters or high precision requirements
- For production applications requiring extreme accuracy or performance, consider specialized libraries
- Limited testing coverage compared to mature packages like SciPy's `special` module

## Contributing

Contributions, bug reports, and suggestions are welcome. Please open an issue or submit a pull request.

## Related Packages

- **Python**: [scipy.special.pro_ang1](https://docs.scipy.org/doc/scipy/reference/special.html) and related functions
- **Mathematica**: Built-in `SpheroidalPS`, `SpheroidalQS` functions
- **Fortran**: Various implementations in SLATEC and elsewhere

## License

MIT License - see LICENSE file for details

