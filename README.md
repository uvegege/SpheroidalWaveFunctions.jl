# SpheroidalWaveFunctions.jl

A Julia implementation for computing angular and radial spheroidal wave functions ($S_{mn}$, $R_{mn}^{(i)}$) and their associated eigenvalues ($\lambda_{mn}$).

## Background

This package was created to support the development of [AnalyticEMModes.jl](https://github.com/uvegege/AnalyticEMModes.jl), where spheroidal wave functions are needed for electromagnetic mode calculations. Since no native Julia implementation was available at the time, this package provides the necessary functionality based on standard references and established numerical methods.

## Features

- **Angular spheroidal wave functions**: Both Legendre expansion and power series methods
- **Radial spheroidal wave functions**: First and second kind
- **Characteristic values (eigenvalues)**: Individual and sequence computation
- **Prolate and oblate cases**: Support for both coordinate systems
- **Pure Julia implementation**: No external dependencies for core computations

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/uvegege/SpheroidalWaveFunctions.jl")
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

## Acknowledgments

Some implementation details (such as matrix sizes and initial values) were informed by examining SciPy's routines and Zhang & Jin's Fortran code. The development benefited from examining SciPy's spheroidal wave function routines and Zhang & Jin's Fortran code for guidance on implementation details and numerical validation.

