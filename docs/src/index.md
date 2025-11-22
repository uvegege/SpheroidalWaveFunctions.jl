# SpheroidalWaveFunctions.jl Documentation

A Julia implementation for computing angular and radial spheroidal wave functions and their associated eigenvalues.

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

## Background

Spheroidal wave functions are solutions to the Helmholtz equation in spheroidal coordinates. They appear in electromagnetic and acoustic scattering problems, among other applications.

This package was developed to support [AnalyticEMModes.jl](https://github.com/uvegege/AnalyticEMModes.jl) and provides a native Julia implementation based on standard mathematical references.

## References

The implementation follows definitions and algorithms from:
- **NIST DLMF** Chapter 30: "Spheroidal Wave Functions"
- **Abramowitz & Stegun**: *Handbook of Mathematical Functions*, Chapter 21
- **Zhang & Jin** (1996): *Computation of Special Functions*
- **Falloon et al.** (2003): "Theory and computation of spheroidal wavefunctions"

See the [API Reference](@ref) for detailed function documentation.
