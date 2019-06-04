# Fast Multicomponent Diffusion

This a header library that can be used as the module to solve for diffusion velocities in a fast manner. The methods that are utilized in this library are detailed in this article.

## Dependencies:

Fast Multicomponent Diffusion makes use of linear algebra routines from the Eigen library. Ensure that the root level of this library is included during the compilation process

## Usage

In order to use the library, include the ``fast_multicomponent_diffusion.hpp`` file in your program. The solver requires various parameters that are required to solve the Maxwell Equation. This solver is capable of solving the following equation.

As we can see above, this means that the inputs such as density, pressure, temperature and properties of the various species are needed for the solver. The solver object gets initialized by passing a vector of ``species`` objects. An instance of the ``species`` class has the following attributes: name of species, mole fraction, mass fraction, thermal diffusivities... :

Density
Molecular Weight of Each Species
Thermal Diffusivities of Each Species
Mole Fraction of Each Species
Body Forces Acting Upon Each Species


## Example

An example for usage along with the benchmarks compared to naive methods can be found under ``example/main.cpp``. This can be compiled by running ``make`` to build the executable using instructions from the associated ``Makefile``.

## License

This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license. The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at <http://mozilla.org/MPL/2.0/.>
