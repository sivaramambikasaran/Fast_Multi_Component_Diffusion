<p align="center">
  <img src="https://github.com/shyams2/Fast_Multi_Component_Diffusion/blob/master/docs/source/images/FMDV.svg" width="250" height="250" alt="Logo of FMDV"/>
</p>

# FMDV

FMDV can be used to solve for multicomponent diffusion velocities in a fast manner. It does so by leveraging the low-rank nature of matrices that arise in the intermediate steps. The methods that are utilized in this library are detailed in [this](https://www.sciencedirect.com/science/article/pii/S1540748916300554) article

## Dependencies:

FMDV makes use of linear algebra routines from the Eigen library. Ensure that the root level of this library is included during the compilation process

## Usage

In order to use FMDV, include the ``FMDV.hpp`` file in your program. The solver requires various parameters that are required to solve the Stefan-Maxwell Equation:

<img src="https://cdn.jsdelivr.net/gh/shyams2/Fast_Multi_Component_Diffusion@master/docs/source/images/equation.svg?invert_in_darkmode" align=middle width=737.11538505pt height=106.33791629999999pt/>

As we can see above, this means that the inputs such as density, pressure, temperature and species properties of the various species are needed for the solver. The main solver object exposes these parameters needed as public attributes. They can be set by the constructor or changed directly through changing the attribute. We hope the examples serve to be explanatory on this front.

## Example

An examples for usage along with the benchmarks compared to naive methods can be found under ``examples``. This can be compiled by running ``make`` to build the executable using instructions from the associated ``Makefile``.

## License

This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license. The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at <http://mozilla.org/MPL/2.0/.>
