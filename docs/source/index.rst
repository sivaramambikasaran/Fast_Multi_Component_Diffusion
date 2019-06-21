.. role:: underline
    :class: underline

Welcome to FMDV's Documentation!
********************************

About :math:`\texttt{FMDV}`:
============================

:math:`\texttt{FMDV}` can be used to solve for multicomponent diffusion velocities in a fast manner. It does so by leveraging the low-rank nature of matrices that arise in the intermediate steps of solving the Stefan-Maxwell equations. The methods that are utilized in this library are detailed in `this <https://www.sciencedirect.com/science/article/pii/S1540748916300554>`_ article
The code is written in C++ and features an easy-to-use interface, where the user provides the  through a ``solver`` object which abstracts data in the inverse binary diffusivities matrix through a member function ``getInverseDiffusionCoefficient(int i, int j)`` which returns the entry at the :math:`i^{\mathrm{th}}` row and :math:`j^{\mathrm{th}}` column of the matrix.

Doc Contents
============
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   tutorial
   benchmarks

Other Links
===========

Learn more about :math:`\texttt{FMDV}` by visiting the

* Code Repository: http://github.com/sivaramambikasaran/Fast_Multi_Component_Diffusion
* Documentation: http://fmdv.rtfd.io
