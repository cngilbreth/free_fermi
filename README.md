# Noninteracting fermions in the canonical ensemble

This project contains several codes for calculating thermodynamic properties of
noninteracting fermions in the *canonical ensemble* (of fixed numbers of
particles) and the grand-canonical ensemble, in arbitrary dimensions.

The code currently is written for isotropic harmonic trapping potentials, but
can easily be changed to support other geometries.

The main file is free_fermi.f90, which contains routines for calculating:

   - The partition function Z
   - The free energy F = E - T S = -T log(Z)
   - The thermal energy E
   - The heat capacity C
   - The single-particle occupations n(k)

It has been tested by doing the calculations in two ways. First, free_fermi.f90
uses particle-number projection (via a Fourier transform) to calculate
canonical-ensemble quantities. Second, an alternate code (in the altsrc
directory) uses arbitrary precision arithmetic and a set for recursive formulas
based on the paper J. Chem. Phys. *98*, 2484 (1993). The two codes give
identical results.

Finally, free_fermi_gc.f90 contains some simple routines for calculating
thermodynamic quantities in the *grand-canonical ensemble*. (They are not as
robust or as extensive as the canonical routines.)



