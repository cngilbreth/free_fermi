# Noninteracting fermions in a harmonic trap

This project contains several codes for calculating thermodynamic properties of
noninteracting fermions in a harmonic trap in the *canonical ensemble* (of fixed
numbers of particles) and the grand-canonical ensemble.

The main file is free_fermi.f90, which utilizes an arbitrary-precision
arithmetic library contained in the files mpfun90.f90 and mpmod90.f90. 
free_fermi.f90 contains routines for calculating:

   - The partition function Z
   - The free energy F = E - T S = -T log(Z)
   - The thermal energy E
   - The heat capacity C

It has been tested by doing the calculations in two ways. First, free_fermi.f90
utilizes arbitrary precision arithmetic and a set for recursive formulas based on
the paper J. Chem. Phys. *98*, 2484 (1993). 

Second, free_fermi_proj.f90 does the same calculations (most of them) using a
particle-number projection method via numerical integration. The results for
reasonable numbers of particles (<= 80 or so) are identical to the recursive
method. (Everything but the heat capacity has been checked in this way.)

Finally, free_fermi_gc.f90 contains some routines for calculaing thermodynamic
quantities in the *grand-canonical ensemble*.



