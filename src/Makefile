#FC = gfortran -cpp -Wall -Wextra
FC = gfortran -cpp -Wall -Wextra -g -fbounds-check

free_fermi_recurse: mpfun90.o mpmod90.o free_fermi_recurse.o
	$(FC) -o $@ $+

free_fermi: free_fermi.o
	$(FC) -o $@ $+

free_fermi_gc: free_fermi_gc.o
	$(FC) -o $@ $+

clean:
	rm -f *.o *.mod

%.o: %.f90
	$(FC) -c $< -o $@
