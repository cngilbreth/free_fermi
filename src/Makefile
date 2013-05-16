#FC = gfortran -cpp -Wall -Wextra
FC = gfortran -cpp -Wall -Wextra -g -fbounds-check

free_fermi: mpfun90.o mpmod90.o free_fermi.o
	$(FC) -o $@ $+

free_fermi_proj: free_fermi_proj.o
	$(FC) -o $@ $+

free_fermi_gc: free_fermi_gc.o
	$(FC) -o $@ $+

clean:
	rm -f *.o *.mod

%.o: %.f90
	$(FC) -c $< -o $@
