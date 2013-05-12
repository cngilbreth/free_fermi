FC = gfortran -cpp -Wall -Wextra -g -fbounds-check

free_fermi: mpfun90.o mpmod90.o main.o
	$(FC) -o $@ $+

free_fermi_proj: main_proj.o
	$(FC) -o $@ $+

free_fermi_gc: main_gc.o
	$(FC) -o $@ $+

clean:
	rm -f *.o *.mod

%.o: %.f90
	$(FC) -c $< -o $@
