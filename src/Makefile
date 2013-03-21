FC = gfortran -Wall -Wextra -g -O0 -fbounds-check

free_fermi: mpfun90.o mpmod90.o main.o
	$(FC) -o $@ $+

clean:
	rm -f *.o *.mod

%.o: %.f90
	$(FC) -c $< -o $@
