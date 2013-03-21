FC = gfortran -Wall -Wextra

free_fermi: mpfun90.o mpmod90.o main.o
	$(FC) -o $@ $+

clean:
	rm -f *.o *.mod

%.o: %.f90
	$(FC) -c $< -o $@
