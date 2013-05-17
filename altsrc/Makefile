#FC = gfortran -cpp -Wall -Wextra
FC = gfortran -cpp -Wall -Wextra -g -fbounds-check

free_fermi_recurse: mpfun90.o mpmod90.o free_fermi_recurse.o
	$(FC) -o $@ $+

test: free_fermi_recurse
	cd ../test && ./test_alt.sh

clean:
	rm -f *.o *.mod

%.o: %.f90
	$(FC) -c $< -o $@
