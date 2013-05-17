#FC = gfortran -cpp -Wall -Wextra
FC = gfortran -cpp -Wall -Wextra -g -fbounds-check

free_fermi: free_fermi.o
	$(FC) -o $@ $+

free_fermi_gc: free_fermi_gc.o
	$(FC) -o $@ $+

test: free_fermi
	cd ../test && ./test.sh

clean:
	rm -f *.o *.mod

%.o: %.f90
	$(FC) -c $< -o $@
