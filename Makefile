EXE = euler
F90 = gfortran
#OPT = -Og -pedantic -g -ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all -Wall
OPT = -O2
OBJ = mod_parameters.o mod_functions.o mod_schemes.o mod_output.o mod_test.o $(EXE).o

$(EXE): $(OBJ)
	$(F90) $(OPT) -o $(EXE) $^

%.o: %.f90
	$(F90) $(OPT) -c $<

%.mod: %.f90
	$(F90) $(OPT) -c $<

clean:
	rm -f *.o *.mod $(EXE) *.vtk sortie.dat

exe: $(EXE)
	rm -f *.vtk
	./$(EXE)

plot_error: error.gp error.dat
	gnuplot -p $<
	
error.dat: $(EXE) convergence.sh
	./convergence.sh $<
