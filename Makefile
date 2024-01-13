EXE = euler
F90 = gfortran
OPT = -Og -pedantic -g -ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all -Wall
#OPT = -O2
OBJ = mod_parameters.o mod_input.o mod_limiters.o mod_fluxes.o  mod_cases.o mod_quadrature.o mod_output.o main.o


PARAMS = parameters.dat

$(EXE): $(OBJ)
	$(F90) $(OPT) -o $(EXE) $^

%.o: %.f90
	$(F90) $(OPT) -c $<

%.mod: %.f90
	$(F90) $(OPT) -c $<

clean:
	make clear
	rm -f $(EXE)
clear:
	rm -f *.o *.mod error.dat output/*

exe: $(EXE)
	rm -f output/*
	mkdir -p output
	./$(EXE) $(PARAMS)

run:
	make exe

plot_error: error.gp error.dat
	gnuplot -p $<
	
error.dat: $(EXE) convergence.sh
	./convergence.sh $<
