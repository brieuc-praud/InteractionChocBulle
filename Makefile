EXE = euler
F90 = gfortran
OPT = -O3 -march=native
#OPT = -Og -pedantic -g -ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all -Wall
OBJ = mod_parameters.o mod_input.o mod_limiters.o mod_fluxes.o mod_cases.o mod_output.o main.o
LDFLAGS = 


PARAMS = parameters.dat

$(EXE): $(OBJ)
	$(F90) $(OPT) -o $(EXE) $^ $(LDFLAGS)

%.o: %.f90
	$(F90) $(OPT) -c $< $(LDFLAGS)

%.mod: %.f90
	$(F90) $(OPT) -c $< $(LDFLAGS)

clean:
	make clear
	rm -f *.o $(EXE)
clear:
	rm -f *.mod error.dat

exe: $(EXE)
	./$(EXE) $(PARAMS)

run:
	make exe

plot_error: error.gp error.dat
	gnuplot -p $<
	
error.dat: $(EXE) convergence.sh
	./convergence.sh $< $(PARAMS)
