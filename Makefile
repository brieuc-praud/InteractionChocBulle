EXE = chaleur
F90 = gfortran
#OPT = -O0 -pedantic -Wall
#OPT = -O2
OPT = -g -ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all -Wall
OBJ = mod_params.o mod_fonctions.o mod_sorties.o $(EXE).o

$(EXE)	: $(OBJ)
	$(F90) $(OPT) -o $(EXE) $^

%.o	: %.f90
	$(F90) $(OPT) -c $<

%.mod	: %.f90
	$(F90) $(OPT) -c $<

clean	:
	rm *.o *.mod $(EXE) *.vtk sortie.dat

exe	: $(EXE)
	./$(EXE)