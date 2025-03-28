COMPILER = g++
FLAGS = -O3 -fopenmp #-fstack-usage
HEADA = latticeeasy.hpp model.hpp parameters.hpp polarization.hpp
SOURCE = ffteasy.o latticeeasy.o evolution.o initialize.o output.o
all: $(HEADA) spectra.hpp element_evo.hpp element_out.hpp $(SOURCE)
	$(COMPILER) $(FLAGS) $(SOURCE) -lm -o latticeeasy

clean:
	rm -f *.o out output.txt latticeeasy *~ *.dat core*
cleaner:
	rm -f *.o out output.txt latticeeasy *~ *.dat *.img

ffteasy.o: ffteasy.cpp
	$(COMPILER) -c $(FLAGS) ffteasy.cpp
latticeeasy.o: latticeeasy.cpp $(HEADA)
	$(COMPILER) -c $(FLAGS) latticeeasy.cpp
evolution.o: evolution.cpp $(HEADA) element_evo.hpp
	$(COMPILER) -c $(FLAGS) evolution.cpp
initialize.o: initialize.cpp $(HEADA)
	$(COMPILER) -c $(FLAGS) initialize.cpp
output.o: output.cpp $(HEADA) element_out.hpp spectra.hpp
	$(COMPILER) -c $(FLAGS) output.cpp