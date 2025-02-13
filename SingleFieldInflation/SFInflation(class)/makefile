COMPILER = g++
FLAGS = -O3 -fopenmp -fstack-usage
HEADA = latticeeasy.hpp model.hpp parameters.hpp
all: $(HEADA) element_out.hpp element_evo.hpp fft_ini.hpp fft_out.hpp latticeeasy.o evolution.o initialize.o output.o
	$(COMPILER) $(FLAGS) latticeeasy.o evolution.o initialize.o output.o -lm -o latticeeasy

clean:
	rm -f *.o out output.txt latticeeasy *~ *.dat core*

cleaner:
	rm -f *.o out output.txt latticeeasy *~ *.dat *.img
	

latticeeasy.o: latticeeasy.cpp $(HEADA)
	$(COMPILER) -c $(FLAGS) latticeeasy.cpp

evolution.o: evolution.cpp latticeeasy.cpp $(HEADA) element_evo.hpp
	$(COMPILER) -c $(FLAGS) evolution.cpp

initialize.o: initialize.cpp latticeeasy.cpp $(HEADA) fft_ini.hpp
	$(COMPILER) -c $(FLAGS) initialize.cpp

output.o: output.cpp latticeeasy.cpp $(HEADA) element_out.hpp fft_out.hpp
	$(COMPILER) -c $(FLAGS) output.cpp


