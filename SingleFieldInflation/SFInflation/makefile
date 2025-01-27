COMPILER = g++
FLAGS = -O3 -fopenmp #-fstack-usage
HEADA = latticeeasy.hpp model.hpp parameters.hpp
all: $(HEADA) element_out.hpp element_evo.hpp ffteasy.o latticeeasy.o evolution.o initialize.o output.o
	$(COMPILER) $(FLAGS) ffteasy.o latticeeasy.o evolution.o initialize.o output.o -lm -o latticeeasy

clean:
	rm -f *.o out output.txt latticeeasy *~ *.dat core*

cleaner:
	rm -f *.o out output.txt latticeeasy *~ *.dat *.img

ffteasy.o: ffteasy.cpp
	$(COMPILER) -c $(FLAGS) ffteasy.cpp

latticeeasy.o: latticeeasy.cpp $(HEADA)
	$(COMPILER) -c $(FLAGS) latticeeasy.cpp

evolution.o: evolution.cpp latticeeasy.cpp $(HEADA) element_evo.hpp
	$(COMPILER) -c $(FLAGS) evolution.cpp

initialize.o: initialize.cpp latticeeasy.cpp $(HEADA)
	$(COMPILER) -c $(FLAGS) initialize.cpp

output.o: output.cpp latticeeasy.cpp $(HEADA) element_out.hpp
	$(COMPILER) -c $(FLAGS) output.cpp


