CXX = g++
CXXFLAGS= -Wall -O
LINKFLAGS = -lpthread
DEBUG=
OBJECTS = ErrorCorrection.o KmerCode.o GetKmers.o

all: lighter

lighter: main.o $(OBJECTS)
	$(CXX) -o $@ $(CXXFLAGS) $(LINKFLAGS) $(OBJECTS) main.o 

main.o: main.cpp utils.h Reads.h Store.h bloom_filter.hpp
#bloom_filter.o: bloom_filter.h 
#Store.o: Store.h bloom_filter.hpp
ErrorCorrection.o: ErrorCorrection.cpp ErrorCorrection.h utils.h
KmerCode.o: KmerCode.cpp KmerCode.h
GetKmers.o: GetKmers.cpp GetKmers.h
#BitTable.o: BitTable.cpp BitTable.h

clean:
	rm -f *.o *.gch lighter
