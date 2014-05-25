CXX = g++
CXXFLAGS= -Wall -O
LINKFLAGS = -lpthread
DEBUG=
OBJECTS = ErrorCorrection.o KmerCode.o GetKmers.o

ifneq (,$(findstring MINGW,$(shell uname)))
	LINKFLAGS = -L. -lpthreadGC2
endif

all: lighter

lighter: main.o $(OBJECTS)
	$(CXX) -o $@ $(CXXFLAGS) $(LINKFLAGS) $(OBJECTS) main.o 

main.o: main.cpp utils.h Reads.h Store.h bloom_filter.hpp
ErrorCorrection.o: ErrorCorrection.cpp ErrorCorrection.h utils.h
KmerCode.o: KmerCode.cpp KmerCode.h
GetKmers.o: GetKmers.cpp GetKmers.h

clean:
	rm -f *.o *.gch lighter
