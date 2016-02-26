CXX = g++
CXXFLAGS= -Wall -O3
LINKFLAGS = -lpthread -lz 
DEBUG=
OBJECTS = ErrorCorrection.o GetKmers.o

# For Windows pthreads library: http://www.sourceware.org/pthreads-win32/
ifneq (,$(findstring MINGW,$(shell uname)))
	LINKFLAGS = -L. -lpthreadGC2
endif

all: lighter

lighter: main.o $(OBJECTS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)

main.o: main.cpp utils.h Reads.h Store.h File.h KmerCode.h bloom_filter.hpp
ErrorCorrection.o: ErrorCorrection.cpp ErrorCorrection.h utils.h
GetKmers.o: GetKmers.cpp GetKmers.h utils.h

clean:
	rm -f *.o *.gch lighter
