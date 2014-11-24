CXX = g++
CXXFLAGS= -Wall -O3
LINKFLAGS = -lpthread -lz
DEBUG=
OBJECTS = ErrorCorrection.o KmerCode.o GetKmers.o

# For Windows pthreads library: http://www.sourceware.org/pthreads-win32/
ifneq (,$(findstring MINGW,$(shell uname)))
	LINKFLAGS = -L. -lpthreadGC2
endif

all: lighter

lighter: main.o $(OBJECTS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)

main.o: main.cpp utils.h Reads.h Store.h File.h bloom_filter.hpp
ErrorCorrection.o: ErrorCorrection.cpp ErrorCorrection.h utils.h
KmerCode.o: KmerCode.cpp KmerCode.h
GetKmers.o: GetKmers.cpp GetKmers.h

clean:
	rm -f *.o *.gch lighter
