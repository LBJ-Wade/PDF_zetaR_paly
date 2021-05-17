MODEL = USRzetaR
SOURCE = source

CXX := g++-10
CXXFLAGS := -std=c++11 -O3 -fopenmp

$(MODEL): $(MODEL).o $(SOURCE)/JacobiPDE.o $(SOURCE)/SRK32.o $(SOURCE)/StocDeltaN.o
	$(CXX) $(CXXFLAGS) -o $(MODEL) $(MODEL).o $(SOURCE)/JacobiPDE.o $(SOURCE)/SRK32.o $(SOURCE)/StocDeltaN.o

$(MODEL).o: $(SOURCE)/JacobiPDE.hpp $(SOURCE)/SRK32.hpp $(SOURCE)/StocDeltaN.hpp
JacobiPDE.o: $(SOURCE)/JacobiPDE.hpp
SRK32.o: $(SOURCE)/SRK32.hpp
StocDeltaN.o: $(SOURCE)/JacobiPDE.hpp $(SOURCE)/SRK32.hpp $(SOURCE)/StocDeltaN.hpp

clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o	
	$(RM) $(MODEL)
