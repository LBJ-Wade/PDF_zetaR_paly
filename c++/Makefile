MODEL = USRzetaR
SOURCE = source

CXX := g++-11
CXXFLAGS := -std=c++11 -O3 -fopenmp

$(MODEL): $(MODEL).o $(SOURCE)/JacobiPDE.o #forchi.o # $(SOURCE)/SRK32.o $(SOURCE)/StocDeltaN.o
	$(CXX) $(CXXFLAGS) -o $(MODEL) $(MODEL).o $(SOURCE)/JacobiPDE.o #forchi.o # $(SOURCE)/SRK32.o $(SOURCE)/StocDeltaN.o

$(MODEL).o: $(SOURCE)/JacobiPDE.hpp #forchi.hpp # $(SOURCE)/SRK32.hpp $(SOURCE)/StocDeltaN.hpp
JacobiPDE.o: $(SOURCE)/JacobiPDE.hpp #forchi.o: $(SOURCE)/JacobiPDEforchi.hpp
#SRK32.o: $(SOURCE)/SRK32.hpp
#StocDeltaN.o: $(SOURCE)/JacobiPDE.hpp $(SOURCE)/SRK32.hpp $(SOURCE)/StocDeltaN.hpp

clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o	
	$(RM) $(MODEL)
