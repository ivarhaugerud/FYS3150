CXXFLAGS = -Wall -O3 -larmadillo -std=c++11
program: main.o SolarSystemSimulator.o SunAtRest.o RelativisticCorrection.o ChangeExponent.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm program *.o
