CXXFLAGS = -Wall -larmadillo -std=c++11
program: ODESolver.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm program *.o
