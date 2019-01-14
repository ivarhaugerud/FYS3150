CXXFLAGS = -Wall -O3 -larmadillo
program: 2_d.o functions.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm program *.o
