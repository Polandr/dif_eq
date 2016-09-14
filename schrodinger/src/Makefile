CXX= g++
HDRS= solver.h solver.hpp parameters.hpp
OPTIONS= -std=c++11

first: test

test: test.cpp $(HDRS)
	$(CXX) test.cpp -o test $(OPTIONS)

clean:
	rm -f test

.PHONY: first clean
