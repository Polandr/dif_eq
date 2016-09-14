CXX= g++
HDRS= solver.h solver.hpp parameters.hpp
OPTIONS= -std=c++0x

first: test

test: test.cpp $(HDRS)
	$(CXX) test.cpp -o test $(OPTIONS)

clean:
	rm -f test

.PHONY: first clean
