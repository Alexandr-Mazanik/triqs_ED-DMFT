CXX = g++
CXXFLAGS = -Ofast -Wall -shared -std=c++11 -fPIC $(shell python3 -m pybind11 --includes) -I$(CONDA_PREFIX)/include -I$(CONDA_PREFIX)/include/eigen3

SOURCES = ./src/solver.cpp ./src/solverLanczos.cpp ./src/helpers.cpp ./src/bindings.cpp
TARGET = ./ed_solver/_ed_solver$(shell python3-config --extension-suffix)

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

clean:
	rm -f ./ed_solver/_ed_solver*.so

.PHONY: all clean
