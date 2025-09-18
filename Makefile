# Contents of this file were written with the help of ChatGPT

# Compiler and flags
CXX := g++ -fopenmp
CXXFLAGS := -O3 -Wall -Wextra
LDFLAGS :=

# Source files and object files
SRCS := main.cpp EquationOfState.cpp Euler.cpp FluxSolver.cpp Mesh.cpp Reconstruction.cpp Solver.cpp STLReader.cpp
OBJS := $(SRCS:.cpp=.o)
TARGET := simple-cfd

# Default target
all: $(TARGET)

# Link the final executable
$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

# Generic compile rule (pattern rule)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean