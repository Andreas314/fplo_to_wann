# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++23 -O2 -Wall -Wextra

# Files
SRC = main.cpp fplo_to_wann.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = fplo_to_wann

# Default target
all: $(TARGET)

# Linking
$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $(TARGET)

# Compilation
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean
