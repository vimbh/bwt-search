# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra

# Source file and output executable
SRC := bwtsearch.cpp
OUT := bwtsearch

all: $(OUT)

$(OUT): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(OUT)

clean:
	rm -f $(OUT)