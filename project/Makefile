# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -O2
LDFLAGS = -lcpgplot -lpgplot -lX11 -lm

# Targets
TARGETS = demo final

# Build targets
all: $(TARGETS)

demo: demo.cpp
	$(CXX) $(CXXFLAGS) demo.cpp -o demo $(LDFLAGS)

final: final.cpp
	$(CXX) $(CXXFLAGS) final.cpp -o final $(LDFLAGS)

# Clean target
clean:
	rm -f $(TARGETS) *.o

.PHONY: all clean

