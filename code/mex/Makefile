# File names
OUTFILE=../mexGMP
SRC=mexGMP.cpp PriorityVector.cpp
OBJ=$(SRC:.cpp=.o)

# Commands
CXX=mex
CXX_FLAGS=-O
LD_FLAGS=-cxx


## COMPILATION

all: $(OBJ)
	$(CXX) $(LD_FLAGS) -output $(OUTFILE) $(OBJ)

%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -c $<

%.o: %.hpp Makefile
mexGMP.o: PriorityVector.hpp


# PHONY RULES

.PHONY: configure clean mrproper color

configure:
	$(CXX) -setup

clean:
	rm -f $(OBJ)

mrproper: clean
	rm -f *~
