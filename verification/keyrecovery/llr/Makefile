#------------------------------------------------------
# Makefile for C/C++ Program
#------------------------------------------------------
# Target: main
# Author: Y. Todo
#------------------------------------------------------

TARGET=main

#------------------------------------------------------
# Default Parameters
#------------------------------------------------------

CC=g++
OPT=-m64 -O2 -funroll-loops -fopenmp -std=c++17
#INC=-I/opt/gurobi701/linux64/include/
#LIB=-lcryptominisat5 -lntl -lgmp

#------------------------------------------------------
# Compile Option
#------------------------------------------------------

-include makefile.opt

#------------------------------------------------------
# Definition
#------------------------------------------------------

.SUFFIXES:.cpp .c .o .h

#---
# Source Files
#---

SRC=$(shell ls *.cpp)
HED=$(shell ls *.h)
OBJ=$(SRC:.cpp=.o)

#------------------------------------------------------
# rules
#------------------------------------------------------

all: $(TARGET)
$(TARGET): $(OBJ)
	$(CC) $(OPT) -w -o $(TARGET) $(OBJ) $(INC) $(LIB) 

.c.o:
	$(CC) $(OPT) -c $< $(INC)
.cpp.o:
	$(CC) $(OPT) -c $< $(INC)

dep:
	g++ -MM -MG $(SRC) >makefile.depend

clean:
	rm -f $(TARGET) $(TARGET).exe
	rm -f *.o *.obj
	rm -f *~ *.~*

tar:
	tar cvzf $(TARGET).tar.gz $(SRC) $(HED) makefile

#--------------------------------------------------
-include makefile.depend
