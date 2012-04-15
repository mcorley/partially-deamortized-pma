# Makefile for testing the pma structure

SHELL = /bin/sh

CXX = g++
CXXPFLAGS = -g -O
LIBS = 

demo: pma_test.o pma.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	-rm -f demo *.o

