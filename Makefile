CXX=g++
INC=-I ./include
CFLAGS=-std=c++11 -W -Wall -pedantic $(INC)

aco_ta: main.cpp ./src/*
	$(CXX) $(CFLAGS) -o $@ $^

clean:
	rm -f aco_ta

.PHONY: clean
