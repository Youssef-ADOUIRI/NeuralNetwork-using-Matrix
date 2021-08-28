CC = g++
CCFLAGS = -g -Wall 
# g++ main.cpp -lX11 -lpthread

all: prog

prog: main.o matrix.o Nn.o
	$(CC) main.o matrix.o Nn.o -Wall -o prog

%.o: %.cpp
	$(CC) -c $(CCFLAGS) $<

clean:
	rm -f core *.o prog