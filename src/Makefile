CC = g++
CCFLAGS = -g -Wall

all: prog

prog: main.o matrix.o Nn.o
	$(CC) main.o matrix.o Nn.o -Wall -o prog

%.o: %.cpp
	$(CC) -c $(CCFLAGS) $<

clean:
	rm -f core *.o prog
