
CC = g++
CCFLAGS = -g -Wall

prog: main.o
	$(CC) main.o matrix.o -Wall -o prog

%.o: %.cpp
	$(CC) -c $(CCFLAGS) $<

clean:
	rm -f core *.o prog
