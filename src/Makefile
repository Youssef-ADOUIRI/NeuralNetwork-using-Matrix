
CC = g++
CCFLAGS = -g -Wall

prog: main.o
	$(CC) -I Eigen main.o -Wall -o prog

%.o: %.cpp
	$(CC) -c $(CCFLAGS) $<

clean:
	rm -f core *.o prog