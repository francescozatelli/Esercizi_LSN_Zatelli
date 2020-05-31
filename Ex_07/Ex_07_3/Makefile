CC = icpc
CFLAGS = -Wall -O3 --std=c++11

main : main.o molecular_dynamics.o histogram.o
	$(CC) molecular_dynamics.o main.o histogram.o -o main
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
molecular_dynamics.o : molecular_dynamics.cpp molecular_dynamics.h histogram.o
	$(CC) -c molecular_dynamics.cpp -o molecular_dynamics.o $(CFLAGS)
histogram.o : histogram.cpp histogram.h
	$(CC) -c histogram.cpp -o histogram.o $(CFLAGS)
reset :
	./clean.sh
clean :
	rm *.o main seed.out