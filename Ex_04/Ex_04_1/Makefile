CC = icpc
CFLAGS = -Wall -O3 --std=c++11

main : main.o molecular_dynamics.o
	$(CC) molecular_dynamics.o main.o -o main
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
molecular_dynamics.o : molecular_dynamics.cpp molecular_dynamics.h
	$(CC) -c molecular_dynamics.cpp -o molecular_dynamics.o $(CFLAGS)
reset :
	./clean.sh && rm old.0 old.final
clean :
	rm *.o main
