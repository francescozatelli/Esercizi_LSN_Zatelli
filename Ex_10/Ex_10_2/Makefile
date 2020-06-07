CC = g++
CFLAGS = -Wall -O3 --std=c++11
MPICOMP = mpicxx

main : main.o random.o functions.o
	$(MPICOMP) main.o random.o functions.o -o main
main.o : main.cpp genetics.h tsp.h
	$(MPICOMP) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)

clean :
	rm *.o main seed.out
