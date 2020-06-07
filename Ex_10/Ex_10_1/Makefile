CC = g++
CFLAGS = -Wall -O3 --std=c++11

main : main.o random.o functions.o annealer.o
	$(CC) main.o random.o functions.o annealer.o -o main
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
annealer.o : annealer.cpp annealer.h genetics.h
	$(CC) -c annealer.cpp -o annealer.o $(CFLAGS)

clean :
	rm *.o main seed.out
