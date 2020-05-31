CC = g++
CFLAGS = -Wall -O3 --std=c++11

main : main.o random.o functions.o ising1D.o
	$(CC) random.o main.o functions.o ising1D.o -o main
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
ising1D.o : ising1D.cpp ising1D.h
	$(CC) -c ising1D.cpp -o ising1D.o $(CFLAGS)

clean :
	rm *.o main seed.out
