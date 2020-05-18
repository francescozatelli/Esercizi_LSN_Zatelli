CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: main_unif main_gauss

unif: main_unif

gauss: main_gauss

main_unif : main_unif.o random.o functions.o metropolis.o
	$(CC) random.o main_unif.o functions.o metropolis.o -o main_unif
main_gauss : main_gauss.o random.o functions.o metropolis.o
	$(CC) random.o main_gauss.o functions.o metropolis.o -o main_gauss
main_unif.o : main_unif.cpp
	$(CC) -c main_unif.cpp -o main_unif.o $(CFLAGS)
main_gauss.o : main_gauss.cpp
	$(CC) -c main_gauss.cpp -o main_gauss.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
metropolis.o : metropolis.cpp metropolis.h
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)

clean :
	rm *.o main seed.out *.points out_*
