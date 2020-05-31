CC = icpc
CFLAGS = -Wall -O3 --std=c++11

all: main_annealing main_sampling

main_annealing : main_annealing.o random.o functions.o annealer.o metropolis.o
	$(CC) random.o main_annealing.o functions.o annealer.o metropolis.o -o main_annealing
main_sampling : main_sampling.o random.o functions.o annealer.o metropolis.o histogram.o
	$(CC) random.o main_sampling.o functions.o annealer.o metropolis.o histogram.o -o main_sampling
main_annealing.o : main_annealing.cpp
	$(CC) -c main_annealing.cpp -o main_annealing.o $(CFLAGS)
main_sampling.o : main_sampling.cpp
	$(CC) -c main_sampling.cpp -o main_sampling.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
annealer.o : annealer.cpp annealer.h metropolis.o
	$(CC) -c annealer.cpp -o annealer.o $(CFLAGS)
metropolis.o : metropolis.cpp metropolis.h
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)
histogram.o : histogram.cpp histogram.h
	$(CC) -c histogram.cpp -o histogram.o $(CFLAGS)

clean :
	rm *.o main_annealing main_sampling seed.out
