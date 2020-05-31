CC = icpc
CFLAGS = -Wall -O3 --std=c++11

main : main.o random.o functions.o MonteCarloNVT.o
	$(CC) random.o main.o functions.o MonteCarloNVT.o -o main
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
MonteCarloNVT.o : MonteCarloNVT.cpp MonteCarloNVT.h
	$(CC) -c MonteCarloNVT.cpp -o MonteCarloNVT.o $(CFLAGS)

clean :
	rm *.o main seed.out
