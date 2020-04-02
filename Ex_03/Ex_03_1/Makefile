CC = g++
CFLAGS = -Wall -O3 --std=c++11

main : main.o random.o functions.o blackscholes.o
	$(CC) random.o main.o functions.o blackscholes.o -o main
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
blackscholes.o : blackscholes.cpp blackscholes.h
	$(CC) -c blackscholes.cpp -o blackscholes.o $(CFLAGS)

clean :
	rm *.o main seed.out
