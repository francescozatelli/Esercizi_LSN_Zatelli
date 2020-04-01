CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: main_a main_b main_c

a: main_a

b: main_b

c: main_c

main_a : main_a.o random.o functions.o
	$(CC) random.o main_a.o functions.o -o main_a
main_b : main_b.o random.o functions.o
	$(CC) random.o main_b.o functions.o -o main_b
main_c : main_c.o random.o functions.o
	$(CC) random.o main_c.o functions.o -o main_c
main_a.o : main_a.cpp
	$(CC) -c main_a.cpp -o main_a.o $(CFLAGS)
main_b.o : main_b.cpp
	$(CC) -c main_b.cpp -o main_b.o $(CFLAGS)
main_c.o : main_c.cpp
	$(CC) -c main_c.cpp -o main_c.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
clean :
	rm *.o main_a main_b main_c seed.out
