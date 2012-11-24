all: cwt

cwt: cwt.o cwtwlets.o main.o
	gcc -s -lm -o cwt cwt.o cwtwlets.o main.o

cwt.o: cwt.c cwt.h cwtwlets.h
	gcc -s -O2 -c cwt.c

cwtwlets.o: cwtwlets.c cwtwlets.h
	gcc -s -O2 -c cwtwlets.c

main.o: main.c cwt.h cwtwlets.h
	gcc -s -O2 -c main.c
