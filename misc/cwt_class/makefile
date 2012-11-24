all: wtclass.exe

wtclass.exe: WTClass.o main.o cwtlib/cwt.o cwtlib/cwtwlets.o
	gpp  -s -lm -o wtclass WTClass.o main.o cwtlib/cwt.o cwtlib/cwtwlets.o

WTClass.o: WTClass.cpp WTClass.h cwtlib/cwt.h cwtlib/cwtwlets.h
	gpp -O2 -c WTClass.cpp

main.o: main.cpp WTClass.h
	gpp -O2 -c main.cpp

cwtlib/cwt.o: cwtlib/cwt.c cwtlib/cwt.h cwtlib/cwtwlets.h
	gpp -O2 -c cwtlib/cwt.c -o cwtlib/cwt.o

cwtlib/cwtwlets.o: cwtlib/cwtwlets.c cwtlib/cwtwlets.h
	gpp -O2 -c cwtlib/cwtwlets.c -o cwtlib/cwtwlets.o
