#    Continuous Wavelet Transform Library
#    Copyright (C) 2004-2009 Stepan V.Karpenko <carp@mail.ru>
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the
#    Free Software Foundation, Inc., 59 Temple Place, Suite 330,
#    Boston, MA  02111-1307  USA

all: cwt

# Compile demo program binary file
cwt: cwt.o cwtwlets.o main.o
	gcc -s -lm -o cwt cwt.o cwtwlets.o main.o

cwt.o: cwt.c cwt.h cwtwlets.h
	gcc -s -O2 -c cwt.c

cwtwlets.o: cwtwlets.c cwtwlets.h
	gcc -s -O2 -c cwtwlets.c

main.o: main.c cwt.h cwtwlets.h
	gcc -s -O2 -c main.c

# Run self-tests
test: cwt testcase/sine.txt testcase/delta.txt testcase/freqmod.txt
	./cwt testcase/sine.txt    wt_sine.txt
	./cwt testcase/delta.txt   wt_delta.txt
	./cwt testcase/freqmod.txt wt_freqmod.txt

# Remove results of compilation and self-tests execution
clean:
	rm *.o
	rm cwt
	rm wt_*.txt
