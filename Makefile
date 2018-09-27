CC = g++
VERSION = -std=c++11
CFLAGS = $(VERSION) -c -Wall

all:	pore

hello:	Source.o
	$(CC) Source.o -o pore

Source.o:	Atom.h Contcar.h Source.cpp
		$(CC) $(CFLAGS) Atom.h Contcar.h Source.cpp

clean:
	rm -rf *o pore
