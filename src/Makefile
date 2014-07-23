.PHONY : all clean test
CC=gcc
CFLAGS= -g -Wall -I/commun/data/packages/hdf5/include
LDFLAGS= -L/commun/data/packages/hdf5/lib
LIBS=-lz -lhdf5

all: ibddb

ibddb: utils.o ibddb.o
	$(CC) -o $@ $^ $(LDFLAGS)  $(LIBS)	

ibddb.o : ibddb.c ibddb.h  utils.h
	$(CC) -c -o $@ $(CFLAGS) $<

utils.o : utils.c utils.h
	$(CC) -c -o $@ $(CFLAGS) $<

clean:
	rm -f *.o ibddb
