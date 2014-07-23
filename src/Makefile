.PHONY : all clean test
HDF5HOME=/commun/data/packages/hdf5
CC=${HDF5HOME}/bin/h5cc
CFLAGS= -noshlib -prefix=${HDF5HOME} -g -Wall `pkg-config --cflags cairo`
LIBS=-lz `pkg-config  --libs cairo`

all: ibddb

ibddb: utils.o ibddb.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)  $(LIBS)	

ibddb.o : ibddb.c ibddb.h  utils.h
	$(CC) $(CFLAGS) -c -o $@  $<

utils.o : utils.c utils.h
	$(CC) -c -o $@ $(CFLAGS) $<

clean:
	rm -f *.o ibddb
