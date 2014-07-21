.PHONY : all clean test
CC=gcc
CFLAGS= -g -Wall -I/commun/data/packages/hdf5/include
LDFLAGS= -L/commun/data/packages/hdf5/lib
LIBS=-lz -lhdf5

all: test ibddb

ibddb: utils.o ibddb.o
	$(CC) -o $@ $^ $(LDFLAGS)  $(LIBS)	

ibddb.o : ibddb.c ibddb.h  utils.h
	$(CC) -c -o $@ $(CFLAGS) $<

utils.o : utils.c utils.h
	$(CC) -c -o $@ $(CFLAGS) $<

test: ibddb
	LD_LIBRARY_PATH=/commun/data/packages/hdf5/lib ./ibddb build \
		-o jeter.h5 \
		--dict /commun/data/pubdb/broadinstitute.org/bundle/1.5/b37/human_g1k_v37.fasta.fai \
		--bed markers.bed.gz \
		--ped LAP.fam \
		--ibd ibd.tx

clean:
	rm -f *.o ibddb
