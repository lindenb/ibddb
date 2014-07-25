.PHONY : all clean test
HDF5HOME=/commun/data/packages/hdf5
CC=gcc
CFLAGS= -fPIC -I${HDF5HOME}/include -g -Wall `pkg-config --cflags cairo` -I/commun/data/packages/R-2.15-3/lib64/R/include
LIBS=-L${HDF5HOME}/lib -L../lib -lz -lm `pkg-config  --libs cairo` -lhdf5


all: ../lib/libibddb.so ../bin/ibddb


../bin/ibddb: ibddb_main.o utils.o ibddb.o hershey.o
	mkdir -p $(dir $@) && \
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)  $(LIBS)	
../lib/libibddb.so: ibddbR.o utils.o ibddb.o hershey.o
	mkdir -p $(dir $@) && \
	$(CC) -fPIC -shared -Wl,-soname,$(basename $@) -o $@ $^
	$(CC) -fPIC -shared -Wl,-soname,ibddb.so -o ../lib/ibddb.so $^

ibddbR.o :  ibddbR.c ibddb.h utils.h
	$(CC) $(CFLAGS) -c -o $@  $<

ibddb_main.o : ibddb_main.c ibddb.h  utils.h
	$(CC) $(CFLAGS) -c -o $@  $<

ibddb.o : ibddb.c ibddb.h  utils.h hershey.h githash.h
	$(CC)  $(CFLAGS) -c -o $@  $<

utils.o : utils.c utils.h githash.h
	$(CC) $(CFLAGS) -c -o $@  $<

hershey.o: hershey.c hershey.h
	$(CC) $(CFLAGS) -c -o $@  $<	

githash.h:
	echo "#ifndef GIT_HASH" > $@
	echo -n '#define GIT_HASH "' >> $@
	-git rev-parse HEAD | tr -d "\n" >> $@
	echo '"' >> $@
	echo '#endif' >> $@


clean:
	rm -f *.o ../bin/ibddb githash.h
