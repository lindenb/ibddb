.PHONY : all clean test
HDF5HOME=/commun/data/packages/hdf5
CC=${HDF5HOME}/bin/h5cc
CFLAGS= -noshlib -prefix=${HDF5HOME} -g -Wall `pkg-config --cflags cairo`
LIBS=-lz -lm `pkg-config  --libs cairo`


all: ../bin/ibddb

../bin/ibddb: utils.o ibddb.o hershey.o
	mkdir -p $(dir $@) && \
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)  $(LIBS)	

ibddb.o : ibddb.c ibddb.h  utils.h hershey.h githash.h
	$(CC) $(CFLAGS) -c -o $@  $<

utils.o : utils.c utils.h githash.h
	$(CC) -c -o $@ $(CFLAGS) $<

hershey.o: hershey.c hershey.h
	$(CC) -c -o $@ $(CFLAGS) $<	

githash.h:
	echo "#ifndef GIT_HASH" > $@
	echo -n '#define GIT_HASH "' >> $@
	-git rev-parse HEAD | tr -d "\n" >> $@
	echo '"' >> $@
	echo '#endif' >> $@


clean:
	rm -f *.o ../bin/ibddb githash.h
