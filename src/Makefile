.PHONY : all clean test
export LD_LIBRARY_PATH

ifeq (${HDF5HOME},)
$(warning variable HDF5HOME is not defined. On kaamelot, could be HDF5HOME=/commun/data/packages/hdf5 )
endif

ifeq (${R_HOME},)
$(warning variable R_HOME is not defined. On kaamelot, could be R_HOME=/commun/data/packages/R/R-3.1.1/R-3.1.1 )
endif


CC=gcc
CFLAGS= -fPIC $(if ${HDF5HOME},-I${HDF5HOME}/include) -g -Wall `pkg-config --cflags cairo` $(if ${R_HOME},-I${R_HOME}/include )
LIBS=$(if ${HDF5HOME},-L${HDF5HOME}/lib) -L../lib -lz -lm `pkg-config  --libs cairo` -lhdf5
#-I${HDF5HOME}/include `pkg-config --cflags cairo`  

all:../bin/ibddb ../lib/libibddb.so
	R --vanilla  < test.R

../bin/ibddb : ibddb_main.o utils.o ibddb.o hershey.o
	mkdir -p $(dir $@) && \
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)  $(LIBS)

../lib/libibddb.so : ibddbR.o utils.o ibddb.o hershey.o
	mkdir -p $(dir $@) && \
	$(CC) -fPIC -shared -Wl,-soname,$(basename $@) -o $@ $^ ${LIBS}

ibddbR.o :  ibddbR.c ibddb.h utils.h githash.h
	$(CC) $(CFLAGS) -c -o $@  $<

ibddb_main.o : ibddb_main.c ibddb.h  utils.h githash.h
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

ifeq ($(realpath ../test),)

test: 

else

test: ../bin/ibddb $(realpath ../test/)
	ls /commun/data/users/simonet/fsimonet/ANALYSES/BRUGADA/famille_belfast/analyse_IBD/belfast_prIBD_*.ibdtxt.gz | grep -E '(Chr22|Chr21)' > ../test/test.list
	rm -f ../test/test.h5
	../bin/ibddb build \
		--out ../test/test.h5 \
		--dict ../test/reference.dict \
		--bed /commun/data/users/simonet/fsimonet/ANALYSES/BRUGADA/famille_belfast/analyse_IBD/belfast_prIBD.bed \
		--pedigree /commun/data/users/simonet/fsimonet/ANALYSES/BRUGADA/famille_belfast/belfast_prIBD.fam \
		--ibd  ../test/test.list

endif



clean:
	rm -f *.o ../bin/ibddb githash.h
