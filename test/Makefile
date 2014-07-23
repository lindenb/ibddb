.PHONY:all test clean ibdexecutable

all: test

test: test.h5
	LD_LIBRARY_PATH=/commun/data/packages/hdf5/lib ../src/ibddb dict $<
	LD_LIBRARY_PATH=/commun/data/packages/hdf5/lib ../src/ibddb markers $< | head
	LD_LIBRARY_PATH=/commun/data/packages/hdf5/lib ../src/ibddb markers -r "22:17202602-17221494" $<

test.h5: ibdexecutable reference.dict markers.bed ibd.txt LAP.fam
	LD_LIBRARY_PATH=/commun/data/packages/hdf5/lib ../src/ibddb build \
		-o $@ \
		--dict reference.dict \
		--bed markers.bed \
		--ped LAP.fam \
		--ibd ibd.txt

ibdexecutable:
	(cd ../src && make ibddb )

clean:
	rm -f test.h5
