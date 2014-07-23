.PHONY:all test clean ibdexecutable

all: test

test: test.h5 manhattan.R
	../src/ibddb dict $< 
	../src/ibddb markers $< | head 
	../src/ibddb markers -r "22:17202602-17221494" $<
	../src/ibddb ped $<
	../src/ibddb pairs $<
	../src/ibddb ibd $<
	../src/ibddb ibd --image jeter.png $<
	
test.h5: ibdexecutable reference.dict markers.bed ibd.txt LAP.fam
	../src/ibddb build \
		-o $@ \
		--dict reference.dict \
		--bed markers.bed \
		--ped LAP.fam \
		--ibd ibd.txt



ibdexecutable:
	(cd ../src && make ibddb )

clean:
	rm -f test.h5
