.PHONY:all test2 test clean ibdexecutable

all: test2

test2: ibdexecutable
	 cat /commun/data/users/simonet/Dina_temp/famille_LAP_pvm/LAP.bim | \
		awk -F ' ' '{printf("%s\t%d\t%d\t%s\n",$$1,$$4,int($$4)-1,$$2);}' |\
		grep -v '^26' |\
		gzip --best > markers.bed.gz
	cat /commun/data/pubdb/broadinstitute.org/bundle/1.5/b37/human_g1k_v37.fasta.fai  |\
		head -n 24 |\
		sed -e 's/^X/23/' -e 's/^Y/24/' > reference.dict
	grep -v LAP26 /commun/data/users/simonet/fsimonet/ANALYSES/PVM/fam_LAP/BDB_20140703/ibdpath.txt > ibd.txt
	../bin/ibddb build \
		-o test.h5 \
		--dict reference.dict \
		--bed markers.bed.gz \
		--ped  /commun/data/users/simonet/Dina_temp/famille_LAP_pvm/LAP.fam  \
		--ibd  ibd.txt
	


test: test.h5 manhattan.R
	../bin/ibddb dict $< 
	../bin/ibddb markers $< | head 
	../bin/ibddb markers -r "22:17202602-17221494" $<
	../bin/ibddb ped $<
	../bin/ibddb pairs $<
	../bin/ibddb ibd $<
	../bin/ibddb ibd --image jeter.png $<
	
test.h5: ibdexecutable reference.dict markers.bed ibd.txt LAP.fam
	../bin/ibddb build \
		-o $@ \
		--dict reference.dict \
		--bed markers.bed \
		--ped LAP.fam \
		--ibd ibd.txt



ibdexecutable:
	(cd ../src && make ibddb )

clean:
	rm -f test.h5
