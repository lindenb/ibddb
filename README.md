IBDDB
#####

# Motivation

Storing & Querying IBD status in a **HDF5** database.

> HDF5 is a data model, library, and file format for storing and managing data.
> It supports an unlimited variety of datatypes, and is designed for flexible and efficient I/O and for high volume and complex data


# Compilation

## Dependencies

* hdf5 library: http://www.hdfgroup.org/HDF5/
* cairo http://cairographics.org/
* GCC compiler and GNU-Make

## Compile

```bash
$ make
```
## error "while loading shared libraries"

If  you get this error message:

```
ibddb: error while loading shared libraries: libhdf5.so.8: cannot open shared object file: No such file or directory
```

you'll have to define the LD_LIBRARY_PATH

```
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/distribution-of/hdf5/lib
```

for my collaborators at @institut_thorax (*.161)

```
export LD_LIBRARY_PATH=/commun/data/packages/hdf5/lib
```

## `build` Building the database

build the HDF5 database. 

### Options


* -o file.h5  : the HDF5 database to write
* --dict (file.dict) : the Reference dictionary. A tab delimited file with at least two columns. First column is the name of the chromosome, the second is the length of the chromosome in 'bp'.
```
1	249250621	52	60	61
2	243199373	253404903	60	61
3	198022430	500657651	60	61
(...)
```
* --bed (marker.bed) : a tab delimited  file containing the position and the name of the markers. Columns: chromosome, start,end,marker-name. the chromosome must be defined in the Reference ictionary.
```
1	754183	754181	rs2121969
1	798959	798958	rs11340777
1	1026959	1026958	rs11579015
1	1040036	1040035	rs6671256
1	1759312	1759313	rs978694
(...)
```
* --pedigree (file.ped) A tab delimited containing a list of individuals. Columns are Family/Name/Father(or '0')/Mother( or '0')/Sex/Status.
```
LDC PA33395 1 10 1 1
LDC PA33141 0 0 2 2
LDC PA33175 L0774 PA09080 2 -9
LDC PA33178 L0775 PA33395 1 -9
LDC PA33303 PA33314 1430 1 -9
LDC PA33304 PA33314 1430 2 -9
LDC PA33334 1410 PA33315 2 -9
LDC PA33341 1 10 1 -9
(...)
```


* --reskin (path.reskin) An optional Reskin files. 

```
LDC	13	13	0	0	0	0	0	0	1	0	0
LDC	13	18	0	0	0	0	0	0	0	0	1
LDC	13	28	0	0	0	0	0	0	0	0	1
LDC	13	34	0	0	0	0	0	0	0	0	1
LDC	13	36	0	0	0	0	0	0	0	0	1
```


* --ibd (path.txt) A text file containing the **file path** to the IBD files.

```
analysis_IBD/BLP10.ibd.gz
analysis_IBD/BLP11.ibd.gz
analysis_IBD/BLP12.ibd.gz
analysis_IBD/BLP13.ibd.gz
analysis_IBD/BLP14.ibd.gz
analysis_IBD/BLP15.ibd.gz
analysis_IBD/BLP16.ibd.gz
analysis_IBD/BLP17.ibd.gz
analysis_IBD/BLP18.ibd.gz
analysis_IBD/BLP19.ibd.gz
```

An IBD file is a tab delimited file.

The first row define the chromosome and the markers:

* number of ibd status ? must be 3
* number of pairs 
* chromosomes
* number of markers
* marker 1
* marker 2
* marker 3
* marker 4
* ....

e.g: 

```
$ gunzip -c file.ibdtxt.gz  | head -n 1 | tr "\t" "\n" | head
3
2556
21
3506
rs574999
rs296537
rs284579
rs578623
rs960546
rs78045
(...)
```

The remaining rows have the following columns:

* Indi1-family
* Indi1-name
* Indi2-family
* Indi2-name
* marker1 IBD0 (float)
* marker1 IBD1
* marker1 IBD2
* marker2 IBD0
* marker2 IBD1
* marker2 IBD2
* marker3 IBD0
* (...)

e.g: 

```
$ gunzip -c ibdtxt.gz  | grep -v rs | head -n 1 | tr "\t" "\n" | head
FAM1
1
FAM1
1
0
0
1
0
0
1
```



All files can be g-zipped.


### Example



```
ibddb build \
	-o out.h5 \
	--dict reference.dict \
	--bed markers.bed.gz \
	--ped  pedigree.fam  \
	--ibd  ibd.txt
```

The content of the h5 file can be checked using `h5dump`: http://www.hdfgroup.org/HDF5/doc/RM/Tools.html#Tools-Dump

```
h5dump -H test.h5  | more 
HDF5 "test.h5" {
GROUP "/" {
   DATASET "dictionary" {
      DATATYPE  H5T_COMPOUND {
         H5T_STRING {
            STRSIZE H5T_VARIABLE;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
         } "name";
         H5T_STD_I32LE "tid";
         H5T_STD_I32LE "length";
      }
(....)
         H5T_STRING {
            STRSIZE H5T_VARIABLE;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
         } "mother";
         H5T_STD_I32LE "sex";
         H5T_STD_I32LE "status";
         H5T_STD_I32LE "index";
      }
      DATASPACE  SIMPLE { ( 50 ) / ( 50 ) }
   }
}
}
```




## `ibd` querying the **h5** database.

### Options:

* -r|--region (chr|chr:start-end) restrict to that region. Optional.
* -i|--individual (fam:name) restrict to that individual. Can be used multiple times.
* -p|--pair (fam1:name1|fam2:name2) restrict to that pair. Can be used multiple times.
* -F|--family (fam) restrict to that family. Can be used multiple times.
* --noselfself ignore all self-self pairs.
* --treshold (float) IBD TRESHOLD default:0.100000 

new in 2016:

* -I|--individualfile (file): tab delimited file containing fam\tname\n to restrict to those individuals.n -p|--pair (fam1:name1|fam2:name2) restrict to that pair. Can be used multiple times.
* -P|--pairfile (file):   tab delimited file containing : fam1\tname1\tfam2\tname2\n to restrict to those pairs.
* -Y|--familyfile (file) read file to restrict to thoses families. Can be used multiple times.
* --reskin 'min/max' if defined and reskin data available, will restrict to pair having reskin ibd0 in this range .



Tabular options

* --noheader don't print data header
* --nopairsinheader don't print pairs in data header


Image Options:

* -g|--image (filename.png) save as PNG picture.
* --width (int) image-width.
* --height (int) image-height.





### Examples:

```
$ ibddb ibd -r 22:50334314-50577522 test.h5 
CHROM POS NAME LDC:PA00313|LDC:PA00313 LDC:L0776|LDC:L0778 LDC:L0777|LDC:L0777 LDC:L0777|LDC:L0778 LDC:L0778|LDC:L0778 COUNT_IBD
22 50523117 rs11568171 0.000000 1.000000 0.000000 0.000000 0.000000 188
22 50530651 rs738190 0.000000 1.000000 0.000000 0.000000 0.000000 188
22 50556331 rs138223 0.000000 1.000000 0.000000 0.000000 0.000000 187
22 50557107 rs138229 0.000000 1.000000 0.000000 0.000000 0.000000 187
22 50559953 rs138233 0.000000 1.000000 0.000000 0.000000 0.000000 186
22 50565311 rs6010138 0.000000 1.000000 0.000000 0.000000 0.000000 186
22 50570852 rs138251 0.000000 1.000000 0.000000 0.000000 0.000000 186
22 50577521 rs1838818 0.000000 1.000000 0.000000 0.000000 0.000000 186
```

```
$ ibddb ibd -r 22:50334314-50577522   --nopairsinheader test.h5
CHROM	POS	NAME	COUNT_IBD
22	50523117	rs11568171	188
22	50530651	rs738190	188
22	50556331	rs138223	187
22	50557107	rs138229	187
22	50559953	rs138233	186
22	50565311	rs6010138	186
22	50570852	rs138251	186
22	50577521	rs1838818	186
```

```
$ ibddb ibd -r 22:50334314-50577522 --nopairsinheader  --treshold 0.5  test.h5 
CHROM	POS	NAME	COUNT_IBD
22	50523117	rs11568171	210
22	50530651	rs738190	210
22	50556331	rs138223	212
22	50557107	rs138229	212
22	50559953	rs138233	212
22	50565311	rs6010138	212
22	50570852	rs138251	212
22	50577521	rs1838818	212
```

```
$ ibddb ibd -r 22:50334314-50577522  --individual LDC:PA00313 --noselfself  test.h5   | head -n 1 | tr "\t" "\n"
CHROM
POS
NAME
LDC:PA00313|LDC:PA00319
LDC:PA00315|LDC:PA00319
LDC:PA00318|LDC:PA00319
LDC:PA00319|LDC:PA00320
LDC:PA00319|LDC:PA00322
LDC:PA00319|LDC:PA00353
LDC:PA00319|LDC:PA00351
LDC:PA00319|LDC:PA00355
LDC:PA00319|LDC:PA00357
LDC:PA00319|LDC:PA00358
LDC:PA00319|LDC:PA00359
LDC:PA00319|LDC:PA00361
LDC:PA00319|LDC:PA00362
LDC:PA00319|LDC:PA00363
LDC:PA00319|LDC:PA00365
LDC:PA00319|LDC:PA00510
LDC:PA00319|LDC:PA06869
LDC:PA00319|LDC:PA09080
LDC:PA00319|LDC:PA11012
LDC:PA00319|LDC:PA11013
LDC:PA00319|LDC:PA11122
LDC:PA00319|LDC:PA11121
LDC:PA00319|LDC:PA11125
LDC:PA00319|LDC:PA11195
LDC:PA00319|LDC:PA11212
LDC:PA00319|LDC:PA11275
LDC:PA00319|LDC:PA11278
LDC:PA00319|LDC:PA11303
LDC:PA00319|LDC:PA11301
LDC:PA00319|LDC:PA11331
LDC:PA00319|LDC:PA11312
LDC:PA00319|LDC:L0316
LDC:PA00319|LDC:L0666
LDC:PA00319|LDC:L0773
LDC:PA00319|LDC:L0771
LDC:PA00319|LDC:L0775
LDC:PA00319|LDC:L0776
LDC:PA00319|LDC:L0777
LDC:PA00319|LDC:L0778
COUNT_IBD
```

```
$ ibddb ibd -r 22:50334314-50577522  --pair "LDC:PA00319|LDC:PA11303" --pair "LDC:L0777|LDC:PA00319" test.h5
CHROM	POS	NAME	LAP:CD00319|LAP:CD11303	LAP:CD00319|LAP:L0777	COUNT_IBD
22	50523117	rs11568171	0.979200	0.966100	0
22	50530651	rs738190	0.978900	0.965500	0
22	50556331	rs138223	0.977300	0.961100	0
22	50557107	rs138229	0.977300	0.961000	0
22	50559953	rs138233	0.977200	0.960500	0
22	50565311	rs6010138	0.976900	0.959800	0
22	50570852	rs138251	0.976500	0.958800	0
22	50577521	rs1838818	0.976000	0.957800	0
```


Generating an image

```
$ ibddb ibd -g out.png test.h5
```


![Screenshot](https://raw.githubusercontent.com/lindenb/ibddb/master/doc/screenshot01.png "Screenshot")



## `markers` dumping the markers

```
$ ibddb markers test.h5 
1	751182	751183	rs3131969
1	798959	798960	rs11210777
1	1036959	1036960	rs11579015
1	1010026	1010027	rs6671356
1	1759213	1759211	rs9786912
1	1776269	1776270	rs1618727
1	1809509	1809510	rs11260621
1	2010991	2010995	rs16821727
1	2292688	2292689	rs2810510
1	2301171	2301172	rs2810536
```

use option `-r` to get a subset or markers:

```
$ ibddb markers -r 22:50565314-50577522 test.h5 
22	50565314	50565315	rs6010138
22	50570852	50570853	rs138251
22	50577521	50577522	rs4838848
```


## `ped` dump the pedigree

```
$ ibddb ped test.h5

LDC	1	1001	1010	1	-9
LDC	10	0	0	2	1
LDC	100	0	0	1	-9
LDC	1001	0	0	1	-9
LDC	1010	0	0	2	-9
LDC	2	1001	1010	1	-9
LDC	20	0	0	2	1
LDC	2110	0	0	2	1
LDC	2120	0	0	1	1
LDC	260	0	0	1	1
```

## `pairs` dump the pairs of individuals

```
$ ibddb pairs test.h5

LDC	PA00313	1	10	2	-9	LDC	PA00313	1	10	2	-9
LDC	PA00313	1	10	2	-9	LDC	PA00315	1	10	2	-9
LDC	PA00313	1	10	2	-9	LDC	PA00318	L0775	PA11195	2	-9
LDC	PA00313	1	10	2	-9	LDC	PA00319	260	PA11312	1	-9
LDC	PA00313	1	10	2	-9	LDC	PA00320	260	PA11312	1	-9
LDC	PA00313	1	10	2	-9	LDC	PA00322	260	PA11312	2	-9
LDC	PA00313	1	10	2	-9	LDC	PA00353	0	0	2	1
LDC	PA00313	1	10	2	-9	LDC	PA00351	1	10	1	-9
LDC	PA00313	1	10	2	-9	LDC	PA00355	PA00351	L0773	1	-9
LDC	PA00313	1	10	2	-9	LDC	PA00357	PA06869	PA00353	1	-9
(...)
```


## `dict` dump the reference dictionary

```
$ ibddb dict test.h5
1	249250621
2	243199373
3	198022430
4	191154276
5	180915260
6	171115067
7	159138663
8	146364022
9	141213431
10	135534747
(...)
```
# R Interface

```
> Opens an IBD-DB HDF5 file
>
> @param filename the HDF5 file
> @keywords HDF5 file
> @return ibd context
> @examples
> ibd.open('file.h5')
ibd.open<-function(filename)
```

```
> Close an ibd context and realease the associated resources
>
> @param ibd the IBD context
> @keywords ibd
> @examples
> ibd.close(ibd)
ibd.close<-function(ibd)
```

```
> returns the number of markers in the IBD context
>
> @param ibd the IBD context
> @keywords ibd
> @return the number of markers
ibd.num.markers<-function(ibd)
```

```
> returns the number of chromosomes in the IBD context
>
> @param ibd the IBD context
> @keywords ibd
> @return the number of chromosomes
ibd.num.chromosomes<-function(ibd)
```

```
> returns the number of pairs in the IBD context
>
> @param ibd the IBD context
> @keywords ibd
> @return the number of pairs
ibd.num.pairs<-function(ibd)
```

```
> returns the number of individuals in the IBD context
>
> @param ibd the IBD context
> @keywords ibd
> @return the number of individuals
ibd.num.individuals<-function(ibd)
```

```
> returns the index-th chromosome in the IBD context. Chromosome are ordered on the original sequence dictionary
> A chromosome is a tuple (name,tid,length)
> @param ibd the IBD context
> @param index 0-based index
> @keywords ibd
> @return the index-th chromosome	
ibd.chromosome<-function(ibd,index)
```

```
> returns the index-th marker in the IBD context. Markers are ordered on the tid,position
> A marker is a tuple (name,chromosome-index,position,self-index)
> @param ibd the IBD context
> @param index 0-based index
> @keywords ibd
> @return the index-th marker	
ibd.marker<-function(ibd,index)
```

```
> returns the index-th pair of individual in the IBD context
> A pair is a tuple (individual1-index,individual2-index,self-index)
> @param ibd the IBD context
> @param index 0-based index
> @keywords ibd
> @return the index-th pair	
ibd.pair<-function(ibd,index)
```

```
> returns the index-th pair of individual in the IBD context
> A pair is a tuple (individual1-index,individual2-index,self-index)
> @param ibd the IBD context
> @param index 0-based index
> @keywords ibd
> @return the index-th pair
ibd.individual<-function(ibd,index)
```

```
> returns the IBD-0 status for the given marker-index,pair-index
> @param ibd the IBD context
> @param marker-index 0-based index of the marker
> @param pair-index 0-based index of the pair
> @keywords ibd 
> @return the IBD-0 (or null if undefined)
ibd.ibd0<-function(ibd,marker_y,pair_x)
```

```
> returns the IBD-1 status for the given marker-index,pair-index
> @param ibd the IBD context
> @param marker-index 0-based index of the marker
> @param pair-index 0-based index of the pair
> @keywords ibd 
> @return the IBD-1 (or null if undefined)
ibd.ibd1<-function(ibd,marker_y,pair_x)
```

```
> returns the IBD-2 status for the given marker-index,pair-index
> @param ibd the IBD context
> @param marker-index 0-based index of the marker
> @param pair-index 0-based index of the pair
> @keywords ibd 
> @return the IBD-2 (or null if undefined)
ibd.ibd2<-function(ibd,marker_y,pair_x)
```

## Example

the following R program print the IBD0 for y=Marker, x=pairs.

```R
source("ibddb.R")

# OPEN IBD Database
ibd <- ibd.open("../test/test.h5");

###  loop over markers
for(y in seq(0,min(50,ibd.num.markers(ibd)-1)) )
	{
	## get y-th marker
	marker <- ibd.marker(ibd,y)
	## get chromosome info for this marker
	chrom <- ibd.chromosome(ibd,marker$tid)
	
	### loop over pairs
	for(x in seq(0,min(10,ibd.num.pairs(ibd)-1)) )
		{
		## get x-th pair
		pair <- ibd.pair(ibd,x)
		## get first individual of this pair
		indi1 <- ibd.individual(ibd,pair$indi1idx)
		## get second individual of this pair
		indi2 <- ibd.individual(ibd,pair$indi2idx)
		
		cat(paste(chrom$name,marker$name,marker$position,
			paste(indi1$family,indi1$name,sep=":"),
			paste(indi2$family,indi2$name,sep=":"),
			sep="\t"))
		cat("\t")
		# print IBD 0
		cat(ibd.ibd0(ibd,y,x))
		cat("\n")
		}
	}



# releases resources associated to IBD database
ibd.close(ibd);
```

# Contribute

- Issue Tracker: http://github.com/lindenb/ibddb/issues`
- Source Code: http://github.com/lindenb/ibddb



# History

* 2016-03: reskin
* 2015-03: new version of ibdld
* 2014 1st version

# Author

Pierre Lindenbaum PhD
@yokofakun


## License

The project is licensed under the MIT license.

