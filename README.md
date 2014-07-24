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

An IBD file is a tab delimited file containing the following columns: Indi1-family,Indi1-name,  Indi2-family,Indi2-name,chromosome,marker-name,IBD0,IBD1,IBD2:

```
LDC	PA00333	LDC	L0336	17	rs31450789	0.0009	0.0050	0.9943	
LDC	PA00333	LDC	L0336	17	rs8078919	0.0007	0.0044	0.9949	
LDC	PA00333	LDC	L0336	17	rs7107537	0.0003	0.0009	0.9993	
LDC	PA00333	LDC	L0336	17	rs9789059	0.0000	0.0007	0.9993	
LDC	PA00333	LDC	L0336	17	rs31453179	0.0000	0.0006	0.9994	
LDC	PA00333	LDC	L0336	17	rs8069970	0.0000	0.0003	0.9997	
LDC	PA00333	LDC	L0336	17	rs31945173	0.0000	0.0001	0.9998	
LDC	PA00333	LDC	L0336	17	rs6565714	0.0000	0.0003	0.9998	
LDC	PA00333	LDC	L0336	17	rs33873157	0.0000	0.0003	0.9999	
LDC	PA00333	LDC	L0336	17	rs9748036	0.0000	0.0003	0.9999	
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

TODO


### Example:

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


# Author

Pierre Lindenbaum PhD
@yokofakun
