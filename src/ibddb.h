/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum PhD.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/
#ifndef IBDDB_H
#define IBDDB_H

#include <hdf5.h>
#include "utils.h"



/** just a flag to map the type as transient (not saved in HDF5 file )*/
#define TRANSIENT(type) type

typedef struct chrom_t
	{
	char* name;
	int tid;
	int length;
	} Chrom,*ChromPtr;


typedef struct marker_t
	{
	char* name;
	int tid;
	int position;
	int index;
	//TRANSIENT(boolean_t) selected;
	} Marker,*MarkerPtr;


typedef struct individual_t
	{
	char* family;
	char* name;
	char* father;
	char* mother;
	int sex;
	int status;
	int index;
	TRANSIENT(boolean_t) selected;
	} Individual,*IndividualPtr;

typedef struct pair_id
	{
	int indi1idx;
	int indi2idx;
	int index;
	TRANSIENT(boolean_t) selected;
	} PairIndi,*PairIndiPtr;

typedef struct region_id
	{
	int tid;
	int start;
	int end;
	} Region,*RegionPtr;

typedef struct context_t
	{
	int argc;
	char** argv;
	char* hdf5_filename;
	char* faidx_filename;
	char* bed_filename;
	char* ped_filename;
	char* ibd_filename;
	hid_t       file_id;   /* file identifier */


	/** chromosomes */
	ChromPtr chromosomes;
	size_t chromosome_count;
	/** markers */
	MarkerPtr markers;
	size_t marker_count;
	/** individuals */
	IndividualPtr individuals;
	size_t individual_count;
	/** pairs **/
	PairIndiPtr pairs;
	size_t pair_count;

	/* default output file */
	FILE* out;

	/** reading flag */
	boolean_t on_read_load_dict;
	boolean_t on_read_load_markers;
	} Context,*ContextPtr;





#endif