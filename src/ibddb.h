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
#include <time.h>
#include <hdf5.h>
#include <cairo.h>
#include "utils.h"
#include "githash.h"


/** just a flag to map the type as transient (not saved in HDF5 file )*/
#define TRANSIENT(type) type

typedef struct chrom_t
	{
	char* name;
	int tid;
	int length;

	TRANSIENT(long) cumulative_start;
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

#define RESKIN_COLUMN_COUNT 9
#define RESKIN_COLUMN_IBD0  ( RESKIN_COLUMN_COUNT -1 )
typedef struct reskin_id
	{
	int pair_id;
	float data[RESKIN_COLUMN_COUNT];
	TRANSIENT(boolean_t) selected;
	} Reskin,*ReskinPtr;


typedef struct context_t
	{
	int argc;
	char** argv;
	char* hdf5_filename;
	char* faidx_filename;
	char* bed_filename;
	char* ped_filename;
	char* ibd_filename;
	char* reskin_filename;
	hid_t       file_id;   /* HDF5 file identifier */

	/** start time */
	time_t startup;

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
	/** reskins **/
	ReskinPtr reskins;
	size_t reskin_count;

	/* default output file */
	FILE* out;

	/** reading flag */
	boolean_t on_read_load_dict;
	boolean_t on_read_load_markers;
	boolean_t on_read_load_pedigree;
	boolean_t on_read_load_pairs;
	boolean_t on_read_load_reskins;
	} Context,*ContextPtr;

/* create a new context from argc/argv */
ContextPtr ContextNew(int,char** );
/* dispose context */
void ContextFree(ContextPtr);
/** open IBD context for reading after config->hdf5_filename and config->on_load_* be assigned */
void ContextOpenForRead(ContextPtr config);

/** placeholder to open/close the '/IBD' dataset, used by standalone and R extension */
typedef struct ibd_dataset_t
	{
	hid_t dataset_id;
	hid_t dataspace_id;
	hid_t  memspace;
	} IbdDataSet,*IbdDataSetPtr;
	
	
/** open the IBD dataset for reading */
IbdDataSetPtr IbdDataSetOpen(ContextPtr config);
/** close the IBD dataset after reading */
void IbdDataSetClose(IbdDataSetPtr ds);

/** graphics Stuff */
typedef struct rectangle2d_t
	{
	double x;
	double y;
	double width;
	double height;
	} Rectangle2D,*Rectangle2DPtr;

typedef struct rectangle_t
	{
	int x;
	int y;
	int width;
	int height;
	} Rectangle,*RectanglePtr;

typedef struct dimension_t
	{
	int width;
	int height;
	} Dimension,*DimensionPtr;

typedef struct dimension2d_t
	{
	double width;
	double height;
	} Dimension2D,*Dimension2DPtr;


typedef struct expdata_t
	{
	size_t marker_index;
	double value;
	} ExpData,*ExpDataPtr;


#define USAGE_PREAMBLE fprintf(stderr,"\n\n%s\nAuthor: Pierre Lindenbaum PhD\nGit-Hash: "GIT_HASH"\nWWW: https://github.com/lindenb/ibddb\nCompilation: %s at %s\n\n",argv[0],__DATE__,__TIME__)

#endif
