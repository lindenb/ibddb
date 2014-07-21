#ifndef IBDDB_H
#define IBDDB_H

#include <hdf5.h>
#include "utils.h"


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
	} Individual,*IndividualPtr;

typedef struct pair_id
	{
	int indi1idx;
	int indi2idx;
	int index;
	} PairIndi,*PairIndiPtr;



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
	} Context,*ContextPtr;




#endif
