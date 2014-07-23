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
#include <unistd.h>
#include <getopt.h>
#include <inttypes.h>
#include <assert.h>
#include "ibddb.h"

#define DATASET_DICTIONARY "/dictionary"
#define DATASET_PAIRS "/pairs"
#define DATASET_IBD "/ibd"
#define DATASET_MARKERS "/markers"
#define DATASET_PEDIGREE "/pedigree"
#define DEFAULT_TRESHOLD_LIMIT 0.1f

static int MarkerCompareByChromName(const void* a,const void *b)
	{
	const MarkerPtr i=(const MarkerPtr)a;
	const MarkerPtr j=(const MarkerPtr)b;
	int d=i->tid - j->tid;
	if(d!=0) return d;
	return strcmp(i->name, j->name);
	}


static int MarkerCompareByLoc(const void* a,const void *b)
	{
	const MarkerPtr i=(const MarkerPtr)a;
	const MarkerPtr j=(const MarkerPtr)b;
	int d=i->tid - j->tid;
	if(d!=0) return d;
	return i->position - j->position;
	}

static int IndividualCompareByFamName(const void* a,const void *b)
	{
	const IndividualPtr i=(const IndividualPtr)a;
	const IndividualPtr j=(const IndividualPtr)b;
	int d=strcmp(i->family,  j->family);
	if(d!=0) return d;
	return strcmp(i->name,  j->name);
	}

static ChromPtr findChromosomeByName(ContextPtr ctx,const char* cname)
	{
	size_t i=0;
	for(i=0;i< ctx->chromosome_count;++i)
		{
		if( strcmp(ctx->chromosomes[i].name,cname)==0)
			{
			return &ctx->chromosomes[i];
			}
		}
	return NULL;
	}

static ChromPtr findExistingChromosomeByName(ContextPtr ctx,const char* cname)
	{
	ChromPtr c=findChromosomeByName(ctx,cname);
	if(c==NULL)
		{
		DIE_FAILURE("Cannot find chromosome \"%s\" in dictionary.\n",cname);
		}
	return c;
	}


static IndividualPtr findIndividualByFamName(ContextPtr ctx,const char* fam,const char* indi)
	{
	IndividualPtr found;
	Individual key;
	key.family=(char*)fam;
	key.name=(char*)indi;
	found=(IndividualPtr)bsearch(
		(const void*)&key,
		ctx->individuals,
		ctx->individual_count,
              	sizeof(Individual),
		IndividualCompareByFamName
		);
	if(found==NULL) DIE_FAILURE("undefined individual %s:%s",fam,indi);
	return found;
	}

static int PairIndiCompare(const void* a,const void *b)
	{
	const PairIndiPtr i=(const PairIndiPtr)a;
	const PairIndiPtr j=(const PairIndiPtr)b;
	int d= i->indi1idx - j->indi1idx;
	if(d!=0) return d;
	return i->indi2idx - j->indi2idx;
	}


static void readIbd(ContextPtr ctx)
	{
	char* line1;
	size_t i,line_len=0UL;
	gzFile in1,in;
	MarkerPtr markersbyname=(MarkerPtr)safeCalloc(ctx->marker_count,sizeof(Marker));
	memcpy((void*)markersbyname,(void*)ctx->markers,ctx->marker_count*sizeof(Marker));
	qsort(
		(void*)markersbyname,
		ctx->marker_count,
		sizeof (Marker),
		MarkerCompareByChromName
		);
	
	if(ctx->ibd_filename==NULL)
		{
		DIE_FAILURE("config->ibd_filename undefined.\n");
		}
	DEBUG("Opening IBD %s",ctx->ibd_filename);
	in1=safeGZOpen(ctx->ibd_filename,"r");
	while((line1=gzReadLine(in1,&line_len))!=NULL)
		{
		char *line;
		DEBUG("Step 1: Opening IBD %s",line1);
		in=safeGZOpen(line1,"r");			
		
		while((line=gzReadLine(in,&line_len))!=NULL)
			{
			IndividualPtr p1;
			IndividualPtr p2;
			PairIndi key;
			PairIndiPtr found;
			char* tokens[9];
			if(strsplit(line,'\t',tokens,9)<9)
				{
				DIE_FAILURE("BOUM IBD in %s",line1);
				}
			p1= findIndividualByFamName(ctx,tokens[0],tokens[1]);
			p2= findIndividualByFamName(ctx,tokens[2],tokens[3]);

			if(p1->index < p2->index)
				{
				key.indi1idx=p1->index;
				key.indi2idx=p2->index;
				}
			else
				{
				key.indi1idx=p2->index;
				key.indi2idx=p1->index;
				}

			found=(PairIndiPtr)bsearch(
				(const void*)&key,
				ctx->pairs,
				ctx->pair_count,
		              	sizeof(PairIndi),
				PairIndiCompare
				);
			if(found==NULL)
				{
				/*DEBUG("insert pair %s %s / %s %s (pairs:%d)",
					tokens[0],tokens[1],
					tokens[2],tokens[3],
					ctx->pair_count
					);*/ 
				ctx->pairs = (PairIndiPtr)safeRealloc(
					ctx->pairs,
					sizeof(PairIndi)*(ctx->pair_count+1)
					);
				memcpy(&ctx->pairs[ctx->pair_count],&key,sizeof(PairIndi));
				ctx->pair_count++;

				qsort(
					(void*)ctx->pairs,
					ctx->pair_count,
					sizeof (PairIndi),
					PairIndiCompare
					);
				}
			free(line);
			}

		gzclose(in);
		free(line1);		
		}
	gzclose(in1);
	for( i=0;i< ctx->pair_count;++i)
		{
		ctx->pairs[i].index=(int)i;
		}
	{
	 hsize_t  dims[1] = {ctx->pair_count};

	int pairtype = H5Tcreate (H5T_COMPOUND, sizeof (PairIndi));
        H5Tinsert(pairtype, "indi1idx", HOFFSET(PairIndi, indi1idx), H5T_NATIVE_INT);
 	H5Tinsert(pairtype, "indi2idx", HOFFSET(PairIndi, indi2idx), H5T_NATIVE_INT);
	H5Tinsert(pairtype, "index", HOFFSET(PairIndi, index), H5T_NATIVE_INT);
	
	

	 /* Create the data space for the dataset. */
	int dataspace_id =  H5Screate_simple (1,dims, NULL);
	/* create dataset */
	int dataset_id = H5Dcreate2(
			ctx->file_id,
			  DATASET_PAIRS,
			  pairtype,
			  dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	 H5Dwrite(
			dataset_id,
			pairtype,
			H5S_ALL, H5S_ALL, H5P_DEFAULT,
			ctx->pairs
			);

	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	H5Tclose(pairtype);
	}	

	
	{

	hsize_t  dims_memory[3]={1,1,3};
	hsize_t  dims[3] = {ctx->marker_count,ctx->pair_count,3};
	hid_t  memspace  = H5Screate_simple(3, dims_memory, NULL); 
	

	int dataspace_id =  H5Screate_simple (3,dims, NULL);
	int dataset_id = H5Dcreate2(
			ctx->file_id,
			  DATASET_IBD,
			  H5T_NATIVE_FLOAT,
			  dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	

	in1=safeGZOpen(ctx->ibd_filename,"r");
	while((line1=gzReadLine(in1,&line_len))!=NULL)
		{
		ChromPtr prev_chrom=NULL;
		char *line;
		DEBUG("Step 2 : Opening IBD %s",line1);
		in=safeGZOpen(line1,"r");			
		
		while((line=gzReadLine(in,&line_len))!=NULL)
			{
			IndividualPtr p1;
			IndividualPtr p2;
			PairIndi key;
			float ibd_values[3]={0.0f,0.0f,0.0f};
			PairIndiPtr found;
			char* tokens[9];
			if(strsplit(line,'\t',tokens,10)<9)
				{
				DIE_FAILURE("BOUM IBD in %s",line1);
				}
			p1= findIndividualByFamName(ctx,tokens[0],tokens[1]);
			p2= findIndividualByFamName(ctx,tokens[2],tokens[3]);

			if(p1->index < p2->index)
				{
				key.indi1idx=p1->index;
				key.indi2idx=p2->index;
				}
			else
				{
				key.indi1idx=p2->index;
				key.indi2idx=p1->index;
				}

			found=(PairIndiPtr)bsearch(
				(const void*)&key,
				ctx->pairs,
				ctx->pair_count,
		              	sizeof(PairIndi),
				PairIndiCompare
				);
			if(found==NULL) 
				{
				DIE_FAILURE("insert pair %s %s / %s %s",
					tokens[0],tokens[1],
					tokens[2],tokens[3]
					); 
				}
				

			
			/* get the chromosome and its' index */
			if(prev_chrom==NULL || strcmp(prev_chrom->name,tokens[4])!=0)
				{
				prev_chrom=findChromosomeByName(ctx,tokens[4]);
				if(prev_chrom==NULL)
					{
					DIE_FAILURE("unknown chromosome %s",tokens[4]);
					}
				}		

			Marker key2;
			key2.tid=prev_chrom->tid;
			key2.name=(char*)tokens[5];
			/* get the marker */
			MarkerPtr marker= (MarkerPtr)bsearch(
				(const void*)&key2,
				markersbyname,
				ctx->marker_count,
		              	sizeof(Marker),
				MarkerCompareByChromName
				);
			if(marker==NULL) DIE_FAILURE("unknown marker %s %s",tokens[4],tokens[5]);

			
			for(i=0;i< 3;++i)
				{
				char* p2;
				ibd_values[i]=(float)strtod(tokens[6+i],&p2);
				if((*p2!=0))
					{
					DIE_FAILURE("bad ibd value column $%zu "PRIuPTR" in %s after '%s'. ",
						(6+i+1),tokens[6+i],p2);
					}
				if(ibd_values[i]<0.0f || ibd_values[i]>1.0f)
					{
					DIE_FAILURE("Column $%zu : bad ibd value =%f in \"%s\".",
						(6+i+1),ibd_values[i],tokens[6+i]);
					}
				}
			
			hsize_t write_start[3] = {marker->index,found->index,0};
			hsize_t write_count[3] = {1,1,3};
			
			assert(write_start[0] < ctx->marker_count);
			assert(write_start[1] < ctx->pair_count);

			
			VERIFY(H5Sselect_hyperslab(
				dataspace_id,
				H5S_SELECT_SET,
				write_start, NULL, 
				write_count, NULL
				));
			
			VERIFY(H5Dwrite(
				dataset_id,
				H5T_NATIVE_FLOAT,
				memspace,
				dataspace_id,
				H5P_DEFAULT,
				ibd_values
				));
		
			free(line);
			}

		gzclose(in);
		free(line1);		
		}
	gzclose(in1);
	
	H5Sclose(memspace);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);

	}



	free(markersbyname);

	}



static void readPed(ContextPtr ctx)
	{
	char* line;
	size_t i,line_len=0UL;
	gzFile in;

	if(ctx->ped_filename==NULL)
		{
		DIE_FAILURE("config->ped_filename undefined.\n");
		}
	DEBUG("Opening PED \"%s\"",ctx->ped_filename);
	in=safeGZOpen(ctx->ped_filename,"r");
	while((line=gzReadLine(in,&line_len))!=NULL)
		{
		IndividualPtr individual=NULL;
		char* tokens[6];
		if(strsplit(line,' ',tokens,6)<6)
			{
			DIE_FAILURE("BOUM PED");
			}
		ctx->individuals = (IndividualPtr)safeRealloc(
			ctx->individuals,
			sizeof(Individual)*(ctx->individual_count+1));
		individual=&ctx->individuals[ ctx->individual_count ];
		individual->family = safeStrDup( tokens[0] );
		individual->name = safeStrDup( tokens[1] );
		individual->father = (strlen(tokens[2])==0 || strcmp(tokens[2],"0")==0 ? NULL: safeStrDup( tokens[2] ));
		individual->mother = (strlen(tokens[3])==0 || strcmp(tokens[3],"0")==0 ? NULL: safeStrDup( tokens[3] ));
		individual->sex=atoi(tokens[4]);
		individual->status=atoi(tokens[5]);
		individual->index=0;

		ctx->individual_count++;
		free(line);		
		}
	gzclose(in);

	qsort(
		(void*)ctx->individuals,
		ctx->individual_count,
		sizeof (Individual),
		IndividualCompareByFamName
		);
	for(i=0;i< ctx->individual_count;++i) ctx->individuals[i].index=i;

	
	{
	 hsize_t  dims[1] = {ctx->individual_count};
	int   strtype = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype, H5T_VARIABLE);
	
	int pedigreetype = H5Tcreate (H5T_COMPOUND, sizeof (Individual));
	H5Tinsert(pedigreetype, "family", HOFFSET(Individual, family), strtype);
	H5Tinsert(pedigreetype, "name", HOFFSET(Individual, name), strtype);
	H5Tinsert(pedigreetype, "father", HOFFSET(Individual, father), strtype);
	H5Tinsert(pedigreetype, "mother", HOFFSET(Individual, mother), strtype);
        H5Tinsert(pedigreetype, "sex", HOFFSET(Individual, sex), H5T_NATIVE_INT);
 	H5Tinsert(pedigreetype, "status", HOFFSET(Individual, status), H5T_NATIVE_INT);
	H5Tinsert(pedigreetype, "index", HOFFSET(Individual, index), H5T_NATIVE_INT);
	
	

	 /* Create the data space for the dataset. */
	int dataspace_id =  H5Screate_simple (1,dims, NULL);
	/* create dataset */
	int dataset_id = H5Dcreate2(
			ctx->file_id,
			  DATASET_PEDIGREE,
			  pedigreetype,
			  dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	 H5Dwrite(
			dataset_id,
			pedigreetype,
			H5S_ALL, H5S_ALL, H5P_DEFAULT,
			ctx->individuals
			);

	H5Tclose (strtype);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	H5Tclose(pedigreetype);
	}

	}



static void readBed(ContextPtr ctx)
	{
	char* line;
	ChromPtr prev_chrom=NULL;
	size_t i,line_len=0UL;
	gzFile in;

	if(ctx->bed_filename==NULL)
		{
		DIE_FAILURE("config->bed_filename undefined.\n");
		}
	DEBUG("Opening BED %s",ctx->bed_filename);
	in=safeGZOpen(ctx->bed_filename,"r");
	while((line=gzReadLine(in,&line_len))!=NULL)
		{
		MarkerPtr marker=NULL;
		char* tokens[4];
		
		if(strsplit(line,'\t',tokens,4)<4)
			{
			DIE_FAILURE("BOUM BED");
			}
		ctx->markers = (MarkerPtr)safeRealloc(
			ctx->markers,
			sizeof(Marker)*(ctx->marker_count+1));
		
		marker=&ctx->markers[ ctx->marker_count ];
		if(prev_chrom==NULL || strcmp(prev_chrom->name,tokens[0])!=0)
			{
			size_t i=0;
			prev_chrom=NULL;
			for(i=0;i< ctx->chromosome_count;++i)
				{
				if( strcmp(ctx->chromosomes[i].name,tokens[0])==0)
					{
					prev_chrom=&ctx->chromosomes[i];
					break;
					}
				}
			if(prev_chrom==NULL)
				{
				DIE_FAILURE("unknown chromosome %s",tokens[0]);
				}
			}		

		marker->tid=prev_chrom->tid;
		marker->index=0;
		marker->name=safeStrDup(tokens[3]);
		marker->position=atoi(tokens[1]);
		if(marker->position<=0)
			{
			DIE_FAILURE("Marker positon<=0 for  %s := %s", tokens[3],tokens[1]);
			}

		ctx->marker_count++;
		free(line);
		}
	gzclose(in);

	qsort(
		(void*)ctx->markers,
		ctx->marker_count,
		sizeof (Marker),
		MarkerCompareByLoc
		);
	for(i=0;i< ctx->marker_count;++i) ctx->markers[i].index=i;
	
	{
	 hsize_t  dims[1] = {ctx->marker_count};
	int   strtype = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype, H5T_VARIABLE);
	
	int markertype = H5Tcreate (H5T_COMPOUND, sizeof (Marker));
	H5Tinsert(markertype, "name", HOFFSET(Marker, name), strtype);
        H5Tinsert(markertype, "tid", HOFFSET(Marker, tid), H5T_NATIVE_INT);
 	H5Tinsert(markertype, "position", HOFFSET(Marker, position), H5T_NATIVE_INT);
	H5Tinsert(markertype, "index", HOFFSET(Marker, index), H5T_NATIVE_INT);
	
	

	 /* Create the data space for the dataset. */
	int dataspace_id =  H5Screate_simple (1,dims, NULL);
	/* create dataset */
	int dataset_id = H5Dcreate2(
			ctx->file_id,
			  DATASET_MARKERS,
			  markertype,
			  dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	 H5Dwrite(
			dataset_id,
			markertype,
			H5S_ALL, H5S_ALL, H5P_DEFAULT,
			ctx->markers
			);

	H5Tclose (strtype);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	H5Tclose(markertype);
	}


	
	

	DEBUG("Closing bed dataset");

	}



static void readFaidx(ContextPtr ctx)
	{
	char* line;

	size_t line_len=0UL;
	gzFile in;
	if(ctx->faidx_filename==NULL)
		{
		DIE_FAILURE("config->faidx_filename undefined.\n");
		}
	DEBUG("Opening FAIDX %s",ctx->faidx_filename);
	in=safeGZOpen(ctx->faidx_filename,"r");
	while((line=gzReadLine(in,&line_len))!=NULL)
		{
		ChromPtr chrom=NULL;
		char* tokens[3];
		
		if(strsplit(line,'\t',tokens,3)<3)
			{
			DIE_FAILURE("BOUM FAIDX");
			}
		ctx->chromosomes = (ChromPtr)safeRealloc(
			ctx->chromosomes,
			sizeof(Chrom)*(ctx->chromosome_count+1));
		chrom=&ctx->chromosomes[ ctx->chromosome_count ];
		

		chrom->tid=(int)ctx->chromosome_count;
		chrom->name=strdup(tokens[0]);
		if( chrom->name == NULL)
			{
			DIE_FAILURE("name == NULL");
			}
		chrom->length=atoi(tokens[1]);
		if(chrom->length<=0)
			{
			DIE_FAILURE("Chromosome length<=0 for  %s := %s", tokens[0],tokens[1]);
			}

		ctx->chromosome_count++;
		free(line);
		}
	gzclose(in);
	
	{
	 hsize_t  dims[1] = {ctx->chromosome_count};
	int   strtype = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype, H5T_VARIABLE);
	
	hid_t chromtype = H5Tcreate (H5T_COMPOUND, sizeof (Chrom));
	H5Tinsert(chromtype, "name", HOFFSET(Chrom, name), strtype);
        H5Tinsert(chromtype, "tid", HOFFSET(Chrom, tid), H5T_NATIVE_INT);
 	H5Tinsert(chromtype, "length", HOFFSET(Chrom, length), H5T_NATIVE_INT);
	


	 /* Create the data space for the dataset. */
	int dataspace_id =  H5Screate_simple (1,dims, NULL);
	/* create dataset */
	int dataset_id = H5Dcreate2(
			ctx->file_id,
			  DATASET_DICTIONARY,
			  chromtype,
			  dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	 H5Dwrite(
			dataset_id,
			chromtype,
			H5S_ALL, H5S_ALL, H5P_DEFAULT,
			ctx->chromosomes
			);

	H5Tclose (strtype);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	H5Tclose(chromtype);
	}


	
	

	DEBUG("Closing faidx dataset");
	//assertGE0(H5Dclose(faidx_dataset_id));*/
	}

static ContextPtr ContextNew(int argc,char** argv)
	{
	ContextPtr config=(ContextPtr)safeCalloc(1,sizeof(Context));
	config->argc = argc;
	config->argv = argv;
	config->out = stdout;
	return config;
	}


static int main_build(int argc,char** argv)
	{
	ContextPtr config = ContextNew(argc,argv);

	
	for(;;)
		{
		struct option long_options[] =
		     {
		      // {"enable-self-self",  no_argument , &config->enable_self_self , 1},
			{"out",    required_argument, 0, 'o'},
 			{"dict",    required_argument, 0, 'D'},
			{"bed",    required_argument, 0, 'b'},
			{"ped",    required_argument, 0, 'p'},
			{"ibd",    required_argument, 0, 'i'},
		       {0, 0, 0, 0}
		     };
		 /* getopt_long stores the option index here. */
		int option_index = 0;
	     	int c = getopt_long (argc, argv, "o:D:b:p:i:",
		                    long_options, &option_index);
		if(c==-1) break;
		switch(c)
			{
			case 'o': config->hdf5_filename=optarg;break;
			case 'D': config->faidx_filename=optarg;break;
			case 'b': config->bed_filename=optarg;break;
			case 'p': config->ped_filename=optarg;break;
			case 'i': config->ibd_filename=optarg;break;
			case 0: break;
			case '?': break;
			default: exit(EXIT_FAILURE); break;
			}
		}
	if(config->hdf5_filename==NULL)
		{
		DIE_FAILURE("config->hdf5_filename undefined.\n");
		}
	if(config->bed_filename==NULL)
		{
		DIE_FAILURE("config->bed_filename undefined.\n");
		}
	if(config->ped_filename==NULL)
		{
		DIE_FAILURE("config->ped_filename undefined.\n");
		}
	if(config->ibd_filename==NULL)
		{
		DIE_FAILURE("config->ibd_filename undefined.\n");
		}
	DEBUG("Opening HDF5 file %s",config->hdf5_filename );
	config->file_id = H5Fcreate(config->hdf5_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if( config->file_id < 0)
		{
		DIE_FAILURE("H5Fcreate failed err=%d.\n", config->file_id);
		}
	readPed(config);
	readFaidx(config);
	readBed(config);
	readIbd(config);
	DEBUG("Closing HDF5 file");
	assertGE0(H5Fclose(config->file_id)); 
	return EXIT_SUCCESS;
	}


#define LOAD_CONFIG_DATASET(DATASETNAME,DATATYPE,ITEM_NAME,ITEM_COUNT) \
		DEBUG("Loading " DATASETNAME); \
		hid_t dataset_id = VERIFY(H5Dopen(config->file_id,DATASETNAME, H5P_DEFAULT)); \
		hid_t dspace = VERIFY(H5Dget_space(dataset_id)); \
		assert(H5Sget_simple_extent_ndims(dspace)==1); \
		int atype  = H5Dget_type(dataset_id);  \
		hsize_t dims[1]; \
		H5Sget_simple_extent_dims(dspace, dims, NULL); \
		config->ITEM_COUNT = dims[0]; \
		config->ITEM_NAME = (DATATYPE*)safeCalloc(config->ITEM_COUNT,sizeof(DATATYPE)); \
		VERIFY(H5Dread(dataset_id, atype, H5S_ALL, H5S_ALL, H5P_DEFAULT, config->ITEM_NAME)); \
		VERIFY(H5Sclose(dspace)); \
		VERIFY(H5Dclose(dataset_id)); \
		DEBUG("End reading " DATASETNAME)

			
void ContextOpenForRead(ContextPtr config)
	{
	DEBUG("Opening HDF5 file %s",config->hdf5_filename );
	config->file_id = H5Fopen(config->hdf5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if( config->file_id < 0)
		{
		DIE_FAILURE("H5Fopen failed err=%d.\n", config->file_id);
		}
	if( config->on_read_load_dict )
		{
		LOAD_CONFIG_DATASET(DATASET_DICTIONARY,Chrom,chromosomes,chromosome_count);
		}

	if( config->on_read_load_markers )
		{
		LOAD_CONFIG_DATASET(DATASET_MARKERS,Marker,markers,marker_count);
		}
	if( config->on_read_load_pedigree )
		{
		LOAD_CONFIG_DATASET(DATASET_PEDIGREE,Individual,individuals,individual_count);
		}
	if( config->on_read_load_pairs )
		{
		LOAD_CONFIG_DATASET(DATASET_PAIRS,PairIndi,pairs,pair_count);
		}
	}

static void ContextFree(ContextPtr config)
	{
	size_t i;

	if(config==NULL) return ;
	if(config->out!=NULL)
		{
		fflush(config->out);
		}
	
	for(i=0;i< config->marker_count;++i)
		{
		free( config->markers[i].name);
		}

	free(config->markers);
	
	
	for(i=0;i< config->chromosome_count;++i)
		{
		free( config->chromosomes[i].name);
		}

	free(config->chromosomes);

	
	if(config->file_id!=0)
		{
		VERIFY(H5Fclose(config->file_id)); 
		}
	free(config);
	}


static int main_dict(int argc,char** argv)
	{
	size_t i;
	ContextPtr config=ContextNew(argc,argv);
	config->on_read_load_dict = 1;
	
	if(optind+1!=argc)
		{
		fprintf(stderr,"Illegal number of arguments.\n");
		return EXIT_FAILURE;
		}	
	config->hdf5_filename=argv[optind];
	ContextOpenForRead(config);
	
	
	for(i=0;i< config->chromosome_count;++i)
		{
		fputs( config->chromosomes[i].name , config->out);
		fputc('\t', config->out);
		fprintf(config->out,"%d", config->chromosomes[i].length );
		fputc('\n', config->out);
		}
	

	ContextFree(config);
	return EXIT_SUCCESS;
	}


void parseRegion(ContextPtr ctx,RegionPtr rgn,const char* s)
	{
	ChromPtr chrom=NULL;
	char* colon=strchr(s,':'); 
	memset((void*)rgn,0,sizeof(Region));
	rgn->tid=-1;
	if(colon==NULL)
		{DEBUG("");
		chrom=findExistingChromosomeByName(ctx,s);
		rgn->tid = chrom->tid;
		rgn->start = 0;
		rgn->end = chrom->length;
		}
	else
		{
		char* tmp;
		char* hyphen=strchr(colon,'-'); 
		
		if(hyphen==NULL) DIE_FAILURE("'-' missing in %s.",s);
		tmp=safeStrNDup(s,colon-s);
		chrom=findExistingChromosomeByName(ctx,tmp);
		free(tmp);
		rgn->tid = chrom->tid;
		
		tmp=safeStrNDup(colon+1,hyphen-(colon+1));
		rgn->start = atoi(tmp);
		free(tmp);

		tmp=safeStrDup(hyphen+1);
		rgn->end = atoi(tmp);
		free(tmp);

		if(rgn->start > rgn->end)
			{
			DIE_FAILURE("bad region: '%s'.",s);
			}

		}
	}

static int main_markers(int argc,char** argv)
	{
	size_t i;
	ContextPtr config=ContextNew(argc,argv);
	config->on_read_load_dict = 1;
	config->on_read_load_markers = 1;
	RegionPtr region=NULL;
	char* rgn_str=NULL;
	for(;;)
		{
		struct option long_options[] =
		     {
		      // {"enable-self-self",  no_argument , &config->enable_self_self , 1},
			{"region",    required_argument, 0, 'r'},
		       {0, 0, 0, 0}
		     };
		 /* getopt_long stores the option index here. */
		int option_index = 0;
	     	int c = getopt_long (argc, argv, "r:",
		                    long_options, &option_index);
		if(c==-1) break;
		switch(c)
			{
			case 'r': rgn_str = optarg ;break;
			
			case 0: break;
			case '?': break;
			default: exit(EXIT_FAILURE); break;
			}
		}	
	if(optind+1!=argc)
		{
		fprintf(stderr,"Illegal number of arguments.\n");
		return EXIT_FAILURE;
		}	
	config->hdf5_filename=argv[optind];
	ContextOpenForRead(config);
	
	if(rgn_str!=NULL)
		{
		region=(RegionPtr)safeCalloc(1,sizeof(Region));
		parseRegion(config,region,rgn_str);
		DEBUG("region: tid=%d:%d-%d",region->tid,region->start,region->end);
		};
	

	for(i=0;i< config->marker_count;++i)
		{
		MarkerPtr marker = &config->markers[i];
		
		if( region!=NULL)
			{
			if( region->tid != marker->tid ) continue;
			if( marker->position < region->start) continue;
			if( marker->position > region->end) continue;
			}
		
		fputs( config->chromosomes[ marker->tid ].name , config->out);
		fputc('\t', config->out);
		fprintf(config->out,"%d", marker->position);
		fputc('\t', config->out);
		fprintf(config->out,"%d", marker->position+1);
		fputc('\t', config->out);
		fputs( marker->name , config->out);
		if( fputc('\n', config->out) < 0) break;
		}
	if(region!=NULL)
		{
		free(region);
		}
	ContextFree(config);
	return EXIT_SUCCESS;
	}

static void printIndividual(const IndividualPtr individual,FILE* out)
	{
	fputs( individual->family , out);
	fputc('\t',out);
	fputs( individual->name , out);
	fputc('\t', out);

	if(individual->father==NULL)
		{
		fputc('0', out);
		}
	else
		{
		fputs( individual->father , out);
		}
	fputc('\t', out);
	if(individual->mother==NULL)
		{
		fputc('0', out);
		}
	else
		{
		fputs( individual->mother , out);
		}
	fputc('\t', out);
	fprintf(out,"%d", individual->sex);
	fputc('\t', out);
	fprintf(out,"%d", individual->status);
	}


static int main_pedigree(int argc,char** argv)
	{
	size_t i;
	ContextPtr config=ContextNew(argc,argv);
	config->on_read_load_pedigree = 1;
	
	for(;;)
		{
		struct option long_options[] =
		     {
		      // {"enable-self-self",  no_argument , &config->enable_self_self , 1},

		       {0, 0, 0, 0}
		     };
		 /* getopt_long stores the option index here. */
		int option_index = 0;
	     	int c = getopt_long (argc, argv, "",
		                    long_options, &option_index);
		if(c==-1) break;
		switch(c)
			{
			case 0: break;
			case '?': break;
			default: exit(EXIT_FAILURE); break;
			}
		}	
	if(optind+1!=argc)
		{
		fprintf(stderr,"Illegal number of arguments.\n");
		return EXIT_FAILURE;
		}	
	config->hdf5_filename=argv[optind];
	ContextOpenForRead(config);
	

	for(i=0;i< config->individual_count;++i)
		{
		IndividualPtr individual = &config->individuals[i];
		
		
		printIndividual(individual,config->out);
		if( fputc('\n', config->out) < 0) break;
		}

	ContextFree(config);
	return EXIT_SUCCESS;
	}

static int main_pairs(int argc,char** argv)
	{
	size_t i;
	ContextPtr config=ContextNew(argc,argv);
	config->on_read_load_pedigree = 1;
	config->on_read_load_pairs = 1;
	
	for(;;)
		{
		struct option long_options[] =
		     {
		      // {"enable-self-self",  no_argument , &config->enable_self_self , 1},

		       {0, 0, 0, 0}
		     };
		 /* getopt_long stores the option index here. */
		int option_index = 0;
	     	int c = getopt_long (argc, argv, "",
		                    long_options, &option_index);
		if(c==-1) break;
		switch(c)
			{
			case 0: break;
			case '?': break;
			default: exit(EXIT_FAILURE); break;
			}
		}	
	if(optind+1!=argc)
		{
		fprintf(stderr,"Illegal number of arguments.\n");
		return EXIT_FAILURE;
		}	
	config->hdf5_filename=argv[optind];
	ContextOpenForRead(config);
	

	for(i=0;i< config->pair_count;++i)
		{
		PairIndiPtr pair = &config->pairs[i];
		printIndividual(&config->individuals[pair->indi1idx],config->out);
		fputc('\t', config->out);
		printIndividual(&config->individuals[pair->indi2idx],config->out);
		if( fputc('\n', config->out) < 0) break;
		}

	ContextFree(config);
	return EXIT_SUCCESS;
	}

static int main_ibd(int argc,char** argv)
	{
	float ibd_values[3];
	float treshold=DEFAULT_TRESHOLD_LIMIT;
	int print_header=TRUE;
	int print_pairs=TRUE;
	size_t i,j;
	ContextPtr config=ContextNew(argc,argv);
	RegionPtr region=NULL;
	char* rgn_str=NULL;
	config->on_read_load_pedigree = 1;
	config->on_read_load_pairs = 1;
	config->on_read_load_dict = 1;
	config->on_read_load_markers = 1;

	/** Graphics2D stuff */
	long genome_size=0L;
	Dimension imageDimension;	
	imageDimension.width=1000;
	imageDimension.height=300;
	ExpData* expData=NULL;
	size_t expData_count=0L;
	double max_y=0.0;
	char* image_filename;
	


	for(;;)
		{
		struct option long_options[] =
		     {
		      // {"enable-self-self",  no_argument , &config->enable_self_self , 1},
		       {"region",    required_argument, 0, 'r'},
		       {"noheader",  no_argument, &print_header, 0},
		       {"image",  required_argument, 0, 'g'},
		       {0, 0, 0, 0}
		     };
		 /* getopt_long stores the option index here. */
		int option_index = 0;
	     	int c = getopt_long (argc, argv, "r:g:",
		                    long_options, &option_index);
		if(c==-1) break;
		switch(c)
			{
			case 'r': rgn_str = optarg ;break;
			case 'g': image_filename = optarg ;break;
			case 0: break;
			case '?': break;
			default: exit(EXIT_FAILURE); break;
			}
		}	
	if(optind+1!=argc)
		{
		fprintf(stderr,"Illegal number of arguments.\n");
		return EXIT_FAILURE;
		}
	
	config->hdf5_filename=argv[optind];
	ContextOpenForRead(config);
	
	if(rgn_str!=NULL)
		{
		region=(RegionPtr)safeCalloc(1,sizeof(Region));
		parseRegion(config,region,rgn_str);
		genome_size=1+(region->end -  region->start);
		}
	else
		{
		for(i=0;i< config->chromosome_count;++i)
			{
			config->chromosomes[i].cumulative_start+=genome_size;
			genome_size+= config->chromosomes[i].length;
			}
		}


	DEBUG("Loading " DATASET_IBD); 
	hid_t dataset_id = VERIFY(H5Dopen(config->file_id,DATASET_IBD, H5P_DEFAULT)); 
	hid_t dataspace_id = VERIFY(H5Dget_space(dataset_id)); 
	hsize_t  dims_memory[3]={1,1,3};
	hid_t  memspace  = H5Screate_simple(3, dims_memory, NULL);
	
	
	


	if(print_header && image_filename==NULL)
		{
		fputs("CHROM\tPOS\tNAME",config->out);
		for(i=0;i< config->pair_count && print_pairs;++i)
			{
			PairIndiPtr pair = &config->pairs[i];
			//if(!pair->selected) continue;
			IndividualPtr indi1=&config->individuals[pair->indi1idx];
			IndividualPtr indi2=&config->individuals[pair->indi2idx];		
			fprintf(config->out,"\t%s:%s|%s:%s",
				indi1->family,indi1->name,
				indi2->family,indi2->name
				);

			}
		fputs("\tCOUNT_IBD",config->out);
		fputc('\n',config->out);
		}

	for(i=0;i< config->marker_count;++i)
		{
		MarkerPtr marker = &config->markers[i];
		int count_pairs=0;
		if( region!=NULL)
			{
			if( region->tid != marker->tid ) continue;
			if( marker->position < region->start) continue;
			if( marker->position > region->end) continue;
			}
		
		if(image_filename==NULL)
			{
			fputs( config->chromosomes[ marker->tid ].name , config->out);
			fputc('\t', config->out);
			fprintf(config->out,"%d", marker->position);
			fputc('\t', config->out);
			fputs( marker->name , config->out);
			}
		for(j=0;j< config->pair_count ;++j)
			{
			PairIndiPtr pair = &config->pairs[j];
			//if(!pair->selected) continue;
			
			

			hsize_t read_start[3] = {marker->index,pair->index,0};
			hsize_t read_count[3] = {1,1,3};
			
			
			VERIFY(H5Sselect_hyperslab(
				dataspace_id,
				H5S_SELECT_SET,
				read_start, NULL, 
				read_count, NULL
				));
			
			VERIFY(H5Dread(
				dataset_id,
				H5T_NATIVE_FLOAT,
				memspace,
				dataspace_id,
				H5P_DEFAULT,
				ibd_values
				));

			if(print_pairs && image_filename==NULL)
				{
				fputc('\t', config->out);
				fprintf(config->out,"%f",ibd_values[0]);
				}
			if( ibd_values[0] < treshold)
				{
				count_pairs++;
				} 
			}
		if(image_filename==NULL)
			{
			fprintf(config->out,"\t%d",count_pairs);
			if( fputc('\n', config->out) < 0) break;
			}
		else
			{
			if(count_pairs>0)
				{
				expData=(ExpDataPtr)safeRealloc(expData,(expData_count+1)*sizeof(ExpData));
				expData[expData_count].marker_index = marker->index;
				expData[expData_count].value = count_pairs;
				expData_count++;
				if( count_pairs > max_y) max_y=count_pairs;
				}
			}
		}
	
	if(image_filename!=NULL)
		{
		Rectangle drawingArea;
 		cairo_surface_t *surface=NULL;
  		cairo_t *cr=NULL;
		surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, imageDimension.width, imageDimension.height);
		if(surface==NULL) DIE_FAILURE("Cannot create image");
  		cr = cairo_create(surface);
		if(cr==NULL) DIE_FAILURE("Cannot create image context");

		drawingArea.x=	100;
		drawingArea.width= imageDimension.width-200;
		drawingArea.y=	50;
		drawingArea.height= imageDimension.height-100;		
		
		cairo_set_font_size(cr, 13);

#define BASE2PIXEL(chrom,POS)  (drawingArea.x+((region==NULL?POS:config->chromosomes[chrom].cumulative_start+POS)/(double)genome_size)*drawingArea.width)

		if(region!=NULL)
			{
			cairo_set_line_width (cr, 0.1);
			for(i=0;i< config->chromosome_count;++i)
				{
				if(i%2==0)
					{
					cairo_set_source_rgb (cr, 0, 0, 0);
					}
				else
					{
					cairo_set_source_rgb (cr, 0, 0, 0);
					}
				cairo_rectangle(cr,
					BASE2PIXEL(i,0),
					drawingArea.y,
					BASE2PIXEL(i,config->chromosomes[i].length)-BASE2PIXEL(i,0),
					drawingArea.height
					);
				cairo_fill(cr);

				cairo_move_to (cr, (BASE2PIXEL(i,0)+BASE2PIXEL(i,config->chromosomes[i].length))/2.0, 12);
				cairo_set_source_rgb (cr, 0.5, 0.5, 0.5);
				cairo_show_text (cr, config->chromosomes[i].name);
				}
			}
		else
			{
			cairo_move_to (cr, drawingArea.x+drawingArea.width/2.0, 12);
			cairo_set_source_rgb (cr, 1, 1, 0.5);
			cairo_show_text (cr, rgn_str);
			}

		for(i=0;i< expData_count;++i)
			{
			MarkerPtr marker=&config->markers[expData[i].marker_index];
			double cx = BASE2PIXEL(marker->tid,marker->position);
			double cy = drawingArea.y + drawingArea.height - (expData[i].value/max_y)*drawingArea.height ;
			
			cairo_move_to (cr, cx-5, cy);
			cairo_line_to (cr, cx+5, cy);
			cairo_move_to (cr, cx, cy-5);
			cairo_line_to (cr, cx, cy+5);
			cairo_stroke (cr);
			}

		//frame
		cairo_set_line_width (cr, 0.1);
		cairo_set_source_rgb (cr,0, 0, 0);
		cairo_rectangle (cr, drawingArea.x, drawingArea.y,drawingArea.width,drawingArea.height);
		cairo_stroke (cr);
		
		DEBUG("saving image as %s",image_filename);
		cairo_surface_write_to_png(surface, image_filename);
  		cairo_destroy(cr);
  		cairo_surface_destroy(surface);
		
		free(expData);
		}


	VERIFY(H5Sclose(memspace));
	VERIFY(H5Sclose(dataspace_id)); 
	VERIFY(H5Dclose(dataset_id)); 

	if(region!=NULL)
		{
		free(region);
		}

	ContextFree(config);
	return EXIT_SUCCESS;
	}


int main(int argc,char** argv)
	{
	int status=EXIT_FAILURE;
	if(argc>1)
		{
		if(strcmp("build",argv[1])==0)
			{
			status=main_build(argc-1,&argv[1]);
			}
		else if(strcmp("dict",argv[1])==0)
			{
			status=main_dict(argc-1,&argv[1]);
			}
		else if(strcmp("markers",argv[1])==0)
			{
			status= main_markers(argc-1,&argv[1]);
			}
		else if(strcmp("ped",argv[1])==0)
			{
			status= main_pedigree(argc-1,&argv[1]);
			}
		else if(strcmp("pairs",argv[1])==0)
			{
			status= main_pairs(argc-1,&argv[1]);
			}
		else if(strcmp("ibd",argv[1])==0)
			{
			status= main_ibd(argc-1,&argv[1]);
			}
		else
			{
			fprintf(stderr,"Unknown sub program %s.\n",argv[1]);
			}
		}
	else
		{
		fprintf(stderr,"sub program missing.\n");
		}
	
	DEBUG("%s: exiting with status=%d",argv[0],status);
	return status;
	}
