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
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "ibddb.h"

/**
 * See also: 
 * http://stackoverflow.com/questions/7032617
 * http://adv-r.had.co.nz/C-interface.html
 * http://stackoverflow.com/questions/8720550
 */



typedef struct ibd_handler_t
	{
	ContextPtr context;
	IbdDataSetPtr ds_param;
	}IbdHandler,*IbdHandlerPtr;

SEXP RIbdDbClose(SEXP handle)
	{	
	void *p = R_ExternalPtrAddr(handle);
	DEBUG("closing H5 file");
	if(p!=NULL)
		{
		IbdHandlerPtr h=(IbdHandlerPtr)p;
		IbdDataSetClose(h->ds_param);
		ContextFree(h->context);
		free(h);
		}
	R_ClearExternalPtr(handle);
	return ScalarLogical(0);
	}

static void _RIbdDbFinalizer(SEXP handle)
	{
	RIbdDbClose(handle);
	}

SEXP RIbdDbOpen(SEXP Rfilename)
	{
	IbdHandlerPtr handler;
	const char* filename= CHAR(STRING_ELT(Rfilename, 0));
	if(filename==NULL) DIE_FAILURE("NULL Filename");
	DEBUG("opening H5 file %s",filename);
	handler = (IbdHandlerPtr)safeCalloc(1,sizeof(IbdHandler));
	handler->context = ContextNew(0,NULL);
	handler->context->on_read_load_pedigree = 1;
	handler->context->on_read_load_pairs = 1;
	handler->context->on_read_load_dict = 1;
	handler->context->on_read_load_markers = 1;
	handler->context->hdf5_filename=(char*)filename;
	ContextOpenForRead(handler->context);
	handler->ds_param = IbdDataSetOpen(handler->context);
	SEXP ext = PROTECT(R_MakeExternalPtr(handler, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(ext,_RIbdDbFinalizer, TRUE);
    	UNPROTECT(1);
	return ext;
	}
#define  METHOD_COUNT_ENTITY(name,field) \
SEXP RIbdDbCount##name(SEXP handle) {\
	SEXP Rlist; int *idim;\
	void *p = R_ExternalPtrAddr(handle);\
	if(p==NULL) return R_NilValue;\
	PROTECT(Rlist = allocVector(INTSXP, 1));\
	idim=INTEGER(Rlist);\
	*idim = (int)((IbdHandlerPtr)p)->context->field;\
	 UNPROTECT(1);\
	return Rlist;\
	}



METHOD_COUNT_ENTITY(Markers,marker_count)
METHOD_COUNT_ENTITY(Individuals,individual_count)
METHOD_COUNT_ENTITY(Pairs,pair_count)
METHOD_COUNT_ENTITY(Chromosomes,chromosome_count)

#define GET_ITEM_AT(ITEM_CLASS,ITEM_NAME,ITEM_COUNT) \
	void *p = R_ExternalPtrAddr(handle);\
	if(p==NULL) return R_NilValue;\
	ContextPtr ctx = ((IbdHandlerPtr)p)->context;\
	if( GET_LENGTH(index_r) !=1) DIE_FAILURE("index length!=1.");\
	int index=(INTEGER(index_r)[0]);\
	if( index<0 || index>= ctx->ITEM_COUNT)\
		{\
		DIE_FAILURE("index out of range.");\
		}\
	ITEM_CLASS item = &(ctx->ITEM_NAME)[index]

SEXP RIbdDbGetMarkerAt(SEXP handle,SEXP index_r)
	{
	GET_ITEM_AT(MarkerPtr,markers,marker_count);
	
	SEXP res = PROTECT(allocVector(VECSXP, 4));
	SET_VECTOR_ELT(res, 0, mkString(item->name));
	SET_VECTOR_ELT(res, 1, ScalarInteger(item->tid));
	SET_VECTOR_ELT(res, 2, ScalarInteger(item->position));
	SET_VECTOR_ELT(res, 3, ScalarInteger(index));

	SEXP sNames = PROTECT(allocVector(STRSXP, 4));
	SET_STRING_ELT(sNames, 0, mkChar("name"));
	SET_STRING_ELT(sNames, 1, mkChar("tid"));
	SET_STRING_ELT(sNames, 2, mkChar("position"));
	SET_STRING_ELT(sNames, 3, mkChar("index"));
	setAttrib(res, R_NamesSymbol, sNames);

	UNPROTECT(2);
	return res;
	}


SEXP RIbdDbGetChromosomeAt(SEXP handle,SEXP index_r)
	{	
	GET_ITEM_AT(ChromPtr,chromosomes,chromosome_count);
	
	SEXP res = PROTECT(allocVector(VECSXP, 3));
	SET_VECTOR_ELT(res, 0, mkString(item->name));
	SET_VECTOR_ELT(res, 1, ScalarInteger(item->tid));
	SET_VECTOR_ELT(res, 2, ScalarInteger(item->length));

	SEXP sNames = PROTECT(allocVector(STRSXP, 3));
	SET_STRING_ELT(sNames, 0, mkChar("name"));
	SET_STRING_ELT(sNames, 1, mkChar("tid"));
	SET_STRING_ELT(sNames, 2, mkChar("length"));
	setAttrib(res, R_NamesSymbol, sNames);

	UNPROTECT(2);
	return res;
	}


SEXP RIbdDbGetPairAt(SEXP handle,SEXP index_r)
	{	
	GET_ITEM_AT(PairIndiPtr,pairs,pair_count);
	
	SEXP res = PROTECT(allocVector(VECSXP, 3));
	SET_VECTOR_ELT(res, 0, ScalarInteger(item->index));
	SET_VECTOR_ELT(res, 1, ScalarInteger(item->indi1idx));
	SET_VECTOR_ELT(res, 2, ScalarInteger(item->indi2idx));

	SEXP sNames = PROTECT(allocVector(STRSXP, 3));
	SET_STRING_ELT(sNames, 0, mkChar("index"));
	SET_STRING_ELT(sNames, 1, mkChar("indi1idx"));
	SET_STRING_ELT(sNames, 2, mkChar("indi2idx"));
	setAttrib(res, R_NamesSymbol, sNames);

	UNPROTECT(2);
	return res;
	}
#define MKSTRING_OR_NULL(s) ((s)==NULL?R_NilValue : mkString(s))

SEXP RIbdDbGetIndividualAt(SEXP handle,SEXP index_r)
	{	
	GET_ITEM_AT(IndividualPtr,individuals,individual_count);
	
	SEXP res = PROTECT(allocVector(VECSXP,7));
	SET_VECTOR_ELT(res, 0, MKSTRING_OR_NULL(item->family));
	SET_VECTOR_ELT(res, 1, MKSTRING_OR_NULL(item->name));
	SET_VECTOR_ELT(res, 2, MKSTRING_OR_NULL(item->father));
	SET_VECTOR_ELT(res, 3, MKSTRING_OR_NULL(item->mother));
	SET_VECTOR_ELT(res, 4, ScalarInteger(item->sex));
	SET_VECTOR_ELT(res, 5, ScalarInteger(item->status));
	SET_VECTOR_ELT(res, 6, ScalarInteger(item->index));


	SEXP sNames = PROTECT(allocVector(STRSXP, 7));
	SET_STRING_ELT(sNames, 0, mkChar("family"));
	SET_STRING_ELT(sNames, 1, mkChar("name"));
	SET_STRING_ELT(sNames, 2, mkChar("father"));
	SET_STRING_ELT(sNames, 3, mkChar("mother"));
	SET_STRING_ELT(sNames, 4, mkChar("sex"));
	SET_STRING_ELT(sNames, 5, mkChar("status"));
	SET_STRING_ELT(sNames, 6, mkChar("index"));
	setAttrib(res, R_NamesSymbol, sNames);

	UNPROTECT(2);
	return res;
	}

SEXP priv_RIbdDbGetIBD(SEXP handle,SEXP marker_y,SEXP pair_y,int ibd_index)
	{	
	void *p = R_ExternalPtrAddr(handle);
	if(p==NULL) return R_NilValue;
	IbdHandlerPtr handler=(IbdHandlerPtr)p;
	float ibd_values[3];
	DEBUG("1 ");
	int marker_index=INTEGER_VALUE(marker_y);//(INTEGER(marker_y)[0]);
	int pair_index=INTEGER_VALUE(pair_y);//((INTEGER(pair_y)[0]);
	DEBUG("2 %d %d",marker_index,pair_index);
	DEBUG("%p",handler);
	if(marker_index<0 || marker_index >=  handler->context->marker_count) 
		{
		DIE_FAILURE("marker index out of range.");
		}
	DEBUG("");
	 pair_index=0;
	if(pair_index<0 || pair_index >=  handler->context->pair_count) 
		{
		DIE_FAILURE("pair index out of range.");
		}	
	DEBUG("");
	hsize_t read_start[3] = {marker_index,pair_index,0};
	hsize_t read_count[3] = {1,1,3};


	
	
	
	VERIFY(H5Sselect_hyperslab(
				handler->ds_param->dataspace_id,
				H5S_SELECT_SET,
				read_start, NULL, 
				read_count, NULL
				));
			
	VERIFY(H5Dread(
		handler->ds_param->dataset_id,
		H5T_NATIVE_FLOAT,
		handler->ds_param->memspace,
		handler->ds_param->dataspace_id,
		H5P_DEFAULT,
		ibd_values
		));
	if( ibd_values[ibd_index] <0.0 || ibd_values[ibd_index]>1.0) return  Rf_ScalarReal(R_NaN);//Rf_ScalarReal(NAN);
	return Rf_ScalarReal(ibd_values[ibd_index]);
	}
#define GET_IBD_AT(INDEX) SEXP RIbdDbGetIBD##INDEX(SEXP handle,SEXP marker_y,SEXP pair_y) \
	{\
	return priv_RIbdDbGetIBD(handle,marker_y,pair_y,INDEX);\
	}
	
GET_IBD_AT(0)
GET_IBD_AT(1)
GET_IBD_AT(2)
