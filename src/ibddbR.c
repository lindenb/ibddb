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

static const char* _argv=__FILENAME__;

typedef struct ibd_handler_t
	{
	ContextPtr context;

	}IbdHandler,*IbdHandlerPtr;

SEXP RIbdDbClose(SEXP handle)
	{	
	void *p = R_ExternalPtrAddr(handle);
	DEBUG("closing H5 file");
	if(p!=NULL)
		{
		FILE* f=(FILE*)p;
		fclose(f);
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
	FILE* in;
	char* filename= CHAR(STRING_ELT(Rfilename, 0));
	if(filename==NULL) DIE_FAILURE("NULL Filename");
	DEBUG("opening H5 file %s",filename);
	in=fopen(filename,"r");
	if(in==NULL) DIE_FAILURE("Cannot open %s",filename);
	SEXP ext = PROTECT(R_MakeExternalPtr(in, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(ext,_RIbdDbFinalizer, TRUE);
    	UNPROTECT(1);
	return ext;
	}
#define  METHOD_COUNT_ENTITY(name,field) \
SEXP RIbdDbCount##name(SEXP handle) {\
	void *p = R_ExternalPtrAddr(handle);\
	if(p==NULL) return R_NilValue;\
	return INTEGER(((IbdHandlerPtr)p)->context->field);\
	}



METHOD_COUNT_ENTITY(Markers,marker_count)
METHOD_COUNT_ENTITY(Individuals,individual_count)
METHOD_COUNT_ENTITY(Pairs,pair_count)
METHOD_COUNT_ENTITY(Chromosome,chromosome_count)

SEXP RIbdDbGetMarkerAt(SEXP handle,SEXP index)
	{	
	void *p = R_ExternalPtrAddr(handle);
	if(p==NULL) return R_NilValue;

	MarkerPtr marker=NULL;
	SEXP res = PROTECT(allocVector(VECSXP, 3));
	SET_VECTOR_ELT(res, 0, mkString(marker->name));
	SET_VECTOR_ELT(res, 1, ScalarInteger(marker->tid));
	SET_VECTOR_ELT(res, 2, ScalarInteger(marker->position));

	SEXP sNames = PROTECT(allocVector(STRSXP, 3));
	SET_STRING_ELT(sNames, 0, mkString("name"));
	SET_STRING_ELT(sNames, 1, mkString("tid"));
	SET_STRING_ELT(sNames, 2, mkString("position"));
	setAttrib(res, R_NamesSymbol, sNames);


	UNPROTECT(2)
	return res;
	}



