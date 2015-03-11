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
#ifndef UTILS_H
#define UTILS_H
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <zlib.h>
#include "githash.h"

typedef signed char boolean_t;

/* gz file */
gzFile safeGZOpen(const char *path, const char *mode);
char* gzReadLine(gzFile in,size_t *str_size);

/* string */
char* ltrim(char* src,size_t* len);
size_t strchcount(const char* src,char c);
size_t strsplit(char* src,char delim,char** offsets,size_t max_tokens);
char* safeStrDup(const char* src);
char* safeStrNDup(const char* src,size_t n);
int strEndsWith(const char* s,const char* suff);
int strStartsWith(const char* s,const char* suff);

/** stdlib */
void* _safeMalloc(const char*,int,size_t);
void* _safeCalloc(const char*,int,size_t,size_t);
void* _safeRealloc(const char*,int,void*,size_t);

#define safeMalloc(n) _safeMalloc(__FILE__,__LINE__,n)
#define safeCalloc(n,m) _safeCalloc(__FILE__,__LINE__,n,m)
#define safeRealloc(p,n) _safeRealloc(__FILE__,__LINE__,p,n)


#define WHERENL fprintf(stderr,"[%s:%d] ",__FILE__,__LINE__)
#define DIE_FAILURE(FormatLiteral,...) do { WHERENL; fprintf (stderr,"Git-Hash:"GIT_HASH". Exiting: " FormatLiteral "\n", ##__VA_ARGS__); exit(EXIT_FAILURE);} while(0)
#define DEBUG(FormatLiteral, ...)  do { fputs("[DEBUG]",stderr); WHERENL; fprintf (stderr,"" FormatLiteral "\n", ##__VA_ARGS__);} while(0)

#ifndef MIN 
#define MIN(a,b) (a<b ? a: b)
#endif
#ifndef MAX 
#define MAX(a,b) (a>b ? a: b)
#endif

#ifndef TRUE 
#define TRUE (1)
#endif
#ifndef FALSE 
#define FALSE (0)
#endif



static inline int _assertGT0(const char* fnamen,int line,int status)
	{
	if(status<=0)
		{
		fprintf(stderr,"%s:%d : Error result status %d<=0\n",fnamen,line,status);
		exit(EXIT_FAILURE);
		}
	return status;
	}
static inline int _assertGE0(const char* fnamen,int line,int status)
	{
	if(status<0)
		{
		fprintf(stderr,"%s:%d : Error result %d<0\n",fnamen,line,status);
		exit(EXIT_FAILURE);
		}
	return status;
	}
#define assertGT0(status) _assertGT0(__FILE__,__LINE__,status)
#define assertGE0(status) _assertGE0(__FILE__,__LINE__,status)
#define VERIFY(status) assertGE0(status)
#endif

