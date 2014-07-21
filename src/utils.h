#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <zlib.h>

/* gz file */
gzFile safeGZOpen(const char *path, const char *mode);
char* gzReadLine(gzFile in,size_t *str_size);

/* string */
size_t strsplit(char* src,char delim,char** offsets,size_t max_tokens);
char* safeStrDup(const char* src);

/** stdlib */
void* _safeMalloc(const char*,int,size_t);
void* _safeCalloc(const char*,int,size_t,size_t);
void* _safeRealloc(const char*,int,void*,size_t);

#define safeMalloc(n) _safeMalloc(__FILE__,__LINE__,n)
#define safeCalloc(n,m) _safeCalloc(__FILE__,__LINE__,n,m)
#define safeRealloc(p,n) _safeRealloc(__FILE__,__LINE__,p,n)


#define WHERENL fprintf(stderr,"[%s]:%d ",__FILE__,__LINE__)
#define DIE_FAILURE(FormatLiteral,...) do { WHERENL; fprintf (stderr,"Exiting: " FormatLiteral "\n", ##__VA_ARGS__); exit(EXIT_FAILURE);} while(0)
#define DEBUG(FormatLiteral, ...)  do { WHERENL; fprintf (stderr,"" FormatLiteral "\n", ##__VA_ARGS__);} while(0)

#ifndef MIN 
#define MIN(a,b) (a<b ? a: b)
#endif
#ifndef MAX 
#define MAX(a,b) (a>b ? a: b)
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
#endif

