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

#include "utils.h"

char* safeStrDup(const char* src)
	{
	char* p=strdup(src);
	if(p==NULL) DIE_FAILURE("strdup failed");
	return p;
	}

char* safeStrNDup(const char* s,size_t n)
	{
	size_t len = strlen(s);
  	char *result;
	if(len<n) n=len;	
	result=(char *)safeMalloc (n + 1);
  	result[n] = '\0';
  	return memcpy ((void*)result,(void*)s, n);
	}

size_t strsplit(char* src,char delim,char** offsets,size_t max_tokens)
	{
	size_t i=0,found=1;

	if(max_tokens<1) DIE_FAILURE("max_limit=0");
	offsets[0]=src;
	while(src[i]!=0 && found+1<=max_tokens)	
		{
		if(src[i]==delim)
			{
			offsets[found++]=&src[i+1];
			src[i]=0;
			}
		++i;
		}
	return found;
	}

char* gzReadLine(gzFile in,size_t *str_size)
	{
	int c;
	char * p=NULL;
	size_t len=0UL;
	size_t capacity=BUFSIZ;
	if(str_size!=NULL) *str_size=0UL;
	if(gzeof(in)) return NULL;
	while((c=gzgetc(in))!=-1 && c!='\n')
		{
		if(p==NULL || len+2>=capacity)
			{
			capacity+=BUFSIZ;
			p=(char*)safeRealloc(p,sizeof(char)*capacity);
			}
		p[len++]=c;
		}
	if(p==NULL) return NULL;
	p[len]=0;
	if(str_size!=NULL) *str_size=len;
	return p;
	}

gzFile safeGZOpen(const char *path, const char *mode)
	{
	gzFile f = gzopen(path,mode);
	if(f==NULL) 
		{
		DIE_FAILURE("Cannot open \"%s\" :%s.\n",
			path,
			strerror(errno));
		}
	return f;
	}


void* _safeMalloc(const char* fname,int line, size_t size)
	{
	void * ptr=malloc(size);
	if(ptr==NULL)
		{
		DIE_FAILURE("Out of memory \"%s\":%d  N=%ld.\n",
			fname,line,size
			);
		}
	return ptr;
	}

void* _safeCalloc(const char* fname,int line, size_t n,size_t size)
	{
	void * ptr=calloc(n,size);
	if(ptr==NULL)
		{
		DIE_FAILURE("Out of memory \"%s\":%d N=%ld*%ld.\n",
			fname,line,n,size
			);
		}
	return ptr;
	}

void* _safeRealloc(const char* fname,int line, void* ptr,size_t size)
	{
	ptr=realloc(ptr,size);
	if(ptr==NULL)
		{
		DIE_FAILURE("Out of memory \"%s\":%d N=%ld\n",
			fname,line,size
			);
		}
	return ptr;
	}
