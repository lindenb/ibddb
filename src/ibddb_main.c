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
#include "ibddb.h"
#define SUBPROG(name) extern int main_##name(int argc,char** argv)
SUBPROG(build);
SUBPROG(ibd);
SUBPROG(dict);
SUBPROG(markers);
SUBPROG(ped);
SUBPROG(pairs);
SUBPROG(pedigree);

static void main_usage(int argc,char** argv)
	{
	USAGE_PREAMBLE;
	fprintf(stderr,"Usage:\n\t%s [subprogram] (options)\n\n",argv[0]);
	fputs("Sub-Programs:\n\n",stderr);
	fputs(" build   : build IBD database.\n",stderr);
	fputs(" ibd     : query ibds.\n",stderr);
	fputs(" dict    : dump reference dictionary.\n",stderr);
	fputs(" ped     : dump pedigree.\n",stderr);
	fputs(" markers : dump markers.\n",stderr);
	fputs(" pairs   : dump pairs of individuals.\n",stderr);
	fputs("\n\n",stderr);
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
		main_usage(argc,argv);
		return EXIT_SUCCESS;
		}
	
	DEBUG("%s: exiting with status=%d",argv[0],status);
	return status;
	}
