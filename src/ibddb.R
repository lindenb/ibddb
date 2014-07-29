
load.dynamic.libraries<-function(libnames)
	{ 
	for(libname in libnames)
		{
		found_file=libname;
        	for(path in unlist(strsplit(Sys.getenv("LD_LIBRARY_PATH"),":",fixed=TRUE))) 
        		{
        	        try_file <- paste0(path,"/",libname);
        	        
        	        if( file.exists(try_file) )
        	                {
        	                
        	                found_file = try_file;
        	                break;
        	                }
        		}
       		write(paste("Loading :", try_file), stderr())
        	dyn.load(found_file);
		}
	
	}

PRIVATE_IBD_LIBRARIES_LOADED=FALSE;

private_load_ibd_libraries<-function()
	{
	if(! PRIVATE_IBD_LIBRARIES_LOADED)
		{
		load.dynamic.libraries(c("libhdf5.so","libibddb.so"))
		}
	PRIVATE_IBD_LIBRARIES_LOADED=TRUE;
	}


ibd.open<-function(filename)
	{
	private_load_ibd_libraries();
	.Call("RIbdDbOpen",filename);
	}

ibd.close<-function(ibd)
	{
	private_load_ibd_libraries();
	.Call("RIbdDbClose",ibd)
	}

ibd.num.markers<-function(ibd)
	{
	.Call("RIbdDbCountMarkers",ibd)
	}
ibd.num.chromosomes<-function(ibd)
	{
	.Call("RIbdDbCountChromosomes",ibd)
	}

ibd.num.pairs<-function(ibd)
	{
	.Call("RIbdDbCountPairs",ibd)
	}

ibd.num.individuals<-function(ibd)
	{
	.Call("RIbdDbCountIndividuals",ibd)
	}
	
ibd.chromosome<-function(ibd,index)
	{
	.Call("RIbdDbGetChromosomeAt",ibd,index);
	}
ibd.marker<-function(ibd,index)
	{
	.Call("RIbdDbGetMarkerAt",ibd,index);
	}
ibd.pair<-function(ibd,index)
	{
	.Call("RIbdDbGetPairAt",ibd,index);
	}
ibd.individual<-function(ibd,index)
	{
	.Call("RIbdDbGetIndividualAt",ibd,index);
	}

ibd.ibd0<-function(ibd,marker_y,pair_x)
	{
	if(!is.integer(marker_y)) stop("not an integer");
	if(!is.integer(pair_x)) stop("not an integer");
	
	.Call("RIbdDbGetIBD0",ibd,marker_y,pair_x);
	}

ibd.ibd1<-function(ibd,marker_y,pair_x)
	{
	if(!is.integer(marker_y)) stop("not an integer");
	if(!is.integer(pair_x)) stop("not an integer");
	.Call("RIbdDbGetIBD1",ibd,marker_y,pair_x);
	}

ibd.ibd2<-function(ibd,marker_y,pair_x)
	{
	if(!is.integer(marker_y)) stop("not an integer");
	if(!is.integer(pair_x)) stop("not an integer");
	.Call("RIbdDbGetIBD2",ibd,marker_y,pair_x);
	}
