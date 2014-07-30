
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

#' Opens an IBD-DB HDF5 file
#'
#' @param filename the HDF5 file
#' @keywords HDF5 file
#' @return ibd context
#' @examples
#' ibd.open('file.h5')
ibd.open<-function(filename)
	{
	private_load_ibd_libraries();
	.Call("RIbdDbOpen",filename);
	}


#' Close an ibd context and realease the associated resources
#'
#' @param ibd the IBD context
#' @keywords ibd
#' @examples
#' ibd.close(ibd)
ibd.close<-function(ibd)
	{
	private_load_ibd_libraries();
	.Call("RIbdDbClose",ibd)
	}

#' returns the number of markers in the IBD context
#'
#' @param ibd the IBD context
#' @keywords ibd
#' @return the number of markers
ibd.num.markers<-function(ibd)
	{
	.Call("RIbdDbCountMarkers",ibd)
	}

#' returns the number of chromosomes in the IBD context
#'
#' @param ibd the IBD context
#' @keywords ibd
#' @return the number of chromosomes
ibd.num.chromosomes<-function(ibd)
	{
	.Call("RIbdDbCountChromosomes",ibd)
	}

#' returns the number of pairs in the IBD context
#'
#' @param ibd the IBD context
#' @keywords ibd
#' @return the number of pairs
ibd.num.pairs<-function(ibd)
	{
	.Call("RIbdDbCountPairs",ibd)
	}

#' returns the number of individuals in the IBD context
#'
#' @param ibd the IBD context
#' @keywords ibd
#' @return the number of individuals
ibd.num.individuals<-function(ibd)
	{
	.Call("RIbdDbCountIndividuals",ibd)
	}

#' returns the index-th chromosome in the IBD context. Chromosome are ordered on the original sequence dictionary
#' A chromosome is a tuple (name,tid,length)
#' @param ibd the IBD context
#' @param index 0-based index
#' @keywords ibd
#' @return the index-th chromosome	
ibd.chromosome<-function(ibd,index)
	{
	.Call("RIbdDbGetChromosomeAt",ibd,index);
	}
	
#' returns the index-th marker in the IBD context. Markers are ordered on the tid,position
#' A marker is a tuple (name,chromosome-index,position,self-index)
#' @param ibd the IBD context
#' @param index 0-based index
#' @keywords ibd
#' @return the index-th marker	
ibd.marker<-function(ibd,index)
	{
	.Call("RIbdDbGetMarkerAt",ibd,index);
	}

#' returns the index-th pair of individual in the IBD context
#' A pair is a tuple (individual1-index,individual2-index,self-index)
#' @param ibd the IBD context
#' @param index 0-based index
#' @keywords ibd
#' @return the index-th pair	
ibd.pair<-function(ibd,index)
	{
	.Call("RIbdDbGetPairAt",ibd,index);
	}
	
	
#' returns the index-th pair of individual in the IBD context
#' A pair is a tuple (individual1-index,individual2-index,self-index)
#' @param ibd the IBD context
#' @param index 0-based index
#' @keywords ibd
#' @return the index-th pair
ibd.individual<-function(ibd,index)
	{
	.Call("RIbdDbGetIndividualAt",ibd,index);
	}


#' returns the IBD-0 status for the given marker-index,pair-index
#' @param ibd the IBD context
#' @param marker-index 0-based index of the marker
#' @param pair-index 0-based index of the pair
#' @keywords ibd 
#' @return the IBD-0 (or null if undefined)
ibd.ibd0<-function(ibd,marker_y,pair_x)
	{
	if(!is.integer(marker_y)) stop("not an integer");
	if(!is.integer(pair_x)) stop("not an integer");
	
	.Call("RIbdDbGetIBD0",ibd,marker_y,pair_x);
	}

#' returns the IBD-1 status for the given marker-index,pair-index
#' @param ibd the IBD context
#' @param marker-index 0-based index of the marker
#' @param pair-index 0-based index of the pair
#' @keywords ibd 
#' @return the IBD-1 (or null if undefined)
ibd.ibd1<-function(ibd,marker_y,pair_x)
	{
	if(!is.integer(marker_y)) stop("not an integer");
	if(!is.integer(pair_x)) stop("not an integer");
	.Call("RIbdDbGetIBD1",ibd,marker_y,pair_x);
	}


#' returns the IBD-2 status for the given marker-index,pair-index
#' @param ibd the IBD context
#' @param marker-index 0-based index of the marker
#' @param pair-index 0-based index of the pair
#' @keywords ibd 
#' @return the IBD-2 (or null if undefined)
ibd.ibd2<-function(ibd,marker_y,pair_x)
	{
	if(!is.integer(marker_y)) stop("not an integer");
	if(!is.integer(pair_x)) stop("not an integer");
	.Call("RIbdDbGetIBD2",ibd,marker_y,pair_x);
	}
