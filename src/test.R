source("ibddb.R")

# OPEN IBD Database
ibd <- ibd.open("../test/test.h5");

print(paste("Number of markers:", ibd.num.markers(ibd),"\n"))
print(paste("Number of chromosomes:", ibd.num.chromosomes(ibd),"\n"))
print(paste("Number of individuals:", ibd.num.individuals(ibd),"\n"))
print(paste("Number of pairs:", ibd.num.pairs(ibd),"\n"))

## print some chromosomes
for(i in seq(0,ibd.num.chromosomes(ibd)-1) )
	{
	print(ibd.chromosome(ibd,i));
	}

## print some markers
for(i in seq(0,min(10,ibd.num.markers(ibd)-1)) )
	{
	print(ibd.marker(ibd,i));
	}



## print some individuals
for(i in seq(0,min(10,ibd.num.individuals(ibd)-1)) )
	{
	print(ibd.individual(ibd,i));
	}
## print some pairs
for(i in seq(0,min(10,ibd.num.pairs(ibd)-1)) )
	{
	print(ibd.pair(ibd,i));
	}
## print some data

###  loop over markers
for(y in seq(0,min(50,ibd.num.markers(ibd)-1)) )
	{
	## get y-th marker
	marker <- ibd.marker(ibd,y)
	## get chromosome info for this marker
	chrom <- ibd.chromosome(ibd,marker$tid)
	
	### loop over pairs
	for(x in seq(0,min(10,ibd.num.pairs(ibd)-1)) )
		{
		## get x-th pair
		pair <- ibd.pair(ibd,x)
		## get first individual of this pair
		indi1 <- ibd.individual(ibd,pair$indi1idx)
		## get second individual of this pair
		indi2 <- ibd.individual(ibd,pair$indi2idx)
		
		cat(paste(chrom$name,marker$name,marker$position,
			paste(indi1$family,indi1$name,sep=":"),
			paste(indi2$family,indi2$name,sep=":"),
			sep="\t"))
		cat("\t")
		# print IBD 0
		cat(ibd.ibd0(ibd,y,x))
		cat("\n")
		}
	}



# releases resources associated to IBD database
ibd.close(ibd);
