# R or Splus code to build a site by spp matrix from an input species data set;
# that is in 'list' form, as described below;

# List form:  each row of the input data set contains ;
#  STRM_ID= unique sample or stream id;
#  SPECCODE = species code;
#  ABUND=abundance of that species in that sample;

# output of this code is a sample by species matrix, ;
# with samples=rows, species=columns;

# Assume "bugref" is the input data frame, with columns defined as above;



########## Begin code;

sample<-unique(bugref$STRM_ID); #vector of sample names;
species<-unique(bugref$SPECCODE); # vector of species names;

nsp<-length(species); #number of species;
nsamp<-length(sample); #number of samples;

#create matrix of 0's for sample by spp;
#Use sample and species names as row and column headers, respectively;

sitspp<-matrix(rep(0,nsp*nsamp),nrow=nsamp,
    dimnames=list(as.character(sample),species));

#loop over samples and fill in the species abundance row vector for each sample;
for(isit in 1:nsamp) {;
onesamp<-bugref[bugref$STRM.ID==sample[isit],]; #extract data containing all species observed at the current site;
matchvec<-match(onesamp$SPECCODE,species,nomatch=0); 
sitspp[isit,matchvec]<-onesamp$ABUND; 
  };

#matrix sitspp is complete;


