#code for harmonizing the taxa columns in test or valdation bug data matrices,
# with the taxa columns in a matrix for reference sites;

# assume that bugref and bugtst are the two matrices. Need to remove all taxa from bugtst that;
# are not in bugref. Also need to add zero columns to bugtst for taxa that are in bugref, but not in bugtst;

#create a new empty matrix for test data having the same taxa (columns) as bugref;
ref.tax<-names(bugref); #taxa names from reference data;
tst.site<-row.names(bugtst); #row (site) names from test data;
new.tst<-data.frame(matrix(rep(0,length(ref.tax)*length(tst.site)),nrow=length(tst.site),
         dimnames=list(tst.site,ref.tax)));

#loop over ref taxa and move nozero data from the old to the new test matrix;
for(curtax in ref.tax){
  if(match(curtax,names(bugtst),nomatch=0)>0) new.tst[,curtax]<-bugtst[,curtax]};

#matrix new.tst should contain the original test bug data, rearranged in columns
# to match the reference bug matrix. The new.tst matrix does not contain any taxa;
# that were found at test sites, but not at reference sites; 