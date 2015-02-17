# R code to build a RIVPACS-type model. ;
# Includes computing dissimilarity matrix, clustering, cluster pruning, and discriminant function analysis (DFA) ;
# options are available for stepwise and all subsets DFA;

#Version 3 , 6/25/07 -- Adds BC index as alternative for model evaluation;

#Version 4, 02/35/10 -- Reorganized, along with model.predict.r, to produce a final model that is more portable;
     #Also includes a function to compare observed versus predicted for all taxa at a single site;

#Version 4.1, )7/15/10 -- More flexible input andd more info in  output for model predictions.
#                       -- Also, option to remove rare taxa before clustering;

# John Van Sickle, US Environmental Protection Agency;

##################################################;

#load required packages;
library(gtools); library(MASS);
library(cluster); library(Hmisc);

############;
# STEP 1 -- INITIAL SETUP -- Input and organize the bug and predictor data;
# The exact code for this step is very particular for each data set, ;
# but similar operations will need to be done on every data set;
# Data setup below is an example for an Oregon DEQ data set which contains both reference and nonreference sites;

###############################;
# Input data are predictor data (all sites) and a (site x taxa) matrix of abundance for all bugs at all sites.;
# The bug matrix is the result of fixed-count subsampling and matrify programs;
# Assumes that predictor data file includes a column to identify the calibration, validation and test sites;

# Step 1a - Read and organize predictor data;

# Input the predictor data, tab delimited. Use the sample/site ID as the row name;
predall<-read.table("enviros27jul04.txt",row.names="Sample",header=T,sep="\t");
#fix error in predall data. See Huff email, 8/31/04;
predall['H111013','Eco1']<-0;
predall['H111013','Mareco']<-1;

head(predall); #look at 1st 5 rows, all columns;
dim(predall); # number of rows and columns;

#specify candidate predictors;
# First, put the (column) names of candidate predictors in a vector;
candvar<-c( "DAYNUM","X_coord","Y_coord","WA_ha_lg","M_Slp_sqt","Mareco",
           "Wcrdeco","Pcp_ann","Tmx_ann","Loolith",
           "Litho1","Litho4","El_m_sqt","Strpw_lg");

#Look at predictor histograms, over all samples, to see if any transformations are desirable;
#first, produce an empty graph window;
windows();
#On window menu bar, click "History", then "recording";
sapply(candvar,function(nam)hist(predall[,nam],main=nam));
#now click on graph window and use PgUp, PgDn to scroll through plots;

## for this example, will not do additional transformations are needed;

###########;
## Step 1b - Input the assemblage data (bug data), as a site-by-taxa matrix;

bugall<-read.table("bugs_matr_27jul.txt",row.names="SAMPLE",header=T,sep="\t");

## Step 1c - Align bug and predictor data, by site/sample;
#check sample(row) alignment of bug and predictor data;
row.names(bugall)==row.names(predall);
#samples are not aligned. Fix by aligning bugs data to predictor data, since latter is sorted by sample type;
bugall<-bugall[row.names(predall),];
#check alignment again -- alignment OK;
row.names(bugall)==row.names(predall);

### Step 1d - Create subsets of calibration and validation data;

#First, create a Presence/absence (1/0) matrix (site by taxa) for the bugs;
bugall.pa<-bugall;
bugall.pa[bugall.pa>0]<-1;

# Extract subsets of bug and predictor data for the calibration ("C") and validation ("V") sites;
predcal<-predall[predall[,'Mdl704']=='C',];  #predictor data -calibration sites;
pred.vld<-predall[substr(as.character(predall[,'Mdl704']),1,1)=='V',];  #predictor data - validation sites;

bugcal<-bugall[predall[,'Mdl704']=='C',]; #Bug Abundance matrix, calibration sites;
bugcal.pa<-bugall.pa[predall[,'Mdl704']=='C',]; #Bug presence/absence matrix, calibration sites;
bug.vld.pa<-bugall.pa[substr(as.character(predall[,'Mdl704']),1,1)=='V',]; #Bug presence/absence matrix, validation sites;

 #Data sets complete and aligned;
 ########################################;
 #STEP 2 -- DISSIMILARITIES AND CLUSTER ANALYSIS;

##########
# 2A. Prepare the bug matrix to be used in calculating site dissimilarities for clustering;

#Next block of script is an option to remove rare taxa before clustering;

#First, get a vector containing the proportion of calibration sites at which each taxon occurs;
psite.occ<-colSums(bugcal.pa)/dim(bugcal.pa)[[1]];

#Now get a vector of names of taxa occuring at greater than 5% of calibration sites;
#Will use only these taxa in clustering;
nonrare.taxa<-names(psite.occ)[psite.occ>0.05];

#Now subset the site-by-taxa matrix to create a new one containing only the nonrare taxa;
#For example, to cluster based on the abundance matrix;
bugcal.nonrare<- bugcal[,nonrare.taxa];
#Alternatively, do the following if you plan to cluster based on P/A data;
bugcal.pa.nonrare<- bugcal.pa[,nonrare.taxa];
#Now use the nonrare version, instead of bugcal or bugcal.pa, in subsequent dissimilarity calculations below;

##########################;
# Bug matrix is ready. Proceed to dissimilarity calculations;

#################################;
# The code below calculates Jaccard, Sorenson or Bray-Curtis dissimilarities;
#Other sources for dissimilarities include vegdist() in vegan package (BC, Jaccard, and others);
#Also, gdist() in package mvpart, and dist() in base R;

# Here, I use the generalized outer product function dapply();
# and choose the desired dissimilarity measure as a called function;
# dapply() output is an (n(n-1)/2) length vector storing ;
# the lower triangle of the site dissimilarity matrix in column major order;
# can be input directly to R clustering functions;

#source dapply;
source("dapply.r");

#execute one of the following dissimilarity blocks;

#####################;
# Option 1 -- Jaccard dissimilarity for P/A data;
# function computes jaccard P/A dissimilarity between one site pair,  siti and sitj;
#input can be P/A or abundance data;

    jaccfun<-function(siti,sitj) {
              shared<-sum((siti>0)&(sitj>0));
              uniquei<-sum((siti>0)&(sitj==0));
              uniquej<-sum((siti==0)&(sitj>0));
              1-(shared/(shared+uniquei+uniquej)); #return Jaccard dissimilarity;
                  } #end of function;
                 
#compute Jaccard dissimilarities among all calibration sites;
dissim<-dapply(bugcal.pa.nonrare,1,bugcal.pa.nonrare,1,jaccfun);
#Proceed to clustering;
#####################;
# Option 2 -- Sorenson dissimilarity for P/A data;
# function computes Sorenson P/A dissimilarity between one site pair,  siti and sitj;
#input can be P/A data or abundance data;

  sornfun<-function(siti,sitj) {
              shared<-sum((siti>0)&(sitj>0));
              uniquei<-sum((siti>0)&(sitj==0));
              uniquej<-sum((siti==0)&(sitj>0));
              1-(2*shared/(2*shared+uniquei+uniquej)); #return Sorenson dissimilarity;
                  } #end of function;
                 
#compute Sorenson dissimilarities among all calibration sites;
dissim<-dapply(bugcal.pa.nonrare,1,bugcal.pa.nonrare,1,sornfun);

#Proceed to clustering;
#########################################;
#Option 3 -- Bray-Curtis (Sorenson) dissimilarity for abundance data;
# in this example, use untransformed relative abundance;

#first compute site by spp matrix of relative abundance;
totabun<-apply(bugcal.nonrare,1,sum); #vector of total abundance, each site;
rel.abun<-sweep(bugcal.nonrare,1,totabun,FUN="/"); #relative abundance matrix;

#function below computes Bray-Curtis dissim between 2 sites;
#siti, sitj are vectors of relative abundance (or absolute abundance) for 2 sites;
#if zero abundance at both sites, then dissimilarity=0;

   bcfun<-function(siti,sitj) {
          bcnum<-sum(abs(siti-sitj));
          bcdenom<-sum(siti+sitj);
          ifelse(bcdenom>0, (bcnum/bcdenom),0); #return BC dissimilarity;
                } #end of function;

#compute Bray-Curtis dissimilarity;
dissim<-dapply(rel.abun,1,rel.abun,1,bcfun);
#Proceed to clustering;
#########################;
####################################;
# Clustering of calibration sites;
# Use flexible-Beta method, with Beta=-0.6;
#See R documentation on agnes() in "cluster" package;
# When using agnes() with Flexible Beta strategy, set Beta=(1-2*Alpha) in Lance-Williams formula;
# A single value for par.method value specifies alpha, so alpha=0.8 gives Beta=-0.6;

clus1<-agnes(x=dissim,diss=T,method="flexible", par.method=0.8,keep.diss=F,keep.data=F);

## Various plots of cluster outcome. Leaf labels are row numbers of dissim matrix;
#that is, the order of sites in the calibration data set;
#plot the dendrogram;
plot(clus1,which.plots=2,labels=row.names(bugcal.pa),cex=.4);

#can also define a more general dendrogram structure. See 'dendrogram' help;
dend1<-as.dendrogram(as.hclust(clus1));
plot(dend1,horiz=T,cex=.2); #Plot dendrogram with horizontal branches;

#######################;
#Pruning the dendrogram to create a small number of groups;
# Level pruning can be done by specifying the number of groups (k parameter);
# Also can prune at a specified height. See cutree help;
# Result is a vector of site group assignments;
# Can repeat this process to generate several candidate groupings from a single dendrogram;

# First example is "level" pruning to afixed number of groups;
grps.6<-cutree(clus1,k=6); #vector of group assignments is in the order of sites in the clustered data;
table(grps.6); #count number of sites in each group;
cbind(row.names(predcal),grps.6); #list calibration sites and their group assignments;
#candidate site groups complete;
#############;
#alternative is non-level pruning.
# Functions below interactively pick out desired clusters;
#First, plot a horizontal dendrogram;
plot(clus1,which.plots=2,labels=row.names(bugcal.pa),cex=.4);
#Next, drag screen corner to spread out the plot horizontally;

#Then execute the following. Move the cursor crosshair to be near a high-level dendrogram branch;
# Then left-click, and function draws a box around all the entire branch, puuting all sites into a group;
# Move horizontally across dendrogram, selecting branches until all sites are in mututally exclusive groups;
# Right-click to stop to process;
ccc<-identify(as.hclust(clus1)); #interactive ID of clusters on dendrogram;
# Cluster selection is complete;

#convert result to a group ID vector;
grp.id<-rep(NA,dim(bugcal)[[1]]);
for(k in 1:length(ccc)) {grp.id[ccc[[k]]]<-k};
table(grp.id); #count number of sites in each group;
cbind(row.names(predcal),grp.id); #list calibration sites and their group assignments;
# The vector "grp.id" now defines the calibration site groups - Use as alternative to grps.6 or similar;
#Cluster analysis is complete;
###########;
#as a check, can plot the geographic  locations of the site clusters;
# Can also use the "maps" package to draw a state or regional boundary;
#This example plots 6 groups;
plot(predcal$X_coord[grps.6==1], predcal$Y_coord[grps.6==1],col='black',
    type='p',xlim=c(-125,-118),ylim=c(42,46), xlab="Longitude",ylab="Latitude");
points(predcal$X_coord[grps.6==2], predcal$Y_coord[grps.6==2],col='red')
points(predcal$X_coord[grps.6==3], predcal$Y_coord[grps.6==3],col='green')
points(predcal$X_coord[grps.6==4], predcal$Y_coord[grps.6==4],col='orange')
points(predcal$X_coord[grps.6==5], predcal$Y_coord[grps.6==5],col='blue')
points(predcal$X_coord[grps.6==5], predcal$Y_coord[grps.6==5],col='purple')
#####################################;

#STEP 3 -- DISCRIMINANT FUNCTION ANALYSIS (DFA);

# Note - Instead of DFA, consider using classification tree model (R packages "tree" or "rpart");
#  or a random forest model (R package "randomForest"):

#Below, I have options for stepwise DFA and also for all-subsets DFA;

########################################;

 ###########################;
 # Option 1 -- Stepwise DFA;
 
#Load the step.dfa function;
source("dfa.step.r");

#then execute one of the two below. First is forward stepwise, second is backwards stepwise;
# See dfa.step documentation. Step.res contains the final chosen model as an lda() object;

step.res<-dfa.step(predcal=predcal,candvar=candvar,grps=grps.6,
                                P.stay.lim=.051,P.enter.lim=.05);
#or could do backwards stepwise;
step.res<-dfa.step(predcal=predcal,candvar=candvar,method='backward',grps=grps.6,
                                 P.stay.lim=.051,P.enter.lim=.05);

################################;

#Option 2 -- All subsets DFA; 

# Feasible for up to about 15 candidate predictors;
# User specifies a small number of best models for selected model orders;
# Wilks lambda, classification accuracy, and statistics of O/E are reported for each best model;
# If user supplies an independent set of validation data (bug data and predictor data), then;
# O/E statistics also computed for validation set;
# In addition (version 4), the CV confusion matrices are reported for each of the best models;
#Also, the null model statistics are available;

#specify a vector describing how many models of each order to keep; 
# The following example is for a case with 14 candidate predictors;
# For example, keep 5 models each for orders 1,2, ...13,
# and also keep the single (saturated) model of order 14,
 nkeep<-c(rep(5,13),1);

#Load the all subsets DFA function;
source("dfa.allsub.v4.r");

#execute the following block of code. dfa.allsub.v3() is surrounded;
#by code that records and prints the execution time;
#Execution may take several minutes;

#In example below, Pc is set to 0.5.
#To retain all taxa in O/E, set Pc to a very small number, like 1.e-14.

start.time=proc.time();
dfm.best<-dfa.allsub.v3(bug.cal=bugcal.pa,bug.vld=bug.vld.pa,pred.cal=predcal,pred.vld=pred.vld,
                    grps=grps.6,candvar=candvar,numkeep=nkeep,Pc=0.5);
elaps<-proc.time()-start.time;
print(c("elapsed time = ",elaps));
                                                   
dfm.best.6grp.0<-dfm.best; #Can store result under a new name, indicating the Pc value used;
dfm.best.6grp.5<-dfm.best; # Can store result of a second run, which had Pc=0.5;
#################;
# Various ideas for exploring the subset of "best" DFA models produced by dfa.allsub;
# dfa.allsub yields a list containing the statistics of O/E: a) predicted by the null model (null.stats),
 #  and b) predicted by several "best" predictive models (subset.stats);   
  
 #A) Performance of the null model;
 dfm.best$null.stats; 
 
# B)Performance of subsets of best predictive models;
bestmods<-dfm.best$subset.stats;
head(bestmods);
#look at model #11;
bestmods[11,];
#look at all best models, sorted by SD(O/E) at CAL sites;
format(bestmods[order(bestmods$SDOE.cal),],digits=3);
#or sort them by overall RMSE of O/E on validation data ;
format(bestmods[order(bestmods$RMSE.vld),],digits=3);
#or else sort them by DF classification accuracy;
format(bestmods[order(bestmods$cls.crct.cv),],digits=3);

# C) Look at crossvalidated error matrix of a selected best model;
dfm.best$CV.error.matrices[11]; #For model #11;

# D) plot a measure of model performance against model size (ie, model order);
    #For example, plot RMSE(O/E) against model order separately for calibration and validation sites;
 par(mfrow=c(1,1));
plot(bestmods$order,bestmods$RMSE.cal,ylim=c(0.15,0.20),type='p',pch='C', col='blue',
      cex=.7,xlab='Model order',ylab='RMSE(O/E)');
points(bestmods$order,bestmods$RMSE.vld,pch='V',cex=.7,col='red');
#put null model RMSE as a baseline, separate for Calibration and validation sites.;
abline(dfm.best$null.stats["RMSE.cal"],0,lty=1);
abline(dfm.best$null.stats["RMSE.vld"],0,lty=2);

# Interactive choice of a model from the above plot;
# identify the C and V points for one model on the plot;
# the Cal point is marked with a solid box, Vld with a solid triangle, and the model is printed;
cc<-identify(bestmods$order,bestmods$RMSE.vld,n=1,plot=F);
points(bestmods$order[cc],bestmods$RMSE.vld[cc],pch=17,cex=1);
points(bestmods$order[cc],bestmods$RMSE.cal[cc],pch=15,cex=1);
print(bestmods$model[cc]); 

#following lines put a title and legend on the plot;
legend(locator(1),legend=c('Calibration sites','Validation sites'),pch=c('C','V'));
title(main=list('ORDEQ models: RMSE(O/E) from 5 best models of each model order',cex=.9));

#Can also experiment with similar plots for BC statistics. "Better" models will have;
# smaller BC90;

###;
# E) Plot the two classification accuracy measures against model order;
 # DFM overfitting starts occurring where the CV accuracy flattens out;
 plot(bestmods$order,bestmods$cls.crct.resub,ylim=c(20,60),type='p',pch='R',
      cex=.5,xlab='Model order',ylab='Percent correct');
points(bestmods$order,bestmods$cls.crct.cv,pch='C',cex=.5);
legend(locator(1),legend=c('Resubstitution','Crossvalidation'),pch=c('R','C'));
title(main=list('ORDEQ models: Classification accuracy of 5 best models of each model order',cex=.9));

#F) PREDICTOR IMPORTANCE. Calculate the percentage of best models that include;
# each of the predictors. Percentage is not weighted by model quality;
 round((100*table(unlist(strsplit(bestmods$model," ")))/dim(bestmods)[[1]]),1);

#G) scatteplot matrix of model size and performance on validation and calibration sites;
pairs(as.matrix(bestmods[,c('order','RMSE.cal','RMSE.vld')]));

##########;
##  Step 3.5 - Based on the stepwise or all-subsets DFA, or on other evaluations,;
## declare your choice of the predictors that will be used in the final model;
# Define the vector "preds.final" to contain the names of the final, chosen predictor variables;

#Option 1 - Choose final predictors from the output of all-subsets run;
  #Option 1A - Choose the predictors from the interactive identification;
    preds.final<-unlist(strsplit(bestmods[cc,'model']," "));
  #Option 1B choose predictors from a specific row in bestmods data frame;
    preds.final<-unlist(strsplit(bestmods[13,'model']," "));

#Option 2 - Choose the final model from stepwise DFA run, which is stored in step.res;
preds.final<-attr(terms(step.res),"term.labels");

# OPTION 3 -- Directly specify the names of chosen predictors;
preds.final<-c("DAYNUM","X_coord", "El_m_sqt");

#############################################################;

#STEP 4 - Finalize calculations of final, chosen predictive model;
# To specify the entire final model, you need to store/export 5 things:
# 4.1) The site-by-taxa matrix of observed presence/absence at calibration sites (bugcal.pa, already available);
# 4.2) Specify the vector of final group membership assignments at calibration sites(grps.final);
      grps.final<-grps.6;
# 4.3) Specify the final selected predictor variables (preds.final) - already done in step 3.5;
# 4.4) Calculate the matrix(grpmns) of mean values of the preds.final variables for the final calibration site groups ;
     datmat<-as.matrix(predcal[,preds.final])
     grpmns<-apply(datmat,2,function(x)tapply(x,grps.final,mean));
# 4.5) Calculate the inverse pooled covariance matrix(covpinv) of the preds.final variables at the calibration sites;
     #First, calculate a list of covariance matrices for each group;
      covlist<-lapply(split.data.frame(datmat,grps.final),cov);
      #pooled cov matrix is weighted average of group matrices, weighted by group size. Johnson & Wichern, 11-64;
      grpsiz<-table(grps.final);
      ngrps<-length(grpsiz);
      npreds<-length(preds.final);
     #zero out an initial matrix for pooled covariance;
     covpool<-matrix(rep(0,npreds*npreds),nrow=npreds,dimnames=dimnames(covlist[[1]]));
     #weighted sum of covariance matrices;
     for(i in 1:ngrps){covpool<-covpool+(grpsiz[i]-1)*covlist[[i]]};
     covpool<-covpool/(sum(grpsiz)-ngrps);#renormalize;
     covpinv<-solve(covpool); #inverse of pooled cov matrix;

#################################;

#Step 5 - Further checks on performance of the final, chosen model;

#Option 5.1 - Make predictions of E and O/E for calibration (reference) sites. Examine O/E statistics and plots;
  # To do this, run the model.predict function, using the calibration data as the 'new' data;
  # See Step 7 below, for more info on making predictions, and also ;
  # see internal documentation of model.predict.v4.1;
  source("model.predict.v4.1.r");
  OE.assess.cal<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,prednew=predcal,bugnew=bugcal.pa,Pc=0.5);
# function prints out mean and SD of O/E, in this case, for the calibration sites;
# If an all-subsets model was selected, then these stats
# should match the values given in dfm.best$bestmods, for the selected model;

#look at other prediction results, for calibration sites;
  names(OE.assess.cal);   #names of 3 components of the prediction results list;
  head(OE.assess.cal$OE.scores); #data frame of O/E scores, 1st 5 rows;
  head(OE.assess.cal$Capture.Probs); #predicted capture probabilties, 1st 5 rows;
  head(OE.assess.cal$Group.Occurrence.Probs); #predicted group occurrence probabilties, 1st 5 rows;
#check distribution of Calibration-site O/E scores. Is it Normal?;
#plot a histogram and a Normal q-q plot;
par(mfrow=c(2,1));
hist(OE.assess.cal$OE.scores$OoverE,xlab="O/E");
qqnorm(OE.assess.cal$OE.scores$OoverE);

#scatterplot of O (on y-axis) vs E (on x-axis). See Pineiro et al. Ecol. Modelling 2008, 316-322, for this choice of axes;
  par(mfrow=c(1,1));
 plot(OE.assess.cal$OE.scores[,c('E','O')],xlab='Expected richness',ylab='Observed richness');
  abline(0,1); #add a 1-1 line;
 ###########;

### Option 5.2 - Repeat Step 5.1, but this time predict for validation data.;
##  Check especially for model bias (mean(O/E) differs from 1.0);
   OE.assess.vld<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,prednew=pred.vld,bugnew=bug.vld.pa,Pc=0.5)  ;
   OE.assess.vld$OE.scores;
  #Can repeat predictions for a different Pc cutoff;
  # This example is for a cutoff of 0 (include all taxa). Need to use a very small positive number for Pc parameter;
   OE.assess.vld.alltax<-model.predict.v4(bugcal.pa,grps.final,preds.final, grpmns,covpinv,prednew=pred.vld,bugnew=bug.vld.pa,Pc=1.e-14)  ;
   OE.assess.vld.alltax$OE.scores;

 #################################;
#Step 6 -- Save/export the final model for use by others;

    #Option 6.1 - Save the model components together in a single .Rdata file.;
    # Any R user can load this file, along with model.predict.v4.r, to make predictions from the model;
    # See Step 7 below;
      save(bugcal.pa, predcal, grps.final, preds.final, grpmns, covpinv,file='MyModel.Version1.Rdata');
    #NOTE - Predcal is not needed to define the model, but is included so that users see the required format for predictor data;
 ###############;

    #Option 6.2 - Export the final model pieces as tab-delimited text files,
    #              in the formats needed for uploading to WCEM website;

      # Export the group means matrix, after transposing so that columns are groups and rows are predictor variables;
       grpmns.2<-t(grpmns);
       grpmns.xport<-data.frame("Group"=row.names(grpmns.2),grpmns.2);
       names(grpmns.xport)<-c("Group",row.names(grpmns));
       write.table(grpmns.xport,file="MyModel.version1.groupmeans.txt", sep = "\t",quote=F,row.names=F);

      #Export the inverse covariance matrix, row and column headers are predictor variables;
       #use full decimal precision;
       covp.xport<-data.frame("Variable"=row.names(covpinv),covpinv);
       write.table(covp.xport,file="MyModel.version1.inv_covariance.txt", sep = "\t",quote=F,row.names=F);
      #export the reference-site bug abundance data for the calibration sites;
      # Need to add the group ID variable as the first column after the sample (site) ID;
       bugcal.xport<-data.frame(sample_name=row.names(bugcal),class_code=grps.final,bugcal);
       write.table(bugcal.xport,file="MyModel.version1.reftaxa.txt", sep = "\t",quote=F,row.names=F);
   ##;
      #Also, if possible. WCEM wants you to export files of reference validation bug data and predictor data for internal validation;
     #For this particular data set, need to first extract the abundance matrix for validation sites (see step 1b, way above)
      bug.vld<-bugall[substr(as.character(predall[,'Mdl704']),1,1)=='V',]; #Bug abundance matrix, validation sites;
     #first, remove all samples from the bug and predictor validation data sets that have a missing ;
     # value for >=1 of the final predictor variables;
        prd.vld2<-pred.vld[complete.cases(pred.vld[,preds.final]),];
        bug.vld2<-bug.vld[complete.cases(pred.vld[,preds.final]),];
      # Export the validation bug abundance matrix. Same format as bugcal.xport, except it lacks the class_code column;
        bugvld.xport<-data.frame(sample_name=row.names(bug.vld2),bug.vld2);
        write.table(bugvld.xport,file="MyModel.version1.reftesttaxa.txt", sep = "\t",quote=F,row.names=F);
      #Now, export the predictor data for validation sites;
      #Keep only the final predictor variables, in same order as in the inv cov matrix;
        prd.vld2<-prd.vld2[,preds.final];
        prdvld.xport<-data.frame(sample_name=row.names(prd.vld2),prd.vld2);
        write.table(prdvld.xport,file="MyModel.version1.refhab.txt", sep = "\t",quote=F,row.names=F);
     #finished with WCEM export files;
################################;

# Step 7 - Making predictions for new (test) data.

# first, source the prediction script and also load the desired model;
   source("model.predict.v4.1.r");
   load('MyModel.Version1.Rdata');

   # User must supply a sample-by-taxa matrix of taxa abundance or else presence/absence(coded as 1 or 0), for all new samples;
   # User must also supply a corresponding file of predictor data for those same samples;
   # These 2 files should have similar formats as the original taxa and predictor data sets used to build the model (see step 1 above);
   # Notes on file formats --
   #   A) The sample ID column in both files should be read into R as a row name (see Step 1 examples).
   #   B) Predictor data set -- Must include columns with the same names, units, etc.,
   #        as the model's predictor variables. All other columns will be ignored;
   #        Column order does not matter;
   #        Predictions and calculations of O/E will be made only for those samples that have;
   #        complete data for all model predictors.;
   #   C)  Sample-by-taxa matrix. Can contain abundance or presence/absence (1 or 0). Missing or empty cells now allowed;
   #       Sample ID's (row names) must match those of predictor data.
   #       Any names for new taxa (column names) are acceptable, in any order;
   #       HOWEVER - Only those new-data taxa names that match the names in the
   #            calibration data can be use to calculate observed richness;
   #            All other taxa (columns) in the new-data bug matrix are ignored;
   #        To see a list of the calibration-taxa names, do:
            names(bugcal.pa)[colSums(bugcal.pa)>0];

##########;
# Example predictions: For nonreference sites in the Oregon DEQ data set that are labeled "N_lc" (see Step 1);

pred.test<-predall[as.character(predall[,'Mdl704'])=='N_lc',];  #predictor data - test sites;
bug.test.pa<-bugall.pa[as.character(predall[,'Mdl704'])=='N_lc',]; #Bug presence/absence matrix, test sites;

#Drop all samples/sites that do not not have complete data for the model predictors;
pred.test<-pred.test[complete.cases(pred.test[,preds.final]),];
bug.test.pa<-bug.test.pa[row.names(pred.test),];

#makes predictions for test data;
OE.assess.test<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,prednew=pred.test,bugnew=bug.test.pa,Pc=0.5);

# look at O/E and BC scores of test-data samples;
 OE.assess.test$OE.scores;
# Look at predicted capture probabilties (1st 5 rows) for all calibration taxa, for test data;
 head(OE.assess.test$Capture.Probs);
# Look at predicted group occurrence probabilties, for all test samples;
  OE.assess.test$Group.Occurrence.Probs;
################################################;

## Assessing an individual sample or site;
source("assess.one.sample.4.1.r")
#This function assesses a single site or sample from a new (test) data set to which
# model.predict.v4.1() has already been applied.
# assess.one.sample() compares observed occurrences with the model-predicted probabilities of occurrence for all taxa;

#Input parameters are:
#       case -- A selected site or sample ID, for which a prediction has already been made using model.predict.v4.1(). ;
# result.prd -- Output of model.predict.v4.1() for new samples that include the chosen case;
# bugnew  -- Sample-by-taxa matrix of new samples that was submitted to model.predict.v4.1().
# Pc -- Cutoff for capture probabilties for inclusion of taxa in O/E;

#The function produces a data frame with one row for each taxon, and the following columns:
 # observed presence(1) or absence(0);
 # predicted capture probability;
# Big.diff = "Yes", if there is a big difference (>=0.5 in magnitude) between observed and predicted;
# In.OtoE = "Yes" if the taxon would be included in the O/E calculation for this sample, given the stated value of Pc;

#By default, the function prints out the results data frame with its rows(taxa) sorted by the magnitude of (observed-predicted),
   # as suggested in Van Sickle, J. (2008), JNABS 27:227-235;
#However, see below for other sorting possibilties;

#Example usage:
 site1.result<-assess.one.sample.4.1(case="01020CSR",result.prd=OE.assess.test, bugnew=bug.test.pa, Pc=0.5);
# Alternative display is to sort the taxa by their predicted occurrence probabilities;
 site1.result[order(site1.result$predicted,decreasing=TRUE),];
# Another alternative is to sort alphabetically by taxon name;
  site1.result[order(row.names(site1.result)),];

## End of model build and prediction examples;










