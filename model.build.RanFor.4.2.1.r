# R code to build a RIVPACS-type model. ;
#This version uses Random Forest (RF), instead of discriminant function analysis, to predict group membership;
# Otherwise, this script is the same as model.build.v4.r;

# J.Van Sickle, 02/25/10;

#Version 4.2.1 (10/6/11) adds calculation of replicate sampling SD for a new model;
##################################################;

#load required packages;

library(cluster); library(Hmisc);
library(randomForest);

############;
# STEP 1 -- INITIAL SETUP -- Input and organize the bug and predictor data;
# The exact code for this step is very particular for each data set, ;
# but similar operations will need to be done on every data set;
# Data setup below is one example that was used for building Oregon DEQ models.

###############################;
# Input data are predictor data (all sites) and a (site x taxa) matrix of abundance for all bugs at all sites.;
# The bug matrix is the result of fixed-count subsampling and matrify programs;
# Assumes that predictor data file includes a column to ID the calibration, validation and test sites;

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
#Now can use the nonrare version, instead of bugcal or bugcal.pa, in subsequent dissimilarity calculations below;

#################################;
# compute dissimilarity matrix for calibration site bugs;
# at this stage, have not excluded any rare taxa;
# The code below cacluates Jaccard, Sorenson and Bray-Curtis dissimilarities;
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
#For this example, use Sorenson on P/A;
#####################;
# Option 1 -- Jaccard dissimilarity;
# function computes jaccard P/A dissimilarity between one site pair,  siti and sitj;
#input can be P/A or abundance data;

    jaccfun<-function(siti,sitj) {
              shared<-sum((siti>0)&(sitj>0));
              uniquei<-sum((siti>0)&(sitj==0));
              uniquej<-sum((siti==0)&(sitj>0));
              1-(shared/(shared+uniquei+uniquej)); #return Jaccard dissimilarity;
                  } #end of function;
                 
#compute Jaccard dissimilarities among all calibration sites;
dissim.jac<-dapply(bugcal.pa,1,bugcal.pa,1,jaccfun);
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
dissim<-dapply(bugcal.pa,1,bugcal.pa,1,sornfun);

#Proceed to clustering;
#########################################;
#Option 3 -- Bray-Curtis (Sorenson) dissimilarity for abundance data;
# in this example, use untransformed relative abundance;

#first compute site by spp matrix of relative abundance;
totabun<-apply(bugcal,1,sum); #vector of total abundance, each site;
rel.abun<-sweep(bugcal,1,totabun,FUN="/"); #relative abundance matrix;

#function below computes Bray-Curtis dissim within dapply();
# Instead, could use gdist() in mvpart package, to do Bray-Curtis;
#siti, sitj are vectors of abundances for 2 sites;
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
#A single value for par.method value specifies alpha, so alpha=0.8 gives Beta=-0.6;

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
# level pruning can be done by specifying the number of groups (k parameter);
#Also can prune at a specified height. See cutree help;
#result is a vector of site group assignments;
#can repeat this process to generate several candidate groupings from a single dendrogram;

grps.6<-cutree(clus1,k=6); #vector of group assignments is in the order of sites in the clustered data;
table(grps.6); #count number of sites in each group;
cbind(row.names(predcal),grps.6); #list calibration sites and their group assignments;
#candidate site groups complete;
#############;
#alternative is non-level pruning. Use the following to interactively;
#pick out desired clusters and store their observation numbers;

plot(clus1,which.plots=2,labels=row.names(bugcal.pa),cex=.4);
#first, spread out the dendrogram plot on screen;
#then run. the following. Right-click to stop;
ccc<-identify(as.hclust(clus1)); #interactive ID of clusters on dendrogram;
#convert result to a group ID vector;
grp.id<-rep(NA,dim(bugcal)[[1]]);
for(k in 1:length(ccc)) {grp.id[ccc[[k]]]<-k};
table(grp.id); #count number of sites in each group;
cbind(row.names(predcal),grp.id); #list calibration sites and their group assignments;
#Cluster analysis is complete;
###########;
#as a check, plot the geographic  locations of the site clusters;
# Can also use the "maps" package to draw a state or regional boundary;
plot(predcal$X_coord[grps.6==1], predcal$Y_coord[grps.6==1],col='black',
    type='p',xlim=c(-125,-118),ylim=c(42,46));
points(predcal$X_coord[grps.6==2], predcal$Y_coord[grps.6==2],col='red')
points(predcal$X_coord[grps.6==3], predcal$Y_coord[grps.6==3],col='green')
points(predcal$X_coord[grps.6==4], predcal$Y_coord[grps.6==4],col='orange')
points(predcal$X_coord[grps.6==5], predcal$Y_coord[grps.6==5],col='blue')
points(predcal$X_coord[grps.6==5], predcal$Y_coord[grps.6==5],col='purple')
#####################################;

# STEP 3 . BUILD RANDOM fOREST MODEL TO PREDICT GROUP MEMBERSHIP;
# First, put the group IDs and predictor variables into a single data frame;
rfdat<-data.frame(predcal[,candvar],Groups=grps.6);
rfdat$Groups<-factor(rfdat$Groups);

#Build RF model;
## NOTE !!!-- I have just used default settings for key fitting parameters, like "mtry". ;
# These choices could use some research, for the RIVPACS context! ;
rf.mod<-randomForest(Groups ~ ., data=rfdat, ntree=500, importance=TRUE, norm.votes=TRUE, keep.forest=TRUE);
 print(rf.mod);   #Displays out-of-bag (OOB) error matrix. Columns are predicted class, rows are true class;;
#various diagnostics;
varImpPlot(rf.mod);  #plots 2 measures of predictor importance;
#Note that all candidate predictors are used in a RF model.;
# Suggest re-running RF with only a few of the most important predictors;
# May give almost-equal performance;

# For key predictors, useful to look at partial dependence plots. ;
#Following example is the role of 1 predictor, "El_m_sqt", in predicting the various site groups; 
windows();
#On graph window menu bar, click "History", then "Recording";
# Then execute;
sapply(unique(rfdat$Groups),function(grp){
    partialPlot(rf.mod,pred.data=rfdat,x.var="El_m_sqt",which.class=grp,main=paste("Group ",grp))});
#now click on graph window and use PgUp, PgDn to scroll through plots;
 #Each plot is for a different site group. X-axis is "El_m_sqt".
 # Y-axis is mean value of logit(p), where p=predicted probability of being in that group,
 # and the mean is taken over all other combinations of the other predictors;
# See Cutler, et al. 2007. Ecology 88:2783-2792 for more info;
 ####

# Once the final RF model is chosen, model development is complete;
#############################################################;

#STEP 4 - Save the final RF predictive model, for future use;
# To specify the entire, final RF model, you need to store 4 things as an R object;
# 4.1) The site-by-taxa matrix of observed presence/absence at calibration sites (bugcal.pa, already available);
# 4.2) Specify the vector of final group membership assignments at calibration sites(grps.final);
      grps.final<-grps.6;
# 4.3) Specify the final predictor variables ;
      preds.final<-candvar;
# 4.4) The final RF model (rfmod);

# Save the model components together in a single .Rdata file.;
    # Any R user can load this file, along with model.predict.RF.r, to make predictions from the model;
     save(bugcal.pa, predcal, grps.final, preds.final, rf.mod, file='My.RF.Model.Version1.Rdata');
    #NOTE - Predcal is not needed to define the model, but is included so that users see the 
    #       required format for predictor data;
#################################;

#Step 5 - Further checks on performance of the final, chosen model;

#Option 5.1 - Make predictions of E and O/E for calibration (reference) sites. Examine O/E statistics and plots;
  # To do this, run the model.predict.RanFor.4.2 function, using the calibration data as the 'new' data;
  # See Step 7 below, for more info on making predictions;
  # Also see internal documentation of model.predict.Ran.For.4.2;
  source("model.predict.RanFor.4.2.r");
  # Two options, for calibration data;
  # Option 1 - Set Cal.OOB=TRUE, for out-of-bag predictions (see Cutler et.al). Gives more realistic (larger) SD(O/E), appropriate for new data; 
  OE.assess.cal<-model.predict.RanFor.4.2(bugcal.pa,grps.final,preds.final, ranfor.mod=rf.mod,prednew=predcal,bugnew=bugcal.pa,Pc=0.5,Cal.OOB=TRUE);
  # Option 2 - Set Cal.OOB=FALSE, for in-bag predictions. Gives optomistically small SD(O/E), because RF models are tightly tuned to in-bag calibration data; 
  OE.assess.cal<-model.predict.RanFor.4.2(bugcal.pa,grps.final,preds.final, ranfor.mod=rf.mod,prednew=predcal,bugnew=bugcal.pa,Pc=0.5,Cal.OOB=FALSE);

#look at other prediction results, for calibration sites;
  names(OE.assess.cal);   #names of 2 components of the prediction results list;
  head(OE.assess.cal$OE.scores); #data frame of O/E scores, 1st 5 rows;
  head(OE.assess.cal$Capture.Probs); #predicted capture probabilties, 1st 5 rows;
  head(OE.assess.cal$Group.Occurrence.Probs); #predicted group occurrence probabilities, 1st 5 rows;
#check distribution of Calibration-site O/E scores. Is it Normal?;
#plot a histogram and a Normal q-q plot;
par(mfrow=c(2,1));
hist(OE.assess.cal$OE.scores$OoverE,xlab="O/E");
qqnorm(OE.assess.cal$OE.scores$OoverE);

#scatterplot of O (on y-axis) vs E (on x-axis). See Pineiro et al. Ecol. Modelling 2008, 316-322, for this choice of axes;
par(mfrow=c(1,1));
plot(OE.assess.cal$OE.scores[,c('E','O')],xlab='Expected richness',ylab='Observed richness');
  abline(0,1); #add a 1-1 line;
#########;

#Option 5.1.5 - Calculate replicate sampling SD of O/E, as a "perfect model" lower bound for SD(O/E) on calibration data;
# reference is Van Sickle et al. null model paper;
#first, compile the following function;
rep.sam.sd<-function(occprb.cal,Pc) {
        cal.cut<-(occprb.cal>=Pc); #TRUE/FALSE matrix denoting which taxa are above predicted P cutoff;
        #use occurrence probabilities only of taxa above the cutoff (cal.cut='TRUE');
        E.cal<-apply(occprb.cal*cal.cut,1,sum); #vector of predicted E ;
        # numerator of site-specific replicate sampling var. Result is a site vector
        RS.cal<-apply(occprb.cal*cal.cut,1,function(x)sum(x*(1-x)));
        SDRS<-sqrt(mean(RS.cal/(E.cal^2))); #replicate sampling SD is sqrt(mean(site-specific replicate sampling variances));
        print(' ',quote=F)
        print(' Replicate sampling SD of O/E: ',quote=F)
        print(SDRS,digits=4);  }; #end of function;
#Then execute the above function, using either in-bag or OOB predicted occurrence probs ('Capture probs') for the calibration data;
       rep.sam.sd(occprb.cal=OE.assess.cal$Capture.Probs,Pc=0.5);
 ###########;

### Option 5.2 - Repeat Step 5.1, but this time use validation data. Check especially for model bias (mean(O/E) differs from 1.0);
   OE.assess.vld<-model.predict.RanFor.4.2(bugcal.pa,grps.final,preds.final, ranfor.mod=rf.mod,prednew=pred.vld,bugnew=bug.vld.pa,Pc=0.5,Cal.OOB=FALSE)  ;
   OE.assess.vld$OE.scores;

 ###############;

## Step 6 --  As of 12/31/10, no facility is available to submit RanFor models on the WCEM website.

################################;

# Step 7 - Making predictions for new data, using random forests;.

# first, source the prediction script and also load the desired model;
   source("model.predict.RanFor.4.2.r");
   load('My.RF.Model.Version1.Rdata');

   # User must supply a sample-by-taxa matrix of taxon abundance or else presence/absence (coded as 1 or 0), for all new samples;
   # User must also supply a corresponding file of predictor data for those same samples;
   # These 2 files should have similar formats as the original taxa and predictor data sets used to build the model (see step 1 above);
   # Notes on format --
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
  OE.assess.test<-model.predict.RanFor.4.2(bugcal.pa,grps.final,preds.final,ranfor.mod=rf.mod, prednew=pred.test,bugnew=bug.test.pa,Pc=0.5,Cal.OOB=FALSE);

# look at O/E scores, for all samples;
 OE.assess.test$OE.scores;

################ ;

## Assessing individual sites;
source("assess.one.sample.4.1.r")
#This function assesses a single site or sample from a new (test) data set to which
# model.predict.RanFor.4.2() has already been applied.
# assess.one.sample() compares observed occurrences with the model-predicted probabilities of occurrence for all taxa;

#Input parameters are:
#       case -- A selected site or sample ID, for which a prediction has already been made using model.predict.v4(). ;
# result.prd -- Output of model.predict.RanFor.4.1() for new samples that include the chosen case;
# bugnew  -- Sample-by-taxa matrix of new samples that was submitted to model.predict.RanFor.4.1.().
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









