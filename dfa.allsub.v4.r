
## Function for all subsets DFA - version 3;
#See dfa.allsub.3.help.tx;
# function output is changed to a list containing bestmods and a vector of null model stats;
#version 3 -- add stats for BC index;
#version 4 - also output the CV confusion matrix for every "best" model;

#function definition;           
dfa.allsub.v3<-function(bug.cal,bug.vld,pred.cal,pred.vld,grps,candvar,numkeep,Pc) {
#####################;
## begin function ops;
 grps<-factor(grps); #declare group factor as factor;
 ngrps<-length(unique(grps)); #number of groups;
 nsite.cal<-dim(pred.cal)[[1]]; #number of calibration sites;
 nu.E<-nsite.cal-ngrps; #within-groups df, null model; 
 nu.H<-ngrps-1; #between groups df;
 flush.console();
 print(" New run of all-subsets DFA",quote=F); 
 grpsiz<-table(grps);
 print("Count of calibration sites in each group",quote=F);     
 print(grpsiz,quote=F);
 
  ncand<-length(candvar); ##number of candidate predictors;
  maxord<-length(numkeep[numkeep>0]); #maximum desired model order;
#############################################;
# BEST SUBSETS DFA.
#Stage 1 -- screen all possible models for wilks lambda;

# 1A -- Construct formulas for all possible model of order 1 to maxord(<=ncand);
#create list containing permutations of all candidates, of all desired orders;
#list components are all groups of models of each order
modls<-sapply(1:maxord, function(xx) {
 combinations(ncand,xx,candvar)});

#count number of models each group, and the total;
mod.cnt<-sapply(modls,dim); #first row = number of models of each order;
totmod<-sum(mod.cnt[1,]); #total number of models (not including null);
print(paste(c(totmod,' models being screened. Please wait ...'),collapse=" "),quote=F);
flush.console();

#curmod<-modls[[4]][841,]; #for testing only;

#compute F-stat of 1-way anova for each order-1 model;
univar<-sapply(modls[[1]],function(curmod) {
cc<-summary(aov(as.matrix(pred.cal[,curmod])~grps));
cc[[1]][1,4]; #return the F-value;
  });

#Wilks lambda and F-stat , compute for all models of order >=2;
#Apply to each model of each component of list;
#THIS COULD TAKE A LOT OF TIME;

#start.time=proc.time();
wilk.all<-lapply(modls[2:maxord],function(xxx)
     {apply(xxx,1,function(curmod) {
      cc<-summary(manova(as.matrix(pred.cal[,curmod])~grps),test='Wilks');
       cc$stats[1,2:3]; #extract & return Wilks lambda and F statistic;
            } ) } );

# Stage 2 -- Extract and explore a small set of best models of each desired order;

#first sort models within each order by their F-statistics, largest to smallest;
#best models are at the top of the list;
#also sort the statistics lists in same fashion; 
#Within a given model order, approximate (multivariate) or exact (univariate) F has same df for all models;
modl.sort<-modls;
wilk.sort<-wilk.all;
modl.sort[[1]]<-modl.sort[[1]][order(-univar)];
univar.sort<-univar[order(-univar)];
for(mord in 2:maxord) {;
  modl.sort[[mord]]<-modl.sort[[mord]][order(-wilk.all[[(mord-1)]]['approx F',]),]
  wilk.sort[[mord-1]]<-wilk.all[[mord-1]][,order(-wilk.all[[(mord-1)]]['approx F',])]
  };
#sorting complete;
print('Model screening complete. Computing O/E and BC stats for subset of best models....',quote=F);
######################################;
#setup for detailed look at best-model subset;

#specify a small subset of best models to keep and explore further;
# for example, keep 5 models for each order, 1-14, and also the single full model;

totmod<-sum(numkeep); #total number to keep;
cumcnt<-cumsum(numkeep);

#before looping over models, compute pieces that are constant for all models;
#compute group occurrence freqs at calibration sites, for all calibration-site taxa;
#matrix of relative occurrences of each spp at sites in each group of reference sites;
grpocc<-apply(bug.cal,2,function(x){tapply(x,grps,function(y){sum(y)/length(y)})});

#compute and print null model results;
   #null model occurrence probs;
    pnull<-apply(bug.cal,2,sum)/dim(bug.cal)[[1]];
    #Compute Expected richness (E) for null model using taxa above Pcutoff.;
    # Note that, with null model, these taxa are fixed for all sites;
   nulltax<-names(pnull[pnull>=Pc]); #subset of taxa with Pnull exceeding cutoff;
   Enull<-sum(pnull[nulltax]);
   O.clb.nul <- apply(bug.cal[,nulltax],1,sum); #observed null richness;
   OE.cal<-O.clb.nul/Enull; #vector of null-model O/E at calibration sites;
   O.vld.nul<- apply(bug.vld[,nulltax],1,sum) ;
   OE.vld<-O.vld.nul/Enull; #vector of null-model O/E at validation sites;

   #compute vectors of null BC, clb & vld;
   sadiff.clb<-apply(bug.cal[,nulltax],1,function(x)sum(abs(x-pnull[nulltax])));
   BC.clb.nul<-sadiff.clb/(O.clb.nul+Enull);
   sadiff.vld<-apply(bug.vld[,nulltax],1,function(x)sum(abs(x-pnull[nulltax])));
   BC.vld.nul<-sadiff.vld/(O.vld.nul+Enull);
     #print null O/E stats;
  print(paste(c('Calibration: Mean, SD, RMSE of Null O/E = ',
            round(mean(OE.cal),3),round(sqrt(var(OE.cal)),3),
            round(sqrt(mean((OE.cal-1)^2)),3)),collapse="  "), quote=F);
  print(paste(c('Calibration: Median, 90%ile of BC = ',
            round(median(BC.clb.nul),3),round(quantile(BC.clb.nul,probs=0.90),3)),collapse="  "),quote=F);
   print(paste(c('Validation: Mean, SD, RMSE of Null O/E = ',
            round(mean(OE.vld),3),round(sqrt(var(OE.vld)),3),
            round(sqrt(mean((OE.vld-1)^2)),3)),collapse="  "), quote=F);
   print(paste(c('Validation: Median, 90%ile of BC = ',
            round(median(BC.vld.nul),3),round(quantile(BC.vld.nul,probs=0.90),3)),collapse="  "),quote=F);
  #VERSION 3 -- save vector of null stats;
   nullstat<-c(mean(OE.cal),sqrt(var(OE.cal)),sqrt(mean((OE.cal-1)^2)),  
                  mean(OE.vld),sqrt(var(OE.vld)),sqrt(mean((OE.vld-1)^2)),
                   median(BC.clb.nul),quantile(BC.clb.nul,probs=0.90),
                   median(BC.vld.nul),quantile(BC.vld.nul,probs=0.90)) ;
    names(nullstat)<-c("Mean.cal","SD.cal","RMSE.cal","Mean.vld","SD.vld","RMSE.vld",
                       "BCMD.cal","BC90.cal","BCMD.vld","BC90.vld");
  

#empty matrix of results for all kept models;
subset.stats<-data.frame(order=rep(NA,totmod),F.stat=rep(NA,totmod), Wilks=rep(NA,totmod), 
     cls.crct.resub=rep(NA,totmod),cls.crct.cv=rep(NA,totmod),
     MNOE.cal=rep(NA,totmod),SDOE.cal=rep(NA,totmod),RMSE.cal=rep(NA,totmod),SDRS.cal=rep(NA,totmod),
     MNOE.vld=rep(NA,totmod),SDOE.vld=rep(NA,totmod),RMSE.vld=rep(NA,totmod),
     BCMD.cal=rep(NA,totmod),BC90.cal=rep(NA,totmod),BCMD.vld=rep(NA,totmod),BC90.vld=rep(NA,totmod),
     model=rep(NA,totmod));
     
#empty list for storing CV confusion matrices;
conmat.list<-vector("list",totmod);
names(conmat.list)<-rep(as.character(" "),totmod);
#loop over model orders and over subset models within each order, fill in results matrix;   
    
    indx<-0; #initialize results index;
for(mord in 1:maxord) {
    for(imod in 1:numkeep[mord]) {
        #mord<-2; #imod<-1;
        indx=indx+1;
        if(mord==1){
           curmod<-modl.sort[[mord]][imod];
           subset.stats[indx,'model']<-curmod;
           subset.stats[indx,'F.stat']<-univar.sort[imod];
           subset.stats[indx,'Wilks']<-1/(1+univar.sort[imod]*nu.H/nu.E)
             }
         else if(mord==ncand){
            curmod<-modl.sort[[mord]];
            subset.stats[indx,'model']<-paste(curmod,collapse=" ");
            subset.stats[indx,'F.stat']<-wilk.sort[[mord-1]]['approx F'];
            subset.stats[indx,'Wilks']<- wilk.sort[[mord-1]]['Wilks'];
             }
         else if(mord>1){
            curmod<-modl.sort[[mord]][imod,];
            subset.stats[indx,'model']<-paste(curmod,collapse=" ");
            subset.stats[indx,'F.stat']<-wilk.sort[[mord-1]]['approx F',imod];
            subset.stats[indx,'Wilks']<- wilk.sort[[mord-1]]['Wilks',imod];
             };
         #flush.console();
         #print(curmod);
        subset.stats[indx,'order']<-mord;
        #compute classification errors for curmod;
        #build model formula;
         ifelse(length(curmod)>1,rhs<-paste(curmod,collapse="+"),rhs<-curmod);
         form<-paste("grps",rhs,sep="~");
        #fit DFA model, get resubs. error;
         qq<-lda(formula(form),data=pred.cal,CV=F); 
         resub<-predict(qq);
         bmat<-table(grps,resub$class);
         subset.stats[indx,'cls.crct.resub']<- 100*sum(diag(bmat))/sum(bmat);
        #repeat above steps to get CV error;
         zz<-lda(formula(form),data=pred.cal,CV=T); 
         cmat<-table(grps,zz$class);
         subset.stats[indx,'cls.crct.cv']<- 100*sum(diag(cmat))/sum(cmat);
        #store the CV confusion matrix, named by the set of model predictors;
         names(conmat.list)[indx]<-subset.stats[indx,'model'];
         conmat.list[[indx]]<-cmat;
        
        #compute O/E and stats for calibration and validation sites;
        #matrices of predicted group membership probs for 2 types of sites;
        calmem.prb<-predict(qq)$posterior; #for calibration sites;
        vldmem.prb<-predict(qq,newdata=pred.vld)$posterior; #for validation sites;
        #matrices of predicted taxa occurrence probs;
        occprb.cal<-calmem.prb%*%grpocc; 
        occprb.vld<-vldmem.prb%*%grpocc; 
        #vectors of E and O;
        cal.cut<-(occprb.cal>=Pc); #TRUE/FALSE matrix denoting which taxa at calibration sites are above predicted P cutoff;
        #for E and O do row sums of predicted prob and observed matrices, only for "TRUE" taxa in cal.cut;
        # achieved by elementwise product of these matrices and cal.cut;
        E.cal<-apply(occprb.cal*cal.cut,1,sum);
        O.cal<-apply(bug.cal*cal.cut,1,sum);
       # numerator of site-specific replicate sampling var;
        RS.cal<-apply(occprb.cal*cal.cut,1,function(x)sum(x*(1-x)));
       # stats of O/E and RS std dev;
        subset.stats[indx,"MNOE.cal"]<-mean(O.cal/E.cal);
        subset.stats[indx,"SDOE.cal"]<-sqrt(var(O.cal/E.cal));
        subset.stats[indx,"RMSE.cal"]<-sqrt(mean((O.cal/E.cal-1)^2));
        subset.stats[indx,"SDRS.cal"]<-sqrt(mean(RS.cal/(E.cal^2)));
       #stats of BC ;
         BC.cal<-apply(abs(bug.cal-occprb.cal)*cal.cut,1,sum)/(O.cal+E.cal);
         subset.stats[indx,"BCMD.cal"] <-median(BC.cal);
         subset.stats[indx,"BC90.cal"] <-quantile(BC.cal,probs=0.90);
       #repeat above steps for validation sites;
        vld.cut<-(occprb.vld>=Pc); 
        E.vld<-apply(occprb.vld*vld.cut,1,sum);
        O.vld<-apply(bug.vld*vld.cut,1,sum);
        # stats of O/E and RS std dev;
        subset.stats[indx,"MNOE.vld"]<-mean(O.vld/E.vld);
        subset.stats[indx,"SDOE.vld"]<-sqrt(var(O.vld/E.vld));
        subset.stats[indx,"RMSE.vld"]<-sqrt(mean((O.vld/E.vld-1)^2));
        #stats of BC for vld;
         BC.vld<-apply(abs(bug.vld-occprb.vld)*vld.cut,1,sum)/(O.vld+E.vld);
         subset.stats[indx,"BCMD.vld"] <-median(BC.vld);
         subset.stats[indx,"BC90.vld"] <-quantile(BC.vld,probs=0.90);
               } #end of model loop;
  } #end of order loop;
  print("All-subsets function finished.");
 return(list(subset.stats=subset.stats,null.stats=nullstat,CV.error.matrices=conmat.list)); #output list contains "best" model summary data frame and vector of null-model stats;
 }; #end of all-subsets function;
