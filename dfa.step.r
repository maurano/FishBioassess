#########    R code for stepwise linear discriminant analysis (DFA.STEP.R);

# John Van Sickle, USEPA/NHEERL/WED, Corvallis OR. ;

#  Version 1.0 (August 5, 2004);

# See dfa.step.help.txt for usage instructions and help;

#############################;
#function for forward or backward stepwise DFA;
dfa.step<-function(predcal,candvar,grps,method='forward',P.enter.lim=0.15,P.stay.lim=0.20) {
#set up additional parameters;
 grps<-factor(grps); #factorize groups vector;
 ncand<-length(candvar);
 ngrps<-length(unique(grps)); #number of groups;
 nu.E<-dim(predcal)[[1]]-ngrps; #within-groups df, null model; 
 nu.H<-ngrps-1; #between-groups df;

if(method=='forward') {
  #initialize;
  looplim<-100;
  curmod=NULL;
  vars.out=candvar;
  step<-0;
  print('Beginning forward stepwise selection',quote=F);

  for(index in 1:looplim) {
  if(length(curmod)==ncand){
     print('All variables entered',quote=F);
     print('',quote=F);
     print('Re-fit final model',quote=F);
     zz<-finalmod(curmod=curmod,grps=grps,dset=predcal);
      return(zz);
        }
  nout<-length(vars.out);
  enter.stats<-data.frame(F.to.enter=rep(NA,nout),P.to.enter=rep(NA,nout))
  row.names(enter.stats)<-vars.out;
  enter.stats<-F.enter(curmod=curmod,vars.out=vars.out,grps=grps,dset=predcal,
         nu.E=nu.E,nu.H=nu.H,es=enter.stats);
  flush.console();
  print('',quote=F);
  print('Variables NOT in current model',quote=F);
  print(round(enter.stats,4));
 
#forward step - check P.to enter and decide on entry variable;
 enter.stats<-enter.stats[order(-enter.stats$F.to.enter),];
 if(enter.stats$P.to.enter[1]>P.enter.lim){
      print('No more significant variables',quote=F);
      print('',quote=F);
      print('Re-fit final model',quote=F);
      print(curmod);
      zz<-finalmod(curmod=curmod,grps=grps,dset=predcal);
      return(zz) ;
        }
 else {
   #add variable with biggest F;
   step<-step+1;
   addvar<-row.names(enter.stats)[1];
   print('',quote=F);
   print(paste(c('Step ',step,' -- Add ',addvar),collapse=""),quote=F);
   vars.out<-vars.out[vars.out!=addvar];
   curmod<-c(curmod,addvar);
     }
 #removal step;
  n.in<-length(curmod);
  stay.stats<-data.frame(F.to.stay=rep(NA,n.in),P.to.stay=rep(NA,n.in))
  row.names(stay.stats)<-curmod;
  stay.stats<-F.stay(curmod=curmod,grps=grps,dset=predcal,
             nu.E=nu.E,nu.H=nu.H,stay=stay.stats);
  flush.console();
  print('',quote=F);
  print('Variables IN current model',quote=F);
  print(round(stay.stats,4));
  stay.stats<-stay.stats[order(stay.stats$F.to.stay),];
  if(stay.stats$P.to.stay[1]>P.stay.lim){
      step<-step+1;
      dropvar<-row.names(stay.stats)[1];
      print('');
      print(paste(c('Step ',step,' -- Remove ',dropvar),collapse=""),quote=F);
      vars.out<-c(vars.out,dropvar);
      curmod<-curmod[curmod!=dropvar];
     }
    } #end of main loop;
    } #end of forward option;

else if(method=='backward') {
  #initialize;
  looplim<-100;
  curmod<-candvar;
  vars.out<-NULL;
  step<-0;
  print('Beginning backward stepwise selection',quote=F);

  for(index in 1:looplim) {
  if(is.null(curmod)){
     print('All variables removed',quote=F);
     print('',quote=F);
     return(curmod);
        }
  #backward step -- check F to stay for all variables in current model;
   n.in<-length(curmod);
  stay.stats<-data.frame(F.to.stay=rep(NA,n.in),P.to.stay=rep(NA,n.in))
  row.names(stay.stats)<-curmod;
  stay.stats<-F.stay(curmod=curmod,grps=grps,dset=predcal,
             nu.E=nu.E,nu.H=nu.H,stay=stay.stats);
  flush.console();
  print('',quote=F);
  print('Variables IN current model',quote=F);
  print(round(stay.stats,4));
  stay.stats<-stay.stats[order(stay.stats$F.to.stay),];
  if(stay.stats$P.to.stay[1]>P.stay.lim){
      step<-step+1;
      dropvar<-row.names(stay.stats)[1];
      print('');
      print(paste(c('Step ',step,' -- Remove ',dropvar),collapse=""),quote=F);
      vars.out<-c(vars.out,dropvar);
      curmod<-curmod[curmod!=dropvar];
         }
  else {
      print('No more variables will be removed',quote=F);
      print('',quote=F);
      print('Re-fit final model',quote=F);
      zz<-finalmod(curmod=curmod,grps=grps,dset=predcal);
      return(zz) ;
        }

 #forward step - compute P.to enter and decide on entry variable; 
  nout<-length(vars.out);
  enter.stats<-data.frame(F.to.enter=rep(NA,nout),P.to.enter=rep(NA,nout))
  row.names(enter.stats)<-vars.out;
  enter.stats<-F.enter(curmod=curmod,vars.out=vars.out,grps=grps,dset=predcal,
         nu.E=nu.E,nu.H=nu.H,es=enter.stats);
  flush.console();
  print('',quote=F);
  print('Variables NOT in current model',quote=F);
  print(round(enter.stats,4));
  enter.stats<-enter.stats[order(-enter.stats$F.to.enter),];
  if(enter.stats$P.to.enter[1]<P.enter.lim){
   #add variable with biggest F;
   step<-step+1;
   addvar<-row.names(enter.stats)[1];
   print('',quote=F);
   print(paste(c('Step ',step,' -- Add ',addvar),collapse=""),quote=F);
   vars.out<-vars.out[vars.out!=addvar];
   curmod<-c(curmod,addvar);
     } #end of entry option;
   
    } #end of main loop;
    } #end of backward option;

   print('',quote=F);
   print('End of program',quote=F);

   } #end of step.dfa function;
###################;

#fit final model using lda(), to determine classification accuracies;
finalmod<-function(curmod,grps,dset) {
  #construct formula;
  ifelse(length(curmod)>1,rhs<-paste(curmod,collapse="+"),rhs<-curmod);  
  bbb<-formula(paste("grps",rhs,sep="~"));
# In-sample classification accuracy ;
  qq<-lda(bbb,data=dset,CV=F); 
  rr<-predict(qq);
print('');
print('Classification accuracy of final model',quote=F);
print("Predicted vs observed classes: Resubstitution",quote=F);
  print("Resubstitution Confusion matrix: Rows = Observed class, Cols. = Predicted class",quote=F);
  bmat<-table(grps,rr$class);
  print(bmat);
  print(c('Overall resubstitution pct. correct = ',round(100*sum(diag(bmat))/sum(bmat),1)),quote=F);

  zz<-lda(bbb,data=dset,CV=T); 
  print('');
  print("Predicted vs observed classes: Leave-one-out crossvalidation (CV)",quote=F);
  print("CV Confusion matrix: Rows = Observed class, Cols. = Predicted class",quote=F);
  cmat<-table(grps,zz$class);
  print(cmat);
  print(c('Overall CV pct. correct = ',round(100*sum(diag(cmat))/sum(cmat),1)),quote=F);
  return(qq) ;#return the lda object w/resubstitution;
         } #end of finalmod function;
  

      ##############;

#############################################;
# function calculates F, P to enter, for all candidate variables not in current model; 

F.enter<-function(curmod,vars.out,grps,dset,nu.E,nu.H,es) {
#if current model is empty, use F-test. Else use Wilks lamda;
if(is.null(curmod)){ 
#compute F-stat of 1-way anova for all vars;
for(curvar in vars.out){
cc<-summary(aov(as.matrix(dset[,curvar])~grps));
es[curvar,'F.to.enter']<-cc[[1]][1,4]; # F-value;
es[curvar,'P.to.enter']<-cc[[1]][1,5]; # P-value;
    } #loop end;
     }

else { 
# Or, current model is not empty. Use Wilks lambda;
# First get Wilks Lambda for current model. ;
#If current model has 1 variable, then Wilks = transformed F from 1-way ANOVA;
    if(length(curmod)==1) {;
       ff<-summary(aov(as.matrix(dset[,curmod])~grps));
       fstat<-ff[[1]][1,4]; #F-statistic for current model;
       wilks.cur<-1/(1+fstat*nu.H/nu.E);
      }
#If current model has >= 2 variables, then get Wilks from MANOVA;
    else {
       man<-summary(manova(as.matrix(dset[,curmod])~grps),test='Wilks');
       wilks.cur<-man$stats[1,2]; #extract & return Wilks lambda;
          }
#compute LR test for newmodel vs curmod, where newmod=(curmod+newvar),;
#for each newvar in the not-in-model variables;
       #print(c("Wilks current =  ",wilks.cur));
        for(curvar in vars.out){
          newmod<-c(curmod,curvar);
          man.new<-summary(manova(as.matrix(dset[,newmod])~grps),test='Wilks');
          wilks.new<-man.new$stats[1,2];
          wilks.part<-wilks.new/wilks.cur;
          df.2<-nu.E-length(curmod);
          F.part<-((1-wilks.part)/wilks.part)*df.2/nu.H;
          es[curvar,'F.to.enter']<-F.part;
          es[curvar,'P.to.enter']<-1-pf(F.part,nu.H,df.2);
                            } #end of LR test;
         } #end of block for current model>=2 vars;
     return(es);
  }; #end of F-to-enter function;

#####################################;

# function calculates F-to-stay and associated P-value, ;
# for all variables in the current model;

F.stay<-function(curmod,grps,dset,nu.E,nu.H,stay) {

n.in<-length(curmod);
#if current model has 1 var, then F-to-stay is from 1-way ANOVA;
if(n.in==1) { 
   cc<-summary(aov(as.matrix(dset[,curmod])~grps));
   stay[curmod,'F.to.stay']<-cc[[1]][1,4]; # F-value;
   stay[curmod,'P.to.stay']<-cc[[1]][1,5]; # P-value;
       }
#current model has >=2 vars;
else {
     #get Wilks.cur for current model, using MANOVA;
     man<-summary(manova(as.matrix(dset[,curmod])~grps),test='Wilks');
     wilks.cur<-man$stats[1,2]; #extract & return Wilks lambda;
     #print(c("Wilks current =  ",wilks.cur));

    # get Wilks for reduced models Use ANOVA or MANOVA, depending on reduced model size;
     for(curvar in curmod) {
       newmod<-curmod[curmod!=curvar];
       n.new<-length(newmod);
       if(n.new==1){
           cc<-summary(aov(as.matrix(dset[,newmod])~grps));
           wilks.new<-1/(1+cc[[1]][1,4]*nu.H/nu.E); }
       else {
           man.new<-summary(manova(as.matrix(dset[,newmod])~grps),test='Wilks');
           wilks.new<-man.new$stats[1,2]; }
        #print(c(curvar,wilks.new));
        
       wilks.part<-wilks.cur/wilks.new;
       df.2<-nu.E-length(newmod);
       F.part<-((1-wilks.part)/wilks.part)*df.2/nu.H;
       stay[curvar,'F.to.stay']<-F.part;
       stay[curvar,'P.to.stay']<-1-pf(F.part,nu.H,df.2);
             } #end of loop over vars in model;
        } #end of option for >=2 vars in model;
      return(stay);
      } #end of function;
   










  
