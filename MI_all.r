# This scripts calculate Mutual Information and use Naive Bayes and Support vector machine to validate the 5mirna.
# Require to load 'spall-dist-ptab.rdata' first.


# Select 88 samples: 22 from Recurrence, 66 from No recurrence
c1recur<-which(dist12$pclst==1&dist12$recur=='Recurrence')
c2nonrecur<-which(dist12$pclst==2&dist12$recur=='No Recurrence')
train88<-c(c1recur,c2nonrecur)
strain88<-spall[train88,]

# Caculate Q1,Q2,Q3,MAX for intervals [Min,Q1),[Q1,Q2),[Q2,Q3),[Q3,Max]
i1<-apply(strain88[,14:860],2,function(x)quantile(x,probs=0.25))
i2<-apply(strain88[,14:860],2,function(x)quantile(x,probs=0.50))
i3<-apply(strain88[,14:860],2,function(x)quantile(x,probs=0.75))
i4<-apply(strain88[,14:860],2,function(x)quantile(x,probs=1))
intv88<-rbind(i1=i1,i2=i2,i3=i3,i4=i4)
write.csv(intv88,'intv88.csv')

# 88 samples fall into four bins.
sint88<-strain88[,14:860]
for( i in 1:88) 
{ for(j in 1:847) 
   { if(sint88[i,j]-intv88[1,j]<0) {sint88[i,j]<-1} 
     else if(sint88[i,j]-intv88[2,j]<0) {sint88[i,j]<-2}
     else if(sint88[i,j]-intv88[3,j]<0) {sint88[i,j]<-3}
     else {sint88[i,j]<-4}
    }
}
sint88$recur<-strain88$recur

# MIfun function: caculate MI value
MIfun<-function(t,value) { N1a<-22
N0a<-66
N<-88
N11<-t['Recurrence',value]
N10<-N1a-N11
N01<-t['No Recurrence',value]
N00<-N0a-N01
Na1<-N01+N11
Na0<-N00+N10
Mi<-0;
if(N11!=0) {Mi<-Mi+(N11/N)*log2(N*N11/(N1a*Na1))  }
if(N01!=0) {Mi<-Mi+(N01/N)*log2(N*N01/(N0a*Na1)) }
if(N10!=0) {Mi<-Mi+(N10/N)*log2(N*N10/(N1a*Na0)) }
if(N00!=0) {Mi<-Mi+(N00/N)*log2(N*N00/(N0a*Na0)) }
return(Mi)
}

# caculate 847*4 MI matrix
MI88<-data.frame(row.names=c('m1','m2','m3','m4'))
for ( i in 1:847) 
 { nc<-colnames(sint88)[i]; 
   t<-table(sint88$recur, sint88[,i]);
   m1<-MIfun(t,'1');
   m2<-MIfun(t,'2');
   m3<-MIfun(t,'3'); 
   m4<-MIfun(t,'4');
   MI88[,nc]<-c(m1,m2,m3,m4)
 }
write.csv(MI88,'MI88.csv')

#top 50 MI values and their miRNA and interval (table S1)

allMI88<-as.vector(as.matrix(MI88))
allMI88sort<-sort(allMI88,decreasing=TRUE)
top50<-unique(allMI88sort[1:50])
mrna<-c();
miv<-c();
interval<-c();
intvalue<-c();

for ( i in 1:length(top50) ) # cacuplate table S1
 {
  tl<-which(MI88==top50[i],arr.ind=TRUE);
  tr<-dim(tl)[1];
  for( j in 1:tr) 
   {
  mrna<-c(mrna, colnames(intv88)[tl[j,2]]);
  miv<-c(miv,top50[i]);
  interval<-c(interval,rownames(tl)[j]);
  intvalue<-c(intvalue, intv88[tl[j,1],tl[j,2]]);
   }
 }
mitop50<-cbind(miRNA=mrna,MI=miv,Interval=interval, Intvalue=intvalue); # matrix for table S1
write.csv(mitop50,'mitop50.csv')

# 176 samples with 847*4 bins ( binary matrix).
sint176<-spall[,14:860]
for( i in 1:176) 
 { for(j in 1:847) 
   { if(sint176[i,j]-intv88[1,j]<0) {sint176[i,j]<-1}
      else if(sint176[i,j]-intv88[2,j]<0) {sint176[i,j]<-2}
      else if(sint176[i,j]-intv88[3,j]<0) {sint176[i,j]<-3}
      else {sint176[i,j]<-4}
    }
  }
write.csv(sint176,'sint176.csv')

# Top 5 miRNAs in [min,Q1] interval
sint5m1<-sint176[,mitop50[1:5,1]]
sint5m1<-ifelse(sint5m1!=1,0,1)

# naive Bayes: using sint5ml binary data: training set: 88 samples; testing data: the left 88 samples
library(e1071)
nb88<-naiveBayes(sint5m1[train88,],as.factor(spall$recur[train88]))
preb88<-predict(nb88,sint5m1)
nbsint5<-cbind(kclst=dist12$kclst,recur=dist12$recur,preb88)
write.csv(nbsint5,'nbsint5.csv')

# support vector machine:  using sint5ml binary data: training set: 88 samples; testing data: the left 88 samples
svm88<-svm(sint5m1[train88,],as.factor(spall$recur[train88]),cost=1000)
psvm88<-predict(svm88,sint5m1)
svmsint5<-cbind(kclst=dist12$kclst,recur=dist12$recur,psvm88)
write.csv(svmsint5,'svmsint5.csv')

#Top 5 miRNA raw data
raw5mRna<-spall[,mitop50[1:5,1]]
# naive Bayes: using 5 miRNAs' rawa data: training set: 88 samples; testing data: the left 88 samples
rawnb88<-naiveBayes(raw5mRna[train88,],as.factor(spall$recur[train88]))
rawpreb88<-predict(rawnb88,raw5mRna)
nbrawmRNA5<-cbind(kclst=dist12$kclst,recur=dist12$recur,rawpreb88)
write.csv(nbrawmRNA5,'nbrawmRNA5.csv')

# support vector machine:  using 5 miRNAs' rawa data: training set: 88 samples; testing data: the left 88 samples
rawsvm88<-svm(raw5mRna[train88,],as.factor(spall$recur[train88]),cost=1000)
rawpsvm88<-predict(rawsvm88,raw5mRna)
svmraw5mRNA<-cbind(kclst=dist12$kclst,recur=dist12$recur,rawpsvm88)
write.csv(svmraw5mRNA,'svmraw5mRNA.csv')

# Kaplan-Meier survival curve for 5 miRNAs
library('survival')
fit122<-survfit(Surv(spall$rfsurv,as.numeric(as.factor(spall$recur)))~as.factor(sint176[,'hsa-miR-122_st']))
fit126<-survfit(Surv(spall$rfsurv,as.numeric(as.factor(spall$recur)))~as.factor(sint176[,'hsa-miR-126_st']))
fit22<-survfit(Surv(spall$rfsurv,as.numeric(as.factor(spall$recur)))~as.factor(sint176[,'hsa-miR-22_st']))
fit15a<-survfit(Surv(spall$rfsurv,as.numeric(as.factor(spall$recur)))~as.factor(sint176[,'hsa-miR-15a_st']))
fit30a<-survfit(Surv(spall$rfsurv,as.numeric(as.factor(spall$recur)))~as.factor(sint176[,'hsa-miR-30a_st']))

pdf('f5miRNA-survival-curvs.pdf')
par(mfrow=c(3,2))
plot(fit122,lty=1:4,col=c('red','black','blue','green'), xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Hsa-miR-122")
text(x=100, y=0.45, paste0("p-value: ",round(pchisq(survdiff(Surv(spall$rfsurv,as.numeric(as.factor(spall$recur)))~as.factor(sint176[,'hsa-miR-122_st']))$chisq,df=3,lower.tail=F),digits=12)))
title(main='A', adj=0)
plot(fit126,lty=1:4,col=c('red','black','blue','green'), xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Hsa-miR-126")
text(x=100, y=0.45, paste0("p-value: ",round(pchisq(survdiff(Surv(spall$rfsurv,as.numeric(as.factor(spall$recur)))~as.factor(sint176[,'hsa-miR-126_st']))$chisq,df=3,lower.tail=F),digits=12)))
title(main='B', adj=0)
plot(fit22,lty=1:4,col=c('red','black','blue','green'), xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Hsa-miR-22")
text(x=100, y=0.45, paste0("p-value: ",round(pchisq(survdiff(Surv(spall$rfsurv,as.numeric(as.factor(spall$recur)))~as.factor(sint176[,'hsa-miR-22_st']))$chisq,df=3,lower.tail=F),digits=12)))
title(main='C', adj=0)
plot(fit30a,lty=1:4,col=c('red','black','blue','green'), xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Hsa-miR-30a")
text(x=100, y=0.4, paste0("p-value: ",round(pchisq(survdiff(Surv(spall$rfsurv,as.numeric(as.factor(spall$recur)))~as.factor(sint176[,'hsa-miR-30a_st']))$chisq,df=3,lower.tail=F),digits=12)))
title(main='D', adj=0)
plot(fit15a,lty=1:4,col=c('red','black','blue','green'), xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Hsa-miR-15a_st")
text(x=100, y=0.4, paste0("p-value: ",round(pchisq(survdiff(Surv(spall$rfsurv,as.numeric(as.factor(spall$recur)))~as.factor(sint176[,'hsa-miR-15a_st']))$chisq,df=3,lower.tail=F),digits=15)))
title(main='F', adj=0)
plot.new()
legend('center',cex=0.8,lty=1:4,col=c('red','black','blue','green'),c('[MIN,Q1)','[Q1,Q2)','[Q2,Q3)','[Q3,MAX]'))
dev.off()