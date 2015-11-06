# This script calculates 60/16/6 min-max features in table 2 of ATJ paper and plot kmeans results.
# Require to load 'spall-dist-ptab.rdata' first.

#60 Min-Max features in table 2 of ATJ paper
hsamin<-read.csv(file='hsaMIN.txt',header=FALSE,colClasses = "character")
hsamax<-read.csv(file="hsaMAX.txt",header=FALSE,colClasses = "character")

# Min-MAX for 847 features
pmax<-as.data.frame(matrix(nrow=89,ncol=0))
for(i in 14:860) 
 { cn<-colnames(spall)[i]
   cp<-tapply(spall[,i],as.factor(spall$pid),max)
   pmax<-cbind(pmax,cp=cp)
   colnames(pmax)[colnames(pmax)=='cp']<-paste0(cn,"_MAX")
 }
 
pmin<-as.data.frame(matrix(nrow=89,ncol=0))
for(i in 14:860) 
 { cn<-colnames(spall)[i]; 
   cp<-tapply(spall[,i],as.factor(spall$pid),min); 
   pmin<-cbind(pmin,cp=cp);
   colnames(pmin)[colnames(pmin)=='cp']<-paste0(cn,"_MIN")
 }

hmin<-pmin[,hsamin$V1]
hmax<-pmax[,hsamax$V1]

# top 60 min-max features + clinical variates
p60new<-cbind(hmin,hmax)
p60new<-merge(ptab,p60new,by.x='pid',by.y='row.names')
write.csv(p60new,file='p60new.csv')

# kmeans for 89 patients with 60 features.
p60km<-kmeans(p60new[,9:68],2)
rownames(p60km$centers)<-c('c1','c2')
pdist60<-rbind(p60new[,9:68],p60km$centers)
pdist60<-as.matrix(dist(pdist60))
pdist60<-as.data.frame(pdist60[1:89,90:91])
pdist60$M<-pdist60$c2-pdist60$c1
pdist60$A<-(pdist60$c2+pdist60$c1)/2
p60new$kclst<-p60km$cluster
p60new<-cbind(p60new,pdist60)
write.csv(p60new,file='p60new.csv')

# Plot the distribution of 89 patients by 60 min-mdax features
pdf('min-max-60f-89p-distribution.pdf')
plot(p60new$A,p60new$M,col=c('black','red')[as.numeric(as.factor(p60new$recur))],pch=c(16,17)[as.numeric(as.factor(p60new$recur))],main="60 Min-Max Features & 89 Patients",xlab="A=(d1+d2)/2",ylab="M=d2-d1")
legend("bottomright",pch=c(16,17),col=c('black','red'),c('No Recurrence','Recurrence'))
dev.off()

# plot the survival curve 
library('survival')
p60fit<-survfit(Surv(p60new$rfsurv,as.numeric(as.factor(p60new$recur))-1)~as.factor(p60new$kclst))
pdf('min-max-60f-89p-survival.pdf')
plot(p60fit,lty=1:2,col=c('black','red'), xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="60 Min-Max Features & 89 Patients")
legend("topright",lty=1:2,col=c('black','red'),c('cluster1','cluster2'))
text(x=90, y=0.5, paste0("p-value: ",round(pchisq(survdiff(Surv(p60new$rfsurv,as.numeric(as.factor(p60new$recur))-1)~as.factor(p60new$kclst))$chisq,df=1,lower.tail=F),digits=8)))
dev.off()

# 16 min-max features
hsa16<-read.csv("hsaHCClit.txt",header=F,stringsAsFactor=F)
phsa16<-p60new[,hsa16$V1]
phsa16km<-kmeans(phsa16,2)
rownames(phsa16km$centers)<-c('c1','c2')
phsa16dist<-rbind(phsa16,phsa16km$centers)
phsa16dist<-as.matrix(dist(phsa16dist))
phsa16dist<-phsa16dist[1:89,90:91]
phsa16<-cbind(phsa16,phsa16dist)
phsa16$M<-phsa16$c2-phsa16$c1
phsa16$A<-(phsa16$c2+phsa16$c1)/2
phsa16$kclst<-ifelse(phsa16$M>0,1,2)
phsa16<-cbind(phsa16,p60new[,1:8])
write.csv(phsa16,file='phsa16.csv')

# plot the distribution of 89 patients by 16 features
pdf('min-max-16f-89p-distribution.pdf')
plot(phsa16$A,phsa16$M,col=c('black','red')[as.numeric(as.factor(phsa16$recur))],pch=c(16,17)[as.numeric(as.factor(phsa16$recur))],main=" 16 Min-Max Features & 89 Patients",xlab="A=(d1+d2)/2",ylab="M=d2-d1")
legend("bottomright",pch=c(16,17),col=c('black','red'),c('Non Recurrence','Recurrence'))
dev.off()
# survival curve
p16fit<-survfit(Surv(phsa16$rfsurv,as.numeric(as.factor(phsa16$recur))-1)~as.factor(phsa16$kclst))
pdf('min-max-16f-89p-survival.pdf')
plot(p16fit,lty=1:2,col=c('black','red'), xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="16 Min-Max Features & 89 Patients")
legend("topright",lty=1:2,col=c('black','red'),c('cluster1','cluster2'))
text(x=90, y=0.5, paste0("p-value: ",round(pchisq(survdiff(Surv(phsa16$rfsurv,as.numeric(as.factor(phsa16$recur))-1)~as.factor(phsa16$kclst))$chisq,df=1,lower.tail=F),digits=8)))
dev.off()

# 6 min-max features
hsa6<-read.csv("hsabold6.txt",header=F,stringsAsFactor=F)
phsa6<-p60new[,hsa6$V1]
phsa6km<-kmeans(phsa6,2)
rownames(phsa6km$centers)<-c('c1','c2')
phsa6dist<-rbind(phsa6,phsa6km$centers)
phsa6dist<-as.matrix(dist(phsa6dist))
phsa6dist<-phsa6dist[1:89,90:91]
phsa6<-cbind(phsa6,phsa6dist)
phsa6$M<-phsa6$c2-phsa6$c1
phsa6$A<-(phsa6$c2+phsa6$c1)/2
phsa6$kclst<-ifelse(phsa6$M>0,1,2)
phsa6<-cbind(phsa6,p60new[,1:8])
write.csv(phsa6,file='phsa6.csv')

# plot the distribution of 89 patients by 6 features
pdf('min-max-6f-89p-distribution.pdf')
plot(phsa6$A,phsa6$M,col=c('black','red')[as.numeric(as.factor(phsa6$recur))],pch=c(16,17)[as.numeric(as.factor(phsa6$recur))],main=" 6 Min-Max Features & 89 Patients",xlab="A=(d1+d2)/2",ylab="M=d2-d1")
legend("right",pch=c(16,17),col=c('black','red'),c('Non Recurrence','Recurrence'))
dev.off()
# survival curve
p6fit<-survfit(Surv(phsa6$rfsurv,as.numeric(as.factor(phsa6$recur))-1)~as.factor(phsa6$kclst))
pdf('min-max-6f-89p-survival.pdf')
plot(p6fit,lty=1:2,col=c('black','red'), xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="6 Min-Max Features & 89 Patients")
legend("topright",lty=1:2,col=c('black','red'),c('cluster1','cluster2'))
text(x=90, y=0.5, paste0("p-value: ",round(pchisq(survdiff(Surv(phsa6$rfsurv,as.numeric(as.factor(phsa6$recur))-1)~as.factor(phsa6$kclst))$chisq,df=1,lower.tail=F),digits=8)))
dev.off()