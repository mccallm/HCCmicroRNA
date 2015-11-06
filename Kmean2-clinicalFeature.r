
# This script applies kmeans (k=20) to 176 samples with 847 features.
# Require to load 'spall-dist-ptab.rdata' first.
# The data frame spall includes all samples info. Colname 14-860 are all miRNA features (raw data)

library(survival)  # survial package

km2<-kmeans(spall[,14:860],2)
rownames(km2$center)<-c('c1','c2')
dist12<-spall[,14:860]
dist12<-rbind(dist12,km2$center)
dist12<-as.matrix(dist(dist12)) # computing the distances between 176 samples and two centers.
dist12<-dist12[1:176,177:178]  # extract the distances from 176 samples to the two cluster centers
dist12<-as.data.frame(dist12)
dist12$M<-dist12$c2-dist12$c1 # the difference between two distances
dist12$A<-(dist12$c2+dist12$c1)/2 # average of two distances
dist12$kclst<-km2$cluster     # cluster assignment 
dist12$recur<-spall$recur       #  clinical feature 
dist12$batchid<-spall$batchid   #  clinical feature
dist12$vasc<-spall$vasc         #  clinical feature
dist12$focality<-spall$focality #  clinical feature
dist12$age<-spall$age           #  clinical feature
dist12$ntumor<-spall$ntumor     #  clinical feature
dist12$milan<-spall$milan       #  clinical feature
dist12$rfsurv<-spall$rfsurv      #  clinical feature
dist12$pid<-spall$pid
dist12$batchid<-spall$batchid   #  clinical feature

# Patients are categorized into three groups:  patients with samples all in cluster 1 
# patients with samples all in cluster 2 Patients with samples both in cluster 1 and 2

p1<-unique(dist12$pid[dist12$kclst==1])
pmix<-which( (dist12$kclst==2) & (dist12$pid %in% p1) )
pmix<-unique(dist12$pid[pmix])
p1<-p1[!p1%in%pmix]
p2<-unique(dist12$pid[dist12$kclst==2])
p2<-p2[!p2%in%pmix]
dist12$pclst<-NA
dist12$pclst[dist12$pid %in% p1]<-1
dist12$pclst[dist12$pid %in% p2]<-2
dist12$pclst[dist12$pid %in% pmix]<-3

#patient-level table
pdist12<-unique(dist12[,c('pid','pclst','rfsurv','recur','focality','ntumor','milan')])
pdist12$miclst<-NA # Milan and cluster features together
pdist12$miclst[which(pdist12$milan=='Outside'& pdist12$pclst==1)]<-'omc1' # ourside milan & cluster 1
pdist12$miclst[which(pdist12$milan=='Outside'& pdist12$pclst==2)]<-'omc2' # ourside milan & cluster 2
pdist12$miclst[which(pdist12$milan=='Outside'& pdist12$pclst==3)]<-'omc3' # ourside milan & mixed
pdist12$miclst[which(pdist12$milan=='Within')]<-'wm'  # within milan

# patient -level surival curves, miclust as groups
pdf('Kaplan-meier-milan-cluster-patients.pdf')
fitmclst<-survfit(Surv(pdist12$rfsurv,as.numeric(as.factor(pdist12$recur))-1)~as.factor(pdist12$miclst))
plot(fitmclst,lty=1:4,col=c('red','black','blue','green'),xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="All Patients (n=89)")
text(x=95, y=0.55, paste0("p-value: ",round(pchisq(survdiff(Surv(pdist12$rfsurv,as.numeric(as.factor(pdist12$recur))-1)~as.factor(pdist12$miclst))$chisq,df=2,lower.tail=F),digits=14)))
legend("topright",cex=0.6,lty=1:4,col=c('red','black','blue','green'),c('Outside Milan & Cluster1','Outside Milan & Cluster2','Outside Milan & Mixed','Within Milan'))
dev.off()

# patient-level surival curves, pclst as groups
poutside<-which(pdist12$milan=="Outside")
pallfit<-survfit(Surv(pdist12$rfsurv,as.numeric(as.factor(pdist12$recur))-1)~as.factor(pdist12$pclst))
moutfit<-survfit(Surv(pdist12$rfsurv[poutside],as.numeric(as.factor(pdist12$recur[poutside]))-1)~as.factor(pdist12$pclst[poutside]))
minfit<-survfit(Surv(pdist12$rfsurv[-poutside],as.numeric(as.factor(pdist12$recur[-poutside]))-1)~as.factor(pdist12$pclst[-poutside]))

pdf('Kaplan-meier-km-milan-paitents.pdf')
par(mfrow=c(1,3))
plot(pallfit,lty=1:3,col=c('red','black','blue'),xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="All Patients ")
legend("topright",cex=0.6,lty=1:3,col=c('red','black','blue'),c('Cluster1','Cluster2','Mixed'))
text(x=95, y=0.55, paste0("p-value: ",round(pchisq(survdiff(Surv(pdist12$rfsurv,as.numeric(as.factor(pdist12$recur))-1)~as.factor(pdist12$pclst))$chisq,df=2,lower.tail=F),digits=8)))
title(main='A', adj=0)
plot(moutfit,lty=1:3,col=c('red','black','blue'),xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Outside Milan Criteria")
legend("topright",cex=0.6,lty=1:3,col=c('red','black','blue'),c('Cluster1','Cluster2','Mixed'))
text(x=95, y=0.55, paste0("p-value: ",round(pchisq(survdiff(Surv(pdist12$rfsurv[poutside],as.numeric(as.factor(pdist12$recur[poutside]))-1)~as.factor(pdist12$pclst[poutside]))$chisq,df=2,lower.tail=F),digits=8)))
title(main='B', adj=0)
plot(minfit,lty=1:3,col=c('red','black','blue'),xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Within Milan Criteria")
legend("topright",cex=0.6,lty=1:2,col=c('red','black'),c('Cluster1','Cluster2'))
text(x=95, y=0.65, paste0("p-value: ",round(pchisq(survdiff(Surv(pdist12$rfsurv[-poutside],as.numeric(as.factor(pdist12$recur[-poutside]))-1)~as.factor(pdist12$pclst[-poutside]))$chisq,df=1,lower.tail=F),digits=8)))
title(main='C', adj=0)
dev.off()


# plot survival curves based on batch id. 
# b12: samples from batch 1 and 2 b34: sample from batch3 and 4
# b34pure: samples from batch 3 and 4 whose patients have not any other samples in batch 1 and 2
# b3: samples from batch 3 b4: samples from batch 4.
# b3pure: only samples in batch 3 whose owners don't have any samples in other batch group.
# b4pure: only samples in batch 4 whose owners don't have any samples in other batch group.

b12<-which(dist12$batchid %in% c(1,2))
b34<-which(dist12$batchid %in% c(3,4))
b3<-which(dist12$batchid==3)
b4<-which(dist12$batchid==4)
bmix<-which(dist12$pid %in% intersect(dist12$pid[b12],dist12$pid[b34]))
b34pure<-b34[!b34%in%bmix]
b3pure<-intersect(b34pure, which(dist12$batchid==3))
b4pure<-b34pure[!b34pure %in% b3]

# survival objects
fitb12<-survfit(Surv(dist12$rfsurv[b12],as.numeric(factor(dist12$recur[b12],levels=c("No Recurrence","Recurrence")))-1)~factor(dist12$pclst[b12],levels=c(1,2,3)))
fitb34<-survfit(Surv(dist12$rfsurv[b34],as.numeric(factor(dist12$recur[b34],levels=c("No Recurrence","Recurrence")))-1)~factor(dist12$pclst[b34],levels=c(1,2,3)))
fitb3<-survfit(Surv(dist12$rfsurv[b3],as.numeric(factor(dist12$recur[b3],levels=c("No Recurrence","Recurrence")))-1)~factor(dist12$pclst[b3],levels=c(1,2,3)))
fitb4<-survfit(Surv(dist12$rfsurv[b4],as.numeric(factor(dist12$recur[b4],levels=c("No Recurrence","Recurrence")))-1)~factor(dist12$pclst[b4],levels=c(1,2,3)))
fitb34pure<-survfit(Surv(dist12$rfsurv[b34pure],as.numeric(factor(dist12$recur[b34pure],levels=c("No Recurrence","Recurrence")))-1)~factor(dist12$pclst[b34pure],levels=c(1,2,3)))
fitb3pure<-survfit(Surv(dist12$rfsurv[b3pure],as.numeric(factor(dist12$recur[b3pure],levels=c("No Recurrence","Recurrence")))-1)~factor(dist12$pclst[b3pure],levels=c(1,2,3)))
fitb4pure<-survfit(Surv(dist12$rfsurv[b4pure],as.numeric(factor(dist12$recur[b4pure],levels=c("No Recurrence","Recurrence")))-1)~factor(dist12$pclst[b4pure],levels=c(1,2,3)))

# plot 1:  Kaplan-meier survival curves: batch 1&2 batch 3&4 (pure), batch 3 batch 4.
pdf('Kaplan-meier-kmclst-batchpure.pdf')
par(mfrow=c(2,2))
plot(fitb12,lty=1:3,col=c('red','black','blue'), xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Batch 1&2 samples")
legend("topright",cex=0.6,lty=1:3,col=c('red','black','blue'),c('cluster1','cluster2','mixed'))
text(x=100, y=0.45, paste0("p-value: ",round(pchisq(survdiff(Surv(dist12$rfsurv[b12],as.numeric(factor(dist12$recur[b12],levels=c("No Recurrence","Recurrence")))-1)~factor(dist12$pclst[b12],levels=c(1,2,3)))$chisq,df=2,lower.tail=F),digits=8)))
title(main='A', adj=0)
plot(fitb34pure,lty=1:3,col=c('red','black','blue'), xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Batch only 3&4 sample")
legend("bottomright",cex=0.6,lty=1:3,col=c('red','black','blue'),c('cluster1','cluster2','mixed'))
text(x=100, y=0.7, paste0("p-value: ",round(pchisq(survdiff(Surv(dist12$rfsurv[b34pure],as.numeric(factor(dist12$recur[b34pure],levels=c("No Recurrence","Recurrence")))-1)~factor(dist12$pclst[b34pure],levels=c(1,2,3)))$chisq,df=2,lower.tail=F),digits=10)))
title(main='B', adj=0)
plot(fitb3pure,conf=F,lty=2,col='black', xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Batch only 3 samples")
legend("bottomright",cex=0.6,lty=1:3,col=c('red','black','blue'),c('cluster1','cluster2','mixed'))
title(main='C', adj=0)
plot(fitb4pure,lty=1:3,col=c('red','black','blue'), xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Batch only 4 samples")
legend("bottomright",cex=0.6,lty=1:3,col=c('red','black','blue'),c('cluster1','cluster2','mixed'))
text(x=100, y=0.7, paste0("p-value: ",round(pchisq(survdiff(Surv(dist12$rfsurv[b4pure],as.numeric(factor(dist12$recur[b4pure],levels=c("No Recurrence","Recurrence")))-1)~factor(dist12$pclst[b4pure],levels=c(1,2,3)))$chisq,df=2,lower.tail=F),digits=10)))
title(main='D', adj=0)
dev.off()

# plot 2: Kaplan-meier survival curves: batch 1&2, batch 3&4,  batch 3 batch 4.
pdf('kaplan-meier-kmclst-batch.pdf')
par(mfrow=c(2,2))
plot(fitb12,lty=1:3,col=c('red','black','blue'), xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Batch 1&2 samples")
legend("topright",cex=0.6,lty=1:3,col=c('red','black','blue'),c('cluster1','cluster2','mixed'))
text(x=100, y=0.45, paste0("p-value: ",round(pchisq(survdiff(Surv(dist12$rfsurv[b12],as.numeric(factor(dist12$recur[b12],levels=c("No Recurrence","Recurrence")))-1)~factor(dist12$pclst[b12],levels=c(1,2,3)))$chisq,df=2,lower.tail=F),digits=8)))
title(main='A', adj=0)
plot(fitb34,lty=1:3,col=c('red','black','blue'), xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Batch 3&4 samples")
legend("bottomright",cex=0.6,lty=1:3,col=c('red','black','blue'),c('cluster1','cluster2','mixed'))
text(x=100, y=0.45, paste0("p-value: ",round(pchisq(survdiff(Surv(dist12$rfsurv[b34],as.numeric(factor(dist12$recur[b34],levels=c("No Recurrence","Recurrence")))-1)~factor(dist12$pclst[b34],levels=c(1,2,3)))$chisq,df=2,lower.tail=F),digits=12)))
title(main='B', adj=0)
plot(fitb3,conf=F,lty=2,col='black', xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Batch 3 samples")
legend("bottomright",cex=0.6,lty=1:3,col=c('red','black','blue'),c('cluster1','cluster2','mixed'))
title(main='C', adj=0)
plot(fitb4,lty=1:3,col=c('red','black','blue'), xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Batch 4 samples")
legend("bottomright",cex=0.6,lty=1:3,col=c('red','black','blue'),c('cluster1','cluster2','mixed'))
text(x=100, y=0.45, paste0("p-value: ",round(pchisq(survdiff(Surv(dist12$rfsurv[b4],as.numeric(factor(dist12$recur[b4],levels=c("No Recurrence","Recurrence")))-1)~factor(dist12$pclst[b4],levels=c(1,2,3)))$chisq,df=2,lower.tail=F),digits=10)))
title(main='D', adj=0)
dev.off()




# plot 176 samples vs recur(color), batchid(symbol)
pdf("km2_nonscale_176_dist_recur.pdf")
plot(dist12$A,dist12$M,pch=c(15,16,17,18)[as.numeric(as.factor(dist12$batchid))],col=c("black","red")[as.numeric(as.factor(dist12$recur))],xlab="A=(d1+d2)/2",ylab="M=d2-d1",main="176 samples distance distribution (KMeans=2)")
legend("topright",c("batch1","batch2","batch3","batch4"), pch=c(15,16,17,18))
legend("bottomright",c("No Recurrence","Recurrence"),col=c("black","red"),pch=16)
abline(h=0,lty=2) # Optional: separate two clusters 
dev.off()

# plot 176 samples vs recur(color), pclst(symbol): patient level cluster: samples from patients whose samples all from cluster 1
# samples from patients whose sample all from cluster 2 samples from patients whose sample both from cluster 1 and 2.
pdf("km2_nonscale_176_dist_cluster.pdf")
plot(dist12$A,dist12$M,pch=c(15,16,17)[as.factor(dist12$pclst)],col=c("black","red")[as.numeric(as.factor(dist12$recur))],xlab="A=(d1+d2)/2",ylab="M=d2-d1",main="176 Samples Distance Distribution (KMeans=2)")
#legend("topright",c("Cluster 1","Cluster 2","Mixed"), pch=c(15,16,17))
legend("topright",c("Cluster 1","Cluster 2","Mixed"), pch=c(0,1,2)) # optional
legend("bottomright",c("No Recurrence","Recurrence"),col=c("black","red"),pch=16)
abline(h=0,lty=2)
dev.off()

#plot 176 samples vs batchid (color), recur(symbol)
pdf("km2_nonscale_176_dist_batch.pdf")
plot(dist12$A,dist12$M,col=c("black","red","blue","green")[as.numeric(as.factor(dist12$batchid))],pch=c(16,17)[as.numeric(as.factor(dist12$recur))],xlab="A=(d1+d2)/2",ylab="M=d2-d1",main="176 samples distance distribution (KMeans=2)")
legend("topright",c("batch1","batch2","batch3","batch4"), col=c("black","red","blue","green"),pch=16)
legend("bottomright",c("No Recurrence","Recurrence"),pch=c(1,2))
abline(h=0,lty=2)
dev.off()

#plot 176 samples vs clinical feature focality
pdf("km2_nonscale_176_focality.pdf")
plot(dist12$A,dist12$M,col=c('black','red','blue')[as.numeric(as.factor(dist12$focality))],pch=c(16,17)[as.numeric(as.factor(dist12$recur))],main="Focality & Recur",xlab="A=(d1+d2)/2",ylab="M=d2-d1")
legend("bottomright",pch=16,col=c('black','red','blue'),c('Bifocal','Multifocal','Unifocal'),title="Focality")
legend("topright",pch=c(1,2),c('Non Recurrence','Recurrence'))
dev.off()

#plot 176 samples vs clinical feature milan
pdf("km2_nonscale_176_milan.pdf")
plot(dist12$A,dist12$M,col=c('black','red')[as.numeric(as.factor(dist12$milan))],pch=c(16,17)[as.numeric(as.factor(dist12$recur))],main="Milan & Recur",xlab="A=(d1+d2)/2",ylab="M=d2-d1")
legend("bottomright",pch=16,col=c('black','red'),c('Outside','Within'),title="Milan")
legend("topright",pch=c(1,2),c('No Recurrence','Recurrence'))
dev.off()

#plot 176 samples vs clinical feature Vasc
pdf("km2_nonscale_176_vasc.pdf")
plot(dist12$A,dist12$M,col=c('black','red')[as.numeric(as.factor(dist12$vasc))],pch=c(16,17)[as.numeric(as.factor(dist12$recur))],main="Vascular & Recur",xlab="A=(d1+d2)/2",ylab="M=d2-d1")
legend("bottomright",pch=16,col=c('black','red'),c('N','Y'),title="Vascular")
legend("topright",pch=c(1,2),c('No Recurrence','Recurrence'))
dev.off()

#plot 176 samples vs clinical feature number tumors
pdf("km2_nonscale_176_ntumor.pdf")
plot(dist12$A,dist12$M,col=c('black','red','blue','green','orange','purple')[as.numeric(as.factor(dist12$ntumor))],pch=c(16,17)[as.numeric(as.factor(dist12$recur))],main="Tumor Number & Recur",xlab="A=(d1+d2)/2",ylab="M=d2-d1")
legend("topright",pch=16,col=c('black','red','blue','green','orange','purple'),c('1','2','3','4','5','Multiple'),title="Tumor Number")
legend("bottomright",pch=c(1,2),c('No Recurrence','Recurrence'))
dev.off()

# plot 176 sample vs four clinical features
pdf('km2-noscale-176-4variates.pdf')
par(mfrow=c(2,2))
plot(dist12$A,dist12$M,col=c('black','red')[as.numeric(as.factor(dist12$vasc))],pch=c(16,17)[as.numeric(as.factor(dist12$recur))],main="Vascularization",xlab="A=(d1+d2)/2",ylab="M=d2-d1");
legend("bottomright",cex=0.6, pch=16,col=c('black','red'),c('N','Y'),title="Vascularization");
legend("topright",cex=0.6,pch=c(1,2),c('No Recurrence','Recurrence'))
title(main='A', adj=0)
plot(dist12$A,dist12$M,col=c('black','red','blue')[as.numeric(as.factor(dist12$focality))],pch=c(16,17)[as.numeric(as.factor(dist12$recur))],main="Focality",xlab="A=(d1+d2)/2",ylab="M=d2-d1");
legend("bottomright",cex=0.6,pch=16,col=c('black','red','blue'),c('Bifocal','Multifocal','Unifocal'),title="Focality");
legend("topright",cex=0.6,pch=c(1,2),c('Non Recurrence','Recurrence'))
title(main='B', adj=0)
plot(dist12$A,dist12$M,col=c('black','red','blue','green','orange','purple')[as.numeric(as.factor(dist12$ntumor))],pch=c(16,17)[as.numeric(as.factor(dist12$recur))],main="Tumor Number",xlab="A=(d1+d2)/2",ylab="M=d2-d1")
legend("topright",cex=0.6,pch=16,col=c('black','red','blue','green','orange','purple'),c('1','2','3','4','5','Multiple'),title="Tumor Number")
legend("bottomright",cex=0.6,pch=c(1,2),c('No Recurrence','Recurrence'))
title(main='C', adj=0)
plot(dist12$A,dist12$M,col=c('black','red')[as.numeric(as.factor(dist12$milan))],pch=c(16,17)[as.numeric(as.factor(dist12$recur))],main="Milan",xlab="A=(d1+d2)/2",ylab="M=d2-d1");
legend("bottomright",cex=0.6,pch=16,col=c('black','red'),c('Outside','Within'),title="Milan");
title(main='D', adj=0)
legend("topright",cex=0.6,pch=c(1,2),c('No Recurrence','Recurrence'));
dev.off()

# patient-level surival curves, pclst as groups: twp figures: outside and within milan
pdf('Kaplan-meier-km-milan-paitents.pdf')
par(mfrow=c(1,2))
plot(moutfit,lty=1:3,col=c('red','black','blue'),xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Outside Milan Criteria")
legend("topright",cex=0.6,lty=1:3,col=c('red','black','blue'),c('Cluster1','Cluster2','Mixed'))
text(x=95, y=0.55, paste0("p-value: ",round(pchisq(survdiff(Surv(pdist12$rfsurv[poutside],as.numeric(as.factor(pdist12$recur[poutside]))-1)~as.factor(pdist12$pclst[poutside]))$chisq,df=2,lower.tail=F),digits=8)))
title(main='A', adj=0)
plot(minfit,lty=1:3,col=c('red','black','blue'),xlab="Recurrence Free Survival (weeks)", ylab="Survival Proportion",main="Within Milan Criteria")
legend("topright",cex=0.6,lty=1:2,col=c('red','black'),c('Cluster1','Cluster2'))
text(x=95, y=0.65, paste0("p-value: ",round(pchisq(survdiff(Surv(pdist12$rfsurv[-poutside],as.numeric(as.factor(pdist12$recur[-poutside]))-1)~as.factor(pdist12$pclst[-poutside]))$chisq,df=1,lower.tail=F),digits=8)))
title(main='B', adj=0)
dev.off()