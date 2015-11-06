# ggplot2 and Rmisc packages are required. Dataset spall and dist12 are required.
library(ggplot2)
library(Rmisc)

# kmeans on 847 features.
km2<-kmeans(spall[,14:860],2)

# set the data in order 
t<-spall   
t$pclst<-dist12$pclst
t<-t[order(t$pclst),]
# set levels in order : good prognosis group, mixed group and poor prognosis group
pidN<-as.numeric(as.factor(t$pid))
lels<-c(unique(pidN[t$pclst==2]),unique(pidN[t$pclst==3]),unique(pidN[t$pclst==1]))
t$pidN<-factor(pidN,lels)
PatientGroup<-factor(t$pclst,labels=c('Poor prognosis','Good prognosis','Mixed'))
# 5 miRNA selected from mutual information.
mir5<-c(
"hsa-miR-122_st",
"hsa-miR-126_st",
"hsa-miR-22_st",
"hsa-miR-15a_st",
"hsa-miR-30a_st"
)
m5<-which(colnames(t) %in% mir5)

### for-loop  each individual miRBNA 
pm<-list()
for( i in 1:5) {

p<-qplot(factor(pidN,lels),t[,m5[i]],color=PatientGroup,xlab='Patients',ylab='miRNA Expression',main=gsub("hsa-|_st","",colnames(t)[m5[i]]))+theme(axis.text=element_text(size=5))+geom_hline(aes(yintercept=c(km2$center[1,m5[i]-13],km2$center[2,m5[i]-13])),color=c('red','black'), linetype='dashed')
p<-p+annotate("text", 5, km2$center[2,m5[i]-13], label = "Good prognosis center",size=3)+annotate("text", 15, km2$center[1,m5[i]-13], label = "Poor prognosis center",size=3)+coord_flip()
pm[[i]]<-p

pdf(paste0(gsub("hsa-|_st","",colnames(t)[m5[i]]),".pdf"))
plot(p)
dev.off()
}

### an alternative way to for-loop
p122<-qplot(factor(pidN,lels),t[,76],color=PatientGroup,xlab='Patients',ylab='miRNA Expression',main="miR-122")+theme(axis.text=element_text(size=5))+geom_hline(aes(yintercept=c(km2$center[1,63],km2$center[2,63])),color=c('red','black'), linetype='dashed')+annotate("text", 5, km2$center[2,63], label = "Good prognosis center",size=3)+annotate("text", 15, km2$center[1,63], label = "Poor prognosis center",size=3)+coord_flip()

p126<-qplot(factor(pidN,lels),t[,119],color=PatientGroup,xlab='Patients',ylab='miRNA Expression',main="miR-126")+theme(axis.text=element_text(size=5))+geom_hline(aes(yintercept=c(km2$center[1,106],km2$center[2,106])),color=c('red','black'), linetype="dashed")+annotate("text", 5, km2$center[2,106], label = "Good prognosis center",size=3)+annotate("text", 15, km2$center[1,106], label = "Poor prognosis center",size=3)+coord_flip()

p15a<-qplot(factor(pidN,lels),t[,226],color=PatientGroup,xlab='Patients',ylab='miRNA Expression',main="miR-15a")+theme(axis.text=element_text(size=5))+geom_hline(aes(yintercept=c(km2$center[1,213],km2$center[2,213])), colour=c("red","black"),linetype="dashed")+annotate("text", 5, km2$center[2,213], label = "Good prognosis center",size=3)+annotate("text", 15, km2$center[1,213], label = "Poor prognosis center",size=3)+coord_flip()

p22<-qplot(factor(pidN,lels),t[,336],color=PatientGroup,xlab='Patients',ylab='miRNA Expression',main="miR-22")+theme(axis.text=element_text(size=5))+geom_hline(aes(yintercept=c(km2$center[1,323],km2$center[2,323])), colour=c("red","black"), linetype="dashed")+annotate("text", 5, km2$center[2,323], label = "Good prognosis center",size=3)+annotate("text", 15, km2$center[1,323], label = "Poor prognosis center",size=3)+coord_flip()

p30a<-qplot(factor(pidN,lels),t[,384],color=PatientGroup,xlab='Patients',ylab='miRNA Expression',main="miR-30a")+theme(axis.text=element_text(size=5))+geom_hline(aes(yintercept=c(km2$center[1,371],km2$center[2,371])), colour=c("red","black" ),linetype="dashed")+annotate("text", 5, km2$center[2,371], label = "Good prognosis center",size=3)+annotate("text", 15, km2$center[1,371], label = "Poor prognosis center",size=3)+coord_flip()


### plot all 5 mRNAs together
pdf('5miRNA-all.pdf')
multiplot(p122,p126,p22,p15a,p30a,cols=2)
dev.off()

pdf('5miRNA-all-test.pdf')
multiplot(pm[[1]],pm[[2]],pm[[3]],pm[[4]],pm[[5]],cols=2)
dev.off()

