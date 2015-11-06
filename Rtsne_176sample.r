## Plot 176 samples distribution by Rtsne package ( Reducing 847 to 2 dimension)
# Require to load 'spall-dist-ptab.rdata' first.

library('Rtsne')
pdf('Rtsne-176sample-recur.pdf')

rtsne_sp<-Rtsne(as.matrix(spall[,14:860]))
plot(rtsne_sp$Y,col=c('black','red')[as.factor(spall$recur)],pch=c(16,17)[as.factor(spall$recur)],main='t-SNE distribution for 176 Samples', xlab='Dimension 1', ylab='Dimension 2')
legend("topright",pch=c(16,17),col=c('black','red'),c('No Recurrence','Recurrence'))
dev.off()