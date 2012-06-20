d=read.table('avalues.txt',sep=" ")
dim(d)<-NULL
d<-as.numeric(d)
setEPS()
postscript("alpha.eps")
hist(d,c(-100,(-40:40)/20,100), xlim=c(-2,2), freq=TRUE)
dev.off()
