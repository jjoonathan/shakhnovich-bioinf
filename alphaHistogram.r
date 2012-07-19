Enterobacter=read.table('Escherichia-Enterobacter/alpha_values.tsv',sep="\t",header=TRUE)
dim(Enterobacter$Alpha)<-NULL
EnterobacterAlpha<-as.numeric(Enterobacter$Alpha)

Shewanella=read.table('Escherichia-Shewanella/alpha_values.tsv',sep="\t",header=TRUE)
dim(Shewanella$Alpha)<-NULL
ShewanellaAlpha<-as.numeric(Shewanella$Alpha)

Photorhabdus=read.table('Escherichia-Photorhabdus/alpha_values.tsv',sep="\t",header=TRUE)
dim(Photorhabdus$Alpha)<-NULL
PhotorhabdusAlpha<-as.numeric(Photorhabdus$Alpha)

Salmonella=read.table('Escherichia-Salmonella/alpha_values.tsv',sep="\t",header=TRUE)
dim(Salmonella$Alpha)<-NULL
SalmonellaAlpha<-as.numeric(Salmonella$Alpha)






#setEPS()
#pdf("alphaNOTC.pdf")
par(mfrow=c(2,2))
hist(EnterobacterAlpha, breaks=c(-100,(-40:40)/20,100), xlim=c(-2,2), freq=TRUE)
hist(ShewanellaAlpha, breaks=c(-100,(-40:40)/20,100), xlim=c(-2,2), freq=TRUE)
hist(PhotorhabdusAlpha, breaks=c(-100,(-40:40)/20,100), xlim=c(-2,2), freq=TRUE)
hist(SalmonellaAlpha, breaks=c(-100,(-40:40)/20,100), xlim=c(-2,2), freq=TRUE)
#dev.off()
