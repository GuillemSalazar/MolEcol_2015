library(spaa)
library(vegan)
library(ape)
library(geiger)
library(indicspecies)
library(picante)
library(MASS)
library(RCurl)
##########################
#### Data preparation ####
##########################

# Load phylogenetic tree
x<-getURL("https://raw.githubusercontent.com/GuillemSalazar/MolEcol_2015/master/FinalTree_Salazar_etal_2015_Molecol.nwk")
tree<-read.tree(text=x)

# Load OTU table and separate taxonomy from counts
x<-getURL("https://raw.githubusercontent.com/GuillemSalazar/MolEcol_2015/master/OTUtable_Salazar_etal_2015_Molecol.txt")
OTUtab.complete.ss<-read.table(text=x,sep="\t",row.names=1,header=TRUE,comment.char="@")
OTUtab.ss<-OTUtab.complete.ss[,1:60]
otus.null<-which(rowSums(OTUtab.ss)>0)
OTUtab.ss<-OTUtab.ss[otus.null,]
OTUtab.complete.ss<-OTUtab.complete.ss[otus.null,]

# Remove from the tree OTUs not present in the table
a.extreure<-which(tree$tip.label %in% rownames(OTUtab.ss)==F)

# Remove from the tree OTUs suspected mitocondrial sequences
a.extreure<-c(a.extreure,which(tree$tip.label %in% c(1499,3618,2090,5140,2389,3602,4420,182)))
tree<-drop.tip(tree,a.extreure)
tree<-multi2di(tree)

# Order table based on tree order
OTUtab.complete.ss<-OTUtab.complete.ss[match(tree$tip.label,rownames(OTUtab.complete.ss)),]
OTUtab.ss<-OTUtab.ss[match(tree$tip.label,rownames(OTUtab.ss)),]
all(rownames(OTUtab.complete.ss)==tree$tip.labels)
all(rownames(OTUtab.ss)==tree$tip.labels)

# Convert RDP taxonomy to character
OTUtab.complete.ss$RDP_Domain<-as.character(OTUtab.complete.ss$RDP_Domain)
OTUtab.complete.ss$RDP_Phylum<-as.character(OTUtab.complete.ss$RDP_Phylum)
OTUtab.complete.ss$RDP_Class<-as.character(OTUtab.complete.ss$RDP_Class)
OTUtab.complete.ss$RDP_Order<-as.character(OTUtab.complete.ss$RDP_Order)
OTUtab.complete.ss$RDP_Family<-as.character(OTUtab.complete.ss$RDP_Family)
OTUtab.complete.ss$RDP_Genus<-as.character(OTUtab.complete.ss$RDP_Genus)

# Label as unclassified assignations with confidence values <90
OTUtab.complete.ss$RDP_Domain[which(OTUtab.complete.ss$RDP_Domain_conf<90)]<-"Unclassified"
OTUtab.complete.ss$RDP_Phylum[which(OTUtab.complete.ss$RDP_Phylum_conf<90)]<-"Unclassified"
OTUtab.complete.ss$RDP_Class[which(OTUtab.complete.ss$RDP_Class_conf<90)]<-"Unclassified"
OTUtab.complete.ss$RDP_Order[which(OTUtab.complete.ss$RDP_Order_conf<90)]<-"Unclassified"
OTUtab.complete.ss$RDP_Family[which(OTUtab.complete.ss$RDP_Family_conf<90)]<-"Unclassified"
OTUtab.complete.ss$RDP_Genus[which(OTUtab.complete.ss$RDP_Genus_conf<90)]<-"Unclassified"

# Convert to factor
OTUtab.complete.ss$RDP_Domain<-factor(OTUtab.complete.ss$RDP_Domain)
OTUtab.complete.ss$RDP_Phylum<-factor(OTUtab.complete.ss$RDP_Phylum)
OTUtab.complete.ss$RDP_Class<-factor(OTUtab.complete.ss$RDP_Class)
OTUtab.complete.ss$RDP_Order<-factor(OTUtab.complete.ss$RDP_Order)
OTUtab.complete.ss$RDP_Family<-factor(OTUtab.complete.ss$RDP_Family)
OTUtab.complete.ss$RDP_Genus<-factor(OTUtab.complete.ss$RDP_Genus)

# Load auxiliary data
x<-getURL("https://raw.githubusercontent.com/GuillemSalazar/MolEcol_2015/master/Metadata_Salazar_etal_2015_Molecol.txt")
metad<-read.table(text=x,sep="\t",header=TRUE)
metad$filtersize<-as.factor(metad$filtersize)
metad08<-metad[which(metad$filtersize==0.8),]
metad02<-metad[which(metad$filtersize==0.2),]
OTUtab.ss<-t(OTUtab.ss)

# Classify Phylum based on SILVA. Divide Proteobacteria by Classes.
Phylums<-NULL
Class<-NULL
Order<-NULL
separat<-strsplit(as.character(OTUtab.complete.ss$SILVA_Taxonomy),";")
for (i in 1:length(separat)){
	if (length(separat[[i]])>1) Phylums<-c(Phylums,separat[[i]][[2]]) else Phylums<-c(Phylums,"Unclassified")
	if (length(separat[[i]])>2) Class<-c(Class,separat[[i]][3]) else Class<-c(Class,"Unclassified")
	if (length(separat[[i]])>3) Order<-c(Order,separat[[i]][4]) else Order<-c(Order,"Unclassified")}

Phylums[which(Phylums=="Proteobacteria")]<-as.character(Class[which(Phylums=="Proteobacteria")])
Phylums<-factor(Phylums)
Phylums<-as.factor(Phylums)
Class<-as.factor(Class)
Order<-as.factor(Order)


#################################
###   RICHNESS / DIVERSITY   ####
#################################

PD<-pd(OTUtab.ss,tree,include.root=F)
otu.estim<-estimateR(OTUtab.ss)
shan<-diversity(OTUtab.ss)
evenness<-shan/log(otu.estim[1,])
MNTD<-mntd(OTUtab.ss,cophenetic(tree))
cbind(t(otu.estim[c(1,2),]),PD=PD$PD,shan)

pdf("Fig1B.pdf",height=6,width=8)
layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
wilcox.test(otu.estim[1,]~metad$filtersize,paired=F)
plot(metad$filtersize,otu.estim[1,],ylab="Richness (Nº of OTUs)",col=c("blue","red"))
wilcox.test(otu.estim[2,]~metad$filtersize,paired=F)
plot(metad$filtersize,otu.estim[2,],ylab="Chao estimator",col=c("blue","red"))
wilcox.test(shan~metad$filtersize,paired=F)
plot(metad$filtersize,shan,ylab="Shannon index",col=c("blue","red"))

wilcox.test(PD$PD~metad$filtersize,paired=F)
plot(metad$filtersize,PD$PD,ylab="Phylogenetic diversity",col=c("blue","red"))
wilcox.test(PD$PD/otu.estim[1,]~metad$filtersize,paired=F)
plot(metad$filtersize,PD$PD/otu.estim[1,],ylab="Phylogenetic diversity / Richness",col=c("blue","red"))
wilcox.test(MNTD~metad$filtersize,paired=F)
plot(metad$filtersize,MNTD,ylab="MNTD",col=c("blue","red"))
dev.off()


#####################
### GENERAL NMDS ####
#####################

# Comparison between MNTDbeta and Bray-Curtis
dist.phylo.MNTD<-comdistnt(OTUtab.ss,cophenetic(tree))
dist.taxo<-vegdist(OTUtab.ss)
dinsfora<-matrix(NA,nrow=dim(metad)[1],ncol=dim(metad)[1])
for (i in 1:dim(dinsfora)[1]){
  for (j in 1:dim(dinsfora)[2]){
    if (metad$filtersize[i]==0.2 & metad$filtersize[j]==0.2) dinsfora[i,j]<-"within FL"
    if (metad$filtersize[i]==0.8 & metad$filtersize[j]==0.8) dinsfora[i,j]<-"within PA"
    if (metad$filtersize[i]==0.2 & metad$filtersize[j]==0.8) dinsfora[i,j]<-"between FL-PA"
    if (metad$filtersize[i]==0.8 & metad$filtersize[j]==0.2) dinsfora[i,j]<-"between FL-PA"
  }
}

mat.taxo<-as.matrix(dist.taxo)
mat.phylo.MNTD<-as.matrix(dist.phylo.MNTD)
diag(mat.taxo)<-NA
diag(mat.phylo.MNTD)<-NA
diag(dinsfora)<-NA

pdf("FigS3.pdf")
layout(matrix(c(2,2,2,4,1,1,1,3,1,1,1,3,1,1,1,3), 4, 4, byrow = TRUE))
par(mar=c(5,5,0,0))
plot(mat.taxo,mat.phylo.MNTD,col=c(rgb(0,0,0,0.2),rgb(0,0,1,0.4),rgb(1,0,0,0.4))[as.factor(dinsfora)],pch=19,cex=0.8,xlab="Bray-Curtis dissimilarity",ylab="betaMNTD dissimilarity")
points(na.omit(mat.taxo[dinsfora=="within FL"]),predict(lm(mat.phylo.MNTD[dinsfora=="within FL"]~mat.taxo[dinsfora=="within FL"])),col="blue",type="l")
points(na.omit(mat.taxo[dinsfora=="within PA"]),predict(lm(mat.phylo.MNTD[dinsfora=="within PA"]~mat.taxo[dinsfora=="within PA"])),col="red",type="l")

mantel(dist.taxo,dist.phylo.MNTD)
par(mar=c(0,5,5,0),xaxt="n")
hist(mat.taxo[dinsfora=="between FL-PA"],xlim=range(mat.taxo,na.rm=T),col=rgb(0,0,0,0.4),nclass=30,main="",plot=T)
hist(mat.taxo[dinsfora=="within PA"],xlim=range(mat.taxo,na.rm=T),col=rgb(1,0,0,0.4),add=T,nclass=30,plot=T)
hist(mat.taxo[dinsfora=="within FL"],xlim=range(mat.taxo,na.rm=T),col=rgb(0,0,1,0.4),add=T,nclass=30,plot=T)

par(mar=c(5,0,0,5),yaxt="n",xaxt="s")
A<-hist(mat.phylo.MNTD[dinsfora=="between FL-PA"],nclass=30,plot=F)
B<-hist(mat.phylo.MNTD[dinsfora=="within PA"],nclass=30,plot=F)
C<-hist(mat.phylo.MNTD[dinsfora=="within FL"],nclass=30,plot=F)
plot(NULL,type="n",xlim=c(0,max(c(A$counts,B$counts,C$counts))),ylim=range(mat.phylo.MNTD,na.rm=T),xlab="Frequency",frame.plot=F,ylab="")
rect(0,A$breaks[1:(length(A$breaks)-1)],A$counts,A$breaks[2:length(A$breaks)],col=rgb(0,0,0,0.4))
rect(0,B$breaks[1:(length(B$breaks)-1)],B$counts,B$breaks[2:length(B$breaks)],col=rgb(1,0,0,0.4))
rect(0,C$breaks[1:(length(C$breaks)-1)],C$counts,C$breaks[2:length(C$breaks)],col=rgb(0,0,1,0.4))
dev.off()

# NMDS
nmds<-metaMDS(OTUtab.ss,autotransform=F,expand=F)

pdf("Fig1A.pdf")
plot(nmds$species,col=metad$filtersize.col,xlab="NMDS 1",ylab="NMDS 2",pch=19,type="n",ylim=c(-1.2,1.2),xlim=c(-1.2,1.2))
ordihull(nmds,metad$filtersize,draw="polygon",col=rgb(0,0,1,0.0001),show.groups=0.2)
ordihull(nmds,metad$filtersize,draw="polygon",col=rgb(1,0,0,0.0001),show.groups=0.8)
points(nmds$points,col=c("blue","red")[metad$filtersize],xlab="NMDS 1",ylab="NMDS 2",pch=as.numeric(metad$filtersize)+14,cex=1.2)
points(nmds$points,col="black",xlab="NMDS 1",ylab="NMDS 2",pch=c(22,21)[as.numeric(metad$filtersize)],cex=1.2)
text(nmds$points, labels = metad$station,cex=0.5,adj=c(0,2))
legend("bottomright",levels(metad$filtersize),col=c("blue","red"),pch=15:16,inset=0.01,cex=1)
dev.off()

# Multivariate Permutational MANOVA
adonis(vegdist(OTUtab.ss)~metad$filtersize,nperm=10000)

# Bray-Curtis partition within station and size-fraction
pdf("FigS2.pdf")
BC<-as.matrix(vegdist(OTUtab.ss))
diag(BC)<-NA
FL<-BC[metad$filtersize==0.2,metad$filtersize==0.2]
PA<-BC[metad$filtersize==0.8,metad$filtersize==0.8]
same.st<-NULL
for (i in unique(metad$station)){
  same.st<-c(same.st,BC[metad$station==i,metad$station==i])}
tot.junt<-c(FL,PA,same.st)
tot.junt.group<-c(rep("within FL",length(FL)),rep("within PA",length(PA)),rep("within station",length(same.st)))
plot(as.factor(tot.junt.group),tot.junt,col=c("blue","red","gray"),ylim=c(0,1),ylab="Bray-Curtis distance")
abline(h=mean(BC,na.rm=T),lty=2)
dev.off()

#######################
## Compute PAN-index ##
#######################
coph.dist<-cophenetic(tree)
diag(coph.dist)<-NA


sf.wm<-NULL
for (i in 1:dim(OTUtab.ss)[2]){
  sf.wm<-c(sf.wm,weighted.mean(decostand(as.numeric(metad$filtersize),"range"),OTUtab.ss[,i]))}

#################################
## Comparison to randomization ##  Remember to run definitely with 1000 iterations (nperm=1000) !!!!
#################################

nperm=10
sf.real<-sf.wm
mat<-OTUtab.ss
sf.wm.null.accum<-NULL
counts.null<-matrix(NA,ncol=3,nrow=nperm)
for (iter in 1:nperm){
  mat.null<-permatswap(mat, "quasiswap",times=1)$perm[[1]]
  rownames(mat.null)<-rownames(mat)
  colnames(mat.null)<-colnames(mat)
  sf.wm.null<-NULL
  for (i in 1:dim(mat.null)[2]){
    sf.wm.null<-c(sf.wm.null,weighted.mean(decostand(as.numeric(metad$filtersize),"range"),mat.null[,i]))
  }
  sf.wm.null.accum<-c(sf.wm.null.accum,sf.wm.null)
  counts.null[iter,]<-c(length(which(sf.wm.null<0.2)),length(which(sf.wm.null>0.8)),length(which(sf.wm.null>=0.2 & sf.wm.null<=0.8)))
}


sf.null.range<-apply(matrix(sf.wm.null.accum,ncol=dim(mat)[2],byrow=T),2,quantile,probs=c(0.025,0.975))
sf.signif<-NULL
state<-NULL
for (i in 1:dim(mat)[2]){
  if (sf.real[i]>sf.null.range[2,i]) state<-"SIGN. PA"
  if (sf.real[i]<sf.null.range[1,i]) state<-"SIGN. FL"
  if (sf.real[i]<=sf.null.range[2,i] & sf.real[i]>=sf.null.range[1,i]) state<-"NON SIGN."
  sf.signif<-c(sf.signif,state)
}

table(sf.signif)
aggregate(colSums(mat),list(sf.signif),sum)



pdf(file="Fig2.pdf",width=6,height=10)
layout(matrix(c(1,2), 2, 1, byrow = TRUE))
hist.null<-hist(sf.wm.null.accum,nclass=30,ylim=c(0,5),freq=F,main="",xlab="Particle-association niche value",col="gray")
points(density(sf.wm.null.accum,bw=1/30),type="l",lwd=1.5)
hist.real<-hist(sf.real,nclass=30,add=T,freq=F,density=30,angle=-45)
points(density(sf.real,bw=1/30),type="l",lwd=1.5,lty=4)
legend("top",c("Null expectation [1000 iterations]","Real data"),lty=c(1,4),col=c("black","black"),cex=0.7)

counts.real<-c(length(which(sf.wm<0.2)),length(which(sf.wm>0.8)),length(which(sf.wm>=0.2 & sf.wm<=0.8)))

ks.test(sf.wm.null.accum,sf.wm)

## Excluding OTUs with less or equal to 10 reads ##
mat<-OTUtab.ss[,which(colSums(OTUtab.ss)>10)]
sf.real.red<-sf.wm[which(colSums(OTUtab.ss)>10)]
dim(mat)
sf.wm.null.accum<-NULL
counts.null.red<-matrix(NA,ncol=3,nrow=nperm)
for (iter in 1:nperm){
  mat.null<-permatswap(mat, "quasiswap",times=1)$perm[[1]]
  rownames(mat.null)<-rownames(mat)
  colnames(mat.null)<-colnames(mat)
  sf.wm.null<-NULL
  for (i in 1:dim(mat.null)[2]){
    sf.wm.null<-c(sf.wm.null,weighted.mean(decostand(as.numeric(metad$filtersize),"range"),mat.null[,i]))
  }
  sf.wm.null.accum<-c(sf.wm.null.accum,sf.wm.null)
  counts.null.red[iter,]<-c(length(which(sf.wm.null<0.2)),length(which(sf.wm.null>0.8)),length(which(sf.wm.null>=0.2 & sf.wm.null<=0.8)))
}


sf.null.range<-apply(matrix(sf.wm.null.accum,ncol=dim(mat)[2],byrow=T),2,quantile,probs=c(0.025,0.975))
sf.signif<-NULL
state<-NULL
for (i in 1:dim(mat)[2]){
  if (sf.real.red[i]>sf.null.range[2,i]) state<-"SIGN. PA"
  if (sf.real.red[i]<sf.null.range[1,i]) state<-"SIGN. FL"
  if (sf.real.red[i]<=sf.null.range[2,i] & sf.real.red[i]>=sf.null.range[1,i]) state<-"NON SIGN."
  sf.signif<-c(sf.signif,state)
}

table(sf.signif)
aggregate(colSums(mat),list(sf.signif),sum)
ks.test(sf.wm.null.accum,sf.real.red)

hist.null<-hist(sf.wm.null.accum,nclass=30,ylim=c(0,5),freq=F,main="",xlab="Particle-association niche value",col="gray")
points(density(sf.wm.null.accum,bw=1/30),type="l",lwd=1.5)
hist.real<-hist(sf.real.red,nclass=30,add=T,freq=F,density=30,angle=-45)
points(density(sf.real.red,bw=1/30),type="l",lwd=1.5,lty=4)
legend("top",c("Null expectation [1000 iterations]","Real data"),lty=c(1,4),col=c("black","black"),cex=0.7)
dev.off()

counts.real.red<-c(length(which(sf.real.red<0.2)),length(which(sf.real.red>0.8)),length(which(sf.real.red>=0.2 & sf.real.red<=0.8)))
counts.complete<-rbind(round(apply(counts.null,2,"range"),3),round(counts.real,3))
colnames(counts.complete)<-c("< 0.2","> 0.8","0.2 - 0.8")
rownames(counts.complete)<-c("Null Max.","Null Min.","Real value")
print(counts.complete)

counts.red<-rbind(round(apply(counts.null.red,2,"range"),3),round(counts.real.red,3))
colnames(counts.red)<-c("< 0.2","> 0.8","0.2 - 0.8")
rownames(counts.red)<-c("Null Max.","Null Min.","Real value")
print(counts.red)
############################################################


# Create reduced dataset (excluding OTUs with less than 10 reads)
tree.reduced<-drop.tip(tree,which(colSums(OTUtab.ss)<=10))
all(tree.reduced$tip.label==colnames(mat))
OTUtab.ss.reduced<-OTUtab.ss[,which(colSums(OTUtab.ss)>10)]
OTUtab.complete.ss.reduced<-OTUtab.complete.ss[which(colSums(OTUtab.ss)>10),]

#########################
## Niche conservatism ##
#########################
coph.dist<-cophenetic(tree.reduced)
diag(coph.dist)<-NA

# Compute mean PAN-values for 0.01 phylogenetic distance bins 
sf.dist<-as.matrix(dist(sf.real.red))
diag(sf.dist)<-NA
breaks=seq(0,max(coph.dist,na.rm=T),by=0.01)
mitjana.sf<-NULL
middle<-NULL
for (i in 1:(length(breaks)-1)){
  middle<-c(middle,mean(breaks[i:(i+1)]))
  a<-which(coph.dist>=breaks[i] & coph.dist<breaks[i+1])
  mitjana.sf<-c(mitjana.sf,mean(sf.dist[a]))
}

# Collect phylogenetic distance within taxa
DISTANCIA<-NULL
TAXA<-NULL

# Domain
distancia<-NULL
taxa<-NULL
for(i in 1:nlevels(OTUtab.complete.ss.reduced$RDP_Domain)){
  quins<-which(OTUtab.complete.ss.reduced$RDP_Domain==levels(OTUtab.complete.ss.reduced$RDP_Domain)[i])
  a<-c(coph.dist[quins,quins])
  distancia<-c(distancia,a)
  b<-rep(levels(OTUtab.complete.ss.reduced$RDP_Domain)[i],length(a))
  taxa<-c(taxa,b)
}
distancia<-distancia[-c(which(taxa=="Unclassified"))]
DISTANCIA<-c(DISTANCIA,distancia)
TAXA<-c(TAXA,rep("Same Domain",length(distancia)))

# Phylum
distancia<-NULL
taxa<-NULL
for(i in 1:nlevels(OTUtab.complete.ss.reduced$RDP_Phylum)){
  quins<-which(OTUtab.complete.ss.reduced$RDP_Phylum==levels(OTUtab.complete.ss.reduced$RDP_Phylum)[i])
  a<-c(coph.dist[quins,quins])
  distancia<-c(distancia,a)
  b<-rep(levels(OTUtab.complete.ss.reduced$RDP_Phylum)[i],length(a))
  taxa<-c(taxa,b)
}
taxa.phylum<-as.factor(taxa)
distancia.phylum<-distancia
distancia<-distancia[-c(which(taxa=="Unclassified"))]
DISTANCIA<-c(DISTANCIA,distancia)
TAXA<-c(TAXA,rep("Same Phylum",length(distancia)))


# Class
distancia<-NULL
taxa<-NULL
for(i in 1:nlevels(OTUtab.complete.ss.reduced$RDP_Class)){
  quins<-which(OTUtab.complete.ss.reduced$RDP_Class==levels(OTUtab.complete.ss.reduced$RDP_Class)[i])
  a<-c(coph.dist[quins,quins])
  distancia<-c(distancia,a)
  b<-rep(levels(OTUtab.complete.ss.reduced$RDP_Class)[i],length(a))
  taxa<-c(taxa,b)
}
distancia<-distancia[-c(which(taxa=="Unclassified"))]
DISTANCIA<-c(DISTANCIA,distancia)
TAXA<-c(TAXA,rep("Same Class",length(distancia)))

# Order
distancia<-NULL
taxa<-NULL
for(i in 1:nlevels(OTUtab.complete.ss.reduced$RDP_Order)){
  quins<-which(OTUtab.complete.ss.reduced$RDP_Order==levels(OTUtab.complete.ss.reduced$RDP_Order)[i])
  a<-c(coph.dist[quins,quins])
  distancia<-c(distancia,a)
  b<-rep(levels(OTUtab.complete.ss.reduced$RDP_Order)[i],length(a))
  taxa<-c(taxa,b)
}
distancia<-distancia[-c(which(taxa=="Unclassified"))]
DISTANCIA<-c(DISTANCIA,distancia)
TAXA<-c(TAXA,rep("Same Order",length(distancia)))

# Family
distancia<-NULL
taxa<-NULL
for(i in 1:nlevels(OTUtab.complete.ss.reduced$RDP_Family)){
  quins<-which(OTUtab.complete.ss.reduced$RDP_Family==levels(OTUtab.complete.ss.reduced$RDP_Family)[i])
  a<-c(coph.dist[quins,quins])
  distancia<-c(distancia,a)
  b<-rep(levels(OTUtab.complete.ss.reduced$RDP_Family)[i],length(a))
  taxa<-c(taxa,b)
}
distancia<-distancia[-c(which(taxa=="Unclassified"))]
DISTANCIA<-c(DISTANCIA,distancia)
TAXA<-c(TAXA,rep("Same Family",length(distancia)))

# Genus
distancia<-NULL
taxa<-NULL
for(i in 1:nlevels(OTUtab.complete.ss.reduced$RDP_Genus)){
  quins<-which(OTUtab.complete.ss.reduced$RDP_Genus==levels(OTUtab.complete.ss.reduced$RDP_Genus)[i])
  a<-c(coph.dist[quins,quins])
  distancia<-c(distancia,a)
  b<-rep(levels(OTUtab.complete.ss.reduced$RDP_Genus)[i],length(a))
  taxa<-c(taxa,b)
}
distancia<-distancia[-c(which(taxa=="Unclassified"))]
DISTANCIA<-c(DISTANCIA,distancia)
TAXA<-c(TAXA,rep("Same Genus",length(distancia)))
TAXA<-factor(TAXA,ordered=T,levels=c("Same Domain","Same Phylum","Same Class","Same Order","Same Family","Same Genus"))

# Make the plot
pdf(file="Fig3.pdf",height=8,width=5)
layout(matrix(c(1,1,1,1,2,3,4,5,6,7,8), 11, 1, byrow = TRUE))
par(mar=c(4,10,4,4))
plot(middle,mitjana.sf,type="b",cex=0.4,ylim=c(0,1),xlab="Between-OTU Phylogenetic Distance",ylab="Between-OTU PANV Difference",pch=19)
abline(v=seq(0,4,by=0.2),lty=2,col="gray")
abline(h=mean(mitjana.sf,na.rm=T),lty=2,col="gray")

par(mar=c(1,10,1,4),cex.axis=0.7,las=1)
hist(DISTANCIA[TAXA==levels(TAXA)[1]],nclass=100,xlim=c(0,max(coph.dist,na.rm=T)),main=NA,col="black",xlab="")
hist(DISTANCIA[TAXA==levels(TAXA)[2]],nclass=100,xlim=c(0,max(coph.dist,na.rm=T)),main=NA,col="black",xlab="")
hist(DISTANCIA[TAXA==levels(TAXA)[3]],nclass=100,xlim=c(0,max(coph.dist,na.rm=T)),main=NA,col="black",xlab="")
hist(DISTANCIA[TAXA==levels(TAXA)[4]],nclass=100,xlim=c(0,max(coph.dist,na.rm=T)),main=NA,col="black",xlab="")
hist(DISTANCIA[TAXA==levels(TAXA)[5]],nclass=100,xlim=c(0,max(coph.dist,na.rm=T)),main=NA,col="black",xlab="")
hist(DISTANCIA[TAXA==levels(TAXA)[6]],nclass=100,xlim=c(0,max(coph.dist,na.rm=T)),main=NA,col="black",xlab="Between-OTU Phylogenetic Distance")
dev.off()


##############################
## Phylogenetic signal test ##
##############################

# Takes a while!
tree.0<-rescale(tree.reduced,"lambda",0)
trait<-sf.real.red
names(trait)<-colnames(OTUtab.ss.reduced)
res<-fitContinuous(tree.reduced,trait,model="lambda",control=list(niter=100))
print(res)

res.null<-fitContinuous(tree.0,trait,model="lambda",control=list(niter=100))
res.null

# Likelihood rate test
LR<-2*(logLik(res)-logLik(res.null))
print(pchisq(LR, df = res$opt$k - res.null$opt$k, lower.tail = FALSE))
print(LR)

#################################
## Indicator OTUs for fraction ##
#################################

indval.filter<-multipatt(data.frame(OTUtab.ss),metad$filtersize,duleg=FALSE,control=permControl(nperm=1000))
indval.filter$sign$p.value<-p.adjust(indval.filter$sign$p.value,method="fdr")

summary(indval.filter,indvalcomp=TRUE,At=0.8,Bt=0.8)

# Filter results (A and B grater than 0.8)
sign.down<-data.frame(OTU.ID=rownames(indval.filter$str)[which(indval.filter$sign$s.0.2==1)],A=indval.filter$A[which(indval.filter$sign$s.0.2==1),1],B=indval.filter$B[which(indval.filter$sign$s.0.2==1),1],stat=indval.filter$sign$stat[which(indval.filter$sign$s.0.2==1)],p.value=indval.filter$sign$p.value[which(indval.filter$sign$s.0.2==1)])
sign.down<-cbind(sign.down,OTUtab.complete.ss$SILVA_Taxonomy[match(sign.down$OTU.ID,paste("X",rownames(OTUtab.complete.ss),sep=""))])
sign.down<-sign.down[which(sign.down$p.value<=0.05 & sign.down$A>=0.8 & sign.down$B>=0.8),]
sign.down<-sign.down[order(sign.down$stat,decreasing=T),]

sign.upp<-data.frame(OTU.ID=rownames(indval.filter$str)[which(indval.filter$sign$s.0.8==1)],A=indval.filter$A[which(indval.filter$sign$s.0.8==1),2],B=indval.filter$B[which(indval.filter$sign$s.0.8==1),2],stat=indval.filter$sign$stat[which(indval.filter$sign$s.0.8==1)],p.value=indval.filter$sign$p.value[which(indval.filter$sign$s.0.8==1)])
sign.upp<-cbind(sign.upp,OTUtab.complete.ss$SILVA_Taxonomy[match(sign.upp$OTU.ID,paste("X",rownames(OTUtab.complete.ss),sep=""))])
sign.upp<-sign.upp[which(sign.upp$p.value<=0.05 & sign.upp$A>=0.8 & sign.upp$B>=0.8),]
sign.upp<-sign.upp[order(sign.upp$stat,decreasing=T),]

write.table(sign.upp,file="IndicatorPA.txt",sep="\t",quote=F,row.names=F,col.names=T)
write.table(sign.down,file="IndicatorFL.txt",sep="\t",quote=F,row.names=F,col.names=T)

#################################
##     Phylogeny with info     ##
#################################

pdf(file="Fig4.pdf",width=10,height=8)
par(mar=c(2,2,7,0))
plot(tree.reduced,type="phylogram",show.tip.label=F,x.lim=13,y.lim=c(-500,2500))
sf.wm.mod<-sf.real.red-0.5
logab<-log10(colSums(OTUtab.ss.reduced))
logab<-logab/max(logab)
exemples<-c("Marine Group II","Marine Group I;","Planctomycetaceae","Phycisphaeraceae","OM190","Arctic97B-4","Opitutae","Myxococcales","Desulfuromonadales","OM27","SAR324","Alteromonadaceae","SAR86","Sinobacteraceae","AEGEAN-169 marine group","DB1-14","SAR116","Sphingomonadales","Saprospiraceae","Flammeovirgaceae","Flavobacteriaceae","SAR406","DA023","BD2-11","AT425-EubC11","OCS155 marine group","Sva0996 marine group","Microcoleus","SHA-109","Anaerolineaceae","SAR202")
# GR-WP33-58 clade because paper:Metagenomic analysis of mesopelagic Antarctic plankton reveals a novel deltaproteobacterial group.
for (i in seq(1,dim(OTUtab.ss.reduced)[2],by=100)){lines(c(4.2,4+(0.2*length(exemples))),c(i,i),col="gray",lty=1)}
j<-0
for (i in 1:length(exemples)){
  j<-j+1
  quins<-grep(as.character(unique(exemples)[i]),OTUtab.complete.ss.reduced$SILVA_Taxonomy)
  coloret<-mean(sf.real.red[quins])
  lines(c(4+(0.2*j),4+(0.2*j)),c(1,length(sf.real.red)),col="gray",lty=1)
  points(rep(4+(0.2*j),length(quins)),quins,cex=0.2,col=rgb(coloret,0,1-coloret,0.5))
  text(4+(0.2*j)-0.18,length(sf.real.red)+100,unique(exemples)[i],cex=0.6,srt=90,pos=4)
}
lines(c(11,11),c(1,length(sf.real.red)),col="gray")
for (i in 1:length(sf.real.red)){lines(c(11,11+sf.wm.mod[i]),c(i,i),col=rgb(sf.real.red[i],0,1-sf.real.red[i],0.5))}
for (i in 1:length(sf.real.red)){lines(c(2.9,2.9+logab[i]),c(i,i),col="gray")}
lines(c(10.5,11.5),c(-40,-40))
lines(c(10.5,10.5),c(-40,-80))
lines(c(11.5,11.5),c(-40,-80))
lines(c(11,11),c(-40,-80))
text(10.5,-150,"0",cex=0.5)
text(11,-150,"0.5",cex=0.5)
text(11.5,-150,"1",cex=0.5)
arrows(10.5,-250,11.5,-250,code=3,length=0.05)
text(10.8,-350,"Free-living",cex=0.7,pos=2,col="blue")
text(11.3,-350,"Particle-attached",cex=0.7,pos=4,col="red")
text(2.9+0.5,-350,"Number of reads (log)",cex=0.7)
logab.un<-log10(colSums(OTUtab.ss.reduced))
lines(c(2.9,2.9+max(logab)),c(-40,-40))
lines(c(2.9,2.9),c(-40,-80))
lines(c(2.9+max(logab),2.9+max(logab)),c(-40,-80))
text(2.9,-150,"0",cex=0.5)
text(2.9+0.5*max(logab),-150,as.character(round(2.9+0.5*max(logab.un)),1),cex=0.5)
text(2.9+max(logab),-150,as.character(round(2.9+max(logab.un)),1),cex=0.5)
pos.upp<-match(substr(sign.upp$OTU.ID,2,10),colnames(OTUtab.ss.reduced))
pos.down<-match(substr(sign.down$OTU.ID,2,10),colnames(OTUtab.ss.reduced))
points(rep(11.7,length(pos.upp)),pos.upp,col="red",cex=0.3,pch=19)
points(rep(11.9,length(pos.down)),pos.down,col="blue",cex=0.3,pch=19)
text(11.7-0.15,length(sf.real.red)+100,"Particle-attached indicator OTUs",pos=4,srt=90,cex=0.6)
text(11.9-0.15,length(sf.real.red)+100,"Free-living indicator OTUs",pos=4,srt=90,cex=0.6)
dev.off()


# Place Abundant Phylums within the NMDS 
Phylums.abundant<-aggregate(t(OTUtab.ss),by=list(Phylums),FUN=sum)
rownames(Phylums.abundant)<-Phylums.abundant$Group.1
Phylums.abundant<-Phylums.abundant[,-c(1)]
Phylums.abundant.red<-Phylums.abundant[which(table(Phylums)>40),]
Phylums.abundant.red



pdf("FigS4.pdf",width=10*(5/3),heigh=10)
layout(matrix(c(1:15), 3, 5, byrow = TRUE))
for (i in rownames(Phylums.abundant.red)){
  par(mar=c(3,3,2,2))
  quin<-i
  z <- kde2d(nmds$species[Phylums==quin,1],nmds$species[Phylums==quin,2], n=100,lims=c(-1.2,1.2,-1.2,1.2))
  plot(nmds$species,col=metad$filtersize.col,xlab="NMDS 1",ylab="NMDS 2",pch=19,type="n",ylim=c(-1.2,1.2),xlim=c(-1.2,1.2),main=i)
  ordihull(nmds,metad$filtersize,draw="polygon",col=rgb(0,0,1,0.0001),show.groups=0.2)
  ordihull(nmds,metad$filtersize,draw="polygon",col=rgb(1,0,0,0.0001),show.groups=0.8)
  points(nmds$points,col=c("blue","red")[metad$filtersize],xlab="NMDS 1",ylab="NMDS 2",pch=as.numeric(metad$filtersize)+14,cex=1.2)
  contour(z, drawlabels=T, nlevels=20,add=TRUE,col="black")
  points(nmds$species[Phylums==quin,1],nmds$species[Phylums==quin,2], xlab="NMDS 1", ylab="NMDS 2", pch=19, cex=log1p(colSums(OTUtab.ss[,Phylums==quin]))/2,col="gray")
  points(nmds$species[Phylums==quin,1],nmds$species[Phylums==quin,2], xlab="NMDS 1", ylab="NMDS 2", pch=1, cex= log1p(colSums(OTUtab.ss[,Phylums==quin]))/2,col="gray30")
}
dev.off()


pdf("FigS5.pdf",width=10,heigh=10*(8/4))
layout(matrix(c(1:32), 8, 4, byrow = TRUE))
for (i in exemples){
  par(mar=c(2,2,2,2))
  quin<-i
  z <- kde2d(nmds$species[grep(quin,OTUtab.complete.ss$SILVA_Taxonomy),1],nmds$species[grep(quin,OTUtab.complete.ss$SILVA_Taxonomy),2], n=100,lims=c(-1.2,1.2,-1.2,1.2))
  plot(nmds$species,col=metad$filtersize.col,xlab="NMDS 1",ylab="NMDS 2",pch=19,type="n",ylim=c(-1.2,1.2),xlim=c(-1.2,1.2),main=i)
  ordihull(nmds,metad$filtersize,draw="polygon",col=rgb(0,0,1,0.0001),show.groups=0.2)
  ordihull(nmds,metad$filtersize,draw="polygon",col=rgb(1,0,0,0.0001),show.groups=0.8)
  points(nmds$points,col=c("blue","red")[metad$filtersize],xlab="NMDS 1",ylab="NMDS 2",pch=as.numeric(metad$filtersize)+14,cex=1.2)
  contour(z, drawlabels=T, nlevels=20,add=TRUE,col="black")
  points(nmds$species[grep(quin,OTUtab.complete.ss$SILVA_Taxonomy),1],nmds$species[grep(quin,OTUtab.complete.ss$SILVA_Taxonomy),2], xlab="NMDS 1", ylab="NMDS 2", pch=19, cex=log1p(colSums(OTUtab.ss[,grep(quin,OTUtab.complete.ss$SILVA_Taxonomy)]))/2,col="gray")
  points(nmds$species[grep(quin,OTUtab.complete.ss$SILVA_Taxonomy),1],nmds$species[grep(quin,OTUtab.complete.ss$SILVA_Taxonomy),2], xlab="NMDS 1", ylab="NMDS 2", pch=1, cex= log1p(colSums(OTUtab.ss[,grep(quin,OTUtab.complete.ss$SILVA_Taxonomy)]))/2,col="gray30")
}
dev.off()

# Compute statistical tests
Wilcoxon.res<-matrix(NA,ncol=4,nrow=length(exemples))
for (i in 1:length(exemples)){
  quin<-exemples[i]
  res<-wilcox.test(sf.wm[grep(quin,OTUtab.complete.ss$SILVA_Taxonomy)],mu=0.5)
  Wilcoxon.res[i,1]<-mean(sf.wm[grep(quin,OTUtab.complete.ss$SILVA_Taxonomy)])
  Wilcoxon.res[i,2]<-sd(sf.wm[grep(quin,OTUtab.complete.ss$SILVA_Taxonomy)])
  Wilcoxon.res[i,3]<-res$statistic
  Wilcoxon.res[i,4]<-res$p.value}
colnames(Wilcoxon.res)<-c("Mean","sd","Statistic","Pvalue")
rownames(Wilcoxon.res)<-exemples
Wilcoxon.res<-as.data.frame(Wilcoxon.res)
Wilcoxon.res$corrPvalue<-p.adjust(Wilcoxon.res$Pvalue,"fdr")
Wilcoxon.res[,4]<-round(Wilcoxon.res[,4],4)
Wilcoxon.res[,5]<-round(Wilcoxon.res[,5],4)

Wilcoxon.res
write.table(Wilcoxon.res,file="Wilcoxon_test_exemples.txt",sep="\t",quote=F)

Wilcoxon.res<-matrix(NA,ncol=4,nrow=length(rownames(Phylums.abundant.red)))
for (i in 1:length(rownames(Phylums.abundant.red))){
  quin<-rownames(Phylums.abundant.red)[i]
  res<-wilcox.test(sf.wm[grep(quin,OTUtab.complete.ss$SILVA_Taxonomy)],mu=0.5)
  Wilcoxon.res[i,1]<-mean(sf.wm[grep(quin,OTUtab.complete.ss$SILVA_Taxonomy)])
  Wilcoxon.res[i,2]<-sd(sf.wm[grep(quin,OTUtab.complete.ss$SILVA_Taxonomy)])
  Wilcoxon.res[i,3]<-res$statistic
  Wilcoxon.res[i,4]<-res$p.value}
colnames(Wilcoxon.res)<-c("Mean","sd","Statistic","Pvalue")
rownames(Wilcoxon.res)<-rownames(Phylums.abundant.red)
Wilcoxon.res<-as.data.frame(Wilcoxon.res)
Wilcoxon.res$corrPvalue<-p.adjust(Wilcoxon.res$Pvalue,"fdr")
Wilcoxon.res[,4]<-round(Wilcoxon.res[,4],4)
Wilcoxon.res[,5]<-round(Wilcoxon.res[,5],4)

Wilcoxon.res

write.table(Wilcoxon.res,file="Wilcoxon_test_phyla.txt",sep="\t",quote=F)

######################################
## Boxplots for Phylum and Examples ##
######################################
pdf(file="Fig5.pdf",width=12,height=8)
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
bons<-which(Phylums %in% rownames(Phylums.abundant.red))
sf.wm.phyl<-sf.wm[bons]
Phylums.phyl<-factor(Phylums[bons])

par(mar=c(5,12,5,5),las=2)
boxplot(sf.wm.phyl~Phylums.phyl,horizontal=T,col="gray",cex=0.2,xlab="Particle-association niche value")
abline(v=0.5,lty=2)

bons<-NULL
clade<-NULL
j<-NULL
for (i in 1:length(exemples)){
  j<-j+1
  quins<-grep(as.character(unique(exemples)[i]),OTUtab.complete.ss$SILVA_Taxonomy)
  bons<-c(bons,quins)
  clade<-c(clade,rep(exemples[i],length(quins)))}

sf.wm.ex<-sf.wm[bons]
clade<-factor(clade,levels=exemples)


par(mar=c(5,12,5,5),las=2)
boxplot(sf.wm.ex~clade,horizontal=T,col="gray",cex=0.2,xlab="Particle-association niche value")
abline(v=0.5,lty=2)
dev.off()

#####################
#### TOP 30 OTUs ####
#####################
ordre<-order(colSums(OTUtab.ss),decreasing=T)
top30<-OTUtab.ss[,ordre[1:30]]
colnames(top30)<-paste("OTU",paste(colnames(top30),OTUtab.complete.ss$SILVA_Taxonomy[ordre[1:30]],sep=": "))


pdf("FigS6.pdf",onefile=T)
for (j in 1:5){
layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
for (i in ((j*6)-5):(j*6)){
	plot(metad$station[metad$filtersize==0.2],top30[metad$filtersize==0.2,i],col="blue",ylim=c(0,max(top30[,i])*1.5),ylab="Reads",xlab="Sampling Station",pch=19,main=colnames(top30)[i],type="n",xlim=c(4,150),cex.main=0.45,xaxt="n")
	axis(side=1,at=unique(metad$station),label=unique(metad$station))
	points(metad$station[metad$filtersize==0.2],top30[metad$filtersize==0.2,i],type="b",col="blue",ylim=c(0,max(top30[,i])*1.5),ylab="Reads",xlab="Sampling Station",pch=19,main=colnames(top30)[i])
	points(metad$station[metad$filtersize==0.8],top30[metad$filtersize==0.8,i],type="b",col="red",pch=19)}}
dev.off()







