#Install packages
packs<-c("DECIPHER","Biostrings",'phangorn')
lapply(packs,require,character.only=T)
#source("https://bioconductor.org/biocLite.R")
#biocLite('ggtree')


#Read in sequences for alignment
fas<-"unaligned/45_seqs.fa"
aa<-readAAStringSet(fas)

#Trim the last 600 basepairs from OleBC fusion sequences (last 6 sequences)
numseqs<-length(aa)
bc<-aa[(numseqs-5):numseqs]
width(bc)
bc.tr<-subseq(bc,start=1,end=width(bc)-600)
width(bc.tr)
aa2<-c(aa[1:(numseqs-6)],bc.tr)

#Align sequences
aa.al <- AlignSeqs(aa2)
BrowseSeqs(aa.al,htmlFile = "45_seqs_aln.html",highlight=0)

#Adjust alignment
aa.adj<-AdjustAlignment(aa.al)
BrowseSeqs(aa.adj,htmlFile="45_seqs_adj.html",highlight=0)
writeXStringSet(aa.adj,"output/45_seqs_aln.fa")

#Stagger alignment
aa.stag<-StaggerAlignment(aa.adj)
writeXStringSet(aa.stag,filepath = "45_seqs_stag.fa",format = "fasta")

#Read in phyDat object for phangorn treeing
phy<-read.phyDat("45_seqs_stag.fa",format="fasta",type="AA")

#Compute a distance matrix
dm<-dist.ml(phy)

#UPGMA
treeUPGMA<-upgma(dm)

#Neighbor-joining
treeNJ<-NJ(dm)

#Parsimony
#treePars<-optim.parsimony(treeUPGMA,phy)
#treeRatchet<-pratchet(phy,trace=0)
#parsimony(c(treePars,treeRatchet),phy)

#Maximum-likelihood (JTT model)
ml<-pml(treeNJ,model="JTT",data=phy)
ml.JTT<-optim.pml(ml,model="JTT")#optimize branch lengths for Jukes-Cantor model
logLik(ml.JTT)

#Plot trees
par(mfrow=c(1,1))
plot(treeUPGMA, main="UPGMA")
plot(treeNJ, main="NJ",mar=0)
plot(ml.JTT$tree,main="ML",mar=0)
treeML<-ml.JTT$tree

#phangorn Bootstrap
#bs<-bootstrap.pml(ml.JTT,bs=100,optNni=TRUE,control=pml.control(trace=0))

#Plot ML with bootstrap
pl<-plotBS(midpoint(ml.JTT$tree),bs,p=50,type="p")

#Fix tree names
pl$tip.label<-gsub("X","",pl$tip.label,fixed=TRUE)
pl$tip.label<-gsub("_"," ",pl$tip.label,fixed=TRUE)
pl$tip.label
#Write tree to newick format
write.tree(pl,file="output/45_seqs_ML.nwk")

#Plot consensus net from bootstrap sample
#cnet<-consensusNet(bs,p=0.2)
#plot(cnet,"2D",show.edge.label=FALSE)


