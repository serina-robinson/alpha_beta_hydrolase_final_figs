#Install packages
packs<-c('seqinr','Biostrings','ape','cluster','fpc','kernlab','mclust','RColorBrewer','ggplot2','bios2mds')
#lapply(packs,install.packages)
lapply(packs,require,character.only=T)



#Using the bio2mds packages
#MDS analysis
#seqs<-read.alignment("output/45_seqs_aln.fa",format="fasta")
seqs<-read.alignment("data/aligned/45_seqs_muscle_trimmed.fa",format="fasta")
phy<-read.phyDat("output/45_seqs_aln.fa",format="fasta",type="AA")
dm<-dist.ml(phy)
mat<-as.matrix(dm)

#Using cmdscale
mds<-cmdscale(dm,eig=TRUE,k=2)
x<-mds$points[,1]
y<-mds$points[,2]
plot(x,y)

#Try mat.dif
#seqs<-import.fasta("output/45_seqs_aln.fa")
#mat2<-mat.dif(seqs,seqs)

#MMDS
mds<-mmds(mat)

#Project
mmds3<-mmds.project(mds,sup=mat)

#2D plot
mmds.2D.plot(mds)

#3D plot
mmds.3D.plot(mds,mmds3)

#Plot my way
x<-mds$coord[,1]
y<-mds$coord[,2]
plot(x,y)

#Find Xanthomonas campestris catalytic nucleophile
xcamp<-seqs$seq[[grep("alpha/beta_hydrolase_Xanthomonas_campestris",seqs$nam)]]

#Find specific sequence motifs
#tri.pos<-words.pos("hdw",xcamp)
#tri.pos<-tri.pos+1
#tri.pos<-words.pos("hyvle",xcamp)
tri.pos<-words.pos("rdi",xcamp)
# oxy.pos<-words.pos("hg",xcamp) 

#Find catalytic nucleophile in all sequences
nuc<-lapply(seqs$seq,function(x) { substr(x,tri.pos,tri.pos) })
nucr<-unlist(nuc)
names(nucr)<-seqs$nam

#Find oxyanion holde in all sequences
nuc2<-lapply(seqs$seq,function(x) { substr(x,oxy.pos,oxy.pos+3) })
nucr<-unlist(nuc2)
names(nucr)<-seqs$nam

#Fix the epoxide hydrolases because they didn't align properly...
#nucr[grep("epoxide",names(nucr),ignore.case=T)]<-"d"
#nucr[grep("thioesterase",names(nucr),ignore.case=T)]<-"s"

####Now try to group into protein subfamilies####
#Set palette
pal<-palette(colorRampPalette(colors=brewer.pal(11,"Spectral"))(10))

clusters<-rep(NA,length(x))
colors<-rep(NA,length(x))

#Epoxide hydrolase
ep.ind<-grep("epoxide",names(nucr),ignore.case=T)
clusters[ep.ind]<-"Epoxide hydrolase"
colors[ep.ind]<-"gray28" 
#colors[ep.ind]<-pal[1]

#Gamma lactamase
gam.ind<-grep("gamma",names(nucr),ignore.case=T)
clusters[gam.ind]<-"Gamma lactamase"
colors[gam.ind]<-"saddlebrown"

#Haloperoxidases
cl.ind<-grep("peroxidase",names(nucr),ignore.case=T)
clusters[cl.ind]<-"Haloperoxidase"
colors[cl.ind]<-"dark blue"

#Peptidase
pep.ind<-grep("peptidase",names(nucr),ignore.case=T)
clusters[pep.ind]<-"Peptidase"
colors[pep.ind]<-"gold4"

#Lipase
lip.ind<-grep("lipase",names(nucr),ignore.case=T)
clusters[lip.ind]<-"Lipase"
colors[lip.ind]<-"orange"

#Thioesterase
thi.ind<-grep("thioesterase",names(nucr),ignore.case=T)
clusters[thi.ind]<-"Thioesterase"
colors[thi.ind]<-"red4"

#Transferase
tr.ind<-grep("transferase",names(nucr),ignore.case=T)
clusters[tr.ind]<-"Transferase"
colors[tr.ind]<-"aquamarine3"

#Dehalogenase
deh.ind<-grep("dehalogenase",names(nucr),ignore.case=T)
clusters[deh.ind]<-"Haloalkane dehalogenase"
colors[deh.ind]<-"deeppink"

#Dienelactone hydrolase
deh.ind<-grep("dienelactone",names(nucr),ignore.case=T)
clusters[deh.ind]<-"Dienelactone hydrolase"
colors[deh.ind]<-"red"

#OleB
oleB.ind<-grep("alpha/beta_hydrolase",names(nucr),ignore.case=T)
clusters[oleB.ind]<-"OleB"
colors[oleB.ind]<-"purple"


allind<-c(ep.ind,gam.ind,cl.ind,pep.ind,lip.ind,thi.ind,tr.ind,deh.ind,oleB.ind)

#OleBC
oleBC.try<-grep(paste(c("hydrolase","acyl_CoA","AMP"),collapse="|"),names(nucr),ignore.case=T)
oleBC.ind<-oleBC.try[!oleBC.try %in% allind]
clusters[oleBC.ind]<-"OleBC fusion"
colors[oleBC.ind]<-"forestgreen"

#For legend
leg<-unique(clusters)
col<-unique(colors)

dat<-data.frame(cbind(x,y))
pdf("output/oxyanion_hole.pdf",width=8,height=8)
par(mar=c(0.01,0.01,0.01,0.01))
mycols<-colorRampPalette(colors=brewer.pal(8,"Dark2"))(10)
ggplot(data = dat,aes(x=x,y=y)) + geom_text(label=toupper(nucr),size=6,
                                              colour=colors) +
  labs(x="Coordinate 1",y="Coordinate 2") +
  theme_bw() +
  #scale_fill_manual(name='Protein family',guide="legend",labels=nam,values=mycols) +
  theme(axis.text = element_text(size = 14),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.position = c(0.14, 0.80),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=element_text(size=20,face="bold",hjust=0)
  )+
  ggtitle("Catalytic nucleophile in diverse \n alpha/beta hydrolase superfamily proteins")
legend("bottomleft",legend=leg,fill=col,border=F,title='Protein family',bty="n")
dev.off()
 
#MMDS # DOES NOT RUN
mds
#mds<-mmds(mat)
matr<-mds$group

newmat<-matrix(nrow = 11,ncol = 2)
newmat[,1]<-leg
newmat[,2]<-col
colnames(newmat)<-c("group","color")
mds$col

newcol<-matrix(nrow=45,ncol=3)
colnames(newcol)<-colnames(mds$col)
newcol[,1]<-rownames(mds$coord)
newcol[,2]<-clusters
newcol[,3]<-colors

mds$group<-newmat
mds$col<-newcol

#2D plot
mmds.2D.plot(mds,active.legend.pos="bottomright",grid = FALSE,active.legend.text = 0.75,outfile.type = "pdf",outfile.name = 
               "output/MDS_plot")

#Screeplot
scree.plot(mds$eigen.perc,lab=TRUE,"Screeplot of metric MDS")

#Project
mds.proj<-mmds.project(mds,sup=mat)
mds.proj$group<-newmat
mds.proj$col<-newcol

mds.proj
#3D plot
mmds.3D.plot(mds,mds.proj,active.type="s",sup.type="p")
