#Install packages
#install.packages("ggtree")
packs<-c("DECIPHER","Biostrings",'phangorn')
lapply(packs,require,character.only=T)
#source("https://bioconductor.org/biocLite.R")
biocLite('bgafun')
biocLite('DECIPHER')

#Make a fancy plot using ggtree
myt<-read.tree("output/45_seqs_ML.nwk")

ggtree(myt) + geom_treescale() + geom_tiplab() + geom_text2(aes(subset=!isTip,label=node),hjust=-.3)


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