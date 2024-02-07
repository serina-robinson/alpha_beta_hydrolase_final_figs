fix_fasta_names<-function(seq_fil,out_fil){
  seqs<-readAAStringSet(seq_fil)
  nams<-names(seqs)
  len<-unique(unlist(lapply(1:length(nams),function(x) { strsplit(nams," ")[[x]][2] })))
  orgs<-unlist(lapply(1:length(nams),function(x) { gsub("_", " ", strsplit(nams," ")[[x]][1]) }))
  names(seqs)<-orgs
  writeXStringSet(seqs,out_fil)
  return(len)
}