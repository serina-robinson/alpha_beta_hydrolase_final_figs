parse_rename_fasta<-function(seq_fil){
  seqs<-readAAStringSet(seq_fil)
  
  nams<-names(seqs)
  nams<-gsub("\\("," ",nams)
  nams<-gsub("\\)"," ",nams)
  nams<-gsub("\\,"," ",nams)
  nams<-gsub("\\:"," ",nams)
  nams<-gsub("-"," ",nams)
  orgs<-unlist(lapply(1:length(nams),function(x) { gsub(" ", "_", gsub("\\]","",strsplit(nams,"\\[")[[x]][2])) }))
  desc<-unlist(lapply(1:length(nams),function(x) { gsub(" ", "_", gsub("\\]","",strsplit(nams,"\\[")[[x]][1])) }))
  #desc<-rep("oleB_",length(names(seqs)))
  
  un_orgs<-orgs[!duplicated(orgs)]
  un_desc<-desc[!duplicated(orgs)]
  un_seqs<-seqs[!duplicated(orgs)]
  
  newnames<-paste0(un_desc,un_orgs)
  names(un_seqs)<-newnames
  
  return(un_seqs)
}
