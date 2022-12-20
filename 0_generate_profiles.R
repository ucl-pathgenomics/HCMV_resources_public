# script to extract alignment per MAregion, Cluster
# this isnt self-contained, and relies o other data oscar has , sorry :/
library(ape)

assignments = read.csv("out/3-clean/gt-assignment.csv")


regions = unique(assignments$region)

#reg = regions[1]


for(reg in regions){
  


  infile = list.files(path = "out/3-clean/", pattern = paste0("^",reg , "-[0-9]{1,6}-[0-9]{1,6}.fasta"),full.names = T)
  msa = read.dna(infile , format = "fasta")
  
  
  t = assignments[assignments$region == reg,]
  clusters = unique(t$genotype)
  for(cluster in clusters){
    
    #cluster = clusters[1]
    
    
    seq_names = t[t$genotype == cluster, ]$sample
    
    # subset alignment
    tmsa = msa[labels(msa) %in% seq_names,]
    
    
    write.FASTA(tmsa, paste0("95_hmmprofiles/msa/",reg, "_", cluster, ".fasta"))
    
    
    
    # --------- generate hmm profile
    
    command = paste0("wsl hmmbuild",
                     " 95_hmmprofiles/hmm/", reg, "_", cluster, ".profile", # hmm name
                     " 95_hmmprofiles/msa/", reg, "_", cluster, ".fasta" )
    system(command)
  }
}