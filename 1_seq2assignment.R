#!/usr/bin/env Rscript
# in fasta file cmv genome
# out assignment of MultiAllelic allele asignment
# requires hmmer in your path
# usage: Rscript 1_seq2assignment.R test.fasta

library(ape)

args = commandArgs(trailingOnly=TRUE)

infile = args[1]

options = list.files("hmm")


region = as.numeric(stringr::str_split(options,"_", simplify = T)[,1])
allele = stringr::str_split(options,"_", simplify = T)[,2]
allele = as.numeric(stringr::str_split(allele,"\\.", simplify = T)[,1])

opts = data.frame(region, allele)
opts = opts[order(region, allele),]

# now we have all the profiles ordered so we can loop through
cat("CMV genotyping by HMM profiles. O Charles & R Goldstein, UCL 2021\n")
cat("Region", "Allele\n")
for(reg in unique(opts$region)){

    best.score = 0
    best.allele = 0
    for(all in unique(opts[opts$region == reg,2])){

        # run hmm

        #command = paste0("nhmmer --noali --tblout=t.t --dna  hmm/", paste0(reg,"_", all, ".profile") , " ", infile)
        command = paste0("nhmmer --noali --dna  hmm/", paste0(reg,"_", all, ".profile") , " ", infile)

        # parse and get score

        t = system(command, intern = T)

        which.line = grep("E-value", t)[1] + 2

        # the best match is always the top, use this to avoid where hmmsearch find other alignment bits
        tscore = read.table( text = t[which.line])[1,2]

        # determine the best scoring allele
        if(tscore > best.score){
            best.score = tscore
            best.allele = all
        }
    }
    cat(paste(reg, best.allele,"\n"))
}