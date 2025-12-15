#!/bin/env Rscript

library(tidyverse)
library(Neoantimon)
#library(future)
#plan("multisession", workers = 20)
#for(r_code in dir("/Users/takanorihasegawa/Documents/GitHub/Neoantimon/R/", ".R")) {
#  print(r_code)
#  source(paste("/Users/takanorihasegawa/Documents/GitHub/Neoantimon/R/", r_code, sep = ""))
#}

path <- scan(file="/lustre/biomed/asilvestris/dati_Cives/b_mandriani_pancreas/hla_paths",
             what='character',
             sep='\t')

nif<- do.call(rbind,
        lapply(path, function(o)
        {
                # Sample
                mem<- regmatches(o, regexpr("ICGC_[0-9]*", o))
                # A1 & A2
                lmao<- read.table(o, header=F, skip=10, sep='\t')$V2
                lmao<- regmatches(lmao, regexpr("[^:]*:[^:]*|[0-9]*$", lmao))
                ayy<- gsub('HLA-|.csv', '', basename(o))
                ifelse(length(lmao)==2,
                        return(c(mem, paste(ayy, lmao, sep='*'), ayy)),
                        return(c(mem, paste(ayy, c('NA', 'NA'), sep='*'), ayy)))
        })
)

fin<- as.data.frame(nif)
nams<- as.data.frame(aggregate(fin, V4~V1, FUN=cbind)$V4)[1,]
Name<- aggregate(fin, V2~V1, FUN=cbind)$V1
one<- as.data.frame(aggregate(fin, V2~V1, FUN=cbind)$V2)
two<- as.data.frame(aggregate(fin, V3~V1, FUN=cbind)$V3)
fin<- cbind(Name, as.data.frame(mapply(cbind, one, two, SIMPLIFY=F)))
yuck<- names(fin)
for (s in 1:ncol(nams))
{
        yuck<- gsub(paste0(colnames(nams)[s], '\\.'), nams[1,s], yuck)
}
names(fin)<- yuck

fin

write.table(fin, file="HLA_types.txt", sep = "\t")


path2 <- scan(file="/lustre/biomed/asilvestris/dati_Cives/b_mandriani_pancreas/gene_expression_path",
             what='character',
             sep='\t')

fin_result_snv <- NULL
fin_result_indel <- NULL
for(file in dir("./", "hg38_multianno.txt")){
  result <- MainSNVClass1(input_annovar_format_file = file, 
                                      refflat_file  = "/lustrehome/asilvestris/exe/miniconda3/envs/Neoantimon/lib/refFlat.grch38.txt", 
                                      refmrna_file = "/lustrehome/asilvestris/exe/miniconda3/envs/Neoantimon/lib/refMrna.grch38.fa", 
                                      netMHCpan_dir = "/lustrehome/asilvestris/exe/miniconda3/envs/Neoantimon/lib/netMHCpan-4.1/netMHCpan", 
                                      hla_types = fin[match(regmatches(file, regexpr('ICGC[^\\.]*', file)), fin[,1]),c(2:7)], 
				      rnaexp_file = path2[match(regmatches(file, regexpr('ICGC[^\\.]*', file)), regmatches(path2, regexpr('ICGC_[^\\/]*', path2)))],
                                      multiple_variants = TRUE, 
                                      depth_tumor_column = 26, 
                                      depth_normal_column = 25)
  if(length(which(!is.na(as.numeric(result[, match("Total_RNA", colnames(result))])))) == 0){
    print("Please add RNA-expression values to calculate Priority scores.")
  } else {
    PS <- CalculatePriorityScores(result = result, useRNAvaf = FALSE)
  }
  fin_result_snv <- rbind(fin_result_snv, c(file, "HLAc1", Export_Summary_SNV(result, Mut_Rank_th = 0.5)))
  
  result <- MainSNVClass2(input_annovar_format_file = file, 
                          refflat_file  = "/lustrehome/asilvestris/exe/miniconda3/envs/Neoantimon/lib/refFlat.grch38.txt", 
                          refmrna_file = "/lustrehome/asilvestris/exe/miniconda3/envs/Neoantimon/lib/refMrna.grch38.fa", 
                          netMHCIIpan_dir = "/lustrehome/asilvestris/exe/miniconda3/envs/Neoantimon/lib/netMHCIIpan-4.3/netMHCIIpan", 
                          hla_types = fin[match(regmatches(file, regexpr('ICGC[^\\.]*', file)), fin[,1]),c(16:27)], 
			  rnaexp_file = path2[match(regmatches(file, regexpr('ICGC[^\\.]*', file)), regmatches(path2, regexpr('ICGC_[^\\/]*', path2)))],
                          multiple_variants = TRUE, 
                          depth_tumor_column = 26, 
                          depth_normal_column = 25)
  if(length(which(!is.na(as.numeric(result[, match("Total_RNA", colnames(result))])))) == 0){
    print("Please add RNA-expression values to calculate Priority scores.")
  } else {
    PS <- CalculatePriorityScores(result = result, useRNAvaf = FALSE)
  }
  fin_result_snv <- rbind(fin_result_snv, c(file, "HLAc2", Export_Summary_SNV(result, Mut_Rank_th = 1)))

}

write.table(x = fin_result_snv, "summary_snv.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

