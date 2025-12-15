# You launch this thing like so
# Rscript alle_gen.R <list_of_sample_names> <folder_with_ALL_files> <list_of_ots_files>
# put pahts with wildcards in quotes e.g. "/lustre/biomed/asilvestris/dati_Cives/b_mandriani_pancreas/result.ID.SNV*"

# gotta provide a sample name list since there's no recurrent pattern

# TESTIN'
#args<- c('mandrianisamps.txt',
#	'/lustre/biomed/asilvestris/dati_Cives/b_mandriani_pancreas/results_no_filtering/result.ID.SNV*',
#	'/lustre/biomed/asilvestris/dati_Cives/b_mandriani_pancreas/denovo_pipeline3_outTable_paths.txt')

library(dplyr)

args<- commandArgs(trailingOnly=T)

nams<- scan(as.character(args[1]), what='character')

snv<- suppressWarnings(system(paste0("find ", as.character(args[2]), " -name '*ALL*' 2>/dev/null"), intern=T))
snv<- snv[sapply(nams, function(s){grep(s, snv)})]

ots<- scan(as.character(args[3]), what='character')
ots<- rep(ots[sapply(nams, function(s){grep(s, ots)})], each=2)

# all(gsub('.*(ICGC_[^\\/]*).*', '\\1', ots)== gsub('.*(ICGC_[^\\.]*).*', '\\1', snv))
# TRUE

mapply(function(to, lel)
{
	ot<- read.table(to, header=T, sep='\t', quote='')
	lal<- read.table(lel, header=T, sep='\t', quote='')
	lal$Frequency<- ot$Frequency[match(lal$Mutation_Position, with(ot, paste(Region, Position, sep='_')))]
	lal<- lal %>% relocate(Frequency, .after=Mutation_Position)
	write.table(lal,
		file= basename(lel),
		row.names=F,
		sep= "\t",
		quote= F)
}, ots, snv)
