

#analysis the potential of phase separation of sars-cov2 proteins


library(biomaRt)
library(magrittr)
library(data.table)
library(Biostrings)

setwd('~/Desktop/projects/phase_separation')
source('./scripts/functions_main.R')


#1. clean up seq data:-------------

ncov_peps <- readAAStringSet('./data/sars-cov2/proteins_sars-cov2.fa')

#only keep peps no shorter than 100:
(width(ncov_peps) >=100) %>% table() # F 6, T 24
ncov_peps_100 <- ncov_peps[width(ncov_peps) >=100]

writeXStringSet(ncov_peps_100, './results/sars-cov2/ncov_peps_100.fa')

#generate separate fa file for each peptide:
for (i in seq(1, length(ncov_peps_100))) {
	gene_name <- names(ncov_peps_100[i])
	file_name <- sprintf('./results/sars-cov2/sep_fa/%s_aa_seq.fa', gene_name)
	writeXStringSet(ncov_peps_100[i], filepath = file_name)
}

# 2 generate iupred2a scores for each gene ----
# 2.1 run iupred2a.py script -----------
#  sh run_iupred2a.sh
# 2.2 collect the iupred2a scores into a data.table -----------------------
gene_score_list <- list.files(path = './results/sars-cov2/iupred2a_results/',
							  full.names = T)
names(gene_score_list) <- gsub('_aa_seq.iupred2a.txt', '', gene_score_list %>% basename)

Score_Calculation(gene_score_list[3])


phase_score_dt <- data.table()
for (i in seq(1, length(gene_score_list))) {
	tmp <- Score_Calculation(gene_score_list[i])
	phase_score_dt <- rbind(phase_score_dt, tmp)
}


#fwrite(phase_score_dt, 
#	   file = './results/sars-cov2/phase_percentag_sars-cov2.txt', sep = '\t')
phase_score_dt <- fread(file = './results/sars-cov2/phase_percentag_sars-cov2.txt',
						sep = '\t')

# 3. calculate the charges by cider ---------------------------------------
# Calculate the FCR and NCPR of each genes by localcider python package.  
# The AA sequences should not contain asterix. 
# python ./scripts/cider/cider_sars-cov2.py

# 4. motif analysis by interproscan ---------------------------------------
# 4.1 get all motifs by interproscan for each protein coding gene ------------
#4.1.1 run_interproscan.sh:
#	sh ~/litt/software/interproscan/interproscan-5.31-70.0/interproscan.sh \
#		-i all_protein_coding_genes_without_asterix_aa_sequence.fa  \
#		-f tsv -f HTML -f GFF3 -f SVG
#4.1.2 clean up motif results:
#	sh process.sh
#	only collect motifs from Pfam, SMART, SUPERFAMILY, MobiDBLite, PRINTS and Coils.	

#motif info:
motifs <-fread(
	'./results/sars-cov2/interproscan/clean_data/clean_motifs_all.txt', sep = '\t')

names(motifs) <- c("gene_symbol", "AA_length", "analysis_methods",
				   "motif", "start","end" )





# 5. merge data together -----------------------------------------
# 5.1 generate clean dfs containing AA, iupred2 scores and cider charge results: ---------
#merge cider results and iupred2a scores by separate genes:
gene_list <- list.files('./results/sars-cov2/cider_results/') %>% gsub('.txt$', '', x = .)
#for (i in gene_list) {
#	print(i)
#	cider_path <- sprintf('./results/sars-cov2/cider_results/%s.txt', i)
#	ps_path <- sprintf('./results/sars-cov2/iupred2a_results/%s_aa_seq.iupred2a.txt', i)
#	cider <- fread(cider_path)
#	ps <- fread(ps_path)
#	names(ps) <- c('location', 'amino_acid', 'iupred2_score', 'achore2_score')
#	cider_ps <- cbind(cider, ps)
#	fwrite(cider_ps, file = sprintf(
#		'./results/sars-cov2/iupred_cider/iupred2a_cider_%s.txt', i), sep = '\t')
#}

#merge separate genes together into one list:
#gene_list <- list.files('./results/sars-cov2/iupred_cider/') %>% 
#	gsub('iupred2a_cider_', '', x = .) %>% gsub('.txt', '', x =.)

#aa_cider_iupred2a_list <- lapply(gene_list, function(x){
#	fread(sprintf('./results/sars-cov2/iupred_cider/iupred2a_cider_%s.txt', x))})
#names(aa_cider_iupred2a_list) <- gene_list
#saveRDS(aa_cider_iupred2a_list, file = './results/sars-cov2/rda/aa_cider_iupred2a_list.rda')
aa_cider_iupred2a_list <- readRDS('./results/sars-cov2/rda/aa_cider_iupred2a_list.rda')
aa_cider_iupred2a_list$M



# 6. visualization --------------------------------------------------------
library(RColorBrewer)
display.brewer.all()
palette(c(brewer.pal(8, 'Set2')))
#color package 1:
library(RSkittleBrewer)
palette(RSkittleBrewer("smarties"))

source('./scripts/functions_main.R')

motifs <- fread('./results/sars-cov2/interproscan/clean_data/clean_motifs_all.txt')
names(motifs) <- c("gene_symbol", "AA_length", "analysis_methods",
				   "motif", "start","end" )

aa_cider_iupred2a_list <- readRDS('./results/sars-cov2/rda/aa_cider_iupred2a_list.rda')

i = 'N'
plot_LCD_Motif_AAs(motif.dt = motifs[gene_symbol == i],
				   lcd.dt = aa_cider_iupred2a_list[[i]], 
				   gene_symbol = i)



names(aa_cider_iupred2a_list)

gene_list <- names(aa_cider_iupred2a_list) #NSP2, NSP6
gene_list <- gene_list[c(-10,-15)]


pdf('./results/sars-cov2/phase_analysis_sars-cov2.pdf', width = 12, height = 10)

for (i in gene_list) {
	print(i)
	plot_LCD_Motif_AAs(motif.dt = motifs[gene_symbol == i],
					   lcd.dt = aa_cider_iupred2a_list[[i]], 
					   gene_symbol = i)
	
}
dev.off()
 

aa_cider_iupred2a_list$NSP2
motifs[gene_symbol == 'NSP2']

