

#pLot the ps scores of ncov2 proteins

library(biomaRt)
library(magrittr)
library(data.table)
library(Biostrings)

setwd('~/Desktop/projects/phase_separation')
source('./scripts/functions_main.R')


#1. clean up seq data:-------------

ncov_peps <- readAAStringSet('./data/sars-cov2/proteins_sars-cov2.fa')



##generate separate fa file for each peptide:
#for (i in seq(1, length(ncov_peps))) {
#	gene_name <- names(ncov_peps[i])
#	file_name <- sprintf('./results/sars-cov2/sep_fa/%s.fa', gene_name)
#	writeXStringSet(ncov_peps[i], filepath = file_name)
#}


# 2 generate iupred2a scores for each gene ----
# 2.1 run iupred2a.py script -----------
#  sh run_iupred2a.sh
# 2.2 collect the iupred2a scores into a data.table -----------------------
gene_score_list <- list.files(path = './results/sars-cov2/iupred2a_results/',
							  full.names = T)
names(gene_score_list) <- gsub('.iupred2a.txt', '', gene_score_list %>% basename)

Score_Calculation(gene_score_list[3])


gene_score_list




phase_score_dt <- data.table()
for (i in seq(1, length(gene_score_list))) {
	tmp <- Score_Calculation(gene_score_list[i])
	phase_score_dt <- rbind(phase_score_dt, tmp)
}



library(tibble)
library(tidyverse)

phase_score <- as.tibble(phase_score_dt)
phase_score %<>% mutate(
	mean_score = (percent_of_iupred2a + percent_of_anchor2)/2) %>%
	arrange(-mean_score)

ps1 <- filter(phase_score,!(gene_symbol %in% c('FUS', 'mEGFP')))
ps2 <- filter(phase_score, gene_symbol %in% c('FUS', 'mEGFP'))
phase_score_order <- rbind(ps1, ps2)

#fwrite(phase_score_order, 
#	   file = './results/sars-cov2/phase_percentag_sars-cov2.txt', sep = '\t')
phase_score_order <- read_tsv( file = './results/sars-cov2/phase_percentag_sars-cov2.txt')


library(pheatmap)


#discard NSP3N, NSP3C

phase_score_order <- phase_score_order[!(phase_score_order$gene_symbol %in% c('NSP3N', 'NSP3C')),]



pheatmap(phase_score_order[,2:3], 
		 color = colorRampPalette(c(gray(98/100), 'springgreen'))(100),
		 lwd =0.4,
		 cluster_rows = F, cluster_cols = F,
		 labels_row = phase_score_order$gene_symbol,
		 labels_col = c('iupred2a', 'anchor2'),
		 gaps_row = c(29),
		 display_numbers = T
)



rgb(0, 158,143, maxColorValue = 255)
pheatmap(phase_score_order[,2:3], 
         color = colorRampPalette(c(gray(98/100), '#009E8F'))(100),
         lwd =0.4,
         cluster_rows = F, cluster_cols = F,
         labels_row = phase_score_order$gene_symbol,
         labels_col = c('iupred2a', 'anchor2'),
         gaps_row = c(29),
         display_numbers = F
)


pheatmap(phase_score_order[,2:3], 
		 color = colorRampPalette(c('seagreen', "white", "firebrick3"))(100),
		 lwd =0.4,
		 cluster_rows = F, cluster_cols = F,
		 labels_row = phase_score_order$gene_symbol,
		 labels_col = c('iupred2a', 'anchor2'),
		 gaps_row = c(31)
)


pheatmap(phase_score_order[,2:3], 
		 lwd =0.4,
		 color = colorRampPalette((brewer.pal(n = 8, name = "Blues")))(100),
		 cluster_rows = F, cluster_cols = F,
		 labels_row = phase_score_order$gene_symbol,
		 labels_col = c('iupred2a', 'anchor2'),
		 gaps_row = c(31)
)


grep('green', colors(), value = T)

barplot(seq(1:length(grep('green', colors(), value = T))),
		col = grep('green', colors(), value = T))










