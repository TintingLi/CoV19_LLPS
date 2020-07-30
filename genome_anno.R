
library(tidyverse)
library(magrittr)
library(Biostrings)

setwd('~/Desktop/projects/ncba/litao/N_phase_separation/')


sars_cov2_genomes <- readDNAStringSet('~/data/gisaid/20200619/sequences_2020-06-18_21-15.fasta')
sars_cov2_genomes_abfreq <- alphabetFrequency(sars_cov2_genomes) %>% as_tibble()

#there are 196 seqs containing '-'
table(sars_cov2_genomes_abfreq$`-` > 0)    #F 48687, T 196

#remove seqs containing '-':
sars_cov2_genomes_removeNs <- sars_cov2_genomes[sars_cov2_genomes_abfreq$`-` == 0]
#writeXStringSet(sars_cov2_genomes_removeNs,
#                '~/data/gisaid/20200619/sequences_2020-06-18_21-15_removeNs.fasta')



# map N protein to genome sequences by exonerate -------------------------------
#pro2genome.sh --> clean_align.txt

n2genome <- read_tsv('./results/exonerate/clean_align.txt',
                     col_names = c('query_range', 'target_genome', 'targe_range'))
#only keep alignments containg complete N protein sequence: 47353 alignments
n2genome <- n2genome[n2genome$query_range == 'Query range: 0 -> 419',]
#write_tsv(n2genome, './results/exonerate/clean_complete_align.txt')

table(n2genome$query_range)
#get N gene locus and subtract nt sequences   --------------
N_locus <- tibble( strain= gsub('Target: ', '', n2genome$target_genome), 
                   N_start = gsub('Target range: ', '', n2genome$targe_range) %>%
                       gsub(' ->.*$', '', x =.) %>% as.integer(),
                   N_end = gsub('Target range: ', '', n2genome$targe_range) %>%
                       gsub('^.*-> ', '', x =.) %>% as.integer())
#left locus: +1, right locus: +3
N_locus$N_start <- N_locus$N_start + 1
N_locus$N_end <- N_locus$N_end + 3
N_locus$width <- N_locus$N_end - N_locus$N_start + 1

table(N_locus$width) 

# 47288 seqs's length is 1260
#only keep gene_length == 1260:
N_locus_1260 <- filter(N_locus, width == 1260)

#subtract all N gene sequences: 
sars_cov2_genomes <- sars_cov2_genomes[names(sars_cov2_genomes) %in% N_locus$strain]
N_seqs <- DNAStringSet() # 47288 sequences
for (i in seq(1, dim(N_locus_1260)[1])) {
    gn <- N_locus_1260$strain[i]
    gs <- sars_cov2_genomes[names(sars_cov2_genomes) == gn]
    tmp <- subseq(gs, 
                  start = N_locus_1260$N_start[i],
                  end = N_locus_1260$N_end[i])
    N_seqs <- append(N_seqs, values = tmp)
}
#writeXStringSet(N_seqs, './results/N_nt_all.fa')

#calculate the freqs for each nt: -----------
N_seqs <- readDNAStringSet('./results/N_nt_all.fa')

width(N_seqs) %>% table()

N_seqs_matrix <- as.matrix(N_seqs)
st_n <- tibble()
st_n <- N_seqs_matrix[, 1] %>% table() %>% as_tibble()
for (i in seq(2, width(N_seqs)[1])) {
    tmp <- N_seqs_matrix[, i] %>% table() %>% as_tibble()
    st_n <- full_join(st_n, tmp, by = '.')
}
names(st_n) <- c('nt', paste0('pos', 1:1260))

st_n_sort_m <- as.matrix(st_n[, 2:1261])
#turn NA into 0:
st_n_sort_m[is.na(st_n_sort_m)] <- 0
rownames(st_n_sort_m) <- st_n$nt
#write_csv(as_tibble(st_n_sort_m, rownames = 'nt'), "./results/Nt_freq_N.csv")
st_n_sort_m <- read_csv('./results/Nt_freq_N.csv')


#get conserved sequences:
conserved_nt <- apply(st_n_sort_m, 2, function(x){
    return(rownames(st_n_sort_m)[which.max(x)])
}) 

N_conserved <- paste(conserved_nt,collapse = '') %>% DNAStringSet()
#writeXStringSet(N_conserved, './results/N_nt_conserved.fa')


#get highest mutation sites:

mutation_freq_N <- apply(st_n_sort_m, 2, function(x){
    tmp <- x[order(x, decreasing = T)]
    names(tmp) <- rownames(st_n_sort_m)[order(x, decreasing = T)]
    num1 <- paste0(tmp[1], " (", names(tmp)[1], ')')
    num2 <- paste0(tmp[2], " (", names(tmp)[2], ')')
    num3 <- paste0(tmp[3], " (", names(tmp)[3], ')')
    return(c(num1, num2, num3))
}) %>% t() %>% as_tibble(rownames = 'position')


colnames(mutation_freq_N) <- c('postion', 'conserved', 'mut1', 'mut2')

mutation_freq_N

mutation_freq_N$mut1_num <- gsub(' (.*)$', '', mutation_freq_N$mut1) %>% as.numeric()
mutation_freq_N$mut1_AA <- gsub('^[0-9].* ', '', mutation_freq_N$mut1)
mutation_freq_N$mut2_num <- gsub(' (.*)$', '', mutation_freq_N$mut2) %>% as.numeric()
mutation_freq_N$mut2_AA <- gsub('^[0-9].* ', '', mutation_freq_N$mut2)
mutation_freq_N %<>% arrange(desc(mut1_num, mut2_num))

write_csv(mutation_freq_N,  './results/mutation_freq_N_nt.csv' )

mutation_freq_N <- read_csv('./results/mutation_freq_N_nt.csv')
#heatmap --------
library(pheatmap)
library(RColorBrewer)

nt_freq_n <- read_csv('./results/Nt_freq_N.csv')

nt_freq_n_m <- as.matrix(nt_freq_n[2:1261])
rownames(nt_freq_n_m) <- nt_freq_n$nt

nt_freq_n_m[nt_freq_n_m > 15000] <- 0

pheatmap(nt_freq_n_m[1:4,], 
         cluster_rows = F, cluster_cols = F, 
         color = colorRampPalette((brewer.pal(n = 8, name = "Blues")))(100),
         labels_col = conserved_nt, cex.axis= 0.8, las = 2,
         show_colnames = F,
         main = 'mutations of N genes of SARS-CoV2',
         ylab ='conserved AA sequences')

pheatmap(nt_freq_n_m[1:4,], 
         cluster_rows = F, cluster_cols = F, 
         color = c("#F8F8F8",
                   colorRampPalette(c( gray(0.90), 
                                       rgb(0, 158, 143, maxColorValue = 255)))(n=599)),
         labels_col = conserved_nt, cex.axis= 0.8, las = 2,
         show_colnames = F,
         main = 'mutations of N genes of SARS-CoV2',
         ylab ='conserved AA sequences')


par(oma = c(4,4,4,4))
tes <- pheatmap(t(colSums(nt_freq_n_m) ), oma = c( 4,4,4,4), mar = c( 4,4,4,4),
         border_color = 'red',
         cluster_rows = F, cluster_cols = F, 
         color = c("#F8F8F8",
                   colorRampPalette(c( gray(0.90), 
                                       rgb(0, 158, 143, maxColorValue = 255)))(n=599)),
          cex.axis= 0.8, las = 2,
         show_colnames = F,
         main = 'mutations of N genes of SARS-CoV2',
         ylab ='conserved AA sequences')





grid.rect(x = c(1, 1260))
box(which = 'figure',lwd = 3, col = 'red')


tes$gtable$widths




redcolors <- grep('red', colors(), value = T)

barplot(rep(1,27), col = redcolors[1:27])



pheatmap(nt_freq_n_m[1:4,], 
         cluster_rows = F, cluster_cols = F, 
         color = colorRampPalette(c(gray(95/100), rgb(0, 158, 143, maxColorValue = 255)))(100),
         labels_col = conserved_nt, cex.axis= 0.8, las = 2,
         show_colnames = F,
         main = 'mutations of N genes of SARS-CoV2',
         ylab ='conserved AA sequences')



nt_freq_n_m

colSums(nt_freq_n_m)
t(colSums(nt_freq_n_m))

image(matrix(colSums(nt_freq_n_m)), add =T)

pheatmap(colSums(nt_freq_n_m) %>% t(),cluster_rows = F, cluster_cols = F)


#add ps score info:
gene_score_list <- list.files(
    path = '~/Desktop/projects/phase_separation/results/sars-cov2/iupred2a_results/',
                              full.names = T)
names(gene_score_list) <- gsub('.iupred2a.txt', '', gene_score_list %>% basename)


source('~/Desktop/projects/phase_separation/scripts/functions_main.R')



library(RSkittleBrewer)
palette(RSkittleBrewer("tropical"))

n_ps <- read_tsv(gene_score_list[names(gene_score_list) == "N"], 
                 skip = 6, 
                 col_names = c('location', 'AA', 'iupred2a_score', 'anchor2_score'))

par(bg =(brewer.pal(n = 8, name = "Blues"))[2], oma= c(4,4,4,4))
plot(n_ps$location, n_ps$iupred2a_score, col =3, ann=F, axes =F,
     type = 'l', lwd = 2, ylim = c(0, 1.1))
lines(n_ps$location, n_ps$anchor2_score,  lwd = 2, col = 4)
#abline(h =0.5, lty = 2)
axis(1, at = c(1, seq(50, 400, by = 50), 419))
axis(4, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0))

legend('topright',c('iupred2a','anchor2'), lwd =2, col = c(3,4), 
       bty = 'n' )
    
    
    
    
    
    
    
    
    
    
    


