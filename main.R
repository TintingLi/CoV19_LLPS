

library(Biostrings)
library(tidyverse)
library(magrittr)

setwd('~/Desktop/projects/ncba/litao/N_phase_separation/')

allproteins <- readAAStringSet('~/data/gisaid/20200619/allprot0609.fasta')

N_proteins <- allproteins[grep('^N\\|', names(allproteins), value = T)]



width(N_proteins) %>% table() 



N_proteins_419 <- N_proteins[width(N_proteins) == 419]


#writeXStringSet(N_proteins_419, './results/N_proteins_419.fa')

N_proteins_419_matrix <- as.matrix(N_proteins_419)
st_n <- tibble()
st_n <- N_proteins_419_matrix[, 1] %>% table() %>% as_tibble()
for (i in seq(2, width(n_aln)[1])) {
    tmp <- N_proteins_419_matrix[, i] %>% table() %>% as_tibble()
    st_n <- full_join(st_n, tmp, by = '.')
}



names(st_n) <- c('AA', paste0('pos', 1:419))
st_n <- st_n[st_n$AA !='X',]


st_n_sort <- arrange_at(st_n, .vars =names(st_n)[2:420], desc)
st_n_sort_m <- as.matrix(st_n_sort[, 2:420])
#turn NA into 0:
st_n_sort_m[is.na(st_n_sort_m)] <- 0
rownames(st_n_sort_m) <- st_n_sort$AA

#write_csv(as_tibble(st_n_sort_m, rownames = 'AA'), "./results/AA_freq_N.csv")


#get conserved sequences:
conserved_aa <- apply(st_n_sort_m, 2, function(x){
    return(rownames(st_n_sort_m)[which.max(x)])
}) 

N_conserved <- paste(conserved_aa,collapse = '') %>% AAStringSet()
writeXStringSet(N_conserved, './results/N_conserved.fa')


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

#write_csv(mutation_freq_N,  './results/mutation_freq_N_pro.csv' )

#heatmap
#heatmap ----
library(pheatmap)

aa_freq_n <- read_csv('./results/AA_freq_N.csv')

aa_freq_n_m <- as.matrix(aa_freq_n[2:420])
rownames(aa_freq_n_m) <- aa_freq_n$AA

aa_freq_n_m[aa_freq_n_m > 5000] <- 0

pheatmap(aa_freq_n_m, 
         cluster_rows = F, cluster_cols = F, 
         labels_col = conserved_aa, cex.axis= 0.8, las = 2,
         show_colnames = F,
         main = 'mutation types of N protein of SARS-CoV2',
         ylab ='conserved AA sequences')

#validation ------------



















