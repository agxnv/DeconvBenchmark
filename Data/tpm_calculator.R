gtf <- rtracklayer::import('Mus_musculus.GRCm38.88.gtf')

length <- as.data.frame(gtf) %>% 
  filter(type == "exon") %>% 
  filter(gene_name %in% genes) %>% #"genes" its a vector that lists the genes in the scRNAseq batch
  group_by(gene_name, gene_version) %>%
  summarise(length = sum(width)) %>% #aggregate length by gene version
  mutate(length_types = sum(length)) %>%
  mutate(unique_types = n_distinct(gene_version)) %>% 
  mutate(length_final = length_types / unique_types) %>% #average the lengths of the present versions of that gene 
  ungroup()  %>%
  distinct(gene_name, .keep_all= TRUE) %>%
  arrange(gene_name) %>% 
  select(gene_name, length)

tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}



