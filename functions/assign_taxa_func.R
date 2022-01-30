# script created 1/30/22

assign.taxa <- function(asv_table, annotations, taxa_col = annotations$order){
  
  asv_taxa_df <- asv_table %>% mutate(class = NA, taxon = NA) 
  identified_hashes <- rep(NA, nrow(asv_taxa_df))
  
  for(i in 1:nrow(asv_taxa_df)){
    taxa_Row <- which(annotations$representative %in% asv_taxa_df$Hash[i])
    
    if(length(taxa_Row) > 0){
      asv_taxa_df$class[i] <- annotations$class[taxa_Row]
      asv_taxa_df$taxon[i] <- taxa_col[taxa_Row]
      identified_hashes[i] <- i
    }
    
    if((i %% 1000) == 0){
      print(i/nrow(asv_taxa_df))
    }
  }
  return(asv_taxa_df)
}
