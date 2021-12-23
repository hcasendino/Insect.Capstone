# Usage eDNAindex(x); x is dataframe with columns: Sample, Hash, nReads, Biological.replicate (and Reach)
# Co-opted from Moncho's code, but simplified syntax (only takes df, not each column as a string)

eDNAindex<- function(df, reach = F){ 

  require(tidyverse)
  
  df %>% 
      group_by(Sample, Hash, Reach, Biological.replicate) %>%
      summarise (sumreads = sum(nReads)) %>%  # This sums technical replicates
      group_by(Sample,Biological.replicate) %>% 
      mutate (Tot = sum(sumreads),
              Row.prop = sumreads / Tot)  %>%         # This creates the proportion on each biological replicate    
      group_by(Sample) %>% 
      mutate (nreps = length(unique(Biological.replicate))) %>% 
      
    {if (reach == TRUE){group_by(., Sample, Reach, Hash)}  
      else{group_by(., Sample, Hash)} }  %>% 
    
      summarise (mean.prop = sum (Row.prop) / max(nreps))   %>% #"Averaging ratios between Biological replicates"
      group_by (Hash) %>%
      mutate (Colmax = max (mean.prop),
              Normalized.reads = mean.prop / Colmax) %>%  # mean Hash proportion across bottles for a given sample, scaled to the max mean Hash proportion across samples
      dplyr::select(-Colmax, -mean.prop) -> output 

  left_join(output, df) -> output
 
  return (output)
}
