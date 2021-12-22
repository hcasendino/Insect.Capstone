# Usage eDNAindex(x) 
# Co-opted from Moncho's code, but simplified syntax (only takes df, not each column as a string)

eDNAindex<- function(df){ 

  require(tidyverse)

      df %>% 
      group_by(sample, Hash, Bottle) %>%
      summarise (sumreads = sum(nReads)) %>%  # This sums technical replicates
      group_by(sample,Bottle) %>% 
      mutate (Tot = sum(sumreads),
              Row.prop = sumreads / Tot)  %>%         # This creates the proportion on each biological replicate    
      group_by(sample) %>% 
      mutate (nreps = length(unique(Bottle))) %>% 
      group_by(sample, Hash) %>% 
      summarise (mean.prop = sum (Row.prop) / max(nreps))   %>% #"Averaging ratios between Biological replicates"
      group_by (Hash) %>%
      mutate (Colmax = max (mean.prop),
              Normalized.reads = mean.prop / Colmax) %>%  # mean Hash proportion across bottles for a given sample, scaled to the max mean Hash proportion across samples
      dplyr::select(-Colmax, -mean.prop) -> output 

  left_join(output, df) -> output
 
  return (output)
}
