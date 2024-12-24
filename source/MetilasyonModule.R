library(Biostrings)
library(ggplot2)

# CpG Adalar?? Bulma Fonksiyonu
findCpGIslands <- function(dna_sequence) {
  # CpG adalar??n??n uzunlu??u ve yo??unlu??u i??in e??ik de??erler
  min_length <- 200
  min_gc_content <- 0.5
  min_obs_exp_ratio <- 0.6
  
  # CpG adalar?? bulma
  cpg_islands <- find_CpG(dna_sequence, min_length, min_gc_content, min_obs_exp_ratio)
  
  if (length(cpg_islands) == 0) {
    return(data.frame(Start = integer(), End = integer(), GC_Content = numeric(), ObsExpRatio = numeric()))
  }
  
  return(cpg_islands)
}

# CpG Adalar??n?? G??rselle??tirme
plotCpGIslands <- function(cpg_islands) {
  if (nrow(cpg_islands) == 0) {
    return(NULL)
  }
  
  ggplot(cpg_islands, aes(x = Start, xend = End, y = 1, yend = 1)) +
    geom_segment(size = 2, color = "blue") +
    labs(title = "CpG Islands", x = "Position", y = "") +
    theme_minimal()
}
