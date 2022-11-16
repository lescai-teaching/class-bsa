so_consequences <- readRDS(url("https://github.com/lescai-teaching/datasets_class/blob/master/reference/annotations/ranked_ensembl_consequences.RData?raw=true", "rb"))

get_from_ann <- function(annotation, element){
  data <- unlist(str_split(annotation, "\\|"))
  return(data[element])
}

get_rank <- function(consequence){
  rank = so_consequences$rank[which(so_consequences$SO_term %in% consequence)]
  return(rank)
}

get_most_severe_index <- function(annotations_list){
  consequences = unlist(lapply(annotations_list, get_from_ann, element = 2))
  ranks = unlist(lapply(consequences, get_rank))
  most_severe = which(min(ranks) %in% ranks)
  return(most_severe[1])
}

get_most_severe_consequence <- function(annotations_list){
  consequence <- tryCatch(
    {
      index <- get_most_severe_index(annotations_list)
      get_from_ann(annotations_list[[index]], element = 2)
    },
    error=function(cond){
      return(NA)
    }
  )
  return(consequence)
}

get_most_severe_gene <- function(annotations_list){
  gene <- tryCatch(
    {
      index <- get_most_severe_index(annotations_list)
      get_from_ann(annotations_list[[index]], element = 4)
    },
    error=function(cond){
      return(NA)
    }
  )
  return(gene)
}

get_most_severe_impact <- function(annotations_list){
  impact <- tryCatch(
    {
      index <- get_most_severe_index(annotations_list)
      get_from_ann(annotations_list[[index]], element = 3)
    },
    error=function(cond){
      return(NA)
    }
  )
  return(impact)
}
