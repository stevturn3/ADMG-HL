library(tidyverse)

###### Utilities

#' Get Bidirected Edges and Undirected Graph with Bidirected Structure
#' 
#' This function identifies bidirected edges in an adjacency matrix and returns 
#' both the list of bidirected edges and the updated adjacency matrix with only 
#' bidirected edges.
#'
#' @param adjmat A square adjacency matrix.
#' @return A list containing a dataframe of bidirected edges and the bidirected 
#'         adjacency matrix.
get_bidirected = function(adjmat) {
  # Initialize a zero matrix for bidirected edges
  bidir = diag(0, nrow(adjmat), ncol(adjmat))
  colnames(bidir) = colnames(adjmat)
  rownames(bidir) = rownames(adjmat)
  
  # Identify bidirected edges 
  bidir[adjmat == 2 & t(adjmat) == 2] = 1
  
  # Create a dataframe of bidirected edges
  edges = (bidir == 1 & t(bidir) == 1) %>% 
    as.table() %>% 
    as.data.frame(., stringsAsFactors=FALSE) %>% 
    filter(Freq == TRUE) %>% 
    mutate("edge" = paste(Var1, ",", Var2, sep = ""))
  
  return(list(edges, bidir))  
}

#' Get Directed Edges and Directed Structure Adjacency Matrix
#' 
#' This function identifies directed edges in an adjacency matrix and returns 
#' both the list of directed edges and the directed adjacency matrix with only 
#' directed edges.
#'
#' @param adjmat A square adjacency matrix.
#' @return A list containing a dataframe of directed edges and the modified 
#'         adjacency matrix.
get_directed = function(adjmat) {
  # Initialize a zero matrix for directed edges
  dir = diag(0, nrow(adjmat), ncol(adjmat))
  colnames(dir) = colnames(adjmat)
  rownames(dir) = rownames(adjmat)
  
  # Identify directed edges 
  dir[(adjmat == 2 & t(adjmat) == 3)] = 1
  
  # Create a dataframe of directed edges
  edges = (dir == 1 & t(dir) == 0) %>% 
    as.table %>% 
    as.data.frame(., stringsAsFactors=FALSE) %>% 
    filter(Freq == TRUE) %>% 
    mutate("edge" = paste(Var1, ",", Var2, sep = ""))
  
  return(list(edges, dir))
}

#' Check for Colliders in Adjacency Matrix
#' 
#' This function checks for colliders (nodes with multiple parents) in an 
#' adjacency matrix.
#'
#' @param adjmat A square adjacency matrix.
#' @return A character vector of colliders.
pid_check_colliders = function(adjmat) {
  # Identify nodes that are colliders
  colliders = names(which(colSums(adjmat == 2) > 1))
  pid_collider = character()
  
  # Check for colliders in the adjacency matrix
  for (c in colliders) {
    v_c = names(which(adjmat[,c] == 2))
    for (i in v_c) {
      for (j in v_c) {
        if (i != j & adjmat[i,j] == 2) {
          pid_collider = c(pid_collider, paste(i, j, sep=","))
        }
      }
    }
  }
  return(pid_collider)
}

#' Transform Adjacency Matrix to GGM Format
#' 
#' This function transforms an adjacency matrix for use in Gaussian graphical 
#' models.
#'
#' @param adj A square adjacency matrix (0,2,3)
#' @return The transformed adjacency matrix (0, 1, 10, 100)
transform_input_ggm = function(adj) {
  adj_ggm = adj
  t_adj_ggm = t(adj_ggm)
  
  # Identify head and tail edges
  tail = adj_ggm == 3
  head = adj_ggm == 2
  t_tail = t_adj_ggm == 3
  t_head = t_adj_ggm == 2
  
  # Transform adjacency matrix according to GGM rules
  adj_ggm[head & t_tail] = 1
  adj_ggm[tail & t_head] = 0
  adj_ggm[tail & t_tail] = 10
  adj_ggm[head & t_head] = 100
  
  return(adj_ggm)
}

#' Transform Adjacency Matrix to pcalg Format
#' 
#' This function transforms an adjacency matrix for use in the PC algorithm.
#'
#' @param adj A square adjacency matrix (0, 1, 10, 100)
#' @return The transformed adjacency matrix (0, 2, 3)
transform_input_pcalg = function(adj) {
  adj_pc = adj
  t_adj_pc = t(adj_pc)
  
  # Identify head and bidirected edges
  head = adj_pc == 1
  bidir = adj_pc == 100
  
  # Transform adjacency matrix according to PC algorithm rules
  adj_pc[head] = 2
  adj_pc[t(head)] = 3
  adj_pc[bidir] = 2
  
  return(adj_pc)
}

#' Add Bidirected Edge
#' 
#' This function adds a bidirected edge between two nodes in an adjacency matrix.
#' Will overwrite existing edges.
#'
#' @param adjmat A square adjacency matrix.
#' @param node_from The starting node.
#' @param node_to The ending node.
#' @return The updated adjacency matrix.
add_bidirected = function(adjmat, 
                          node_from, 
                          node_to) {
  aj = adjmat
  aj[node_from, node_to] = 2
  aj[node_to, node_from] = 2
  return(aj)
}

#' Add Directed Edge
#' 
#' This function adds a directed edge between two nodes in an adjacency matrix.
#' Will overwrite existing edges.
#'
#' @param adjmat A square adjacency matrix.
#' @param node_from The starting node.
#' @param node_to The ending node.
#' @return The updated adjacency matrix.
add_directed = function(adjmat, 
                        node_from, 
                        node_to) {
  aj = adjmat 
  aj[node_from, node_to] = 2
  aj[node_to, node_from] = 3
  return(aj)
}

#' Remove Edge
#' 
#' This function removes an edge between two nodes in an adjacency matrix.
#'
#' @param adjmat A square adjacency matrix.
#' @param node_from The starting node.
#' @param node_to The ending node.
#' @return A list containing the updated adjacency matrix and the type of edge 
#'         removed.
remove_edge = function(adjmat, 
                       node_from, 
                       node_to) {
  aj = adjmat 
  edge_transfer = NA
  
  # Check if the edge is bidirected
  bidir = aj[node_from, node_to] == 2 & aj[node_to, node_from] == 2
  if (bidir) {
    edge_transfer = c("bidir", "none")
  } else {
    edge_transfer = c("dir", "none")
  }
  
  # Remove the edge
  aj[node_from, node_to] = 0
  aj[node_to, node_from] = 0
  
  return(list(aj, edge_transfer))
}

#' Reverse Edge
#' 
#' This function reverses the direction of an edge represented as a string.
#'
#' @param edge_string A string representing an edge (e.g., "A,B").
#' @return A string representing the reversed edge (e.g., "B,A").
reverse_edge = function(edge_string) {
  ed = strsplit(edge_string, split = ",")[[1]]
  rev_ed = paste0(ed[2], ",", ed[1])
  return(rev_ed)
}

#' Get Blanket Edges
#' 
#' This function gets the blanket edges for each node in an adjacency matrix.
#'
#' @param adjmat A square adjacency matrix.
#' @return A list of blanket edges for each node.
get_blanket_edges = function(adjmat) {
  V = colnames(adjmat)
  
  # Get the directed graph from the adjacency matrix
  graph = get_directed(adjmat)[[2]]
  graph = graph.adjacency(graph, mode = "directed")
  
  total_paths = list()
  
  # Find all simple paths for each node
  for (v in V) {
    paths = all_simple_paths(graph, v)
    p_v = names(sapply(paths, function(x) x[length(x)]))
    if (length(p_v) != 0) {
      total_paths[[v]] = p_v
    } else {
      total_paths[[v]] = numeric(0)
    }
  }
}
