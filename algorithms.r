library(pcalg)
library(igraph)
library(ggm)

source("util.r")
source("ricf.r")

#' Run FCI Algorithm
#'
#' This function runs the Fast Causal Inference (FCI) algorithm to estimate the
#' Partial Ancestral Graph (PAG) using the given sufficient statistics and 
#' independence test.
#'
#' @param suffStat A list containing the sufficient statistics.
#' @param indepTest A function for testing conditional independence.
#' @param alpha Significance level for the independence tests.
#' @param labels Labels for the variables.
#' @param verbose Logical indicating whether to print detailed output.
#' @param max_condition_set The maximum size of the conditioning set (default: 4).
#' @return A PAG object resulting from the FCI algorithm.
run_fci <- function(suffStat, 
                    indepTest, 
                    alpha, 
                    labels, 
                    verbose,
                    max_condition_set = 4) {
  #Run fci from pcalg package
  fci.pag = fci(suffStat, indepTest, alpha = alpha, labels = labels,
                verbose = verbose, type = "anytime", m.max = max_condition_set)
  return(fci.pag)
}

#' Run RICF for Gaussian Graphical Models
#'
#' This function runs the Recursive Iterative Conditional Fitting (RICF) 
#' algorithm for Gaussian Graphical Models.
#'
#' @param amat Adjacency matrix.
#' @param suffStat A list containing the sufficient statistics.
#' @param tol Tolerance for convergence (default: 1e-04).
#' @param beta_0 Initial beta parameter (default: NULL).
#' @param omega_0 Initial omega parameter (default: NULL).
#' @param max_iter Maximum number of iterations (default: 10000).
#' @return A list containing Beta, Sigma, BIC, log-likelihood, number of 
#'         parameters, Ohat, and Bhat.
ricf_ggm <- function(amat, 
                     suffStat, 
                     tol = 1e-04, 
                     beta_0 = NULL, 
                     omega_0 = NULL, 
                     max_iter = 10) {
  #Get parameters for RICF
  S = suffStat$C
  n = suffStat$n
  p = ncol(S)
  V = colnames(S)
  adj_ggm = transform_input_ggm(amat)
  py <- length(which(colSums(adj_ggm) != 0))
  n_par = sum(adj_ggm == 1) + sum(adj_ggm == 100) / 2 + py
  out = NULL
  #Run Ricf
  tryCatch({
    out = fitAncestralGraph(adj_ggm, S, n, beta_0, omega_0, tol, max_iter)
  }, error = function(e) e)
  if (is.null(out)) {
    return(list(NULL, NULL, Inf, NULL, NULL, NULL, NULL))
  }
  #Get parameters and BIC
  Beta = -1 * out$Bhat + diag(1, nrow(out$Bhat), ncol(out$Bhat))
  Sigma = out$Shat[V, V]
  ll_Sigma = -n / 2 * (log(det(Sigma)) + sum(diag(solve(Sigma) %*% S))) - (n * p) / 2 * log(2 * pi)
  BIC = -2 * ll_Sigma + log(n) * n_par
  return(list(Beta, Sigma, BIC, ll_Sigma, n_par, out$Ohat, out$Bhat))
}

#' Calculate True BIC
#'
#' This function calculates the Bayesian Information Criterion (BIC) for the 
#' true model.
#'
#' @param amat Adjacency matrix.
#' @param suffStat A list containing the sufficient statistics.
#' @param Beta The Beta matrix.
#' @param Sigma The Sigma matrix.
#' @return A list containing Beta, Sigma, BIC, log-likelihood, and 
#'         number of parameters.
true_bic <- function(amat, suffStat, Beta, Sigma) {
  S = suffStat$C
  n = suffStat$n
  p = ncol(S)
  V = colnames(S)
  adj_ggm = transform_input_ggm(amat)
  py <- length(which(colSums(adj_ggm) != 0))
  n_par = sum(adj_ggm == 1) + sum(adj_ggm == 100) / 2 + py
  ll_Sigma = -n / 2 * (log(det(Sigma)) + sum(diag(solve(Sigma) %*% S))) - (n * p) / 2 * log(2 * pi)
  BIC = -2 * ll_Sigma + 4 * log(n) * n_par
  return(list(Beta, Sigma, BIC, ll_Sigma, n_par, NULL, NULL))
}

#' Build PBI
#'
#' This function builds the possible bidirected edges (PBI) based on the 
#' adjacency matrix.
#'
#' @param adjmat Adjacency matrix.
#' @return A character vector of possible bidirected edges.
build_pbi <- function(adjmat) {
  #Get directed edges
  res = get_directed(adjmat)
  adjmat_directed = res[[2]]
  tmp_mat = adjmat_directed
  pbi = character()
  #iterate through lengths, l, in {2, ..., p}
  for (k in seq(nrow(adjmat))) {
    tmp_mat = tmp_mat %*% adjmat_directed
    #find walks of length l
    walks = (tmp_mat == 1) %>%
      as.table() %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      filter(Freq == TRUE) %>%
      mutate("edge" = paste(Var1, ",", Var2, sep = ""))
    if (nrow(walks) == 0) {
      break
    }
    #check for pbi
    for (j in seq(nrow(walks))) {
      from = walks[j, 1]
      to = walks[j, 2]
      if (res[[2]][from, to] == 1) {
        pbi = c(pbi, paste(from, to, sep = ","))
      }
    }
  }
  return(pbi)
}

#' Build PID
#'
#' This function builds the possible directed edges (PID) based on the 
#' adjacency matrix.
#'
#' @param adjmat Adjacency matrix.
#' @param max_path_length The maximum length of paths to consider (default: 20).
#' @return A character vector of possible directed edges.
build_pid <- function(adjmat, max_path_length = 20) {
  #Build undirected graph from bidirected structure
  V = rownames(adjmat)
  res = get_bidirected(adjmat)
  bimat = res[[2]]
  pid = pid_check_colliders(adjmat)
  tmp_mat = bimat
  graph = graph.adjacency(bimat, mode = "undirected")
  for (k in seq(2, max_path_length)) {
    tmp_mat = tmp_mat %*% bimat
    #ignore paths with the same start and end node
    diag(tmp_mat) = 0
    tmp_mat[tmp_mat != 0] = 1
    #Get walks
    walks = (tmp_mat >= 1) %>%
      as.table() %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      filter(Freq == TRUE) %>%
      mutate("edge" = paste(Var1, ",", Var2, sep = ""))
    if (nrow(walks) == 0) {
      break
    }
    #find edges to add to pid
    for (j in seq(nrow(walks))) {
      from = walks[j, 1]
      to = walks[j, 2]
      if (adjmat[from, to] %in% c(2, 3)) {
        next
      }
      v_from = names(which(adjmat[V[V != to], from] == 2))
      v_to = names(which(adjmat[V[V != to], to] == 2))
      for (i in v_from) {
        for (j in v_to) {
          if (i != j & adjmat[i, j] == 2) {
            pid = c(pid, paste(i, j, sep = ","))
          }
        }
      }
    }
  }
  return(unique(pid))
}

#' Identify possible directed and bidirected edges in a graph
#' 
#' This function identifies potential directed and bidirected edges based on the 
#' given adjacency matrix and node set.
#' 
#' @param adjmat The adjacency matrix representing the graph
#' @param V The set of nodes in the graph
#' @return A matrix containing potential directed and bidirected edges
possible_edges = function(adjmat, 
                          V)
{
  edges = matrix(NA, nrow = 0, ncol = 3)
  for(v in V)
  {
    avoid_nodes = c()
    de_v = pcalg::searchAM(adjmat, v, type = "de")
    sps_v = pcalg::searchAM(adjmat, v, type = "sp")
    avoid_nodes = c(avoid_nodes, de_v)
    pos_parents = V[!(V %in% avoid_nodes)]
    for (i in pos_parents)
    {
      sps = pcalg::searchAM(adjmat, i, type = "sp")
      an_i = pcalg::searchAM(adjmat, i, type = "an")
      if(adjmat[i, v] == 0 & !any(sps %in% de_v) & !(any(sps_v %in% an_i)))
      {
        edges = rbind(edges, c(i, v, "dir"))
        edges = rbind(edges, c(i, v, "bi"))
      }
    }
  }
  return(edges)
}

#' Score potential edges based on Bayesian Information Criterion (BIC)
#' 
#' This function scores potential edges using the Bayesian Information Criterion 
#' (BIC) based on the given adjacency matrix, edges, sufficient statistics, 
#' and optional parameters for initial conditions.
#' 
#' @param adjmat The adjacency matrix representing the graph
#' @param edges A matrix containing potential directed and bidirected edges
#' @param suffStat Sufficient statistics for the graph
#' @param omega_0 Initial value for the omega matrix
#' @param beta_0 Initial value for the regression coefficients
#' @param hot_start Boolean indicating whether to use hot start initialization
#' @return A list containing scores for edges, next values of omega, 
#' and next values of beta
score_edges = function(adjmat, 
                       edges, 
                       suffStat, 
                       omega_0 = NULL,
                       beta_0 = NULL, 
                       hot_start = F)
{
  adjmat_pros = adjmat
  scores = numeric(nrow(edges))
  min_score = Inf
  omega_next = NA
  beta_next = NA
  for(i in seq_along(scores))
  {
    e = edges[i,]
    from = as.numeric(e[1])
    to = as.numeric(e[2])
    type = e[3]
    if (type == "dir")
    {
      adjmat_prop = add_directed(adjmat, from, to)
    }
    else
    {
      adjmat_prop = add_bidirected(adjmat, from, to)
    }
    if(!suppressWarnings(isADMG(transform_input_ggm(adjmat_prop))))
    {
      scores[i] = Inf
      next
    }
    if(hot_start)
    {
      ricf_res = ricf_ggm(adjmat_prop, suffStat, omega_0 = omega_0, beta_0 = beta_0, tol = 1, max_iter = 10)
    }
    else
    {
      ricf_res = ricf_ggm(adjmat_prop, suffStat, tol =1, max_iter = 10)
    }
    scores[i] = ricf_res[[3]]
    if(min_score > scores[i])
    {
      min_score = scores[i]
      omega_next = ricf_res[[6]]
      beta_next = ricf_res[[7]]
    }
  }
  return(list(scores, omega_next, beta_next))
}


#' Perform greedy edge addition to improve graph structure
#' 
#' This function performs greedy edge addition to improve the structure of the 
#' graph based on the given adjacency matrix, node set, current BIC score, 
#' sufficient statistics, and optional parameters for initial conditions 
#' and control.
#' 
#' @param adjmat The adjacency matrix representing the graph
#' @param V The set of nodes in the graph
#' @param current_bic Current Bayesian Information Criterion (BIC) score
#' @param suffStat Sufficient statistics for the graph
#' @param omega_0 Initial value for the concentration matrix
#' @param beta_0 Initial value for the regression coefficients
#' @param hot_start Boolean indicating whether to use hot start initialization
#' @param passes Number of passes for the greedy algorithm
#' @param tol Tolerance for convergence criteria
#' @param max_row Maximum number of rows to consider in the edge addition process
#' @return The updated adjacency matrix after greedy edge addition
greedy_add = function(adjmat, 
                      V, 
                      current_bic, 
                      suffStat, 
                      omega_0 = NULL, 
                      beta_0 = NULL, 
                      hot_start = F, 
                      passes = 1, 
                      tol=0, 
                      max_row = Inf)
{
  #T = when there is an edge offering a decrease
  omega_t = omega_0
  beta_t = beta_0
  add = T
  n = suffStat$n
  adjmat_prop = adjmat
  last_bic = NULL
  omega_t = omega_0
  beta_t = beta_0
  #Identify possible edges and BIC with added edge
  edges = possible_edges(adjmat_prop, 1:nrow(suffStat[[1]]))
  if (nrow(edges) == 0)
  {
    return (adjmat_prop)
  }
  #score edges
  for(j in seq(passes))
  {
    print(paste0("pass ", j))
    if (hot_start)
    {
      diag(omega_t) = diag(suffStat[[1]])
    }
    res_sc = score_edges(adjmat_prop, edges, suffStat, omega_t, beta_t, hot_start = hot_start)
    edges = cbind(edges, res_sc[[1]])
    bic_dec = which(edges[,4] < current_bic) 
    edges_to_add = edges[bic_dec,, drop=FALSE]
    if(nrow(edges_to_add) == 0)
    {
      return (adjmat_prop)
    }
    #Select edges
    edges_to_add = edges_to_add[order(edges_to_add[,4]),]
    added = 0
    #Add edges
    for (i in seq(nrow(edges_to_add)))
    {
      if(added >= max_row)
      {
        break
      }
      e = edges_to_add[i,]
      from = as.numeric(e[1])
      to = as.numeric(e[2])
      type = e[3]
      aj = adjmat_prop
      if (type == "dir")
      {
        aj = add_directed(adjmat_prop, from, to)
      }
      else
      {
        aj = add_bidirected(adjmat_prop, from, to)
      }
      if (hot_start)
      {
        diag(omega_t) = diag(suffStat[[1]])
      }
      if(!isADMG(transform_input_ggm(aj)))
      {
        next
      }
      if(hot_start)
      {
        omega_t[from, to] = 0
        omega_t[to, from] = 0
        beta_t[from, to] = 0
        beta_t[to, from] = 0
        ricf_res = ricf_ggm(aj, suffStat, omega_0 = omega_t, beta_0 = beta_t, tol = 1, max_iter = 10)
      }
      else
      {
        ricf_res = ricf_ggm(aj, suffStat,  tol = 1, max_iter = 10)
      }
      if (current_bic > ricf_res[[3]])
      {
        if(!is.null(last_bic))
        {
          if( (abs(current_bic - last_bic)/n < tol))
          {
            next
          }
        }
        adjmat_prop = aj
        last_bic = current_bic
        current_bic = ricf_res[[3]]
        omega_t = ricf_res[[6]]
        beta_t = ricf_res[[7]]
        added = added + 1
      }
    }
  }
  return(adjmat_prop)
}

#' Perform greedy edge removal to improve graph structure
#' 
#' This function performs greedy edge removal to improve the 
#' structure of the graph based on the given adjacency matrix, node set, 
#' current BIC score, sufficient statistics, and optional parameters for control.
#' 
#' @param adjmat The adjacency matrix representing the graph
#' @param V The set of nodes in the graph
#' @param current_bic Current Bayesian Information Criterion (BIC) score
#' @param suffStat Sufficient statistics for the graph
#' @param tol Tolerance for convergence criteria
#' @return The updated adjacency matrix after greedy edge removal
greedy_remove = function(adjmat, V, current_bic, suffStat, tol = 0)
{
  dir = get_directed(adjmat)[[1]]
  bidir = get_bidirected(adjmat)[[1]]
  edges = rbind(dir, bidir)
  adjmat_prop = adjmat
  for(i in seq(nrow(edges)))
  {
    aj = adjmat_prop
    e_i = edges[i,]
    from = unlist(e_i[1])
    to = unlist(e_i[2])
    aj[from, to] = 0
    aj[to, from] = 0
    ricf_res = ricf_ggm(aj, suffStat,  tol = 1, max_iter = 10)
    if (current_bic > ricf_res[[3]])
    {
      adjmat_prop = aj
      last_bic = current_bic
      current_bic = ricf_res[[3]]
      omega_t = ricf_res[[6]]
      beta_t = ricf_res[[7]]
    }
  }
  return(adjmat_prop)
}


