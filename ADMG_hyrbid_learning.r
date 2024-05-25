source("proposals.r")

#' ADMG Hybrid Learning Algorithm
#'
#' This function performs hybrid learning for Ancestral Directed Mixed Graphs 
#' (ADMGs) using a combination of constraint-based and score-based methods. It 
#' integrates Gaussian Conditional Independence tests and iterates through 
#' possible models to find the best structure according to the 
#' Bayesian Information Criterion (BIC).
#'
#' @param data A dataframe or matrix where rows are observations and columns 
#'                  are variables.
#' @param CItest A string specifying the type of conditional independence 
#'                  test to use ("gaussian").
#' @param alpha_fci A numeric value for the significance level of the FCI 
#'                  algorithm.
#' @param max_iter An integer specifying the maximum number of iterations for 
#'                  the search algorithm.
#' @param temperture A numeric value for the initial temperature in the 
#'                  simulated annealing process.
#' @param cooldown A numeric value for the cooldown rate in the simulated 
#'                  annealing process.
#' @param verbose A logical value indicating whether to print detailed messages 
#'                  during execution.
#' @param hot_start A logical value indicating whether to use a hot start for 
#'                  optimization procedures.
#' @param alpha_adjust A logical value indicating whether to adjust alpha based 
#'                  on the sample size.
#' @param max_greedyadd An integer specifying the maximum number of greedy 
#'                  additions.
#' @param final.MLE.iter An integer specifying the maximum number of iterations 
#'                  for the final RICF CALL. Results in more specific estimates 
#' @return A list containing the final BIC score, the final MAG, the acceptance 
#'         rates, the final likelihood, the final number of parameters, 
#'         the final sigma matrix, and the final beta matrix.
ADMG_hybrid_learning = function(data, 
                                CItest = "gaussian", 
                                alpha_fci = 0.2, 
                                max_iter = 10,
                                temperture = 0,
                                cooldown = 0,
                                verbose = FALSE,
                                hot_start = TRUE,
                                alpha_adjust = FALSE,
                                max_greedyadd = 500,
                                final.MLE.iter = 100) {
  # Initialize final results with default values
  final_bic = Inf
  final_likeli = -1
  final_npar = -1
  final_mag = NULL
  indepTest = NA
  
  # Select the appropriate conditional independence test
  if (CItest == "gaussian") {
    indepTest = gaussCItest
  } else {
    stop("Not Implemented")
  }
  
  # Prepare the data and sufficient statistics
  n = nrow(data)
  V = colnames(data)
  cov_dat = cov(data)
  cor_dat = cov2cor(cov_dat)
  suffStat = list("C" = cor_dat, "n" = n)
  
  # Adjust alpha if required
  if (alpha_adjust) {
    alpha = alpha_fci / sqrt(n)
  } else {
    alpha = alpha_fci
  }
  
  # Run FCI to establish CI constraints
  fci.pag = run_fci(suffStat, indepTest, alpha = alpha, labels = V, verbose = verbose)
  
  # Build search space using MAGs
  mags = list()
  nodes = colnames(data)
  p_amat = attributes(fci.pag)$amat
  i = 1
  for (n_i in nodes) {
    mags[[i]] = pag2magAM(p_amat, n_i)
    i = i + 1
  }
  search_space = unique(mags)
  
  # Iterate through search space and run the proposal method on each MAG
  suffStat$C = cov_dat
  for (m in search_space) {
    res_i = ricf_ggm(m, suffStat)
    m = greedy_add(m, 1:nrow(suffStat), res_i[[3]], suffStat, omega_0 = res_i[[6]], beta_0 = res_i[[7]], hot_start = hot_start, max_row = max_greedyadd)
    if (verbose) print("GREEDY DONE")
    
    fci.mag = fci.pag
    attributes(fci.mag)$amat = m
    if (verbose) print("MAG SEARCH")
    
    res_m = mag_search(fci.mag, V, suffStat, max_iter = max_iter, temperture = temperture, cooldown = cooldown, hot_start = hot_start)
    if (verbose) print("MAG SEARCH DONE")
    
    if (res_m[[2]] < final_bic) {
      final_bic = res_m[[2]]
      final_mag = res_m[[1]]
      final_accept = res_m[[3]]
      final_likeli = res_m[[4]]
      final_npar = res_m[[5]]
    }
  }
  
  if (verbose) print("GREEDY REM")
  
  # Final greedy removal and MLE optimization
  m = greedy_remove(attributes(final_mag)$amat, colnames(suffStat$C), final_bic, suffStat)
  
  if (verbose) print("RICF Start")
  
  res_f = ricf_ggm(m, suffStat, max_iter = final.MLE.iter)
  
  if (verbose) print("RICF END")
  
  final_bic = res_f[[3]]
  final_likeli = res_f[[4]]
  final_npar = res_f[[5]]
  attributes(final_mag)$amat = m
  final_sigma = res_f[[2]]
  final_beta = res_f[[1]]
  
  # Return the final results
  return(list(final_bic, 
              final_mag, 
              final_accept, 
              final_likeli, 
              final_npar,
              final_sigma,
              final_beta))
}
