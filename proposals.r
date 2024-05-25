source("algorithms.r")
#' Perform a phase in the MAG search algorithm
#' 
#' This function performs a phase in the MAG search algorithm based on the given 
#' proposal edges, current MAG, adjacency matrix,sufficient statistics, 
#' type of phase, temperature, and optional parameters for initial conditions 
#' and control.
#' 
#' @param proposal_edges Proposal edges for the phase
#' @param mag_t Current MAG (Magnitude Adjusted Graph)
#' @param adjmat The adjacency matrix representing the graph
#' @param suffStat Sufficient statistics for the graph
#' @param type Type of phase (1: PBI, 2: PDI, 3: Reverse PBI, 4: Reverse PDI)
#' @param temperture Temperature for the Metropolis-Hastings algorithm
#' @param omega_0 Initial value for the concentration matrix
#' @param beta_0 Initial value for the regression coefficients
#' @param edge_d0 Initial set of directed edges
#' @param edge_b0 Initial set of bidirected edges
#' @param hot_start Boolean indicating whether to use hot start initialization
#' @return A list containing information about the accepted phase, including MAG, 
#'      edge transfer, and force accept status
phase <- function(proposal_edges, 
                  mag_t, 
                  adjmat, 
                  suffStat, 
                  type, 
                  temperture,
                  omega_0 = NULL,
                  beta_0 = NULL,
                  edge_d0 = NULL, 
                  edge_b0 = NULL,
                  hot_start = FALSE)
{
    force_accept = FALSE  
    edge_accepted = NA
    ricf_res_accepted = NA
    bic_accepted = Inf
    mag_accepted = mag_t
    edge_transfer_accepted = NA
    unif_r = runif(1, 0, 1)
    if(unif_r < temperture)
    {
      force_accept = TRUE
      proposal_edges = sample(proposal_edges, 1)
    }
    for (edge in proposal_edges)
    {
        beta_e = beta_0
        omega_e = omega_0
        diag(omega_e) = diag(suffStat[[1]])
        edge_transfer_proposal = NA
        adjmat_proposal = NA
        ed = str_split(edge, ",")[[1]]
        if (type == 1)
        {
            #PBI: Directed -> Bidirected
            adjmat_proposal = add_bidirected(adjmat, ed[1], ed[2])
            edge_transfer_proposal = c("dir", "bidir")
            if(hot_start)
            {
              omega_e[ed[1], ed[2]] = 0
              omega_e[ed[2], ed[1]] = 0
              beta_e[ed[1], ed[2]] = 0
              beta_e[ed[2], ed[1]] = 0
            }
        }
        else if(type == 2)
        {
            #PDI: Edge -> No Edge
            res = remove_edge(adjmat, ed[1], ed[2])
            adjmat_proposal = res[[1]]
            edge_transfer_proposal = res[[2]]
            if(hot_start)
            {
              beta_e[ed[1], ed[2]] = 0
              beta_e[ed[2], ed[1]] = 0
              omega_e[ed[1], ed[2]] = 0
              omega_e[ed[2], ed[1]] = 0
            }
        }
        else if(type == 3)
        {
            #Reverse PBI : Bidir -> dir
            adjmat_proposal = add_directed(adjmat, ed[1], ed[2])
            edge_transfer_proposal = c("bidir", "dir")
            if(hot_start)
            {
              omega_e[ed[1], ed[2]] = 0
              omega_e[ed[2], ed[1]] = 0
            }
        }
        else if(type == 4)
        {
            #Reverse PID : No Edge -> Edge
            if(is.null(edge_d0) | is.null(edge_b0))
            {
                stop("Needs Initial Edge Sets")
            }
            if (edge %in% edge_d0)
            {
                adjmat_proposal = add_directed(adjmat, ed[1], ed[2])
                edge_transfer_proposal = c("none", "dir")
            }
            else
            {       
                adjmat_proposal = add_bidirected(adjmat, ed[1], ed[2])
                edge_transfer_proposal = c("none", "bidir")
            }
        }

        if(hot_start)
        {
          ricf_proposal = ricf_ggm(adjmat_proposal, 
                                   suffStat, 
                                   omega_0 = omega_e,
                                   beta_0 = beta_e, 
                                   tol = 1, 
                                   max_iter = 10)
        }
        else
        {
          ricf_proposal = ricf_ggm(adjmat_proposal, 
                                   suffStat, 
                                   tol = 1, 
                                   max_iter = 10)
        }
        bic_proposal = ricf_proposal[[3]]
        if(is.na(bic_proposal))
        {
            next
        }
        else if (bic_proposal <= bic_accepted)
        {
            bic_accepted = bic_proposal
            ricf_res_accepted = ricf_proposal
            edge_accepted = edge
            attributes(mag_accepted)$amat = adjmat_proposal
            edge_transfer_accepted  = edge_transfer_proposal
        }
    }
    return(list(ricf_res_accepted, 
                mag_accepted, 
                edge_accepted,
                edge_transfer_accepted,
                force_accept))
}

#' Update MAG and related variables based on results
#' 
#' This function updates the MAG and related variables based on the results 
#' of a phase in the MAG search algorithm.
#' 
#' @param results Results from the phase including accepted information
#' @param edge_d Current set of directed edges
#' @param edge_b Current set of bidirected edges
#' @param current_bic Current Bayesian Information Criterion (BIC) score
#' @return A list containing updated MAG and related variables if an update 
#'         is performed, NULL otherwise
update_mag = function(results, 
                      edge_d, 
                      edge_b, 
                      current_bic)
{
    ricf_res_accepted = results[[1]]
    bic_update = ricf_res_accepted[[3]]
    mag_accepted = results[[2]]
    edge_accepted = results[[3]]
    edge_transfer_accepted = results[[4]]
    force_accept = results[[5]]
    update = FALSE
    if ( ((bic_update >= current_bic) | is.na(edge_accepted)) & !force_accept)
    {
        return(list("update" = update, "edge_d" = NULL, "edge_b" = NULL, 
                "beta" = NULL, "omega" = NULL,
                "bic" = NULL, "mag" = NULL))    
    }
    else 
    {
        update = TRUE
        edge_dnew = edge_d
        edge_bnew = edge_b
        remove_edge = edge_transfer_accepted[1]
        add_edge = edge_transfer_accepted[2]

        if(add_edge == "dir")
        {
            edge_dnew = c(edge_d, edge_accepted)
        }
        else if(add_edge == "bidir")
        {
            rev_ed = reverse_edge(edge_accepted)
            edge_bnew = c(edge_b, edge_accepted, rev_ed)
        }

        if(remove_edge == "dir")
        {
            edge_dnew = edge_d[edge_d != edge_accepted]
        }
        else if(remove_edge == "bidir")
        {
            rev_ed = reverse_edge(edge_accepted)
            edge_bnew = edge_b[edge_b != edge_accepted & edge_b != rev_ed]
        }
        return(list("update" = update, "edge_d" = edge_dnew, "edge_b" = edge_bnew, 
                "beta" = ricf_res_accepted[[7]], "omega" = ricf_res_accepted[[6]],
                "bic" = ricf_res_accepted[[3]], "mag" = mag_accepted, 
                "likeli" = ricf_res_accepted[[4]], "npar" = ricf_res_accepted[[5]]))
    }

}

#' Perform the MAG search algorithm
#' 
#' This function performs the MAG search algorithm to find the optimal graph 
#' structure based on the given initial MAG, node set, sufficient statistics, 
#' and optional parameters for control and tuning.
#' 
#' @param fci.mag Initial MAG (Magnitude Adjusted Graph) to start the search
#' @param V The set of nodes in the graph
#' @param suffStat Sufficient statistics for the graph
#' @param max_iter Maximum number of iterations for the search algorithm
#' @param temperture Temperature for the Metropolis-Hastings algorithm
#' @param cooldown Cooldown factor for decreasing temperature
#' @param test Boolean indicating whether to run a test search 
#' @param hot_start Boolean indicating whether to use hot start initialization 
#'      in the search
#' @return A list containing the updated MAG, Bayesian Information Criterion 
#'   (BIC) score, accepted phases, likelihood, and number of parameters
mag_search = function(fci.mag, 
                      V, 
                      suffStat, 
                      max_iter = 10, 
                      temperture = 0,
                      cooldown = 0,
                      test = FALSE,
                      hot_start = FALSE)
{  

  #Init Graphical Variables 
  mag_accepted= fci.mag
  adjmat = attributes(fci.mag)$amat
  edge_d = get_directed(adjmat)[[1]]$edge
  edge_b = get_bidirected(adjmat)[[1]]$edge
  edge_d0 = edge_d
  edge_b0 = edge_b
  #Initalize SEM Params
  pbi = build_pbi(adjmat)
  print("BUILT PDI")
  pid = build_pid(adjmat)
  print("BUILT PID")
  #Get inital estimates
  res_init = ricf_ggm(adjmat, suffStat)
  if(is.na(res_init[[3]]))
  {
    return(list(NA, Inf))
  }
  beta_accepted = res_init[[7]]
  omega_accepted = res_init[[6]]
  bic_accepted = res_init[[3]]
  likeli_accepted = res_init[[4]]
  par_accepted = res_init[[5]]
  #Set up iterator
  it = 1
  accepted_props = numeric(0)
  while(it <= max_iter)
  {
    changes = 0
    #phase 1
    #Find proposal edges 
    pbi_prop = edge_d[edge_d %in% pbi] 
    if(length(pbi_prop) > 0)
    {
      #Run phase
      out = phase(pbi_prop, 
                  mag_accepted, 
                  adjmat, 
                  suffStat, 
                  1, 
                  temperture,
                  omega_0 = omega_accepted, 
                  beta_0 = beta_accepted, 
                  hot_start = hot_start)
      updates = update_mag(out, edge_d, edge_b, bic_accepted)
      if(updates[["update"]])
      {
        print("INNER UPDATE")
        edge_d = updates[["edge_d"]]
        edge_b = updates[["edge_b"]]
        omega_accepted = updates[["omega"]]
        beta_accepted = updates[["beta"]]
        bic_accepted = updates[["bic"]]
        print(bic_accepted)
        mag_accepted = updates[["mag"]]
        likeli_accepted = updates[["likeli"]]
        par_accepted = updates[["npar"]]
        adjmat = attributes(mag_accepted)$amat
        changes = changes + 1
        accepted_props = c(accepted_props, 1)
        if(test)
        {
          return(list(mag_accepted, bic_accepted))
        }
      }
    }
    #phase 2
    #Find proposal edges 
    pid_prop = c(edge_d, edge_b)[c(edge_d, edge_b) %in% pid]
    if(length(pid_prop) > 0)
    {
      #Run phase 2
      out = phase(pid_prop, 
                  mag_accepted, 
                  adjmat, 
                  suffStat, 
                  2,
                  temperture,  
                  omega_0 = omega_accepted,
                  beta_0 = beta_accepted,
                  hot_start = hot_start)
      updates = update_mag(out, edge_d, edge_b, bic_accepted)
      if(updates[["update"]])
      {
        print("INNER UPDATE")
        edge_d = updates[["edge_d"]]
        edge_b = updates[["edge_b"]]
        omega_accepted = updates[["omega"]]
        beta_accepted = updates[["beta"]]
        bic_accepted = updates[["bic"]]
        print(bic_accepted)
        mag_accepted = updates[["mag"]]
        likeli_accepted = updates[["likeli"]]
        par_accepted = updates[["npar"]]
        adjmat = attributes(mag_accepted)$amat
        accepted_props = c(accepted_props, 2)
        changes = changes + 1
        if(test)
        {
          return(list(mag_accepted, bic_accepted))
        }        
      }
    }
    #phase 3
    #Get proposal edges
    pbi_rev_prob = edge_b[edge_b %in% pbi]
    if(length(pbi_rev_prob) > 0)
    {
      #Run phase 3
      out = phase(pbi_rev_prob, 
                  mag_accepted, 
                  adjmat, 
                  suffStat, 
                  3, 
                  temperture,
                  omega_0 = omega_accepted,
                  beta_0 = beta_accepted,
                  hot_start = hot_start)
      updates = update_mag(out, edge_d, edge_b, bic_accepted)
      if(updates[["update"]])
      {
        print("INNER UPDATE")
        edge_d = updates[["edge_d"]]
        edge_b = updates[["edge_b"]]
        omega_accepted = updates[["omega"]]
        beta_accepted = updates[["beta"]]
        bic_accepted = updates[["bic"]]
        print(bic_accepted)
        mag_accepted = updates[["mag"]]
        likeli_accepted = updates[["likeli"]]
        par_accepted = updates[["npar"]]
        adjmat = attributes(mag_accepted)$amat
        accepted_props = c(accepted_props, 3)
        changes = changes + 1
        if(test)
        {
          return(list(mag_accepted, bic_accepted))
        }        
      }
    }
    #phase 4
    #Get proposal edges
    pid_rev_add = pid[!(pid %in% c(edge_d, edge_b))]
    if(length(pid_rev_add) > 0)
    {
      #Run phase 4
      out = phase(pid_rev_add, 
                  mag_accepted, 
                  adjmat, 
                  suffStat,
                  4, 
                  temperture,
                  edge_d0,
                  edge_b0, 
                  omega_0 = omega_accepted,
                  beta_0 = beta_accepted,
                  hot_start = hot_start)
      updates = update_mag(out, edge_d, edge_b, bic_accepted)
      if(updates[["update"]])
      {
        print("INNER UPDATE")
        edge_d = updates[["edge_d"]]
        edge_b = updates[["edge_b"]]
        omega_accepted = updates[["omega"]]
        beta_accepted = updates[["beta"]]
        bic_accepted = updates[["bic"]]
        print(bic_accepted)
        mag_accepted = updates[["mag"]]
        likeli_accepted = updates[["likeli"]]
        par_accepted = updates[["npar"]]
        adjmat = attributes(mag_accepted)$amat
        accepted_props = c(accepted_props, 4)
        changes = changes + 1
        if(test)
        {
          return(list(mag_accepted, bic_accepted))
        }        
      }
    }
    
    it = it + 1
    temperture = temperture * cooldown
    if (changes == 0)
    {
      break
    }
  }
  return(list(mag_accepted, 
              bic_accepted, 
              accepted_props, 
              likeli_accepted, 
              par_accepted))
}
