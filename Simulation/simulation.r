rm(list = ls())

library(greedyBAPs)
library(lrpsadmm)
library(mvtnorm)
library(bnlearn)
library(mgcv)
source("../ADMG_hyrbid_learning.r", chdir = TRUE)
source("../algorithms.r", chdir = TRUE)

source("./utils_lrps/generate_data_for_GES.R")
source("./utils_lrps/run_GES.R")
source("./utils_lrps/compute_metrics.R")
source("Nowzohour_sim_and_equivclass.R")

calculate_f1_score = function(TP, FP, FN, R)
{
  return( (2 * TP) / (2 * TP + FP + FN + 2 * R))
}

decode_fromto = function(i,j, x)
{
  res = c()
  if(x[i, j] == 100 & x[j, i] == 100)
  {
    res = c(res, "<->")
  }
  else if(x[i, j] == 1 & x[j, i] == 0)
  {
    res = c(res, "->")
  }
  else if(x[i, j] == 0 & x[j, i] == 1)
  {
    res = c(res, "<-")
  }
  else if(x[i, j] == 0 & x[j, i] == 0)
  {
    res = c(res, "N")
  }
  return(res)
}

compare_equivclass = function(mg, true_mg, p, n, covMat, data, eps = 1e-5, maxIter = 1000)
{
  #FP = (Skeleton Comparison) True equivalence class denotes that there should be an edge between (i,j), not present in estimated equivalence class
  #FN = (Skeleton Comparison) Estimated equivalence denotes that there should be an edge between (i, j), not present in the true equivalence class
  #R = (Orientation Comparison) Estimated equivalence class include an orientation of the present edge (i,j) that is not present in true equivalence class
  #TP = For an existin edge (i,j), all orientations present in the estimated equivalence class are present in the true equivalence class
  depth.max <- p*(p-1)/2
  tmp <- findEquivalentModels(list(mg=mg), c(), list(),
                              NA, eps, depth.max=depth.max,
                              n=n, maxIter=maxIter, covMat=covMat,
                              faithful.eps=0)
  scores <- tmp$scores
  mg_equiv <- unique(c(lapply(tmp$res, function(x) x[[1]]), list(mg)))
  tmp <- findEquivalentModels(list(mg=true_mg), c(), scores, NA, eps,
                              depth.max=depth.max, time.max=Inf,
                              n=n, maxIter=maxIter, covMat=covMat, faithful.eps=0)
  true_mg_equiv = unique(c(lapply(tmp$res, function(x) x[[1]]), list(true_mg)))
  rev_true = list()
  rev_estim = list()
  for(i in seq(p))
  {
    if(i + 1 > p)
    {
      next
    }
    for(j in seq(i+1,p))
    {
      true_options = unique(sapply(true_mg_equiv, function(x) decode_fromto(i, j, x)))
      estim_options = unique(sapply(mg_equiv, function(x) decode_fromto(i, j, x)))
      rev_true[[paste0(i, j)]] = true_options
      rev_estim[[paste0(i, j)]] = estim_options
    }
  }
  edges = c("->", "<-", "<->")
  R = 0
  FP = 0
  TP = 0
  FN = 0
  TN = 0
  for (i in names(rev_true))
  {
    rev_estim[[i]] = sort(unlist(rev_estim[[i]]))
    rev_true[[i]] = sort(unlist(rev_true[[i]]))
    if (identical(rev_estim[[i]], c("N")) & identical(rev_true[[i]], c("N")))
    {
      TN = TN + 1
    } else if (identical(rev_true[[i]], rev_estim[[i]])) 
    {
      TP = TP + 1
    } else if (identical(rev_estim[[i]], c("N")) & !(identical(rev_true[[i]], c("N")))
    ) {
      FN = FN + 1
    } else if ((identical(rev_true[[i]], c("N"))) & !(identical(rev_estim[[i]], c("N"))))
    {
      FP = FP + 1
    } else {
      R = R + 1
    }
  }
  f1 = calculate_f1_score(TP, FP, FN, R)
  return(list(f1, TP, TN, FP, FN, R, list(mg_equiv), list(true_mg_equiv), rev_estim, rev_true))
}

metrics = function(method_bics, true_bics, 
                   method_likeli, true_likeli,
                   method_test_likeli, true_test_likeli,
                   method_npar, true_npar, f1,times, n = 100000, eps = 10e-6)
{
  methods = length(method_bics)
  MAD_bic = numeric(methods)
  MSE_bic = numeric(methods)
  exp_diff_bic = numeric(methods)
  MAD_test_likeli = numeric(methods)
  test_likeli = numeric(methods)
  MAD_likeli = numeric(methods)
  MSE_likeli = numeric(methods)
  exp_diff_likeli = numeric(methods)
  p_likeli = numeric(methods)
  p_bic = numeric(methods)
  avg_bic = numeric(methods)
  avg_likeli = numeric(methods)
  var_test_likeli = numeric(methods)
  npar_diff = numeric(methods)
  var_diff_bic = numeric(methods)
  auc = numeric(methods)
  f1_m = numeric(methods)
  f1_var = numeric(methods)
  avg_time = numeric(methods)
  var_time = numeric(methods)
  for (i in seq(methods))
  {
    avg_bic[i] = mean(method_bics[[i]]/n)
    avg_likeli[i] = mean(method_likeli[[i]]/n)
    MAD_bic[i] = mean(abs(method_bics[[i]]/n - true_bics/n))
    MSE_bic[i] = mean((method_bics[[i]]/n - true_bics/n)^2)
    exp_diff_bic[i] = mean((method_bics[[i]]/n - true_bics/n))
    var_diff_bic[i] = var((method_bics[[i]]/n - true_bics/n))
    MAD_test_likeli[i] =  mean(abs(method_test_likeli[[i]] - true_test_likeli))
    test_likeli[i] = mean(method_test_likeli[[i]])
    var_test_likeli[i] = var(method_test_likeli[[i]])
    MAD_likeli[i]  = mean(abs(method_likeli[[i]] - true_likeli))
    MSE_likeli[i]  = mean((method_likeli[[i]] - true_likeli)^2)
    p_likeli[i]  = mean(method_likeli[[i]] - true_likeli > log(eps))
    p_bic[i]  = mean(method_bics[[i]]/n - true_bics/n < eps)
    npar_diff[i] = mean(method_npar[[i]] - true_npar)
    f1_m[i] = mean(f1[[i]])
    f1_var[i] = var(f1[[i]])
    avg_time[i] = mean(times[[i]])
    var_time[i] = var(times[[i]])
  }
  return(list("avg_bic" = avg_bic,
              "avg_likeli" = avg_likeli,
              "MAD_bic" = MAD_bic, 
              "MSE_bic" = MSE_bic, 
              "MAD_test_likeli" = MAD_test_likeli,
              "test_likeli" = test_likeli,
              "MAD_likeli" = MAD_likeli,
              "MSE_likeli" = MSE_likeli,
              "p_eq_likeli" = p_likeli,
              "p_eq_bic" = p_bic,
              "npar_diff" = npar_diff,
              "exp_diff_bic" = exp_diff_bic,
              "avg_test_likeli" = test_likeli,
              "var_diff_bic" = var_diff_bic,
              "var_test_likeli" = var_test_likeli,
              "f1" = f1_m,
              "f1_v" = f1_var,
              "avg_time" = avg_time,
              "var_time" = var_time))
}

gather_metrics = function(metrics_list)
{
  settings = names(metrics_list)[order(as.numeric(names(metrics_list)))]
  avg_bic = list()
  avg_likeli = list()
  avg_test_likeli = list()
  MAD_bic = list()
  MSE_bic = list()
  MAD_test_likeli = list()
  MSE_test_likeli = list()
  MAD_likeli = list()
  MSE_likeli = list()
  p_eq_bic = list()
  p_eq_likeli = list()
  npar_diff = list()
  expec_diff = list()
  auc = list()
  var_expec_diff = list()
  var_test_likeli = list()
  f1 = list()
  f1_var = list()
  avg_time = list()
  var_time = list()
  for(j in seq(length(metrics_list[[1]]$MAD_bic)))
  {
    avg_bic[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[1]][j]))
    avg_likeli[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[2]][j]))
    MAD_bic[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[3]][j]))
    MSE_bic[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[4]][j]))
    MAD_test_likeli[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[5]][j]))
    MSE_test_likeli[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[6]][j]))
    MAD_likeli[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[7]][j]))
    MSE_likeli[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[8]][j]))
    p_eq_bic[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[9]][j]))
    p_eq_likeli[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[10]][j]))
    npar_diff[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[11]][j]))
    npar_diff[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[11]][j]))
    expec_diff[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[12]][j]))
    avg_test_likeli[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[13]][j]))
    var_expec_diff[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[14]][j]))
    var_test_likeli[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[15]][j]))
    f1[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[16]][j]))
    f1_var[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[17]][j]))
    avg_time[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[18]][j]))
    var_time[[j]] = sapply(settings, function(x) as.numeric(metrics_list[[x]][[19]][j]))
  }
  return(list("avg_bic" = avg_bic,
              "avg_likeli" = avg_likeli,
              "avg_test_likeli" = avg_test_likeli,
              "MAD_bic" = MAD_bic, 
              "MSE_bic" = MSE_bic, 
              "MAD_likeli" = MAD_likeli,
              "MSE_likeli" = MSE_likeli,
              "MAD_test_likeli" = MAD_test_likeli,
              "MSE_test_likeli" = MSE_test_likeli,
              "p_eq_likeli" = p_eq_likeli,
              "p_eq_bic" = p_eq_bic,
              "npar_diff" = npar_diff,
              "epec_bic_diff" = expec_diff,
              "var_bic_diff" = var_expec_diff,
              "var_test_likeli" = var_test_likeli,
              "f1" = f1,
              "fl_var" = f1_var,
              "avg_time" = avg_time,
              "var_time" = var_time))
}

getCausalEffects <- function(B) {
  diag(B) <- 0
  return(t(solve(diag(ncol(B))-t(B))))
}


calc_data_likeli = function(test_cov, Sigma, Beta, n, p)
{
  ll_Sigma = -n/2 * (log(det(Sigma)) + sum(diag(solve(Sigma) %*% test_cov))) - (n * p)/2 * log(2 * pi)
  return(ll_Sigma)
}


sim = function(num_tests, 
               n = 100000, 
               e_p = 0.4, 
               p = 10, 
               temperture = 0,
               cooldown = 0,
               max_iter = 10,
               random_graphs = 10,
               alternative_methods = TRUE,
               sparsity = .3,
               store_mags = F,
               hot_start = F,
               equiv = T,
               verbose = F)
{
  hybrid_bics = numeric(num_tests)
  true_bics = numeric(num_tests)
  greepybap_bics = numeric(num_tests)
  lrps_bics = numeric(num_tests)
  
  hybrid_likeli = numeric(num_tests)
  true_likeli = numeric(num_tests)
  greepybap_likeli = numeric(num_tests)
  lrps_likeli = numeric(num_tests)
  
  hybrid_likeli_test = numeric(num_tests)
  true_likeli_test = numeric(num_tests)
  greepybap_likeli_test = numeric(num_tests)
  lrps_likeli_test = numeric(num_tests)
  
  hybrid_npar = numeric(num_tests)
  true_npar = numeric(num_tests)
  greepybap_npar = numeric(num_tests)
  lrps_npar = numeric(num_tests)
  
  hybrid_equiv = list()
  greedy_equiv = list()
  lrps_equiv = list()
  
  hybrid_f1 = numeric(num_tests)
  greedy_f1 = numeric(num_tests)
  lrps_f1 = numeric(num_tests)  
  
  time_hybrid = numeric(num_tests)
  time_greedy = numeric(num_tests)
  time_lrps = numeric(num_tests)
  
  hybrid_mags = list()
  true_mags = list()
  props = list()
  
  for(i in seq(num_tests))
  {
    if(verbose) print(paste0("START SIM: ", i))
    res = GenerateGroundTruth(p, max.in.degree = ceiling(p * sparsity))
    dat = scale(GenerateData(n = n, res$params), center = T, scale = F)
    test_dat = scale(GenerateData(n = 10000, res$params), center = T, scale = F)
    suffStat = list("C" = cov(dat), "n" = nrow(dat))
    suffStat_test = list("C" = cov(test_dat), "n" = nrow(test_dat))
    ricf_true_mag = true_bic(res$mg, suffStat, res$params$B, res$covMat)
    indepTest = gaussCItest
    V = colnames(suffStat$C)
    #RUN HYBRID LEARN
    
    start = Sys.time()
    hybrid_learned = ADMG_hybrid_learning(dat, temperture = temperture, cooldown = cooldown, max_iter = max_iter, hot_start = hot_start, max_greedyadd = 700)
    end = Sys.time()
    
    time_hybrid[i] = end - start
    bic = hybrid_learned[[1]]
    equivclass = compare_equivclass(transform_input_ggm(attributes(hybrid_learned[[2]])$amat), res$mg, p, n, suffStat[[1]])
    hybrid_equiv[[i]]=equivclass
    hybrid_f1[i] = equivclass[[1]]
    if(verbose) print("HYBRID LEARN DONE")
    
    ricf_bap = NULL
    ricf_lrps = NULL
    if(alternative_methods)
    {
      if(verbose) print("RUNNING BAP")
      #RUN greepybap https://arxiv.org/abs/1508.01717
      
      start = Sys.time()
      bap_learn = greedyBAPs::greedySearch(suffStat[[1]], suffStat[[2]], n.restarts = 10)
      end = Sys.time()
      ricf_bap = ricf_ggm(bap_learn$final.bap, suffStat, max_iter = 250000)
      
      time_greedy[i] = end - start
      equivclass = compare_equivclass(bap_learn$final.bap, res$mg, p, n, suffStat[[1]])
      greedy_equiv[[i]] = equivclass
      greedy_f1[i] = equivclass[[1]]
      
      if(verbose) print("ENDING BAP")
      #RUN lrpsadmm + GES (https://academic.oup.com/jrsssb/article/81/3/459/7048382#396338676)
      if(verbose) print("LRPS RUNNING")
      start = Sys.time()
      gammas <- c(0.05, 0.07, 0.1, 0.12, 0.15, 0.17, 0.2)
      xval_path = lrpsadmm.cv(X = dat, gammas = gammas, covariance.estimator = cor, n.folds = 5, 
                             verbose = FALSE, n.lambdas = 40, lambda.ratio=1e-04, backend='RcppEigen')
      if(xval_path$best.fit$number.of.edges == 0)
      {
        lrps_adj = matrix(0, nrow = nrow(suffStat[[1]]), ncol = nrow(suffStat[[1]]))
      }
      else
      {
        selected_S = xval_path$best.fit$fit$S
        fake_data = generate.data.for.GES(Sest = selected_S, n=n, p=ncol(selected_S))
        lrps_ges_output = run.GES.and.select.with.BIC(obs.data = fake_data, nv = ncol(dat), sim.data = NA)
        lrps_adj = lrps_ges_output$best.essgraph
        lrps_adj[!upper.tri(lrps_adj)] = 0
      }
      end = Sys.time()
      
      time_lrps[i] = end - start
      dimnames(lrps_adj) = dimnames(suffStat[[1]])
      equivclass = compare_equivclass(lrps_adj, res$mg, p, n, suffStat[[1]])
      lrps_equiv[[i]] = equivclass
      lrps_f1[i] = equivclass[[1]]
      
      ricf_lrps = ricf_ggm(lrps_adj, suffStat, max_iter = 250000)
      if(verbose) print("LRPS ENDING")
    }
    #METRICS
    hybrid_bics[i] = bic
    hybrid_likeli[i] = hybrid_learned[[4]]
    hybrid_npar[i] = hybrid_learned[[5]]
    hybrid_likeli_test[i] = calc_data_likeli(suffStat_test[[1]], hybrid_learned[[6]], hybrid_learned[[7]], suffStat_test[[2]], p)
    if(store_mags)
    {
      hybrid_mags[[i]] = attributes(hybrid_learned[[2]])$amat
      true_mags[[i]] = adjmat_mag
    }
    true_bics[i] = ricf_true_mag[[3]]
    true_likeli[i] = ricf_true_mag[[4]]
    true_npar[i] = ricf_true_mag[[5]]
    true_likeli_test[i] = calc_data_likeli(suffStat_test[[1]], ricf_true_mag[[2]], ricf_true_mag[[1]], suffStat_test[[2]], p)
    if(alternative_methods)
    {
      if(length(ricf_bap[[4]]) == 0)
      {
        next
      }
      greepybap_bics[i] = ricf_bap[[3]]
      greepybap_likeli[i] = ricf_bap[[4]]
      greepybap_npar[i] = ricf_bap[[5]]
      greepybap_likeli_test[i] = calc_data_likeli(suffStat_test[[1]], ricf_bap[[2]], ricf_bap[[1]], suffStat_test[[2]], p)
      lrps_bics[i] = ricf_lrps[[3]]
      lrps_likeli[i] = ricf_lrps[[4]]
      lrps_npar[i] = ricf_lrps[[5]] 
      lrps_likeli_test[i] = calc_data_likeli(suffStat_test[[1]], ricf_lrps[[2]], ricf_lrps[[1]], suffStat_test[[2]], p)
    }
    props[[i]] = hybrid_learned[[3]]
    if(verbose) print(paste0("Simulation ", i, " DONE"))
    if(verbose) print("######SCORE#########")
    if(verbose) print(hybrid_bics[i])
    if(verbose) print(greepybap_bics[i])
    if(verbose) print(lrps_bics[i])
    if(verbose) print(ricf_true_mag[[3]])
    if(verbose) print("###################")
    i = i + 1
  }
  
  bics_methods = list("Hybrid Learn BIC" = hybrid_bics,
              "greepyBAP BIC" = greepybap_bics,
              "lrpsGES BIC" = lrps_bics)
  likelis_methods = list("Hybrid Learn likeli" = hybrid_likeli,
              "greepyBAP likeli" = greepybap_likeli,
              "lrpsGES likeli" = lrps_likeli)
  likelis_test_methods = list("Hybrid Learn likeli" = hybrid_likeli_test,
                         "greepyBAP likeli" = greepybap_likeli_test,
                         "lrpsGES likeli" = lrps_likeli_test)
  npars_methods = list("Hybrid Learn npar" = hybrid_npar,
               "greepyBAP npar" = greepybap_npar,
               "lrpsGES npar" = lrps_npar)
  hybrid_ADMG_results = list("Hybrid Learned ADMG" = hybrid_mags,
                             "Accepted Proposals" = props)
  f1 = list("Hybrid Learned f1" = hybrid_f1,
            "GreedyBAP f1" = greedy_f1,
            "lrpsGES F1" = lrps_f1)
  times = list("Hybrid Learned" = time_hybrid,
               "GreedyBAP" = time_greedy,
               "lrpsGES" = time_lrps)
  eq = list(hybrid_equiv, greedy_equiv, lrps_equiv)
  return(list("bics" = bics_methods,
              "likelis" = likelis_methods,
              "likelis_test" = likelis_test_methods,
              "npars" = npars_methods,
              "ADMG Hybrid Learn results" = hybrid_ADMG_results,
              "True ADMGs" = true_mags,
              "True ADMG BIC" = true_bics,
              "True ADMG test likeli" = true_likeli_test,
              "True ADMG likeli" = true_likeli,
              "True ADMG npar" = true_npar,
              "F1 Score" = f1,
              "times" = times,
              "eq" = eq))
}

sim_across_param = function(sim_range, 
                        num_tests, 
                        sim_over = "p",
                        n = 10000,
                        e_p = 0.4, 
                        p = 10, 
                        temperture = 0, 
                        cooldown = 0, 
                        random_graphs = 10,
                        sparsity =.3,
                        save = "",
                        max_iter = 10,
                        equiv = T, 
                        verbose = F)
{
  sim_results = list()
  args = list("num_tests" = num_tests,
              "n" = n, 
              "e_p" = e_p, 
              "p" = p, 
              "temperture" = temperture,
              "cooldown" = cooldown,
              "random_graphs" = random_graphs,
              "alternative_methods" = T,
              "sparsity" = sparsity,
              "hot_start" = T,
              "equiv"= equiv,
              "verbose" = verbose)
  for(s in sim_range)
  {
    args[[sim_over]] = s
    sim_s = do.call(sim, args)
    if(sim_over == "n")
    {
      n = s
    }
    true_bics = sim_s[["True ADMG BIC"]]
    true_likeli = sim_s[["True ADMG likeli"]]
    true_test_likeli = sim_s[["True ADMG test likeli"]]
    true_npar = sim_s[["True ADMG npar"]]
    bics = sim_s[["bics"]]
    likeli = sim_s[["likelis"]]
    likeli_test = sim_s[["likelis_test"]]
    f1 = sim_s[["F1 Score"]]
    npar = sim_s[["npars"]]
    times = sim_s[["times"]]
    save(sim_s, file = paste0("res_", save, s))
    sim_results[[paste0(s)]] = metrics(bics, true_bics,
                                               likeli, true_likeli,
                                               likeli_test, true_test_likeli,
                                               npar, true_npar, f1, n = n, times = times)
  }
  metrics = gather_metrics(sim_results)
  return(metrics)
}

res_p = sim_across_param(c(5), 2, save = "tol_", verbose = T)

res_p = sim_across_param(c(5,10, 15, 20, 30), 50, save = "tol_", verbose = T)
res_p_small_n = sim_across_param(c(5,10, 15, 20, 30), 50, n = 100, save = "smalln_")
res_n = sim_across_param(c(50,100, 500, 1000, 5000, 10000), 50, max_iter = 50, p = 10, sim_over = "n", save = "acrossn_")
