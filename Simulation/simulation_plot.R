library(ggplot2)
library(reshape2) 
library(tidyverse)
gather_mag = function(x, type)
{
  if(typeof(x) == "list")
  {
    return(x[[type]])
  }
  else
  {
    return(x)
  }
}

load_res = function(dir, save_name = 'res_acrossn', sim_over = "n")
{
  fs = dir(dir)
  fs = fs[grepl(save_name, fs)]
  sim_results = list()
  sim_metrics = list()
  for(i in fs)
  {
    var_i = str_split(i, "_")[[1]][[3]]
    load(paste0(dir,i))
    sim_results[[var_i]] = sim_s
    if(sim_over == "n")
    {
      n = as.numeric(var_i)
    }
    else
    {
      n = 10000
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
    sim_metrics[[var_i]] = metrics(bics, true_bics,
                                   likeli, true_likeli,
                                   likeli_test, true_test_likeli,
                                   npar, true_npar, f1, n = n, times = times)
    rm(sim_s)
  }
  sim_metrics = gather_metrics(sim_metrics)
  return(list(sim_results, sim_metrics))
}

display_one_sim = function(fs, metric_name, sim_range, n = NULL, title = "",  xlab = "metric", metric_title = "", no_diff = T, log_v = T)
{
  plot_dat = data.frame("ADMG Hybrid Learn" = numeric(length(fs) * 50), "GreedyBAP Learning" = numeric(length(fs) * 50), "LRPS+GES Estimator" = numeric(length(fs) * 50), "True ADMG" = numeric(length(fs) * 50), "sim" = numeric(length(fs) * 50))
  for (i in seq(length(fs)))
  {
    print(fs[i])
  f = fs[i]
  load(f)
  metric_map = c("True ADMG BIC", "True ADMG test likeli","True ADMG likeli", "True ADMG npar")
  names(metric_map) = c("bics", "likelis_test", "likelis", "npars")
  metric_list = sim_s[[metric_name]]
  metric_true = sim_s[[metric_map[metric_name]]]
  d = NULL
  if(!no_diff)
  {
    if (!is.null(n))
    {
      a = metric_list[[1]]/n - metric_true/n
      b = metric_list[[3]]/n- metric_true/n
      c = metric_list[[4]]/n- metric_true/n
    }
    else
    {
      a = metric_list[[1]] - metric_true
      b = metric_list[[3]] - metric_true
      c = metric_list[[4]] - metric_true  
    }
    plot_dat[seq(1 + (i - 1) * 50,50 +(i - 1) * 50),1:5] = data.frame("ADMG Hybrid Learn" = a, "GreedyBAP Learning" = b, "LRPS+GES Estimator" = c, "True ADMG" = 0, "sim" = sim_range[i])
  }
  else
  {
    if (!is.null(n))
    {
      if(length(n) > 1)
      {
        print(n[i])
      a = metric_list[[1]]/n[i]
      b = metric_list[[3]]/n[i]
      c = metric_list[[4]]/n[i]
      }
      else
      {
        a = metric_list[[1]]/n
        b = metric_list[[3]]/n
        c = metric_list[[4]]/n
      }
      if(is.null(metric_true))
      {
        d = 0
      }
      else
      {
        if(length(n) > 1)
        {
          d = metric_true/n[i]
          
        }
        else
        {
          d = metric_true/n
          
        }
      }
    }
    else
    {
      a = metric_list[[1]] 
      b = metric_list[[3]]
      c = metric_list[[4]] 
      if(is.null(metric_true))
      {
        d = 0
      }
      else
      {
        d = metric_true
      }    
    }
    plot_dat[seq(1 + (i - 1) * 50,50 +(i - 1) * 50), 1:5]  = data.frame("ADMG Hybrid Learn" = a, "GreedyBAP Learning" = b, "LRPS+GES Estimator" = c, "True ADMG" = d, "sim" = sim_range[i])
  }
  }
  if(is.null(metric_true))
  {
    plot_dat = plot_dat %>% select(-`True.ADMG`)
  }
  plot_dat = reshape2::melt(plot_dat, id = "sim")
  if(log_v)
  {
    plot_dat["value"] = log(plot_dat["value"])
  }
  print(metric_title)
  p = ggplot(plot_dat, aes(x = variable, y = value, fill = variable)) + geom_violin() + 
    labs(x = "method", y = metric_title, title = title, fill = "method")  + 
    theme_bw(base_size = 12, base_family = "Helvetica") + 
    facet_wrap(~ sim, ncol=2,scales = "free") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.x = element_text(size=18, face="bold"),
          axis.title.y = element_text(size=18),
          axis.text.y = element_text(size=14),
          legend.title=element_text(size=18),
          legend.text=element_text(size=14), legend.position = c(1, 0),
  legend.justification = c(1, 0))
  return(p)
}

gg_display = function(metric, metric_sd, seq, sim_over = "log(n)",
                      metric_title = "Expected BIC Difference", title = "",
                      methods = c("Hybrid ADMG Learning", "GreedyBAP Learning", "LRPS + GES Estimator"),bar = F)
{
  metric = data.frame(metric)
  metric_sd = data.frame(metric_sd)
  l = seq(1, length(metric))
  names(metric) = methods
  names(metric_sd) = methods
  metric['seq'] = as.numeric(rownames(metric))
  metric_sd['seq'] = as.numeric(rownames(metric_sd))
  metric = reshape2::melt(metric, id = c("seq"))
  metric_sd = reshape2::melt(metric_sd, id = c("seq"))
  metric = metric %>% inner_join(metric_sd, by = c("seq", "variable"))
  if(bar)
  {
    gfg_plot <- ggplot(metric,             
                       aes(x = variable,
                           y = value.x,
                           ymin=ifelse(value.x-value.y > 0, value.x-value.y, 0) , 
                           ymax=value.x+value.y,
                           fill = variable,
                           label = round(value.x, 2))) + 
      geom_bar(stat = "identity") + 
      facet_wrap(~ seq, ncol=3, scales = "free") +
      geom_errorbar() + 
      labs(x = sim_over, y = metric_title, title = title, fill = "method")  + 
      theme_bw(base_size = 12, base_family = "Helvetica") + 
      theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    strip.text.x = element_text(size=18, face="bold"),
                    axis.text.y = element_text(size=14),
                    legend.title=element_text(size=18),
                    legend.text=element_text(size=14),
                    legend.position = c(1, 0),
                    legend.justification = c(1, 0))
  }
  else
  {
  gfg_plot <- ggplot(metric,             
                     aes(x = seq, 
                         y = value.x, 
                         color = variable)) + geom_line() +   
              geom_errorbar(aes(ymin=value.x-value.y, ymax=value.x+value.y), width=.2) + 
              scale_x_continuous(breaks=c(5, 10, 15, 20, 30)) + 
              labs(x = sim_over, y = metric_title, title = title, color = "method")  + 
              theme_bw(base_size = 12, base_family = "Helvetica") + 
              theme(axis.text=element_text(size=12),
                axis.title=element_text(size=18), 
                legend.title=element_text(size=14),
                legend.text=element_text(size=12))
  }
  print(gfg_plot)
}




acrossn_res = load_res("./")
p_res = load_res("./", save_name = "res_tol", sim_over = "p")
smalln_res = load_res("./", save_name = "res_smalln", sim_over = "p")
smalln_res = load_res("./", save_name = "res_acrossn", sim_over = "n")

gg_display(acrossn_res[[2]]$epec_bic_diff[c(1,3,4)], lapply(acrossn_res[[2]]$var_bic_diff[c(1,3,4)], function(x) sqrt(x)), c(50,100,500,1000,5000,10000))
gg_display(acrossn_res$f1[c(1,3,4)], lapply(acrossn_res$fl_var[c(1,3,4)], function(x) sqrt(x)), c(50,100,500,1000,5000,10000), bar = T, metric_title = "F1 Score")
gg_display(p_res[[2]]$epec_bic_diff[c(1,3,4)], lapply(p_res[[2]]$var_bic_diff[c(1,3,4)], function(x) sqrt(x)), c(5, 10, 15, 20, 30), sim_over = "p")
gg_display(p_res[[2]]$epec_bic_diff[c(1,3)], lapply(p_res[[2]]$var_bic_diff[c(1,3)], function(x) sqrt(x)), c(5, 10, 15, 20, 30), sim_over = "p", methods = c("Hybrid ADMG Learning", "GreedyBAP Learning"))

gg_display(p_res[[2]]$f1[c(1,3,4)], lapply(p_res[[2]]$fl_var[c(1,3,4)], function(x) sqrt(x)), c(5, 10, 15, 20), bar = T, metric_title = "F1 Score")
gg_display(smalln_res[[2]]$epec_bic_diff[c(1,3,4)], lapply(smalln_res[[2]]$var_bic_diff[c(1,3,4)], function(x) sqrt(x)), c(5, 10, 15, 20), sim_over = "p")
gg_display(smalln_res[[2]]$f1[c(1,3,4)], lapply(smalln_res[[2]]$fl_var[c(1,3,4)], function(x) sqrt(x)), c(5, 10, 15, 20), bar = T, metric_title = "F1 Score")
gg_display(p_res[[2]]$avg_test_likeli[c(1,3,4)], lapply(p_res[[2]]$var_test_likeli[c(1,3,4)], function(x) sqrt(x)), c(5, 10, 15, 20), bar = F, metric_title = "Test-Likelihood", sim_over = "p")



display_one_sim(c('./res_tol_5', './res_tol_10', './res_tol_15', './res_tol_20', './res_tol_30'), "times", c(5, 10, 15, 20, 30), n = NULL, title = "",  xlab = "metric", metric_title = "time(log seconds) ", no_diff = T, log_v = T)
display_one_sim(c('./res_tol_5', './res_tol_10', './res_tol_15', './res_tol_20', './res_tol_30'), "likelis_test", sim_range=c(5, 10, 15, 20, 30),n = 10000, title = "",  xlab = "metric", metric_title = "Test Log Likelihood (normalized by n)", no_diff = T, log_v = F)

