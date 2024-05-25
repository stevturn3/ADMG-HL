source("simulation.R")
source("simulation_plot.R")

setting_acrossp <- function()
{
  p_res = sim_across_param(c(5,10, 15, 20, 30), 50, save = "tol_", verbose = T)
  plot1 = gg_display(p_res$epec_bic_diff, lapply(p_res$var_bic_diff, function(x) sqrt(x)), c(5, 10, 15, 20, 30), sim_over = "p")
  plot2 = gg_display(p_res$f1, lapply(p_res$fl_var, function(x) sqrt(x)), c(5, 10, 15, 20), bar = T, metric_title = "F1 Score")
  plot3 = display_one_sim(c('./res_tol_5', './res_tol_10', './res_tol_15', './res_tol_20', './res_tol_30'), "times", c(5, 10, 15, 20, 30), n = NULL, title = "",  xlab = "metric", metric_title = "time(log seconds) ", no_diff = T, log_v = T)
  plot4 = display_one_sim(c('./res_tol_5', './res_tol_10', './res_tol_15', './res_tol_20', './res_tol_30'),"likelis_test", c(5, 10, 15, 20, 30), n = 10000, title = "",  xlab = "metric", metric_title = "Test Likelihood (normalized by 10,000)", no_diff = T, log_v = F)
  plots = list(plot1, plot2, plot3, plot4)
  return(list("results" = p_res,
              "plots" = plots))
}

setting_acrossn <- function()
{
#res_n = sim_across_param(c(50,100, 500, 1000, 5000, 10000), 50, max_iter = 50, p = 10, sim_over = "n", save = "acrossn_")
  plot1 = gg_display(res_n$epec_bic_diff, lapply(res_n$var_bic_diff, function(x) sqrt(x)), seq = log(c(50,100,500,1000,5000,10000)), sim_over = "n")
  plot2 = gg_display(res_n$f1, lapply(res_n$fl_var, function(x) sqrt(x)), c(50,100,500,1000,5000,10000), bar = T, metric_title = "F1 Score")
  plot3 = display_one_sim(c('./res_acrossn_50', './res_acrossn_100', './res_acrossn_500', './res_acrossn_1000', './res_acrossn_5000', 'res_acrossn_10000'), "times", c(50,100,500,1000,5000,10000), n = NULL, title = "",  xlab = "metric", metric_title = "time(log seconds) ", no_diff = T, log_v = T)
  plot4 = display_one_sim(c('./res_acrossn_50', './res_acrossn_100', './res_acrossn_500', './res_acrossn_1000', './res_acrossn_5000', 'res_acrossn_10000'), "likelis_test", c(50,100,500,1000,5000,10000), n = NULL, title = "",  xlab = "metric", metric_title = "Test Likelihood (normalized by 10,000)", no_diff = T, log_v = F)
  plots = list(plot1, plot2, plot3, plot4)
  return(list("results" = res_n,
              "plots" = plots))
}