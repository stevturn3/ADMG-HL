#rm(list = ls())
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#setwd("..")
source("../algorithms.r", chdir = TRUE)

test_run_fci = function()
{
    #Function implemented by pcalg
    print("PASSED")
}

test_build_pbi_positive = function()
{
    adjmat = matrix(c(0,3,0,3,2,0,3,2,0,2,0,3,2,2,2,0), ncol = 4, byrow = TRUE)
    colnames(adjmat) = c("A", "B", "C", "D")
    rownames(adjmat) = c("A", "B", "C", "D")
    pbi = build_pbi(adjmat)
    stopifnot(pbi == c("D,A"))
    print("PASSED")
}

test_build_pbi_negative = function()
{
    adjmat = matrix(c(0,3,0,0,2,0,3,2,0,2,0,3,0,2,2,0), ncol = 4, byrow = TRUE)
    colnames(adjmat) = c("A", "B", "C", "D")
    rownames(adjmat) = c("A", "B", "C", "D")
    pbi = build_pbi(adjmat)
    stopifnot(pbi == character())
    print("PASSED")
}

test_build_pid_positive_longbidir = function()
{
    adjmat = matrix(c(0, 2, 0,0,2, 2, 0, 2,0,0, 0,2,0,2,0, 0,0,2,0,3, 3,0,0, 2, 0), byrow = TRUE, ncol = 5)
    colnames(adjmat) = 1:5
    rownames(adjmat) = 1:5
    pid = build_pid(adjmat)
    stopifnot(pid == c("1,5", "5,1"))
    print("PASSED")
}

test_build_pid_positive_length1bidir = function()
{
    adjmat = matrix(c(0, 2,0,2, 3,0,2,2,0,3,0,2,3,2,3,0), byrow = TRUE, ncol = 4)
    colnames(adjmat) = 1:4
    rownames(adjmat) = 1:4
    pid = build_pid(adjmat)
    stopifnot(pid == c("1,4", "1,2", "2,3"))
    print("PASSED")
}


test_build_pid_negative = function()
{

}

test_ricf_ggm  = function()
{
    #Mostly to test BIC since this function is driven by fitAncestralGraph
    data(gmD)
    adjmat = as(gmD$g, "matrix")   
    x = gmD$x
    n = nrow(x)
    p = ncol(x)
    adjmat[adjmat == 1] = 2
    adjmat[t(adjmat) == 2] = 3
    colnames(adjmat) = colnames(x)
    rownames(adjmat) = colnames(x)
    adjmat_ggm = transform_input_ggm(adjmat)
    py = length(which(colSums(adjmat_ggm) != 0))
    n_par = sum(adjmat_ggm == 1) + sum(adjmat_ggm == 100)/2 + py
    suffStat = list("C" = cov(x), "n" = nrow(x))
    res = fitAncestralGraph(adjmat_ggm, suffStat$C, suffStat$n)
    ll_s= -n/2 * (log(det(suffStat$C)) + sum(diag(solve(suffStat$C) %*% suffStat$C))) - (n * p)/2 * log(2 * pi)
    ll_sigma = - res$dev / 2 + ll_s
    BIC_1 = -2 * ll_sigma + n_par * log(suffStat$n)
    res2 = ricf_ggm(adjmat, suffStat)
    stopifnot(all.equal(res2[[3]], BIC_1))
    print("PASSED")
}

test_build_pbi_negative()
test_build_pbi_positive()
test_ricf_ggm()
test_run_fci()
test_build_pid_positive_longbidir()
