#rm(list = ls())
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#setwd("..")
source("../proposals.r", chdir = TRUE)

##TO DO:
##IMPLEMENT UNIT TESTS

test_phase_type1 = function()
{
    set.seed(0)
    data(gmD)
    suffStat = list("C" = cor(gmD[[1]]), "n" = nrow(gmD[[1]]))
    indepTest = gaussCItest
    alpha = 0.05
    labels = names(gmD[[1]])
    fci.pag = run_fci(suffStat, indepTest, alpha, labels)
    fci.mag = fci.pag
    adjmat = pag2magAM(attributes(fci.pag)$amat, 1)
    adjmat["X1", "X2"] = 2
    adjmat["X2", "X1"] = 3
    adjmat["X3", "X2"] = 2
    adjmat["X2", "X3"] = 3
    attributes(fci.mag)$amat = adjmat
    pbi = build_pbi(attributes(fci.mag)$amat)
    res = phase(pbi, fci.mag, adjmat, suffStat, 1)
    adjmat_proposal = adjmat
    adjmat_proposal["X1", "X2"] = 2
    adjmat_proposal["X2", "X1"] = 2
    stopifnot(round(res[[1]][[3]],1) == 140192.2)
    stopifnot(attributes(res[[2]])$amat == adjmat_proposal)
    stopifnot(res[[3]] == "X1,X2")
    stopifnot(res[[4]] == c("dir", "bidir"))
    print("PASSED")
}

test_phase_type2 = function()
{
    set.seed(0)
    data(gmD)
    suffStat = list("C" = cor(gmD[[1]]), "n" = nrow(gmD[[1]]))
    indepTest = gaussCItest
    alpha = 0.05
    labels = names(gmD[[1]])
    fci.pag = run_fci(suffStat, indepTest, alpha, labels)
    fci.mag = fci.pag
    adjmat = pag2magAM(attributes(fci.pag)$amat, 1)
    adjmat["X3", "X4"] = 2
    adjmat["X4", "X3"] = 3
    adjmat_proposal = adjmat
    adjmat_proposal["X3", "X4"] = 0
    adjmat_proposal["X4", "X3"] = 0
    attributes(fci.mag)$amat = adjmat
    pid = build_pid(adjmat)
    res = phase(pid, fci.mag, adjmat, suffStat, 2)
    stopifnot(attributes(res[[2]])$amat == adjmat_proposal)
    stopifnot(round(res[[1]][[3]],1) == 140175.2)
    stopifnot(res[[3]] == "X3,X4")
    stopifnot(res[[4]] == c("dir", "none"))
    adjmat["X3", "X4"] = 2
    adjmat["X4", "X3"] = 2
    attributes(fci.mag)$amat = adjmat
    pid = build_pid(adjmat)
    res = phase(pid, fci.mag, adjmat, suffStat, 2)
    stopifnot(attributes(res[[2]])$amat == adjmat_proposal)
    stopifnot(round(res[[1]][[3]],1) == 140175.2)
    stopifnot(res[[3]] == "X3,X4" | res[[3]] == "X4,X3")
    stopifnot(res[[4]] == c("bidir", "none"))
    print("PASSED")
}

test_phase_type3 = function()
{
    set.seed(0)
    data(gmD)
    suffStat = list("C" = cor(gmD[[1]]), "n" = nrow(gmD[[1]]))
    indepTest = gaussCItest
    alpha = 0.05
    labels = names(gmD[[1]])
    fci.pag = run_fci(suffStat, indepTest, alpha, labels)
    fci.mag = fci.pag
    adjmat = pag2magAM(attributes(fci.pag)$amat, 1)
    adjmat["X1", "X2"] = 2
    adjmat["X2", "X1"] = 3
    adjmat["X3", "X2"] = 2
    adjmat["X2", "X3"] = 3
    attributes(fci.mag)$amat = adjmat
    pbi = build_pbi(attributes(fci.mag)$amat)
    res = phase(pbi, fci.mag, adjmat, suffStat, 1)
    pbi_reversal = pbi
    res2 = phase(pbi_reversal, res[[2]], attributes(res[[2]])$amat, suffStat, 3)
    stopifnot(round(res2[[1]][[3]],0) == 140183)
    stopifnot(attributes(res2[[2]])$amat == adjmat)
    stopifnot(res2[[3]] == "X1,X2")
    stopifnot(res2[[4]] == c("bidir", "dir"))
    print("PASSED")
}

test_phase_type4 = function()
{
    set.seed(0)
    data(gmD)
    suffStat = list("C" = cor(gmD[[1]]), "n" = nrow(gmD[[1]]))
    indepTest = gaussCItest
    alpha = 0.05
    labels = names(gmD[[1]])
    fci.pag = run_fci(suffStat, indepTest, alpha, labels)
    fci.mag = fci.pag
    adjmat = pag2magAM(attributes(fci.pag)$amat, 1)
    adjmat["X3", "X4"] = 2
    adjmat["X4", "X3"] = 3
    attributes(fci.mag)$amat = adjmat
    E_d0 = get_directed(adjmat)[[1]]$edge
    E_b0 = get_bidirected(adjmat)[[1]]$edge
    pid = build_pid(adjmat)
    res = phase(pid, fci.mag, adjmat, suffStat, 2)

    pid_reversal = res[[3]]
    res2 = phase(pid_reversal, res[[2]], attributes(res[[2]])$amat, suffStat, 4, E_d0, E_b0)
    stopifnot(round(res2[[1]][[3]],1) == 140184.4)
    stopifnot(attributes(res2[[2]])$amat == adjmat)
    stopifnot(res2[[3]] == "X3,X4")
    stopifnot(res2[[4]] == c("none", "dir"))
    print("PASSED")
}

test_update_mag_positive = function()
{
    set.seed(0)
    data(gmD)
    suffStat = list("C" = cor(gmD[[1]]), "n" = nrow(gmD[[1]]))
    indepTest = gaussCItest
    alpha = 0.05
    labels = names(gmD[[1]])
    fci.pag = run_fci(suffStat, indepTest, alpha, labels)
    fci.mag = fci.pag
    adjmat = pag2magAM(attributes(fci.pag)$amat, 1)
    adjmat["X1", "X2"] = 2
    adjmat["X2", "X1"] = 3
    adjmat["X3", "X2"] = 2
    adjmat["X2", "X3"] = 3
    attributes(fci.mag)$amat = adjmat
    E_d0 = get_directed(adjmat)[[1]]$edge
    E_b0 = get_bidirected(adjmat)[[1]]$edge
    pbi = build_pbi(adjmat)
    ricf_orig = ricf_ggm(adjmat, suffStat)
    res = phase(pbi, fci.mag, adjmat, suffStat, 1)
    updates = update_mag(res, E_d0, E_b0, Inf)
    stopifnot(updates[[1]])
    stopifnot(updates[[2]] == c("X3,X2", "X1,X3", "X5,X4"))
    stopifnot(updates[[3]] == c("X4,X2", "X2,X4", "X1,X2", "X2,X1"))
    stopifnot(any(updates[[4]] != ricf_orig[[1]]))
    stopifnot(any(updates[[5]] != ricf_orig[[2]]))
    stopifnot(any(updates[[6]] != ricf_orig[[3]]))
    pbi_reversal = pbi
    res2 = phase(pbi_reversal, res[[2]], attributes(res[[2]])$amat, suffStat, 3)
    updates = update_mag(res2, updates$edge_d, updates$edge_b, Inf)
    stopifnot(updates[[1]])
    stopifnot(sort(updates[[2]]) == sort(E_d0))
    stopifnot(sort(updates[[3]]) == sort(E_b0))
    stopifnot(all(updates[[4]] == ricf_orig[[1]]))
    stopifnot(all(updates[[5]] == ricf_orig[[2]]))
    stopifnot(all(updates[[6]] == ricf_orig[[3]]))

    adjmat = pag2magAM(attributes(fci.pag)$amat, 1)
    adjmat["X3", "X4"] = 2
    adjmat["X4", "X3"] = 3
    attributes(fci.mag)$amat = adjmat
    E_d0 = get_directed(adjmat)[[1]]$edge
    E_b0 = get_bidirected(adjmat)[[1]]$edge
    pid = build_pid(adjmat)
    ricf_orig = ricf_ggm(adjmat, suffStat)
    res = phase(pid, fci.mag, adjmat, suffStat, 2)
    updates = update_mag(res, E_d0, E_b0, Inf)
    stopifnot(updates[[1]])
    stopifnot(updates[[2]] == c("X1,X3", "X5,X4"))
    stopifnot(updates[[3]] == c("X3,X2", "X4,X2", "X2,X3", "X2,X4"))
    stopifnot(any(updates[[4]] != ricf_orig[[1]]))
    stopifnot(any(updates[[5]] != ricf_orig[[2]]))
    stopifnot(any(updates[[6]] != ricf_orig[[3]]))
    pid_reversal = "X3,X4"
    res2 = phase(pid_reversal, res[[2]], attributes(res[[2]])$amat, suffStat, 4, E_d0, E_b0)
    updates = update_mag(res2, updates$edge_d, updates$edge_b, Inf)
    stopifnot(updates[[1]])
    stopifnot(sort(updates[[2]]) == sort(E_d0))
    stopifnot(sort(updates[[3]]) == sort(E_b0))
    stopifnot(all(updates[[4]] == ricf_orig[[1]]))
    stopifnot(all(updates[[5]] == ricf_orig[[2]]))
    stopifnot(all(updates[[6]] == ricf_orig[[3]]))
    print("PASSED")
}

test_update_mag_negative = function()
{
    set.seed(0)
    data(gmD)
    suffStat = list("C" = cor(gmD[[1]]), "n" = nrow(gmD[[1]]))
    indepTest = gaussCItest
    alpha = 0.05
    labels = names(gmD[[1]])
    fci.pag = run_fci(suffStat, indepTest, alpha, labels)
    fci.mag = fci.pag
    adjmat = pag2magAM(attributes(fci.pag)$amat, 1)
    adjmat["X3", "X4"] = 2
    adjmat["X4", "X3"] = 3
    attributes(fci.mag)$amat = adjmat
    E_d0 = get_directed(adjmat)[[1]]$edge
    E_b0 = get_bidirected(adjmat)[[1]]$edge
    pid = build_pid(adjmat)
    res = phase(pid, fci.mag, adjmat, suffStat, 2)
    updates = update_mag(res,E_d0, E_b0, -Inf)
    stopifnot(!updates[[1]])
    print("PASSED")
}

test_mag_search_pid = function()
{
    set.seed(0)
    data(gmD)
    suffStat = list("C" = cor(gmD[[1]]), "n" = nrow(gmD[[1]]))
    indepTest = gaussCItest
    alpha = 0.05
    labels = names(gmD[[1]])
    fci.pag = run_fci(suffStat, indepTest, alpha, labels)
    fci.mag = fci.pag
    adjmat = pag2magAM(attributes(fci.pag)$amat, 1)
    adjmat["X3", "X4"] = 2
    adjmat["X4", "X3"] = 3
    attributes(fci.mag)$amat = adjmat

    adjmat_proposal = adjmat 
    adjmat_proposal["X3", "X4"] = 0
    adjmat_proposal["X4", "X3"] = 0
    res_init = ricf_ggm(adjmat,suffStat)
    res = mag_search(fci.mag, labels, suffStat, test = TRUE)
    #Correct Calcualtion
    stopifnot(round(res[[2]], 1) == 140175.2)
    #Change in BIC
    stopifnot(round(res[[2]], 1) != res_init[[3]])
    stopifnot(attributes(res[[1]])$amat == adjmat_proposal)
    print("PASSED")
}

test_mag_search_pbi = function()
{
    set.seed(0)
    n = 10000
    X = rnorm(n)
    U = rnorm(n)
    Y = X + U + rnorm(n)
    S = Y + rnorm(n)
    Z = S + U + rnorm(n)
    W = Z + rnorm(n)
    data = cbind(X, Y, Z, W, S)
    colnames(data) = c("X", "Y","Z","W", "S")
    adjmat_given = matrix(c(0, 2, 0, 0, 0, 3, 0, 2, 0, 2,  0, 3, 0, 2, 3, 0, 0, 3, 0, 0, 0, 3, 2, 0, 0), byrow = T, ncol =5)
    colnames(adjmat_given) = c("X", "Y","Z","W", "S")
    rownames(adjmat_given) = c("X", "Y","Z","W", "S")
    
    suffStat = list("C" = cor(data), "n" = n)
    indepTest = gaussCItest
    alpha = 0.05
    labels = colnames(suffStat$C)
    fci.pag = run_fci(suffStat, indepTest, alpha, labels)
    fci.mag = fci.pag
    #Test Case 
    attributes(fci.mag)$amat = adjmat_given
    suffStat = list("C" = cov(data), "n" = n)
    
    res_init = ricf_ggm(adjmat_given,suffStat)
    res = mag_search(fci.mag, labels, suffStat, test = TRUE)
    #Correct Calcualtion
    stopifnot(round(res[[2]], 1) == 152719.2)
    #Change in BIC
    stopifnot(round(res[[2]], 1) < res_init[[3]])
    print("PASSED")
}




test_phase_type1()
test_phase_type2()
test_phase_type3()
test_phase_type4()

test_update_mag_positive()
test_update_mag_negative()

test_mag_search_pid()
test_mag_search_pbi()
