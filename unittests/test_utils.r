#rm(list = ls())
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#setwd("..")
source("../util.R", chdir = TRUE)

test_get_bidirected = function()
{
    adjmat = matrix(c(0,2,0,0,2,0,2,2,0,2,0,2,0,2,2,0), ncol = 4)
    colnames(adjmat) = c("A", "B", "C", "D")
    rownames(adjmat) = c("A", "B", "C", "D")
    res = adjmat
    res[adjmat == 2] = 1
    stopifnot(all(get_bidirected(adjmat)[[2]] == res))
    print("PASSED")
}

test_get_directed = function()
{
    adjmat = matrix(c(0,2,0,0,3,0,2,2,0,3,0,2,0,2,3,0), ncol = 4)
    colnames(adjmat) = c("A", "B", "C", "D")
    rownames(adjmat) = c("A", "B", "C", "D") 
    res = matrix(c(0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0), ncol = 4)
    colnames(res) = c("A", "B", "C", "D")
    rownames(res) = c("A", "B", "C", "D")
    stopifnot(all(get_directed(adjmat)[[2]] == res))   
    print("PASSED") 
}

test_pid_check_colliders = function()
{
    adjmat = matrix(c(0,2,0,0,3,0,2,2,0,3,0,2,0,2,3,0), ncol = 4)
    colnames(adjmat) = c("A", "B", "C", "D")
    rownames(adjmat) = c("A", "B", "C", "D")
    pid_collider = pid_check_colliders(adjmat)
    stopifnot(pid_collider == c("D,C"))
    print("PASSED")
}

test_transform_input_ggm = function()
{
    res = matrix(c(0,0,0,0,1,0,0,100,0,1,0,0,0,100,1,0), byrow = TRUE, ncol = 4)
    adjmat = matrix(c(0,2,0,0,3,0,2,2,0,3,0,2,0,2,3,0), ncol = 4)
    colnames(adjmat) = c("A", "B", "C", "D")
    rownames(adjmat) = c("A", "B", "C", "D")
    stopifnot(transform_input_ggm(adjmat) == res)
    print("PASSED")
}

test_add_bidirected = function()
{
    adjmat = matrix(c(0,2,0,0,3,0,2,2,0,3,0,2,0,2,3,0), ncol = 4)
    adjmat_res = matrix(c(0,2,0,2,3,0,2,2,0,3,0,2,2,2,3,0), ncol = 4)
    colnames(adjmat) = c("A", "B", "C", "D")
    rownames(adjmat) = c("A", "B", "C", "D")
    adjmat_new = add_bidirected(adjmat, "A", "D")
    stopifnot(adjmat_new == adjmat_res)
    print("PASSED")
}

test_add_directed = function()
{
    adjmat = matrix(c(0,2,0,0,3,0,2,2,0,3,0,2,0,2,3,0), ncol = 4)
    adjmat_res = matrix(c(0,2,0,3,3,0,2,2,0,3,0,2,2,2,3,0), ncol = 4)
    colnames(adjmat) = c("A", "B", "C", "D")
    rownames(adjmat) = c("A", "B", "C", "D")
    adjmat_new = add_directed(adjmat, "A", "D")
    stopifnot(adjmat_new == adjmat_res)
    print("PASSED") 
}

test_remove_edge = function()
{
    adjmat = matrix(c(0,2,0,3,3,0,2,2,0,3,0,2,2,2,3,0), ncol = 4)
    adjmat_res = matrix(c(0,2,0,0,3,0,2,2,0,3,0,2,0,2,3,0), ncol = 4)
    colnames(adjmat) = c("A", "B", "C", "D")
    rownames(adjmat) = c("A", "B", "C", "D")
    res = remove_edge(adjmat, "A", "D")
    adjmat_new = res[[1]]
    type = res[[2]]
    stopifnot(adjmat_new == adjmat_res)
    stopifnot(type == c("dir", "none"))
    adjmat_res = matrix(c(0,2,0,3,3,0,2,0,0,3,0,2,2,0,3,0), ncol = 4)
    colnames(adjmat) = c("A", "B", "C", "D")
    rownames(adjmat) = c("A", "B", "C", "D")
    res = remove_edge(adjmat, "B", "D")
    adjmat_new = res[[1]]
    type = res[[2]]
    stopifnot(adjmat_new == adjmat_res)
    stopifnot(type == c("bidir", "none"))
    print("PASSED") 
}

test_reverse_edge = function()
{
    edge = "B,A"
    rev = reverse_edge(edge)
    stopifnot(rev == "A,B")
    print("PASSED")
}



test_get_bidirected()
test_get_directed()
test_pid_check_colliders()
test_transform_input_ggm()
test_add_bidirected()
test_add_directed()
test_remove_edge()
test_reverse_edge()

