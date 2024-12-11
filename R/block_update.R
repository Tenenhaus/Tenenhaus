#' @importFrom expm sqrtm 
#TODO remove last import if I manage to use fct sqrt_matrix from RGCCA\tests\test_rgcca.R

block_update <- function(x, grad) {
  UseMethod("block_update")
}

#' @export
block_update.block <- function(x, grad) {
  x$a <- pm(t(x$x), grad, na.rm = x$na.rm)
  return(block_project(x))
}

#' @export
block_update.dual_block <- function(x, grad) {
  x$alpha <- grad
  return(block_project(x))
}

#' @export
block_update.ac_block <- function(x, grad) {
  mu <- x$penalty_coef * 
    eigen(x = x$B, symmetric = T, only.values = T)$values[1] *
    eigen(x = x$M, symmetric = T, only.values = T)$values[1]
  M <- x$tau * diag(x$p) + (1 - x$tau) * pm(t(x$x), x$x, na.rm = x$na.rm) / x$N
  x$f <- 1/(2*mu) * pm(x$M, 
                       (pm(sqrtm(x$M), 
                           pm(t(x$x), grad, na.rm = x$na.rm), na.rm = x$na.rm) 
                        - 2 * x$penalty_coef * pm(x$B, 
                                                  pm(sqrtm(M), x$a, na.rm = x$na.rm), na.rm = x$na.rm)
                        ), na.rm = x$na.rm)
  return(block_project(x))
}