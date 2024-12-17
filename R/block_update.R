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
<<<<<<< HEAD
  x$f <- 1/(2 * x$mu) * (pm(x$sqrt_M_inv,
                            pm(t(x$x),
                               grad,
                               na.rm = x$na.rm),
                            na.rm = x$na.rm)
                         - 2 * x$penalty_coef * pm(x$B,
                                                   pm(x$sqrt_M,
                                                      x$a,
                                                      na.rm = x$na.rm),
                                                   na.rm = x$na.rm))
=======
  x$f <- 1/(2 * x$mu) * pm(x$M_inv,
                           (pm(x$sqrt_M_inv,
                               pm(t(x$x),
                                  grad,
                                  na.rm = x$na.rm),
                               na.rm = x$na.rm)
                            - 2 * x$penalty_coef * pm(x$B,
                                                      pm(x$sqrt_M,
                                                         x$a,
                                                         na.rm = x$na.rm),
                                                      na.rm = x$na.rm)),
                           na.rm = x$na.rm)
>>>>>>> 9c9aadad608e602f237b0c0a94005d9896d929bb
  return(block_project(x))
}