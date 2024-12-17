#' @importFrom MASS ginv

block_init <- function(x, init = "svd") {
  print("block_init")
  print(class(x))
  UseMethod("block_init")
}

#' @export
block_init.block <- function(x, init = "svd") {
  if (init == "svd") {
    x$a <- initsvd(x$x, dual = FALSE)
  } else {
    x$a <- rnorm(x$p)
  }
  print("block_init.block is being used")
  return(block_project(x))
}

#' @export
block_init.dual_block <- function(x, init = "svd") {
  if (init == "svd") {
    x$alpha <- initsvd(x$x, dual = TRUE)
  } else {
    x$alpha <- rnorm(x$n)
  }
  print("block_init.dual_block is being used")
  return(block_project(x))
}

#' @export
block_init.primal_regularized_block <- function(x, init = "svd") {
  x$M <- ginv(
    x$tau * diag(x$p) + (1 - x$tau) * pm(t(x$x), x$x, na.rm = x$na.rm) / x$N
  )
  print("block_init.primal_regularized_block is being used")
  NextMethod()
}

#' @export
block_init.dual_regularized_block <- function(x, init = "svd") {
  x$M <- ginv(x$tau * diag(x$n) + (1 - x$tau) * x$K / x$N)
  print("block_init.dual_regularized_block is being used")
  NextMethod()
}

#' @export
block_init.ac_block <- function(x, init = "svd") {
  print("block_init.ac_block is being used")
  x$M <- x$tau * diag(x$p) + (1 - x$tau) * pm(t(x$x), x$x, na.rm = x$na.rm) / x$N
  x$M_inv <- ginv(x$M)
  x$sqrt_M <- sqrt_matrix(x$M)
  x$sqrt_M_inv <- sqrt_matrix(x$M_inv) #TODO check if result is the same with sqrt_matrix(x$M, inv = T)

  XM_tmp <- pm(x$x, x$sqrt_M_inv, na.rm = x$na.rm)
  x$B <- pm(
    t(XM_tmp), pm(
      x$confounders, XM_tmp, na.rm = x$na.rm), 
    na.rm = x$na.rm) #TODO should I use crossprod instead of pm?
  
  x$mu <- x$penalty_coef * 
<<<<<<< HEAD
    eigen(x = x$B, symmetric = T, only.values = T)$values[1]
=======
    eigen(x = x$B, symmetric = T, only.values = T)$values[1] *
    eigen(x = x$M_inv, symmetric = T, only.values = T)$values[1]
>>>>>>> 9c9aadad608e602f237b0c0a94005d9896d929bb
  
  #NextMethod()
  
  if (init == "svd") {
    x$a <- initsvd(x$x, dual = FALSE)
  } else {
    x$a <- rnorm(x$p)
  }
  return(block_project(x))
}
