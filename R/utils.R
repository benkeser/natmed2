#' Helper function to truncate values away from 0 and 1
#' 
#' @param x Vector of values to truncate
#' @param tol Truncation tolerance
#' @param bound_away_from A vector of 0 and/or 1.
#' @return Vector with truncated values replaced by \code{tol}.
g_truncate <- function(g, tol, bound_away_from = c(0, 1)){
  if(0 %in% bound_away_from){
    g[g < tol] <- tol
  }
  if(1 %in% bound_away_from){
    g[g > (1 - tol)] <- 1 - tol
  }
  return(g)
}

#' Helper function to format risk and effect estimates
#' 
#' @param psi11n Estimate of E[Y(1, S(1))]
#' @param psi00n Estimate of E[Y(0, S(0))]
#' @param psi10n Estimate of E[Y(1, S(0))]
#' @param psi01n Estimate of E[Y(0, S(1))]
#' @param cov_matrix Covariance matrix
#' @param cov_matrix_cv Cross-validated covariance matrix
#' @return A \code{list} with named objects \code{risk}, \code{eff}, \code{eff2}

get_risks_and_effects <- function(psi11n, psi00n, psi10n, psi01n,
                                  cov_matrix, cov_matrix_cv){
  # risks
  one_step_est <- c(psi11n, psi00n, psi10n, psi01n)
  cil <- one_step_est - 1.96 * sqrt(diag(cov_matrix))
  ciu <- one_step_est + 1.96 * sqrt(diag(cov_matrix))
  cil_cv <- one_step_est - 1.96 * sqrt(diag(cov_matrix_cv))
  ciu_cv <- one_step_est + 1.96 * sqrt(diag(cov_matrix_cv))
  risk <- data.frame(one_step = one_step_est, cil = cil, ciu = ciu, cil_cv = cil_cv, ciu_cv = ciu_cv)
  row.names(risk) <- c("E[Y(1,S(1))]", "E[Y(0,S(0))]", "E[Y(1,S(0))]", "E[Y(0,S(1))]")

  # effects 
  total_eff <- psi11n / psi00n
  grad_log_total_eff <- matrix(c(1 / psi11n, -1 / psi00n, 0, 0), ncol = 1)
  ci_total_eff <- get_ci(total_eff, grad_log_total_eff, cov_matrix)
  ci_total_eff_cv <- get_ci(total_eff, grad_log_total_eff, cov_matrix_cv)

  indirect_eff <- psi11n / psi10n
  grad_log_indirect_eff <- matrix(c(1 / psi11n, 0, -1 / psi10n, 0), ncol = 1)
  ci_indirect_eff <- get_ci(indirect_eff, grad_log_indirect_eff, cov_matrix)
  ci_indirect_eff_cv <- get_ci(indirect_eff, grad_log_indirect_eff, cov_matrix_cv)

  direct_eff <- psi10n / psi00n
  grad_log_direct_eff <- matrix(c(0, -1 / psi00n, 1 / psi10n, 0), ncol = 1)
  ci_direct_eff <- get_ci(direct_eff, grad_log_direct_eff, cov_matrix)
  ci_direct_eff_cv <- get_ci(direct_eff, grad_log_direct_eff, cov_matrix_cv)

  direct_eff2 <- psi11n / psi01n
  grad_log_direct_eff2 <- matrix(c(1 / psi11n, 0, 0, -1 / psi01n), ncol = 1)
  ci_direct_eff2 <- get_ci(direct_eff2, grad_log_direct_eff2, cov_matrix)
  ci_direct_eff2_cv <- get_ci(direct_eff2, grad_log_direct_eff2, cov_matrix_cv)
  
  indirect_eff2 <- psi01n / psi00n
  grad_log_indirect_eff2 <- matrix(c(0, -1 / psi00n, 0, 1 / psi01n), ncol = 1)
  ci_indirect_eff2 <- get_ci(indirect_eff2, grad_log_indirect_eff2, cov_matrix)
  ci_indirect_eff2_cv <- get_ci(indirect_eff2, grad_log_indirect_eff2, cov_matrix_cv)


  eff <- data.frame(effect = c("E[Y(1,S(1))] / E[Y(0,S(0))]", "E[Y(1,S(0))] / E[Y(0,S(0))]", "E[Y(1,S(1))] / E[Y(1,S(0))]"),
                    one_step_est = c(total_eff, direct_eff, indirect_eff),
                    cil = c(ci_total_eff[1], ci_direct_eff[1], ci_indirect_eff[1]),
                    ciu = c(ci_total_eff[2], ci_direct_eff[2], ci_indirect_eff[2]),
                    cil_cv = c(ci_total_eff_cv[1], ci_direct_eff_cv[1], ci_indirect_eff_cv[1]),
                    ciu_cv = c(ci_total_eff_cv[2], ci_direct_eff_cv[2], ci_indirect_eff_cv[2]))
  row.names(eff) <- c("Total", "Direct", "Indirect")

  eff2 <- data.frame(effect = c("E[Y(1,S(1))] / E[Y(0,S(0))]", "E[Y(1,S(1))] / E[Y(0,S(1))]", "E[Y(0,S(1))] / E[Y(0,S(0))]"),
                    one_step_est = c(total_eff, direct_eff2, indirect_eff2),
                    cil = c(ci_total_eff[1], ci_direct_eff2[1], ci_indirect_eff2[1]),
                    ciu = c(ci_total_eff[2], ci_direct_eff2[2], ci_indirect_eff2[2]),
                    cil_cv = c(ci_total_eff_cv[1], ci_direct_eff2_cv[1], ci_indirect_eff2_cv[1]),
                    ciu_cv = c(ci_total_eff_cv[2], ci_direct_eff2_cv[2], ci_indirect_eff2_cv[2]))
  row.names(eff2) <- c("Total", "Direct", "Indirect")
  return(list(risk = risk, eff = eff, eff2 = eff2))
}

#' Helper function to compute 95 percent confidence intervals on log-scale
#' 
#' @param est The effect estimate (on original scale)
#' @param grad The gradient of the effect transformation
#' @param cov_matrix The estimated covariance matrix of the counterfactual means 
#' 
#' @return A vector of confidence interval limits. 
get_ci <- function(est, grad, cov_matrix){
  se_eff <- sqrt(t(grad) %*% cov_matrix %*% grad)      
  ci <- exp(log(est) + c(-1.96, 1.96) * rep(se_eff, 2))
  return(ci)
}

#' Helper function to replace missing values
#' @param x The vector with missing values to be replaced
#' @param val The value to replace with
#' @return A vector of the same length as \code{x}, but with \code{NA}s replaced. 
replace_nas <- function(x, val = -999){
  x[is.na(x)] <- val
  return(x)
}

#' Helper function make outcome of extra nuisance regression
make_D1a <- function(a, A, Y, C, gA, gC, gAS, QY_WAS){
  H1 <- as.numeric(A == a & C == 1) / ( gA * gC )
  H2 <- ( gA / ( 1 - gA ) ) * ( 1 - gAS ) / gAS
  D1 <- H1 * H2 * (Y - QY_WAS)
  return(as.numeric(D1))
}



#' Temporary fix for convex combination method mean squared error
#' Relative to existing implementation, we reduce the tolerance at which
#' we declare predictions from a given algorithm the same as another
tmp_method.CC_LS <- function() {
  computeCoef <- function(Z, Y, libraryNames, verbose, obsWeights,
                          errorsInLibrary = NULL, ...) {
    cvRisk <- apply(Z, 2, function(x) {
      mean(obsWeights * (x -
        Y)^2)
    })
    names(cvRisk) <- libraryNames
    compute <- function(x, y, wt = rep(1, length(y))) {
      wX <- sqrt(wt) * x
      wY <- sqrt(wt) * y
      D <- crossprod(wX)
      d <- crossprod(wX, wY)
      A <- cbind(rep(1, ncol(wX)), diag(ncol(wX)))
      bvec <- c(1, rep(0, ncol(wX)))
      fit <- tryCatch(
        {
          quadprog::solve.QP(
            Dmat = D, dvec = d, Amat = A,
            bvec = bvec, meq = 1
          )
        },
        error = function(e) {
          out <- list()
          class(out) <- "error"
          out
        }
      )
      invisible(fit)
    }
    modZ <- Z
    naCols <- which(apply(Z, 2, function(z) {
      all(z == 0)
    }))
    anyNACols <- length(naCols) > 0
    if (anyNACols) {
      warning(paste0(
        paste0(libraryNames[naCols], collapse = ", "),
        " have NAs.", "Removing from super learner."
      ))
    }
    tol <- 4
    dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
    anyDupCols <- length(dupCols) > 0
    # if (anyDupCols) {
    #   warning(paste0(
    #     paste0(libraryNames[dupCols], collapse = ", "),
    #     " are duplicates of previous learners.", " Removing from super learner."
    #   ))
    # }
    if (anyDupCols | anyNACols) {
      rmCols <- unique(c(naCols, dupCols))
      modZ <- Z[, -rmCols, drop = FALSE]
    }
    fit <- compute(x = modZ, y = Y, wt = obsWeights)
    if (class(fit) != "error") {
      coef <- fit$solution
    } else {
      coef <- rep(0, ncol(Z))
      coef[which.min(cvRisk)] <- 1
    }
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] <- 0
    }
    if (class(fit) != "error") {
      if (anyDupCols | anyNACols) {
        ind <- c(seq_along(coef), rmCols - 0.5)
        coef <- c(coef, rep(0, length(rmCols)))
        coef <- coef[order(ind)]
      }
      coef[coef < 1e-04] <- 0
      coef <- coef / sum(coef)
    }
    if (!sum(coef) > 0) {
      warning("All algorithms have zero weight", call. = FALSE)
    }
    list(cvRisk = cvRisk, coef = coef, optimizer = fit)
  }
  computePred <- function(predY, coef, ...) {
    predY %*% matrix(coef)
  }
  out <- list(
    require = "quadprog", computeCoef = computeCoef,
    computePred = computePred
  )
  invisible(out)
}


#' Temporary fix for convex combination method negative log-likelihood loss
#' Relative to existing implementation, we reduce the tolerance at which
#' we declare predictions from a given algorithm the same as another.
#' Note that because of the way \code{SuperLearner} is structure, one needs to
#' install the optimization software separately.
tmp_method.CC_nloglik <- function() {
  computePred <- function(predY, coef, control, ...) {
    if (sum(coef != 0) == 0) {
      stop("All metalearner coefficients are zero, cannot compute prediction.")
    }
    stats::plogis(trimLogit(predY[, coef != 0], trim = control$trimLogit) %*%
      matrix(coef[coef != 0]))
  }
  computeCoef <- function(Z, Y, libraryNames, obsWeights, control,
                          verbose, ...) {
    tol <- 4
    dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
    anyDupCols <- length(dupCols) > 0
    modZ <- Z
    if (anyDupCols) {
      # warning(paste0(
      #   paste0(libraryNames[dupCols], collapse = ", "),
      #   " are duplicates of previous learners.", " Removing from super learner."
      # ))
      modZ <- modZ[, -dupCols, drop = FALSE]
    }
    modlogitZ <- trimLogit(modZ, control$trimLogit)
    logitZ <- trimLogit(Z, control$trimLogit)
    cvRisk <- apply(logitZ, 2, function(x) {
      -sum(2 * obsWeights *
        ifelse(Y, stats::plogis(x, log.p = TRUE), stats::plogis(x,
          log.p = TRUE,
          lower.tail = FALSE
        )))
    })
    names(cvRisk) <- libraryNames
    obj_and_grad <- function(y, x, w = NULL) {
      y <- y
      x <- x
      function(beta) {
        xB <- x %*% cbind(beta)
        loglik <- y * stats::plogis(xB, log.p = TRUE) + (1 -
          y) * stats::plogis(xB, log.p = TRUE, lower.tail = FALSE)
        if (!is.null(w)) {
          loglik <- loglik * w
        }
        obj <- -2 * sum(loglik)
        p <- stats::plogis(xB)
        grad <- if (is.null(w)) {
          2 * crossprod(x, cbind(p - y))
        } else {
          2 * crossprod(x, w * cbind(p - y))
        }
        list(objective = obj, gradient = grad)
      }
    }
    lower_bounds <- rep(0, ncol(modZ))
    upper_bounds <- rep(1, ncol(modZ))
    if (anyNA(cvRisk)) {
      upper_bounds[is.na(cvRisk)] <- 0
    }
    r <- tryCatch(
      {
        nloptr::nloptr(
          x0 = rep(1 / ncol(modZ), ncol(modZ)),
          eval_f = obj_and_grad(Y, modlogitZ, w = obsWeights), 
          lb = lower_bounds,
          ub = upper_bounds, eval_g_eq = function(beta) {
            (sum(beta) -
              1)
          }, eval_jac_g_eq = function(beta) rep(1, length(beta)),
          opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_abs = 1e-08)
        )
      },
      error = function(e) {
        out <- list()
        class(out) <- "error"
        out
      }
    )
    if (class(r) != "error") {
      if(r$status == 0){
        coef <- r$solution
      } else {
        coef <- rep(0, ncol(modZ))
        coef[which.min(cvRisk)[1]] <- 1
      }
    } else {
      coef <- rep(0, ncol(modZ))
      coef[which.min(cvRisk)[1]] <- 1
    }
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] <- 0
    }
    if (anyDupCols) {
      ind <- c(seq_along(coef), dupCols - 0.5)
      coef <- c(coef, rep(0, length(dupCols)))
      coef <- coef[order(ind)]
    }
    coef[coef < 1e-04] <- 0
    coef <- coef / sum(coef)
    out <- list(cvRisk = cvRisk, coef = coef, optimizer = r)
    return(out)
  }
  list(require = "nloptr", computeCoef = computeCoef, computePred = computePred)
}
