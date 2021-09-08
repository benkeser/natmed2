#' Natural in/direct effects with two-phase sampling of the mediator
#' 
#' @param W A \code{data.frame} of baseline covariates
#' @param B A \code{vector} of 0/1 indicating whether obs. is in bridging population
#' @param A A numeric vector of treatment assignment
#' @param T A numeric vector of trial
#' @param R A numeric indicator of having mediator measured in phase 3 study
#' @param S A \code{data.frame} of measured mediators
#' @param Y The binary outcome of interest
#' @param glm_gR A glm formula for the RHS of the estimation of 
#' P(R = 1 | W, A, T, Y). The variable-names that can be used in the formula are:
#' \code{colnames(W)}, \code{A}, \code{T}, and \code{Y}.
#' @param gRn A vector of known or external estimate of P(R = 1 | W, A, T, Y). If supplied then
#' \code{glm_gR} is ignored. 
#' @param glm_gA A glm formula for the RHS of the estimation of P(A = 1 | W). Defaults
#' to a main terms regression.
#' @param SL_gA A \code{SuperLearner} library for estimation of P(A = 1 | W).
#' @param glm_gAS A glm formula for the RHS of the estimation of P(A = 1 | W, S)
#' and P(A = 2 | A != 1, W, S). Defaults to a main terms regression. 
#' @param SL_gAS A \code{SuperLearner} library for estimation of P(A = 1 | W, S)
#' and P(A = 2 | A != 1, W, S).
#' @param glm_QY_WAS A glm formula for the RHS of the estimation of P(Y = 1 | T = 2, A, W, S). Defaults
#' to a main terms regression. This regression is fit using inverse weights \code{R}/\code{gRn}, which may trigger
#' warnings from \code{glm} about the outcome not being between 0 and 1. These warnings are 
#' safe to ignore.
#' @param SL_QY_WAS A \code{SuperLearner} library for estimation of the outcome regression 
#' described above. This regression is fit using inverse weights \code{R}/\code{gRn}. 

#' @param glm_QD A glm formula for the RHS of the estimation of 
#' E[( 1 / (gA2n) ) ( gA2n / (gA1n) ) ( (gA1S) / gA2S ) (Y - QY_WAS) | T = 2, A = 2, W, Y]. 
#' The available names for the regression formula are: \code{colnames(W)}, \code{Y}.
#' @param SL_QD A \code{SuperLearner} library for estimation of a relevant piece of the projection
#' of the full data EIF onto the tangent space of gR. 
#' @param glm_QY_W A \code{glm} formula for the RHS of the estimation of 
#' E[QY_WAS | T = 1, A = 1, W].
#' @param SL_QY_W A \code{SuperLearner} library for estimation of 
#' E[QY_WAS | T = 1, A = 1, W].
#' @param glm_QY_WA A \code{glm} formula for the outcome regression E[Y | A, W], which is used
#' in estimation of E[Y(a, S(a))]. Available variable names are: \code{colnames(W)}, \code{A}. 
#' @param SL_QY_WA A \code{SuperLearner} library for estimation of the outcome regression used to 
#' compute estimates of E[Y(a, S(a))]. 
#' @param seed The seed set before each \code{SuperLearner} fit.
#' @param ... Other options (not currently used)
#' 
#' @export
#' @importFrom stats predict glm
#' @importFrom SuperLearner SuperLearner
#' 
#' @return A list with named entries \describe{
#' \item{blah}{blah blah}
#' }
#' 
#' @examples
#' n <- 500
#' W1 <- rbinom(n, 1, 0.5)
#' W2 <- rnorm(n, 0, 1)
#' A <- rbinom(n, 1, 0.5)
#' S <- W1 / 4 - W2 / 3 + A + rnorm(n)
#' Y <- rbinom(n, 1, plogis(-2 + A + W1 / 2 - S / 2))
#' 
#' # add censoring
#' C <- rbinom(n, 1, plogis(2 + W1 / 2 - W2 / 3))
#' # arbitrary fill in
#' Y[C == 0] <- -999
#' 
#' R <- rep(0, n)
#' # case-cohort sampling
#' R <- rbinom(n, 1, 0.25)
#' R[Y == 1] <- 1
#' 
#' fit <- natmed2(
#'   W = data.frame(W1 = W1, W2 = W2), 
#'   A = A, R = R, S = S, C = C, Y = Y
#' )
natmed2_fx <- function(
  W, B, A, T, R, S, Y,
  glm_gR = ".^2",
  gRn = NULL,
  glm_gA = ".", 
  SL_gA,
  glm_gAS = paste0(paste0(colnames(W), collapse = " + "), " + ", paste0(colnames(S), collapse = " + ")), 
  SL_gAS,
  glm_QY_WAS = paste0(paste0(colnames(W), collapse = " + "), " + A + ", paste0(colnames(S), collapse = " + ")), 
  SL_QY_WAS, # Y | R = 1, T = 2, W, A, S
  glm_QD = ".", 
  SL_QD, # First piece of EIF | R = 1, T = 1, A = 1, W, Y
  glm_QY_W = ".", 
  SL_QY_W, # QY_WAS | A = a_2, W
  glm_QY_WA = ".", 
  SL_QY_WA, # Y | W, A
  seed = 1, 
  tol_gA = 1 / sqrt(length(A)),
  tol_gAS = 1 / sqrt(sum(R == 1)),
  ...
){
  n <- length(Y)

  if(!is.null(glm_gR) & is.null(gRn)){
    gRn_1 <- rep(1, n)

    gR_fit <- stats::glm(paste0("R ~ ", glm_gR), family = binomial(),
                  data = data.frame(R = R, W, Y = Y, A = A)[T == 2,])
    gRn_1[T == 2] <- stats::predict(gR_fit, type = "response")
  }else if(!is.null(gRn)){
    gRn_1 <- gRn
  }else{
    stop("specify glm formula for glm_gR or gRn")
  }

  if(!is.null(glm_gA)){
    # P(A = 1 | W)
    gA1_fit <- stats::glm(paste0("A1 ~ ", glm_gA), family = binomial(), 
                  data = data.frame(A1 = as.numeric(A == 1), W))
    # P(A = 2 | A != 1, W)
    gA2_fit <- stats::glm(paste0("A2 ~ ", glm_gA), family = binomial(), 
                  data = data.frame(A2 = as.numeric(A == 2), W)[A != 1, ])

    gAn_1 <- stats::predict(gA1_fit, type = "response")
    gAn_2 <- stats::predict(gA2_fit, type = "response", newdata = data.frame(A2 = A, W)) * 
                (1 - gAn_1)
    gAn_0 <- 1 - (gAn_1 + gAn_2)
  }else{
    # P(A = 1 | W)
    gA1_fit <- SuperLearner::SuperLearner(Y = as.numeric(A == 1), X = W,
                                         family = binomial(), 
                                         SL.library = SL_gA, 
                                         method = tmp_method.CC_nloglik())
    # P(A = 2 | A != 1, W)
    gA2_fit <- SuperLearner::SuperLearner(Y = as.numeric(A[A != 1] == 2), X = W[A != 1, ,drop = FALSE],
                                         family = binomial(), 
                                         SL.library = SL_gA, 
                                         method = tmp_method.CC_nloglik())
    gAn_1 <- stats::predict(gA1_fit, type = "response")[[1]]
    gAn_2 <- stats::predict(gA2_fit, type = "response", newdata = data.frame(A = A, W))[[1]] * 
                (1 - gAn_1)
    gAn_0 <- 1 - (gAn_1 + gAn_2)
  }
  gAn_1 <- g_truncate(gAn_1, tol = tol_gA)
  gAn_2 <- g_truncate(gAn_2, tol = tol_gA)
  gAn_0 <- g_truncate(gAn_0, tol = tol_gA)

  if(!is.null(glm_gAS)){
    # P(A = 1 | W, S)
    gA1S_fit <- stats::glm(paste0("A1 ~ ", glm_gAS), family = binomial(),
                   data = data.frame(A1 = as.numeric(A == 1), S, W, wt = R / gRn_1)[R == 1,],
                   weights = wt)
    # P(A = 2 | A != 1, W, S)
    gA2S_fit <- stats::glm(paste0("A2 ~ ", glm_gAS), family = binomial(),
                   data = data.frame(A2 = as.numeric(A == 2), S, W, wt = R / gRn_1)[R == 1 & A != 1,],
                   weights = wt)

    gASn_1 <- rep(NA, n)
    gASn_1[R == 1] <- stats::predict(gA1S_fit, type = "response", 
                              newdata = data.frame(S, W)[R == 1,])
    gASn_2 <- rep(NA, n)
    gASn_2[R == 1] <- stats::predict(gA2S_fit, type = "response", 
                              newdata = data.frame(S, W)[R == 1,]) * (1 - gASn_1[R == 1])
  }else{
    set.seed(seed)
    # P(A = 1 | W, S)
    gA1S_fit <- SuperLearner::SuperLearner(Y = as.numeric(A[R == 1] == 1), 
                            X = data.frame(S, W)[R == 1, ],
                            obsWeights = (R / gRn_1)[R == 1],
                            family = binomial(), 
                            SL.library = SL_gAS,
                            method = tmp_method.CC_nloglik())
    # P(A = 2 | A!= 1, W, S)
    gA2S_fit <- SuperLearner::SuperLearner(Y = as.numeric(A[R == 1 & A != 1] == 2), 
                            X = data.frame(S, W)[R == 1 & A != 1, ],
                            obsWeights = (R / gRn_1)[R == 1 & A != 1],
                            family = binomial(), 
                            SL.library = SL_gAS,
                            method = tmp_method.CC_nloglik())
    gASn_1 <- rep(NA, n)    
    gASn_1[R == 1] <- gA1S_fit$SL.predict
    gASn_2 <- rep(NA, n)    
    gASn_2[R == 1] <- predict(gA2S_fit, newdata = data.frame(S, W)[R == 1,])[[1]] * 
                          (1 - gASn_1[R == 1])    
  }
  gASn_1 <- g_truncate(gASn_1, tol = tol_gAS)
  gASn_2 <- g_truncate(gASn_2, tol = tol_gAS)

  if(!is.null(glm_QY_WAS)){
    QY_WAS_fit <- stats::glm(paste0("Y ~ ", glm_QY_WAS), family = binomial(),
                      data = data.frame(Y = Y, A = A, S, W, wt = (R / gRn_1))[R == 1 & T == 2, ],
                      weights = wt)
    QY_WASn_A2 <- rep(NA, n)
    QY_WASn_A2[R == 1] <- stats::predict(QY_WAS_fit, type = "response", 
                                  newdata = data.frame(Y = Y, A = 2, S, W)[R == 1, ])
  }else{
    set.seed(seed)
    QY_WAS_fit <- SuperLearner::SuperLearner(Y = Y[R == 1 & T == 2], 
                               X = data.frame(A = A, S, W)[R == 1 & T == 2, ],
                               obsWeights = (R / gRn_1)[R == 1 & T == 2],
                               family = binomial(), 
                               SL.library = SL_QY_WAS,
                               method = tmp_method.CC_nloglik())
    QY_WASn_A2 <- rep(NA, n)
    QY_WASn_A2[R == 1] <- stats::predict(QY_WAS_fit, newdata = data.frame(A = 2, S, W)[R == 1,])[[1]]    
  }

  gBn_1 = rep(mean(B == 1), n)

  DY = as.numeric(A == 2 & T == 2 & B == 1) / ( gAn_2 * gBn_1 ) * 
        ( gAn_2 / gAn_1 ) * ( gASn_1 / gASn_2 ) * ( Y - QY_WASn_A2 )
  
  if(!is.null(glm_QD)){
    QD_fit <- stats::glm(paste0("DY ~ ", glm_QD),
       data = data.frame(DY = DY, W, Y = Y)[B == 1 & A == 2 & T == 2 & R == 1, ])
    E_QD <- rep(NA, n)
    E_QD[!(B == 1 & A == 2 & T == 2)] <- 0
    E_QD[B == 1 & A == 2 & T == 2] <- stats::predict(QD_fit, type = "response", 
          newdata = data.frame(DY = DY, W, Y = Y)[B == 1 & A == 2 & T == 2, ])
  }else{
    set.seed(seed)
    QD_fit <- SuperLearner::SuperLearner(Y = DY[B == 1 & A == 2 & T == 2 & R == 1],
                                         X = data.frame(W, Y = Y)[B == 1 & A == 2 & T == 2 & R == 1, ],
                                         family = stats::gaussian(),
                                         SL.library = SL_QD,
                                         method = tmp_method.CC_nloglik())
    E_QD <- rep(NA, n)
    E_QD[!(B == 1 & A == 2 & T == 2)] <- 0
    E_QD[B == 1 & A == 2 & T == 2] <- stats::predict(QD_fit, type = "response", 
              newdata = data.frame(W, Y = Y)[B == 1 & A == 2 & T == 2, ])[[1]]
  }
  if(!is.null(glm_QY_W)){
    QY_W_fit <- stats::glm(paste0("QY_WASn_A2 ~ ", glm_QY_W), family = binomial(),
                          data = data.frame(QY_WASn_A2 = QY_WASn_A2, W)[T == 1 & A == 1, ])
    QY_W <- rep(NA, n)
    QY_W <- predict(QY_W_fit, newdata = W, type = "response")
  }else{
    set.seed(seed)
    QY_W_fit <- SuperLearner::SuperLearner(Y = QY_WASn_A2[T == 1 & A == 1], 
                                   X = W[T == 1 & A == 1, ],
                                   family = gaussian(), 
                                   SL.library = SL_QY_W,
                                   method = tmp_method.CC_LS())
    QY_W <- as.numeric(
      stats::predict(QY_W_fit, newdata = W)[[1]]
    )
    QY_W[QY_W <= 0] <- 0
    QY_W[QY_W >= 1] <- 1
  }

  plug_in <- mean(QY_W[B == 1])
  eif1 <- R / gRn_1 * ( (B == 1) / gBn_1 ) * ( (A == 2 & T == 2) / (gAn_1) ) * ( gASn_1 / gASn_2 ) * ( Y - QY_WASn_A2 )
  eif1[is.na(eif1)] <- 0
  eif2 <- ( (B == 1) / gBn_1 ) * ( (A == 1 & T == 1) / gAn_1 ) * ( QY_WASn_A2 -  QY_W )
  eif2[is.na(eif2)] <- 0
  eif3 <- ( (B == 1) / gBn_1 ) * ( QY_W - plug_in )
  eif4 <- E_QD / gRn_1 * ( R - gRn_1 )
  
  eif <- eif1 + eif2 + eif3 - eif4

  one_step <- plug_in + mean(eif)
  se <- sqrt( var(eif) / n )

  return(list(plug_in = plug_in, one_step = one_step, se = se))
}