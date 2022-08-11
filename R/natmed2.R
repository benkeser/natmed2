#' Natural in/direct effects with two-phase sampling of the mediator
#' 
#' @param W A \code{data.frame} of baseline covariates
#' @param A A numeric vector of binary treatment assignment
#' @param R A numeric indicator of having mediator measured
#' @param S A data.frame containing mediator variable/s
#' @param C A numeric indicator of having the outcome \code{Y} measured, 
#' i.e., a measure of being not right-censored. 
#' @param Y The binary outcome of interest
#' @param glm_gR A glm formula for the RHS of the estimation of 
#' P(R = 1 | W, A, C, CY). The variable-names that can be used in the formula are:
#' \code{colnames(W)}, \code{A}, \code{CY11} (indicator that \code{C}=1 and \code{Y}=1),
#' and \code{CY10} (indicator that \code{C}=1 and \code{Y}=0). 
#' @param gRn A vector of known or external estimate of P(R = 1 | W, A, C, CY). If supplied then
#' \code{glm_gR} is ignored. 
#' @param glm_gC A glm formula for the RHS of the estimation of 
#' P(C = 1 | A, W). Defaults to a main terms regression.
#' @param SL_gC A \code{SuperLearner} library for estimation of censoring probability
#' described above. 
#' @param glm_gAS A glm formula for the RHS of the estimation of P(A = 1 | W, S). Defaults
#' to a main terms regression. 
#' @param SL_gAS A \code{SuperLearner} library for estimation of the mediator-conditional
#' propensity score described above.
#' @param glm_QY_WAS A glm formula for the RHS of the estimation of P(Y = 1 | A, W, S). Defaults
#' to a main terms regression. This regression is fit using inverse weights \code{R}/\code{gRn}, which may trigger
#' warnings from \code{glm} about the outcome not being between 0 and 1. These warnings are 
#' safe to ignore.
#' @param SL_QY_WAS A \code{SuperLearner} library for estimation of the outcome regression 
#' described above. This regression is fit using inverse weights \code{R}/\code{gRn}. 
#' @param glm_QY_WACY A glm formula for the RHS of the estimation of E[P(Y = 1 | W, A = a1, S) | W, A = a2, C, CY].
#' The available names for the regression formula are: \code{colnames(W)}, \code{CY11}, \code{CY10}.
#' @param SL_QY_WACY A \code{SuperLearner} library for estimation of the sequential outcome regression
#' described above. 
#' @param glm_QD_WACY A glm formula for the RHS of the estimation of 
#' E[I(A = a1, C = 1) / (gAn gCn) gAn / (1 - gAN) (1 - gAS) / gAS (Y - QY_WAS) | R = 1, W, A = a1, C, CY]. 
#' The available names for the regression formula are: \code{colnames(W)}, \code{CY11}, \code{CY10}.
#' @param SL_QD_WACY A \code{SuperLearner} library for estimation of a relevant piece of the projection
#' of the full data EIF onto the tangent space of gR. 
#' @param glm_QY_W A \code{glm} formula for the RHS of the estimation of 
#' E[QY_WACY | A = a2, W]. 
#' @param SL_QY_W A \code{SuperLearner} library for estimation of the final sequential regression. 
#' @param glm_QY_WA A \code{glm} formula for the outcome regression E[Y | A, W], which is used
#' in estimation of E[Y(a, S(a))]. Available variable names are: \code{colnames(W)}, \code{A}. 
#' @param SL_QY_WA A \code{SuperLearner} library for estimation of the outcome regression used to 
#' compute estimates of E[Y(a, S(a))]. 
#' @param seed The seed set before each \code{SuperLearner} fit.
#' @param lazy If \code{TRUE}, then the non-simplified version of the EIF is used, in lieu of
#' the extra sequential regression. 
#' @param ... Other options (not currently used)
#' 
#' @export
#' @importFrom stats predict glm
#' @importFrom SuperLearner SuperLearner
#' 
#' @return A list with named entries \describe{
#' \item{risk}{A \code{data.frame} that contains point estimates and confidence intervals 
#' of several counterfactual probabilities. Two confidence intervals are shown. \code{cil} 
#' and \code{ciu} indicate the limits of a 95% confidence interval based on the empirical 
#' variance of the estimated influence function. \code{cil_cv} and \code{ciu_cv} indicate
#' the limits of a 95% confidence interval based on a cross-validated estimate of the variance 
#' of the estimated influence function. If \code{SuperLearner} was not used for any of the 
#' nuisance regressions, then the two confidence intervals will coincide. }
#' \item{eff}{A \code{data.frame} that shows the total, direct, and indirect effects along 
#' with confidence intervals. These intervals are created on the log-scale using the delta-method.}
#' \item{eff2}{A \code{data.frame} that shows the total, direct, and indirect effects based
#' on an alternative decomposition.}
#' \item{cov}{The estimated covariance matrix of the estimates}
#' \item{cov_cv}{The cross-validated estimate of the covariance matrix of the estimates}
#' }
#' 
#' @examples
#' n <- 500
#' W1 <- rbinom(n, 1, 0.5)
#' W2 <- rnorm(n, 0, 1)
#' A <- rbinom(n, 1, 0.5)
#' S <- data.frame(S = W1 / 4 - W2 / 3 + A + rnorm(n))
#' Y <- rbinom(n, 1, plogis(-2 + A + W1 / 2 - S$S / 2))
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
natmed2 <- function(
  W, A, R, S, C, Y,
  glm_gR = ".^2",
  gRn = NULL,
  glm_gC = ".", 
  SL_gC = NULL,
  glm_gA = ".", 
  glm_gAS = paste0(paste0(colnames(W), collapse = " + "), " + ",
                   paste0(colnames(S), collapse = " + ")), 
  SL_gAS,
  glm_QY_WAS = paste0(paste0(colnames(W), collapse = " + "), "+ A + ",
                      paste0(colnames(S), collapse = " + ")), 
  SL_QY_WAS, # Y | R = 1, C = 1, W, A, S
  glm_QY_WACY = ".", 
  SL_QY_WACY, # QY_WAS | R = 1, W, A, C, CY 
  glm_QD_WACY = ".", 
  SL_QD_WACY, # First piece of EIF | R = 1, W, A = a, C, CY 
  glm_QD_WACY_lazy = ".", 
  SL_QD_WACY_lazy, # EIF | R = 1, W, A, C, CY 
  glm_QY_W = ".", 
  SL_QY_W, # QY_WACY | A = a_2, W
  glm_QY_WA = ".", 
  SL_QY_WA, # Y | C = 1, W, A
  lazy = FALSE, 
  seed = 1, 
  tol_gA = 1 / sqrt(length(A)),
  tol_gAS = 1 / sqrt(sum(R == 1)),
  tol_gC = 1 / sqrt(length(A)), 
  ...
){
  n <- length(Y)

  # collapse C and Y into dummy variables
  tmpY <- Y; tmpY[is.na(tmpY)] <- 0
  CY10 = as.numeric(tmpY == 0 & C == 1)
  CY11 = as.numeric(tmpY == 1 & C == 1)

  if(!is.null(glm_gR) & is.null(gRn)){
    # glm_gR a formula of colnames(W), A, CY10, CY11
    gR_fit <- stats::glm(paste0("R ~ ", glm_gR), family = binomial(),
                  data = data.frame(R = R, W, CY10 = CY10, CY11 = CY11, A = A))
    gRn_1 <- stats::predict(gR_fit, type = "response")
  }else if(!is.null(gRn)){
    gRn_1 <- gRn
  }else{
    stop("specify glm formula for glm_gR or gRn")
  }

  if(!is.null(glm_gA)){
    gA_fit <- stats::glm(paste0("A ~ ", glm_gA), family = binomial(), 
                  data = data.frame(A = A, W))
    gAn_1 <- stats::predict(gA_fit, type = "response")
  }else{
    stop("specify glm formula for glm_gA")
  }
  gAn_1 <- g_truncate(gAn_1, tol = tol_gA)

  if(!is.null(glm_gC)){
    gC_fit <- stats::glm(paste0("C ~ ", glm_gC), family = binomial(),
                  data = data.frame(C = C, A = A, W))
    gCn_1_A0 <- stats::predict(gC_fit, type = "response", 
                        newdata = data.frame(C = C, A = 0, W))
    gCn_1_A0_cv <- gCn_1_A0
    gCn_1_A1 <- stats::predict(gC_fit, type = "response", 
                        newdata = data.frame(C = C, A = 1, W))
    gCn_1_A1_cv <- gCn_1_A1
  }else{
    set.seed(seed)
    gC_fit <- SuperLearner::SuperLearner(Y = C, X = data.frame(A = A, W),
                           family = binomial(), 
                           SL.library = SL_gC,
                           # method = tmp_method.CC_nloglik(),
                           control = list(saveCVFitLibrary = TRUE))
    gCn_1_A0 <- stats::predict(gC_fit, type = "response", 
                        newdata = data.frame(C = C, A = 0, W))[[1]]
    gCn_1_A1 <- stats::predict(gC_fit, type = "response", 
                        newdata = data.frame(C = C, A = 1, W))[[1]]

    # get partially cross-validated predictions
    gCn_1_A0_cv <- partial_cv_preds(gC_fit, a_0 = 0, W = W, include = rep(TRUE, n), family = stats::binomial())
    gCn_1_A1_cv <- partial_cv_preds(gC_fit, a_0 = 1, W = W, include = rep(TRUE, n), family = stats::binomial())
  }
  gCn_1_A0 <- g_truncate(gCn_1_A0, tol = tol_gC, bound_away_from = 0)
  gCn_1_A0_cv <- g_truncate(gCn_1_A0_cv, tol = tol_gC, bound_away_from = 0)
  gCn_1_A1 <- g_truncate(gCn_1_A1, tol = tol_gC, bound_away_from = 0)
  gCn_1_A1_cv <- g_truncate(gCn_1_A1_cv, tol = tol_gC, bound_away_from = 0)

  if(!is.null(glm_gAS)){
    gAS_fit <- stats::glm(paste0("A ~ ", glm_gAS), family = binomial(),
                   data = data.frame(A = A, S, W, wt = R / gRn_1)[R == 1,],
                   weights = wt)
    gASn_1 <- rep(NA, n)
    gASn_1[R == 1] <- stats::predict(gAS_fit, type = "response", 
                              newdata = data.frame(C = C, S, W)[R == 1,])
    gASn_1_cv <- gASn_1
  }else{
    set.seed(seed)
    gAS_fit <- SuperLearner::SuperLearner(Y = A[R == 1], X = data.frame(S, W)[R == 1, ],
                            obsWeights = (R / gRn_1)[R == 1],
                            family = binomial(), 
                            SL.library = SL_gAS,
                            cvControl = list(stratify = TRUE),
                            # method = tmp_method.CC_nloglik(),
                            control = list(saveCVFitLibrary = TRUE))
    gASn_1 <- rep(NA, n)
    gASn_1[R == 1] <- gAS_fit$SL.predict

    # get partially cross-validated predictions
    gASn_1_cv <- rep(NA, n)
    gASn_1_cv[R == 1] <- partial_cv_preds(gAS_fit, easy = TRUE, family = stats::binomial())
  }
  gASn_1 <- g_truncate(gASn_1, tol = tol_gAS)
  gASn_1_cv <- g_truncate(gASn_1_cv, tol = tol_gAS)

  if(!is.null(glm_QY_WAS)){
    QY_WAS_fit <- stats::glm(paste0("Y ~ ", glm_QY_WAS), family = binomial(),
                      data = data.frame(Y = Y, A = A, S, W, wt = (R / gRn_1))[R == 1 & C == 1, ],
                      weights = wt)
    QY_WASn_A1 <- rep(NA, n)
    QY_WASn_A0 <- rep(NA, n)
    QY_WASn_A1[R == 1] <- stats::predict(QY_WAS_fit, type = "response", 
                                  newdata = data.frame(Y = Y, A = 1, S, W)[R == 1, ])
    QY_WASn_A1_cv <- QY_WASn_A1
    QY_WASn_A0[R == 1] <- stats::predict(QY_WAS_fit, type = "response", 
                                  newdata = data.frame(Y = Y, A = 0, S, W)[R == 1, ])    
    QY_WASn_A0_cv <- QY_WASn_A0
  }else{
    set.seed(seed)
    QY_WAS_fit <- SuperLearner::SuperLearner(Y = Y[R == 1 & C == 1], 
                               X = data.frame(A = A, S, W)[R == 1 & C == 1, ],
                               obsWeights = (R / gRn_1)[R == 1 & C == 1],
                               family = binomial(), 
                               SL.library = SL_QY_WAS,
                               # method = tmp_method.CC_nloglik(),
                               control = list(saveCVFitLibrary = TRUE))
    QY_WASn_A1 <- rep(NA, n)
    QY_WASn_A0 <- rep(NA, n)
    QY_WASn_A1[R == 1] <- stats::predict(QY_WAS_fit, newdata = data.frame(A = 1, S, W)[R == 1,])[[1]]
    QY_WASn_A0[R == 1] <- stats::predict(QY_WAS_fit, newdata = data.frame(A = 0, S, W)[R == 1,])[[1]]
    
    # get partially cross-validated predictions
    QY_WASn_A0_cv <- partial_cv_preds_QY_WASn(QY_WAS_fit, newdata = data.frame(A = 0, S, W), R = R, C = C, family = stats::binomial())
    QY_WASn_A1_cv <- partial_cv_preds_QY_WASn(QY_WAS_fit, newdata = data.frame(A = 1, S, W), R = R, C = C, family = stats::binomial())
  }

  # compute outcome of extra nuisance regression
  DY_A1 <- make_D1a(a = 1, A = A, Y = Y, C = C, 
                    gA = gAn_1, gAS = gASn_1, gC = gCn_1_A1, QY_WAS = QY_WASn_A1)
  DY_A0 <- make_D1a(a = 0, A = A, Y = Y, C = C, 
                    gA = 1 - gAn_1, gAS = 1 - gASn_1, gC = gCn_1_A0, QY_WAS = QY_WASn_A0)

  if(!is.null(glm_QD_WACY)){
    QD_WACY_fit_A0 <- stats::glm(paste0("DY_A0 ~ ", glm_QD_WACY),
                          data = data.frame(DY_A0 = DY_A0, W, CY11 = CY11, CY10 = CY10)[A == 0 & C == 1 & R == 1, ])
    QD_WACYn_A0 <- rep(NA, n)
    QD_WACYn_A0[A == 1 | C == 0] <- 0
    QD_WACYn_A0[A == 0 & C == 1] <- stats::predict(QD_WACY_fit_A0, type = "response", 
                                            newdata = data.frame(DY_A0 = DY_A0, W, CY11 = CY11, CY10 = CY10)[A == 0 & C == 1, ])
    QD_WACYn_A0_cv <- QD_WACYn_A0

    QD_WACY_fit_A1 <- stats::glm(paste0("DY_A1 ~ ", glm_QD_WACY),
                          data = data.frame(DY_A1 = DY_A1, W, CY11 = CY11, CY10 = CY10)[A == 1 & C == 1 & R == 1, ])
    QD_WACYn_A1 <- rep(NA, n)
    QD_WACYn_A1[A == 0 | C == 0] <- 0
    QD_WACYn_A1[A == 1 & C == 1] <- stats::predict(QD_WACY_fit_A1, type = "response", 
                                            newdata = data.frame(DY_A1 = DY_A1, W, CY11 = CY11, CY10 = CY10)[A == 1 & C == 1, ])  
    QD_WACYn_A1_cv <- QD_WACYn_A1
  }else{
    set.seed(seed)
    if(length(unique(DY_A0[A == 0 & C == 1 & R == 1])) < 3){
      SL_QD_WACY <- "SL.mean"
    }
    QD_WACY_fit_A0 <- SuperLearner::SuperLearner(Y = DY_A0[A == 0 & C == 1 & R == 1], 
                                   X = data.frame(W, CY11 = CY11, CY10 = CY10)[A == 0 & C == 1 & R == 1, ], 
                                   family = gaussian(), 
                                   SL.library = SL_QD_WACY,
                                   # method = tmp_method.CC_LS(),
                                   control = list(saveCVFitLibrary = TRUE))
    QD_WACYn_A0 <- rep(NA, n)
    QD_WACYn_A0[A == 1 | C == 0] <- 0
    QD_WACYn_A0[A == 0 & C == 1] <- as.numeric(stats::predict(QD_WACY_fit_A0, newdata = data.frame(W, CY11 = CY11, CY10 = CY10)[A == 0 & C == 1, ])[[1]])

    set.seed(seed)
    if(length(unique(DY_A1[A == 1 & C == 1 & R == 1])) < 3){
      SL_QD_WACY <- "SL.mean"
    }
    QD_WACY_fit_A1 <- SuperLearner::SuperLearner(Y = DY_A1[A == 1 & C == 1 & R == 1], 
                                   X = data.frame(W, CY11 = CY11, CY10 = CY10)[A == 1 & C == 1 & R == 1, ], 
                                   family = gaussian(), 
                                   SL.library = SL_QD_WACY,
                                   # method = tmp_method.CC_LS(),
                                   control = list(saveCVFitLibrary = TRUE))
    QD_WACYn_A1 <- rep(NA, n)
    QD_WACYn_A1[A == 0 | C == 0] <- 0
    QD_WACYn_A1[A == 1 & C == 1] <- as.numeric(stats::predict(QD_WACY_fit_A1, newdata = data.frame(W, CY11 = CY11, CY10 = CY10)[A == 1 & C == 1, ])[[1]])


    # get partially cross-validated predictions
    QD_WACYn_A0_cv <- partial_cv_preds_QD_WACYn(QD_WACY_fit_A0, a = 0, newdata = data.frame(W, CY11 = CY11, CY10 = CY10), A = A, R = R, C = C, family = stats::gaussian())
    QD_WACYn_A1_cv <- partial_cv_preds_QD_WACYn(QD_WACY_fit_A1, a = 1, newdata = data.frame(W, CY11 = CY11, CY10 = CY10), A = A, R = R, C = C, family = stats::gaussian())
  }

  if(!is.null(glm_QY_WACY)){
    QY_WACY_fit_A0_A1 <- stats::glm(paste0("QY_WASn_A0 ~ ", glm_QY_WACY), family = binomial(),
                          data = data.frame(QY_WASn_A0 = QY_WASn_A0, 
                                            W, CY11 = CY11, CY10 = CY10)[A == 1 & R == 1, ])
    QY_WACYn_A0_A1 <- stats::predict(QY_WACY_fit_A0_A1, type = "response", 
                              newdata = data.frame(QY_WASn_A0 = QY_WASn_A0, 
                                                   W, CY11 = CY11, CY10 = CY10))
    QY_WACYn_A0_A1_cv <- QY_WACYn_A0_A1

    QY_WACY_fit_A1_A0 <- stats::glm(paste0("QY_WASn_A1 ~ ", glm_QY_WACY), family = binomial(),
                          data = data.frame(QY_WASn_A1 = QY_WASn_A1, 
                                            W, CY11 = CY11, CY10 = CY10)[A == 0 & R == 1, ])
    QY_WACYn_A1_A0 <- stats::predict(QY_WACY_fit_A1_A0, type = "response", 
                              newdata = data.frame(QY_WASn_A1 = QY_WASn_A1, 
                                                   W, CY11 = CY11, CY10 = CY10))
    QY_WACYn_A1_A0_cv <- QY_WACYn_A1_A0
  }else{
    set.seed(seed)
    if(length(unique(QY_WASn_A0[A == 1 & R == 1])) < 3){
      SL_QY_WACY <- "SL.mean"
    }
    QY_WACY_fit_A0_A1 <- SuperLearner::SuperLearner(Y = QY_WASn_A0[A == 1 & R == 1], 
                                      X = data.frame(W, CY11 = CY11, CY10 = CY10)[A == 1 & R == 1, ], 
                                      family = gaussian(), 
                                      SL.library = SL_QY_WACY,
                                      # method = tmp_method.CC_LS(),
                                      control = list(saveCVFitLibrary = TRUE))
    QY_WACYn_A0_A1 <- as.numeric(
      stats::predict(QY_WACY_fit_A0_A1, 
              newdata = data.frame(W, CY11 = CY11, CY10 = CY10))[[1]]
    )
    QY_WACYn_A0_A1[QY_WACYn_A0_A1 <= 0] <- 0
    QY_WACYn_A0_A1[QY_WACYn_A0_A1 >= 1] <- 1

    set.seed(seed)
    if(length(unique(QY_WASn_A1[A == 0 & R == 1])) < 3){
      SL_QY_WACY <- "SL.mean"
    }
    QY_WACY_fit_A1_A0 <- SuperLearner::SuperLearner(Y = QY_WASn_A1[A == 0 & R == 1], 
                                      X = data.frame(W, CY11 = CY11, CY10 = CY10)[A == 0 & R == 1, ], 
                                      family = gaussian(), 
                                      SL.library = SL_QY_WACY,
                                      # method = tmp_method.CC_LS(),
                                      control = list(saveCVFitLibrary = TRUE))
    QY_WACYn_A1_A0 <- as.numeric(
      stats::predict(QY_WACY_fit_A1_A0, 
              newdata = data.frame(W, CY11 = CY11, CY10 = CY10))[[1]]
    )
    QY_WACYn_A1_A0[QY_WACYn_A1_A0 <= 0] <- 0
    QY_WACYn_A1_A0[QY_WACYn_A1_A0 >= 1] <- 1

    # get partially cross-validated predictions
    QY_WACYn_A1_A0_cv <- partial_cv_preds_QY_WACYn(QY_WACY_fit_A1_A0, a = 0, newdata = data.frame(W, CY11 = CY11, CY10 = CY10), A = A, R = R, family = stats::gaussian())
    QY_WACYn_A1_A0_cv[QY_WACYn_A1_A0_cv <= 0] <- 0
    QY_WACYn_A1_A0_cv[QY_WACYn_A1_A0_cv >= 1] <- 1
    
    QY_WACYn_A0_A1_cv <- partial_cv_preds_QY_WACYn(QY_WACY_fit_A0_A1, a = 1, newdata = data.frame(W, CY11 = CY11, CY10 = CY10), A = A, R = R, family = stats::gaussian())
    QY_WACYn_A0_A1_cv[QY_WACYn_A0_A1_cv <= 0] <- 0
    QY_WACYn_A0_A1_cv[QY_WACYn_A0_A1_cv >= 1] <- 1
  }

  if(!is.null(glm_QY_W)){
    QY_W_fit_A0_A1 <- stats::glm(paste0("QY_WACYn_A0_A1 ~ ", glm_QY_W), family = binomial(),
                          data = data.frame(QY_WACYn_A0_A1 = QY_WACYn_A0_A1, W)[A == 1, ])
    QY_Wn_A0_A1 <- stats::predict(QY_W_fit_A0_A1, type = "response", 
                              newdata = data.frame(QY_WACYn_A0_A1 = QY_WACYn_A0_A1, W))
    QY_Wn_A0_A1_cv <- QY_Wn_A0_A1

    QY_W_fit_A1_A0 <- stats::glm(paste0("QY_WACYn_A1_A0 ~ ", glm_QY_W), family = binomial(),
                          data = data.frame(QY_WACYn_A1_A0 = QY_WACYn_A1_A0, W)[A == 0, ])
    QY_Wn_A1_A0 <- stats::predict(QY_W_fit_A1_A0, type = "response", 
                           newdata = data.frame(QY_WACYn_A1_A0 = QY_WACYn_A1_A0, W))
    QY_Wn_A1_A0_cv <- QY_Wn_A1_A0
  }else{
    set.seed(seed)
    if(length(unique(QY_WACYn_A0_A1[A == 1])) < 3){
      SL_QY_W <- "SL.mean"
    }
    QY_W_fit_A0_A1 <- SuperLearner::SuperLearner(Y = QY_WACYn_A0_A1[A == 1], 
                                   X = W[A == 1, ], 
                                   family = gaussian(), 
                                   SL.library = SL_QY_W,
                                   # method = tmp_method.CC_LS(),
                                   control = list(saveCVFitLibrary = TRUE))
    QY_Wn_A0_A1 <- as.numeric(
      stats::predict(QY_W_fit_A0_A1, newdata = W)[[1]]
    )
    QY_Wn_A0_A1[QY_Wn_A0_A1 <= 0] <- 0
    QY_Wn_A0_A1[QY_Wn_A0_A1 >= 1] <- 1

    set.seed(seed)
    if(length(unique(QY_WACYn_A1_A0[A == 0])) < 3){
      SL_QY_W <- "SL.mean"
    }
    QY_W_fit_A1_A0 <- SuperLearner::SuperLearner(Y = QY_WACYn_A1_A0[A == 0], 
                                   X = W[A == 0, ], 
                                   family = gaussian(), 
                                   SL.library = SL_QY_W,
                                   # method = tmp_method.CC_LS(),
                                   control = list(saveCVFitLibrary = TRUE))
    QY_Wn_A1_A0 <- as.numeric(
      stats::predict(QY_W_fit_A1_A0, newdata = W)[[1]]
    )
    QY_Wn_A0_A1[QY_Wn_A0_A1 <= 0] <- 0
    QY_Wn_A0_A1[QY_Wn_A0_A1 >= 1] <- 1

    QY_Wn_A1_A0_cv <- partial_cv_preds_QY_Wn(QY_W_fit_A1_A0, a = 0, newdata = W, A = A, family = stats::gaussian())
    QY_Wn_A1_A0_cv[QY_Wn_A1_A0_cv <= 0] <- 0
    QY_Wn_A1_A0_cv[QY_Wn_A1_A0_cv >= 1] <- 1
    
    QY_Wn_A0_A1_cv <- partial_cv_preds_QY_Wn(QY_W_fit_A0_A1, a = 1, newdata = W, A = A, family = stats::gaussian())
    QY_Wn_A0_A1_cv[QY_Wn_A0_A1_cv <= 0] <- 0
    QY_Wn_A0_A1_cv[QY_Wn_A0_A1_cv >= 1] <- 1
  }


  # for the lazy evaluation of projection
  # E_X[QY_WAS | A = a, W] using inverse weights
  if(!is.null(glm_QY_W)){
    QY_W_fit_A0_A1_lazy <- stats::glm(paste0("QY_WASn_A0 ~ ", glm_QY_W), family = binomial(),
                          data = data.frame(QY_WASn_A0 = QY_WASn_A0, W, wt = R / gRn_1)[A == 1 & R == 1, ],
                          weights = wt)
    QY_Wn_A0_A1_lazy <- stats::predict(QY_W_fit_A0_A1_lazy, type = "response", 
                              newdata = data.frame(QY_WASn_A0 = QY_WASn_A0, W, wt = R / gRn_1))
    QY_Wn_A0_A1_lazy_cv <- QY_Wn_A0_A1_lazy

    QY_W_fit_A1_A0_lazy <- stats::glm(paste0("QY_WASn_A1 ~ ", glm_QY_W), family = binomial(),
                          data = data.frame(QY_WASn_A1 = QY_WASn_A1, W, wt = R / gRn_1)[A == 0 & R == 1, ],
                          weights = wt)
    QY_Wn_A1_A0_lazy <- stats::predict(QY_W_fit_A1_A0_lazy, type = "response", 
                           newdata = data.frame(QY_WASn_A1 = QY_WASn_A1, W, wt = R / gRn_1))
    QY_Wn_A1_A0_lazy_cv <- QY_Wn_A1_A0_lazy
  }else{
    set.seed(seed)
    if(length(unique(QY_WASn_A0[A == 1 & R == 1])) < 3){
      SL_QY_W <- "SL.mean"
    }
    QY_W_fit_A0_A1_lazy <- SuperLearner::SuperLearner(Y = QY_WASn_A0[A == 1 & R == 1], 
                                   X = W[A == 1 & R == 1, ], 
                                   obsWeights = (R / gRn_1)[A == 1 & R == 1],
                                   family = gaussian(), 
                                   SL.library = SL_QY_W,
                                   # method = tmp_method.CC_LS(),
                                   control = list(saveCVFitLibrary = TRUE))
    QY_Wn_A0_A1_lazy <- as.numeric(
      stats::predict(QY_W_fit_A0_A1_lazy, newdata = W)[[1]]
    )
    QY_Wn_A0_A1_lazy[QY_Wn_A0_A1_lazy <= 0] <- 0
    QY_Wn_A0_A1_lazy[QY_Wn_A0_A1_lazy >= 1] <- 1

    set.seed(seed)
    if(length(unique(QY_WASn_A1[A == 0 & R == 1])) < 3){
      SL_QY_W <- "SL.mean"
    }
    QY_W_fit_A1_A0_lazy <- SuperLearner::SuperLearner(Y = QY_WASn_A1[A == 0 & R == 1], 
                                   X = W[A == 0 & R == 1, ], 
                                   obsWeights = (R / gRn_1)[A == 0 & R == 1],
                                   family = gaussian(), 
                                   SL.library = SL_QY_W,
                                   # method = tmp_method.CC_LS(),
                                   control = list(saveCVFitLibrary = TRUE))
    QY_Wn_A1_A0_lazy <- as.numeric(
      stats::predict(QY_W_fit_A1_A0_lazy, newdata = W)[[1]]
    )
    QY_Wn_A0_A1_lazy[QY_Wn_A0_A1_lazy <= 0] <- 0
    QY_Wn_A0_A1_lazy[QY_Wn_A0_A1_lazy >= 1] <- 1

    QY_Wn_A1_A0_lazy_cv <- partial_cv_preds_QY_Wn_lazy(QY_W_fit_A1_A0_lazy, a = 0, A = A, newdata = W, family = stats::gaussian())
    QY_Wn_A1_A0_lazy_cv[QY_Wn_A1_A0_lazy_cv <= 0] <- 0
    QY_Wn_A1_A0_lazy_cv[QY_Wn_A1_A0_lazy_cv >= 1] <- 1
    
    QY_Wn_A0_A1_lazy_cv <- partial_cv_preds_QY_Wn_lazy(QY_W_fit_A0_A1_lazy, a = 1, A = A, newdata = W, family = stats::gaussian())
    QY_Wn_A0_A1_lazy_cv[QY_Wn_A0_A1_lazy_cv <= 0] <- 0
    QY_Wn_A0_A1_lazy_cv[QY_Wn_A0_A1_lazy_cv >= 1] <- 1
  }

  # E[ D(P_X)(O) | R = 1, W, A, C, CY]
  # need eif for the outcome
  D_A1_A0 <- make_full_data_eif(a1 = 1, a2 = 0, A = A, C = C, 
                                gA = gAn_1, gC = gCn_1_A1, 
                                gAS = gASn_1, Y = Y, QY_WAS = QY_WASn_A1, 
                                QY_W = QY_Wn_A1_A0)
  D_A0_A1 <- make_full_data_eif(a1 = 0, a2 = 1, A = A, C = C, 
                                gA = 1 - gAn_1, gC = gCn_1_A1,
                                gAS = 1 - gASn_1, Y = Y, QY_WAS = QY_WASn_A0, 
                                QY_W = QY_Wn_A0_A1)


  if(!is.null(glm_QD_WACY_lazy)){
    QD_WACY_fit_A1_A0 <- stats::glm(paste0("D_A1_A0 ~ ", glm_QD_WACY_lazy), 
                                    data = data.frame(D_A1_A0 = D_A1_A0, A = A, W, 
                                                      CY11 = CY11, CY10 = CY10)[R == 1, ])
    QD_WACYn_A1_A0 <- stats::predict(QD_WACY_fit_A1_A0, type = "response", 
                                     newdata = data.frame(D_A1_A0 = D_A1_A0, 
                                                          A = A, W, CY11 = CY11, 
                                                          CY10 = CY10))
    QD_WACYn_A1_A0_cv <- QD_WACYn_A1_A0

    QD_WACY_fit_A0_A1 <- stats::glm(paste0("D_A0_A1 ~ ", glm_QD_WACY_lazy),
                                    data = data.frame(D_A0_A1 = D_A0_A1, A = A, W, 
                                                      CY11 = CY11, CY10 = CY10)[R == 1, ])
    QD_WACYn_A0_A1 <- stats::predict(QD_WACY_fit_A0_A1, type = "response", 
                                     newdata = data.frame(D_A0_A1 = D_A0_A1, 
                                                          A = A, W, CY11 = CY11, 
                                                          CY10 = CY10))
    QD_WACYn_A0_A1_cv <- QD_WACYn_A0_A1
  }else{
    set.seed(seed)
    if(length(unique(D_A1_A0[R == 1])) < 3){
      SL_QD_WACY_lazy <- "SL.mean"
    }
    QD_WACY_fit_A1_A0 <- SuperLearner::SuperLearner(Y = D_A1_A0[R == 1], 
                                   X = data.frame(W, A = A, CY11 = CY11, CY10 = CY10)[R == 1, ], 
                                   family = gaussian(), 
                                   SL.library = SL_QD_WACY_lazy,
                                   # method = tmp_method.CC_LS(),
                                   control = list(saveCVFitLibrary = TRUE))
    QD_WACYn_A1_A0 <- as.numeric(
      stats::predict(QD_WACY_fit_A1_A0, newdata = data.frame(W, A = A, CY11 = CY11, CY10 = CY10))[[1]]
    )

    set.seed(seed)
    if(length(unique(D_A0_A1[R == 1])) < 3){
      SL_QD_WACY_lazy <- "SL.mean"
    }
    QD_WACY_fit_A0_A1 <- SuperLearner::SuperLearner(Y = D_A0_A1[R == 1], 
                                   X = data.frame(W, A = A, CY11 = CY11, CY10 = CY10)[R == 1, ], 
                                   family = gaussian(), 
                                   SL.library = SL_QD_WACY_lazy,
                                   # method = tmp_method.CC_LS(),
                                   control = list(saveCVFitLibrary = TRUE))
    QD_WACYn_A0_A1 <- as.numeric(
      stats::predict(QD_WACY_fit_A0_A1, newdata = data.frame(W, A = A, CY11 = CY11, CY10 = CY10))[[1]]
    )

    # get partially cross-validated predictions
    QD_WACYn_A1_A0_cv <- partial_cv_preds_QD_WACYn_lazy(QD_WACY_fit_A1_A0, newdata = data.frame(W, A = A, CY11 = CY11, CY10 = CY10), R = R, family = stats::gaussian())
    QD_WACYn_A0_A1_cv <- partial_cv_preds_QD_WACYn_lazy(QD_WACY_fit_A0_A1, newdata = data.frame(W, A = A, CY11 = CY11, CY10 = CY10), R = R, family = stats::gaussian())
  }

  # E[Y(1, S(1))] and E[Y(0, S(0))]
  if(!is.null(glm_QY_WA)){
    QY_WA_fit <- stats::glm(paste0("Y ~ ", glm_QY_WA), family = binomial(),
                      data = data.frame(Y = Y, A = A, W)[C == 1,])
    QY_WAn_A1 <- rep(NA, n)
    QY_WAn_A0 <- rep(NA, n)
    QY_WAn_A1 <- stats::predict(QY_WA_fit, type = "response", 
                         newdata = data.frame(Y = Y, A = 1, W))
    QY_WAn_A0 <- stats::predict(QY_WA_fit, type = "response", 
                         newdata = data.frame(Y = Y, A = 0, W)) 
    QY_WAn_A1_cv <- QY_WAn_A1
    QY_WAn_A0_cv <- QY_WAn_A0
  }else{
    set.seed(seed)
    QY_WA_fit <- SuperLearner::SuperLearner(Y = Y[C == 1], 
                              X = data.frame(A = A, W)[C == 1, ],
                              family = binomial(), 
                              SL.library = SL_QY_WA,
                              # method = tmp_method.CC_nloglik(),
                              control = list(saveCVFitLibrary = TRUE))
    QY_WAn_A1 <- as.numeric(stats::predict(QY_WA_fit, newdata = data.frame(A = 1, W))[[1]])
    QY_WAn_A0 <- as.numeric(stats::predict(QY_WA_fit, newdata = data.frame(A = 0, W))[[1]])
    
    # get partially cross-validated predictions
    QY_WAn_A0_cv <- partial_cv_preds_QY_WAn(QY_WA_fit, newdata = data.frame(A = 0, W), C = C, family = stats::binomial())
    QY_WAn_A1_cv <- partial_cv_preds_QY_WAn(QY_WA_fit, newdata = data.frame(A = 1, W), C = C, family = stats::binomial())
  }


  eif_psi10 <- make_eif_ya1_sa2(a1 = 1, a2 = 0, R = R, 
                                gR = gRn_1, A = A, C = C, 
                                gA = gAn_1, gC = gCn_1_A1, 
                                gAS = gASn_1, Y = Y, 
                                QY_WAS = QY_WASn_A1, 
                                QD_WACY = QD_WACYn_A1, 
                                QY_WACY = QY_WACYn_A1_A0, 
                                QY_W = QY_Wn_A1_A0)
  eif_psi10_cv <- make_eif_ya1_sa2(a1 = 1, a2 = 0, R = R, 
                                   gR = gRn_1, A = A, C = C, 
                                   gA = gAn_1, gC = gCn_1_A1_cv, 
                                   gAS = gASn_1_cv, Y = Y, 
                                   QY_WAS = QY_WASn_A1_cv, 
                                   QD_WACY = QD_WACYn_A1_cv, 
                                   QY_WACY = QY_WACYn_A1_A0_cv, 
                                   QY_W = QY_Wn_A1_A0_cv)

  eif_psi01 <- make_eif_ya1_sa2(a1 = 0, a2 = 1, R = R, 
                                gR = gRn_1, A = A, C = C, 
                                gA = 1 - gAn_1, gC = gCn_1_A0, 
                                gAS = 1 - gASn_1, Y = Y, 
                                QY_WAS = QY_WASn_A0, 
                                QD_WACY = QD_WACYn_A0, 
                                QY_WACY = QY_WACYn_A0_A1, 
                                QY_W = QY_Wn_A0_A1)
  eif_psi01_cv <- make_eif_ya1_sa2(a1 = 0, a2 = 1, R = R, 
                                   gR = gRn_1, A = A, C = C, 
                                   gA = 1 - gAn_1, gC = gCn_1_A0, 
                                   gAS = 1 - gASn_1_cv, Y = Y, 
                                   QY_WAS = QY_WASn_A0_cv, 
                                   QD_WACY = QD_WACYn_A0_cv, 
                                   QY_WACY = QY_WACYn_A0_A1_cv, 
                                   QY_W = QY_Wn_A0_A1_cv)

  psi10n <- mean(QY_Wn_A1_A0) + mean(eif_psi10)
  if(psi10n <= 0){
    psi10n <- mean(QY_Wn_A1_A0)
  }
  psi01n <- mean(QY_Wn_A0_A1) + mean(eif_psi01)
  if(psi01n <= 0){
    psi01n <- mean(QY_Wn_A0_A1)
  }
  

  eif_psi10_lazy <- make_eif_ya1_sa2_lazy(a1 = 1, a2 = 0, R = R, 
                                          gR = gRn_1, A = A, C = C, 
                                          gA = gAn_1, gC = gCn_1_A1, 
                                          gAS = gASn_1, Y = Y, 
                                          QY_WAS = QY_WASn_A1, 
                                          QD_WACY = QD_WACYn_A1_A0, 
                                          QY_W = QY_Wn_A1_A0_lazy)

  eif_psi10_lazy_cv <- make_eif_ya1_sa2_lazy(a1 = 1, a2 = 0, R = R, 
                                          gR = gRn_1, A = A, C = C, 
                                          gA = gAn_1, gC = gCn_1_A1, 
                                          gAS = gASn_1_cv, Y = Y, 
                                          QY_WAS = QY_WASn_A1_cv, 
                                          QD_WACY = QD_WACYn_A1_A0_cv, 
                                          QY_W = QY_Wn_A1_A0_lazy_cv)

  eif_psi01_lazy <- make_eif_ya1_sa2_lazy(a1 = 0, a2 = 1, R = R, 
                                          gR = gRn_1, A = A, C = C, 
                                          gA = 1 - gAn_1, gC = gCn_1_A1, 
                                          gAS = 1 - gASn_1, Y = Y, 
                                          QY_WAS = QY_WASn_A0, 
                                          QD_WACY = QD_WACYn_A0_A1, 
                                          QY_W = QY_Wn_A0_A1_lazy)

  eif_psi01_lazy_cv <- make_eif_ya1_sa2_lazy(a1 = 0, a2 = 1, R = R, 
                                          gR = gRn_1, A = A, C = C, 
                                          gA = 1 - gAn_1, gC = gCn_1_A1, 
                                          gAS = 1 - gASn_1_cv, Y = Y, 
                                          QY_WAS = QY_WASn_A0_cv, 
                                          QD_WACY = QD_WACYn_A0_A1_cv, 
                                          QY_W = QY_Wn_A0_A1_lazy_cv)

  psi10n_lazy <- mean(QY_Wn_A1_A0_lazy) + mean(eif_psi10_lazy)
  psi01n_lazy <- mean(QY_Wn_A0_A1_lazy) + mean(eif_psi01_lazy)
  

  eif_psi11 <- make_eif_ya_sa(a = 1, A = A, C = C, gA = gAn_1, gC = gCn_1_A1, Y = Y, QY_WA = QY_WAn_A1)
  eif_psi00 <- make_eif_ya_sa(a = 0, A = A, C = C, gA = 1 - gAn_1, gC = gCn_1_A0, Y = Y, QY_WA = QY_WAn_A0)
  eif_psi11_cv <- make_eif_ya_sa(a = 1, A = A, C = C, gA = gAn_1, gC = gCn_1_A1_cv, Y = Y, QY_WA = QY_WAn_A1_cv)
  eif_psi00_cv <- make_eif_ya_sa(a = 0, A = A, C = C, gA = 1 - gAn_1, gC = gCn_1_A0_cv, Y = Y, QY_WA = QY_WAn_A0_cv)
  
  psi11n <- mean(QY_WAn_A1) + mean(eif_psi11)
  if(psi11n <= 0){
    psi11n <- mean(QY_WAn_A1)
  }
  psi00n <- mean(QY_WAn_A0) + mean(eif_psi00)
  if(psi00n <= 0){
    psi00n <- mean(QY_WAn_A0)
  }

  eif_matrix <- cbind(eif_psi11, eif_psi00, eif_psi10, eif_psi01)
  cov_matrix <- cov(eif_matrix) / n

  eif_matrix_cv <- cbind(eif_psi11_cv, eif_psi00_cv, eif_psi10_cv, eif_psi01_cv)
  cov_matrix_cv <- cov(eif_matrix_cv) / n

  eif_matrix_lazy <- cbind(eif_psi11, eif_psi00, eif_psi10_lazy, eif_psi01_lazy)
  cov_matrix_lazy <- cov(eif_matrix_lazy) / n

  eif_matrix_lazy_cv <- cbind(eif_psi11_cv, eif_psi00_cv, eif_psi10_lazy_cv, eif_psi01_lazy_cv)
  cov_matrix_lazy_cv <- cov(eif_matrix_lazy_cv) / n

  risks_and_effects <- get_risks_and_effects(psi11n = psi11n, psi00n = psi00n,
                                             psi10n = psi10n, psi01n = psi01n,
                                             cov_matrix = cov_matrix,
                                             cov_matrix_cv = cov_matrix_cv)
  risks_and_effects_lazy <- get_risks_and_effects(psi11n = psi11n, psi00n = psi00n,
                                                  psi10n = psi10n_lazy, psi01n = psi01n_lazy,
                                                  cov_matrix = cov_matrix_lazy,
                                                  cov_matrix_cv = cov_matrix_lazy_cv)
  
  # TO DO: Add fitted nuisance models to output?
  out <- list(risk = risks_and_effects$risk, 
              eff = risks_and_effects$eff, 
              eff2 = risks_and_effects$eff2, 
              cov = cov_matrix, cov_cv = cov_matrix_cv,
              risk_lazy = risks_and_effects_lazy$risk, 
              eff_lazy = risks_and_effects_lazy$eff, 
              eff2_lazy = risks_and_effects_lazy$eff2, 
              cov_lazy = cov_matrix_lazy, cov_lazy_cv = cov_matrix_lazy_cv)
  return(out)
}
