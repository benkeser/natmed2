# make newdata data.frame(A = a, W, S)
partial_cv_preds_QY_WAn <- function(fit_sl, newdata, C, family){
  n <- length(C)
  n_algo <- length(fit_sl$cvRisk)
  n_folds <- length(fit_sl$validRows)
  alpha_hat <- matrix(fit_sl$coef, nrow = n_algo)
  rslt_list <- vector(mode = "list", length = n_folds)
  for(v in seq_len(n_folds)){
    # who's in the validation fold of the included folks
    foldv_ids <- fit_sl$validRows[[v]]
    # these are the people who were not used to fit these models
    foldv_models <- fit_sl$cvFitLibrary[[v]]
    rslt <- matrix(NA, nrow = length(foldv_ids), ncol = n_algo)
    for(k in seq_len(n_algo)){
      rslt[ , k] <- stats::predict(foldv_models[[k]], 
                                   newdata = newdata[C == 1,][foldv_ids,],
                                   family = family)
    }
    rslt_list[[v]] <- rslt
  }
  # combine using weights from full super learner
  pred_list <- vector(mode = "list", length = n_folds)
  for(v in seq_len(n_folds)){
    pred_list[[v]] <- rslt_list[[v]] %*% alpha_hat
  }
  reorder_preds <- rep(NA, n)
  # fill in observations in regression with cross-validated prediction
  reorder_preds[C == 1][unlist(fit_sl$validRows)] <- unlist(pred_list, use.names = FALSE)
  # all others fill in with prediction from super learner 
  if(any(C == 0)){
    reorder_preds[C == 0] <- stats::predict(fit_sl, newdata = newdata[C == 0, , drop = FALSE])[[1]]
  }
  
  return(reorder_preds)
}

# make newdata data.frame(A = a, W, S)
partial_cv_preds_QY_Wn_lazy <- function(fit_sl, newdata, a, A, family){
  n <- length(A)
  n_algo <- length(fit_sl$cvRisk)
  n_folds <- length(fit_sl$validRows)
  alpha_hat <- matrix(fit_sl$coef, nrow = n_algo)
  rslt_list <- vector(mode = "list", length = n_folds)
  for(v in seq_len(n_folds)){
    # who's in the validation fold of the included folks
    foldv_ids <- fit_sl$validRows[[v]]
    # these are the people who were not used to fit these models
    foldv_models <- fit_sl$cvFitLibrary[[v]]
    rslt <- matrix(NA, nrow = length(foldv_ids), ncol = n_algo)
    for(k in seq_len(n_algo)){
      rslt[ , k] <- stats::predict(foldv_models[[k]], 
                                   newdata = newdata[A == a, , drop = FALSE][foldv_ids,],
                                   family = family)
    }
    rslt_list[[v]] <- rslt
  }
  # combine using weights from full super learner
  pred_list <- vector(mode = "list", length = n_folds)
  for(v in seq_len(n_folds)){
    pred_list[[v]] <- rslt_list[[v]] %*% alpha_hat
  }
  reorder_preds <- rep(NA, n)
  # fill in observations in regression with cross-validated prediction
  reorder_preds[A == a][unlist(fit_sl$validRows)] <- unlist(pred_list, use.names = FALSE)
  # all others fill in with prediction from super learner 
  reorder_preds[A != a] <- stats::predict(fit_sl, newdata = newdata[A != a, ])[[1]]
  
  return(reorder_preds)
}

# make newdata data.frame(A = a, W, S)
partial_cv_preds_QY_Wn <- function(fit_sl, newdata, a, A, family){
  n <- length(A)
  n_algo <- length(fit_sl$cvRisk)
  n_folds <- length(fit_sl$validRows)
  alpha_hat <- matrix(fit_sl$coef, nrow = n_algo)
  rslt_list <- vector(mode = "list", length = n_folds)
  for(v in seq_len(n_folds)){
    # who's in the validation fold of the included folks
    foldv_ids <- fit_sl$validRows[[v]]
    # these are the people who were not used to fit these models
    foldv_models <- fit_sl$cvFitLibrary[[v]]
    rslt <- matrix(NA, nrow = length(foldv_ids), ncol = n_algo)
    for(k in seq_len(n_algo)){
      rslt[ , k] <- stats::predict(foldv_models[[k]], newdata = newdata[A == a, , drop = FALSE][foldv_ids,],
                                   family = family)
    }
    rslt_list[[v]] <- rslt
  }
  # combine using weights from full super learner
  pred_list <- vector(mode = "list", length = n_folds)
  for(v in seq_len(n_folds)){
    pred_list[[v]] <- rslt_list[[v]] %*% alpha_hat
  }
  reorder_preds <- rep(NA, n)
  # fill in observations in regression with cross-validated prediction
  reorder_preds[A == a][unlist(fit_sl$validRows)] <- unlist(pred_list, use.names = FALSE)
  # all others fill in with prediction from super learner 
  reorder_preds[A != a] <- stats::predict(fit_sl, newdata = newdata[A != a, ])[[1]]
  
  return(reorder_preds)
}



# make newdata data.frame(A = a, W, S)
partial_cv_preds_QY_WACYn <- function(fit_sl, newdata, a, A, R, family){
  n <- length(R)
  n_algo <- length(fit_sl$cvRisk)
  n_folds <- length(fit_sl$validRows)
  alpha_hat <- matrix(fit_sl$coef, nrow = n_algo)
  rslt_list <- vector(mode = "list", length = n_folds)
  for(v in seq_len(n_folds)){
    # who's in the validation fold of the included folks
    foldv_ids <- fit_sl$validRows[[v]]
    # these are the people who were not used to fit these models
    foldv_models <- fit_sl$cvFitLibrary[[v]]
    rslt <- matrix(NA, nrow = length(foldv_ids), ncol = n_algo)
    for(k in seq_len(n_algo)){
      rslt[ , k] <- stats::predict(foldv_models[[k]], 
                                   newdata = newdata[A == a & R == 1,][foldv_ids,],
                                   family = family)
    }
    rslt_list[[v]] <- rslt
  }
  # combine using weights from full super learner
  pred_list <- vector(mode = "list", length = n_folds)
  for(v in seq_len(n_folds)){
    pred_list[[v]] <- rslt_list[[v]] %*% alpha_hat
  }
  reorder_preds <- rep(NA, n)
  # fill in observations in regression with cross-validated prediction
  reorder_preds[A == a & R == 1][unlist(fit_sl$validRows)] <- unlist(pred_list, use.names = FALSE)
  # all others fill in with prediction from super learner 
  reorder_preds[A != a | R == 0] <- stats::predict(fit_sl, newdata = newdata[A != a | R == 0, ])[[1]]
  
  return(reorder_preds)
}

  


# make newdata data.frame(A = a, W, S)
partial_cv_preds_QD_WACYn <- function(fit_sl, newdata, a, A, R, C,
                                      family){
  n <- length(C)
  n_algo <- length(fit_sl$cvRisk)
  n_folds <- length(fit_sl$validRows)
  alpha_hat <- matrix(fit_sl$coef, nrow = n_algo)
  rslt_list <- vector(mode = "list", length = n_folds)
  for(v in seq_len(n_folds)){
    # who's in the validation fold of the included folks
    foldv_ids <- fit_sl$validRows[[v]]
    # these are the people who were not used to fit these models
    foldv_models <- fit_sl$cvFitLibrary[[v]]
    rslt <- matrix(NA, nrow = length(foldv_ids), ncol = n_algo)
    for(k in seq_len(n_algo)){
      rslt[ , k] <- stats::predict(foldv_models[[k]], 
                                   newdata = newdata[A == a & C == 1 & R == 1,][foldv_ids,],
                                   family = family)
    }
    rslt_list[[v]] <- rslt
  }
  # combine using weights from full super learner
  pred_list <- vector(mode = "list", length = n_folds)
  for(v in seq_len(n_folds)){
    pred_list[[v]] <- rslt_list[[v]] %*% alpha_hat
  }
  reorder_preds <- rep(0, n)
  # fill in observations in regression with cross-validated prediction
  reorder_preds[A == a & C == 1 & R == 1][unlist(fit_sl$validRows)] <- unlist(pred_list, use.names = FALSE)
  # all others fill in with prediction from super learner 
  if(sum(A == a & C == 1 & R == 0) > 1){
    reorder_preds[A == a & C == 1 & R == 0] <- stats::predict(fit_sl, newdata = newdata[A == a & C == 1 & R == 0, ])[[1]]
  }
  return(reorder_preds)
}

# make newdata data.frame(A = a, W, S)
partial_cv_preds_QD_WACYn_lazy <- function(fit_sl, newdata, R, family){
  n <- length(R)
  n_algo <- length(fit_sl$cvRisk)
  n_folds <- length(fit_sl$validRows)
  alpha_hat <- matrix(fit_sl$coef, nrow = n_algo)
  rslt_list <- vector(mode = "list", length = n_folds)
  for(v in seq_len(n_folds)){
    # who's in the validation fold of the included folks
    foldv_ids <- fit_sl$validRows[[v]]
    # these are the people who were not used to fit these models
    foldv_models <- fit_sl$cvFitLibrary[[v]]
    rslt <- matrix(NA, nrow = length(foldv_ids), ncol = n_algo)
    for(k in seq_len(n_algo)){
      rslt[ , k] <- stats::predict(foldv_models[[k]], 
                                   newdata = newdata[R == 1,][foldv_ids,],
                                   family = family)
    }
    rslt_list[[v]] <- rslt
  }
  # combine using weights from full super learner
  pred_list <- vector(mode = "list", length = n_folds)
  for(v in seq_len(n_folds)){
    pred_list[[v]] <- rslt_list[[v]] %*% alpha_hat
  }
  reorder_preds <- rep(0, n)
  # fill in observations in regression with cross-validated prediction
  reorder_preds[R == 1][unlist(fit_sl$validRows)] <- unlist(pred_list, use.names = FALSE)
  # all others fill in with prediction from super learner 
  if(any(R == 0)){
    reorder_preds[R == 0] <- stats::predict(fit_sl, newdata = newdata[R == 0, ])[[1]]
  }
  return(reorder_preds)
}

# make newdata data.frame(A = a, W, S)
partial_cv_preds_QY_WASn <- function(fit_sl, newdata, R, C, family){
  n <- length(C)
  n_algo <- length(fit_sl$cvRisk)
  n_folds <- length(fit_sl$validRows)
  alpha_hat <- matrix(fit_sl$coef, nrow = n_algo)
  rslt_list <- vector(mode = "list", length = n_folds)
  for(v in seq_len(n_folds)){
    # who's in the validation fold of the included folks
    foldv_ids <- fit_sl$validRows[[v]]
    # these are the people who were not used to fit these models
    foldv_models <- fit_sl$cvFitLibrary[[v]]
    rslt <- matrix(NA, nrow = length(foldv_ids), ncol = n_algo)
    for(k in seq_len(n_algo)){
      rslt[ , k] <- stats::predict(foldv_models[[k]], 
                                   newdata = newdata[R == 1 & C == 1,][foldv_ids,],
                                   family = family)
    }
    rslt_list[[v]] <- rslt
  }
  # combine using weights from full super learner
  pred_list <- vector(mode = "list", length = n_folds)
  for(v in seq_len(n_folds)){
    pred_list[[v]] <- rslt_list[[v]] %*% alpha_hat
  }
  reorder_preds <- rep(NA, n)
  # fill in observations in regression with cross-validated prediction
  reorder_preds[R == 1 & C == 1][unlist(fit_sl$validRows)] <- unlist(pred_list, use.names = FALSE)
  # all others fill in with prediction from super learner 
  reorder_preds[R == 1] <- stats::predict(fit_sl, newdata = newdata[R == 1,,drop = FALSE])[[1]]
  
  return(reorder_preds)
}




#' Helper function to properly format partially cross-validated predictions 
#' from a fitted super learner.
#' 
#' @param fit_sl A fitted \code{SuperLearner} object with 
#' \code{control$saveCVFitLibrary = TRUE}
#' @param a_0 Treatment level to set. If \code{NULL}, assume this function
#' is being used to get partially cross-validated propensity score predictions.
#' @param W A \code{data.frame} of named covariates.
#' @param include A boolean vector indicating which observations were actually
#' used to fit the regression. 
#' @param easy A boolean indicating whether the predictions can be 
#' computed the "easy" way, i.e., based just on the Z matrix from SuperLearner.
#' This is possible for propensity score models when no missing data AND no 
#' stratification. 
partial_cv_preds <- function(fit_sl, a_0, W = NULL, subset = TRUE, 
                             include = NULL, easy = FALSE, family){
  n_algo <- length(fit_sl$cvRisk)
  n_folds <- length(fit_sl$validRows)

  if(!easy & all(subset)){ 
    n <- length(W[,1])
  }else if(!easy & !all(subset)){
    n <- sum(subset)
  }else{ # if used in easy scenario, then fit_sl will have been 
         # fit using all observations and no W will enter, so we'll check
         # the Z matrix that is used to generate predictions below
    n <- length(fit_sl$Z[,1])
  }
  alpha_hat <- matrix(fit_sl$coef, nrow = n_algo)
  if(!easy){
    # rslt_list will eventually hold cross-validated predictions for 
    # all observations that were actually used to fit the regression
    rslt_list <- vector(mode = "list", length = n_folds)
    for(v in seq_len(n_folds)){
      # who's in the validation fold of the included folks
      foldv_ids <- fit_sl$validRows[[v]]
      # these are the people who were not used to fit these models
      foldv_models <- fit_sl$cvFitLibrary[[v]]
      rslt <- matrix(NA, nrow = length(foldv_ids), ncol = n_algo)
      for(k in seq_len(n_algo)){
        # predict under a_0
        if(!is.null(a_0)){
          rslt[ , k] <- predict(foldv_models[[k]], 
                                newdata = data.frame(A = a_0, W)[include,][foldv_ids,],
                                family = family)
        }else{
          rslt[ , k] <- predict(foldv_models[[k]], 
                                newdata = W[include,][foldv_ids,],
                                family = family)
        }
      }
      rslt_list[[v]] <- rslt
    }
    # combine using weights from full super learner
    pred_list <- vector(mode = "list", length = n_folds)
    for(v in seq_len(n_folds)){
      pred_list[[v]] <- rslt_list[[v]] %*% alpha_hat
    }
    reorder_preds <- rep(NA, n)
    # fill in observations in regression with cross-validated prediction
    reorder_preds[include][unlist(fit_sl$validRows)] <- unlist(pred_list, use.names = FALSE)
    # all others fill in with prediction from super learner 
    if(any(!include)){
      if(!is.null(a_0)){
        reorder_preds[!include] <- predict(fit_sl, newdata = data.frame(A = a_0, W)[!include,,drop = FALSE])[[1]]
      }else{
        reorder_preds[!include] <- predict(fit_sl, newdata = W[!include,, drop = FALSE])[[1]]
      }
    }
  }else{ # when a_0 is NULL
    # in this case, we're operating on a propensity score model 
    # in which case, fit_sl$Z already has cross-validated predictions
    reorder_preds <- as.numeric(fit_sl$Z %*% alpha_hat)
  }
  return(reorder_preds)
}

