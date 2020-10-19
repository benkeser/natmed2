#' Evaluate EIF of E[Y(a1, S(a2))] on data
#' 
#' @param a1 Either 0 or 1
#' @param a2 Should be \code{1 - a1}
#' @param R A numeric indicator of having mediator measured
#' @param A A numeric vector of binary treatment assignment
#' @param C A numeric indicator of having the outcome \code{Y} measured, 
#' i.e., a measure of being not right-censored. 
#' @param gA Estimate of P(A = a1 | W_i), i = 1,...,n
#' @param gR Estimate of P(R = 1 | W_i, A_i, C_i, C_iY_i)
#' @param gAS Estimate of P(A = a1 | S_i, W_i)
#' @param gC Estimate of P(C = 1 | A = a1, W_i)
#' @param QY_WAS Estimate of E[Y | A = a1, W_i, S_i]
#' @param QD_WACY Estimate of extra nuisance regression
#' @param QY_WACY Estimate of E[QY_WAS | A = a2, W, C, CY]
#' @param QY_W Estimate of E[QY_WACY | A = a2, W]
#' 
#' @return An n-length vector of the estimated EIF evaluated on the observations
make_eif_ya1_sa2 <- function(a1, a2, R, gR, A, C, gA, gC, gAS, Y, QY_WAS,
                             QD_WACY, QY_WACY, QY_W){
  Y <- replace_nas(Y)
  QY_WAS <- replace_nas(QY_WAS)
  gAS <- replace_nas(gAS)
  p1 <- R / gR * as.numeric(A == a1 & C == 1) / (gA * gC) * (1 - gAS) / gAS * (Y - QY_WAS)
  p2 <- QD_WACY / gR * (R - gR)
  p3 <- R / gR * as.numeric(A == a2) / (1 - gA) * (QY_WAS - QY_WACY)
  p4 <- as.numeric(A == a2) / (1 - gA) * (QY_WACY - QY_W)
  p5 <- QY_W - mean(QY_W)
  return(p1 - p2 + p3 + p4 + p5)
}

#' Evaluate full data of E[Y(a1, S(a2))] on data
#' 
#' @param a1 Either 0 or 1
#' @param a2 Should be \code{1 - a1}
#' @param R A numeric indicator of having mediator measured
#' @param A A numeric vector of binary treatment assignment
#' @param C A numeric indicator of having the outcome \code{Y} measured, 
#' i.e., a measure of being not right-censored. 
#' @param gA Estimate of P(A = a1 | W_i), i = 1,...,n
#' @param gR Estimate of P(R = 1 | W_i, A_i, C_i, C_iY_i)
#' @param gAS Estimate of P(A = a1 | S_i, W_i)
#' @param gC Estimate of P(C = 1 | A = a1, W_i)
#' @param QY_WAS Estimate of E[Y | A = a1, W_i, S_i]
#' @param QY_W Estimate of E[QY_WACY | A = a2, W]
#' 
#' @return An n-length vector of the estimated EIF evaluated on the observations
make_full_data_eif <- function(a1, a2, R, gR, A, C, gA, gC, gAS, Y, QY_WAS, QY_W){
  Y <- replace_nas(Y)
  QY_WAS <- replace_nas(QY_WAS)
  gAS <- replace_nas(gAS)
  p1 <- as.numeric(A == a1 & C == 1) / (gA * gC) * (1 - gAS) / gAS * (Y - QY_WAS)
  p2 <- as.numeric(A == a2) / (1 - gA) * (QY_WAS - QY_W)
  p3 <- QY_W - mean(QY_W)
  return(p1 + p2 + p3)
}


#' evaluate EIF of E[Y(a, S(a))] on data
#' @param gA P(A = a | ...)
#' @param gC P(C = 1 | A = a, ...)
#' @param a value a 
#' @param QY_WA E[Y | A = a, W]
make_eif_ya_sa <- function(a, A, C, gA, gC, Y, QY_WA){
  Y <- replace_nas(Y)
  p1 <- as.numeric(A == a & C == 1) / (gA * gC) * (Y - QY_WA)
  p2 <- QY_WA - mean(QY_WA)
  return(p1 + p2)
}