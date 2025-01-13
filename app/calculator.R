
prob_rule_out2 <- function(
    ref_tx = 0.06,
    ref_cx = 0.03,
    fold_increase = 3,
    ped_n_cx = 100,
    ped_n_tx = 100,
    prior_tx_a = 1,
    prior_tx_b = 1,
    prior_cx_a = 1,
    prior_cx_b = 1,
    ped_cx_event = 0,
    ped_cx_n = 0,
    ped_tx_event = 0,
    ped_tx_n = 0) {
  num_ouputs <- 6
  m <- fold_increase # for length of interval of treatment responses
  trt_nseq <- ped_n_tx
  pbo_nseq <- ped_n_cx
  pbo_pcseq_low <- ref_cx # kept at 0.03 for common ADR; usually 0.02 to 0.05
  trt_pcseq_low <- ref_tx
  trt_pcseq <- seq(trt_pcseq_low, pbo_pcseq_low + m * (trt_pcseq_low - pbo_pcseq_low), length = num_ouputs)
  pbo_pcseq <- rep(pbo_pcseq_low, length(trt_pcseq))
  delta_p <- trt_pcseq - pbo_pcseq
  delta_p_low <- trt_pcseq_low - pbo_pcseq_low
  delta_pmult <- delta_p / delta_p_low
  # browser()


  alpha1 <- prior_tx_a + ped_tx_event
  beta1 <- prior_tx_b + (ped_tx_n - ped_tx_event)
  alpha2 <- prior_cx_a + ped_cx_event
  beta2 <- prior_cx_b + (ped_cx_n - ped_cx_event)



  MatrixA <- matrix(nrow = num_ouputs, ncol = 3)

  for (j in seq_len(num_ouputs)) {
    a1 <- trt_nseq * trt_pcseq[j] + alpha1
    b1 <- trt_nseq - trt_nseq * trt_pcseq[j] + beta1
    a2 <- pbo_nseq * pbo_pcseq[j] + alpha2
    b2 <- pbo_nseq - pbo_nseq * pbo_pcseq[j] + beta2
    x <- seq(delta_p[j], 1, 0.001) # interval for integration
    y <- sapply(x, diff_beta_plus, a1 = a1, a2 = a2, b1 = b1, b2 = b2)
    xy <- cbind(x, y)
    xy.narm <- xy[complete.cases(xy), ]
    Pr.less.x <- 1 - pracma::trapz(xy.narm[, 1], xy.narm[, 2])
    MatrixA[j, ] <- c(delta_p[j], delta_pmult[j], Pr.less.x)
  }


  Safety.mat <- data.frame(MatrixA)
  colnames(Safety.mat) <- c("Incidence Difference", "Incidence Difference Factor", "Probability")
  return(Safety.mat)
}

diff_beta_plus <- function(a1, a2, b1, b2, x) {
  A <- beta(a1, b1) * beta(a2, b2)
  p1.density <- try(beta(a2, b1) * x^(b1 + b2 - 1) * (1 - x)^(a2 + b1 - 1) * F1(b1, a1 + a2 + b1 + b2 - 2, 1 - a1, b1 + a2, 1 - x, 1 - x^2) / A)
  if (class(p1.density) != "try-error") {
    p1.density
  } else {
    p1.density <- "NA"
  }
}

diff_beta_plus_log <- function(a1, a2, b1, b2, x) {
  logA <- lbeta(a1, b1) + lbeta(a2, b2)
  lbeta(a2, b1) +
    (b1 + b2 - 1) * log(x) +
    (a2 + b1 - 1) * log(1 - x) +
    log(F1(b1, a1 + a2 + b1 + b2 - 2, 1 - a1, b1 + a2, 1 - x, 1 - x^2)) -
    logA
}

F1 <- function(a, b, b.prime, c, x, y, ...) {
  A1.f1 <- try({
    integrate(
      f = A1.simple,
      lower = 0,
      upper = 1,
      a = a,
      b = b,
      b.prime = b.prime,
      c = c,
      x = x,
      y = y,
      subdivisions = 5000,
      rel.tol = .Machine$double.eps^0.5,
      stop.on.error = FALSE, ...)$value
  },
  silent = TRUE
  )
  if (class(A1.f1) != "try-error") {
    A1.f2 <- gamma(c) / (gamma(a) * gamma(c - a)) * A1.f1
  } else {
    A1.f2 <- NA
  }
}


A1.simple <- function(u, a, b, b.prime, c, x, y) {
  u^(a - 1) * (1 - u)^(c - a - 1) * (1 - u * x)^(-b) * (1 - u * y)^(-b.prime)
}
