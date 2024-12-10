qprob <- function(corr) {
  0.5 - mets::pmvn(upper = cbind(0, 0), sigma = corr, cor = TRUE)
}

prob.fct <- function(x, alpha, corr) {
  q <- qprob(corr)
  0.5 * pchisq(x, 1, lower.tail = FALSE) +
    q * pchisq(x, 2, lower.tail = FALSE) - alpha
}

q.fct <- function(alpha, corr) {
  uniroot(prob.fct,
          alpha = alpha,
          corr = corr,
          interval = c(0, 10)
  )$root
}

###############################################################
## Calculating test statistics
## Calculating p-values
###############################################################
##' @title Signed intersection Wald test
##' @param thetahat1 parameter estimate 1
##' @param se1 standard error of parameter estimate 1
##' @param thetahat2 parameter estimate 2
##' @param se2 standard error of parameter estimate 2
##' @param noninf1 non-inferiority margin for paramter 1
##' @param noninf2 non-inferiority margin for paramter 2
##' @param corr correlation between parameter 1 and 2
##' @param alpha nominal level
##' @author
##' Christian Bressen Pipper,
##' Klaus Kähler Holst
##' @return list with Wald
##' @author Klaus Kähler Holst
test_intersectsignedwald <- function(thetahat1,
                                     se1,
                                     thetahat2,
                                     se2,
                                     noninf1,
                                     noninf2,
                                     corr,
                                     alpha) {
  z1 <- (thetahat1 - noninf1) / se1
  z2 <- (thetahat2 - noninf2) / se2
  zmin <- min(z1, z2)
  zmax <- max(z1, z2)
  SignWald.intersect <- ifelse(zmax >= 0 & zmin <= (corr * zmax), 1, 0) *
    zmax * zmax + ifelse(zmax >= 0 & zmin > (corr * zmax),
      (zmax * zmax + zmin * zmin - 2 * corr * zmax * zmin) / (1 - corr * corr),
      0
    )
  SignWald1 <- ifelse(z1 >= 0, 1, 0) * z1^2
  SignWald2 <- ifelse(z2 >= 0, 1, 0) * z2^2
  critval.intersect <- q.fct(alpha, corr)
  pval.intersect <- ifelse(SignWald.intersect > 0,
    ## prob.fct(SignWald.intersect, alpha, corr) + alpha, 1
    prob.fct(SignWald.intersect, 0, corr), 1
    )
  pval1 <- ifelse(SignWald1 > 0,
    0.5 * pchisq(SignWald1, 1, lower.tail = FALSE), 1
  )
  pval2 <- ifelse(SignWald2 > 0,
    0.5 * pchisq(SignWald2, 1, lower.tail = FALSE), 1
  )
  test.int <- structure(list(
    data.name = "H1 ^ H2",
    statistic = c("Q" = unname(SignWald.intersect)),
    parameter = NULL,
    method = "Signed Wald Intersection Test",
    ## null.value = "one",
    ## alternative = "one.sided",
    p.value = pval.intersect
    ## estimate = 1
  ), class = "htest")
  test.1 <- structure(list(
    data.name = sprintf("H1: b1 ≤ %g", noninf1),
    statistic = c("Q" = unname(SignWald1)),
    estimate = c("b1" = unname(thetahat1)),
    parameter = NULL,
    method = "Signed Wald Test",
    ## null.value = noninf1,
    alternative = sprintf("HA1: b1 > %g", noninf1),
    p.value = pval1
  ), class = "htest")
  test.2 <- structure(list(
    data.name = sprintf("H2: b2 ≤ %g", noninf2),
    statistic = c("Q" = unname(SignWald2)),
    estimate = c("b2" = unname(thetahat2)),
    parameter = NULL,
    method = "Signed Wald Test",
    alternative = sprintf("HA2: b2 > %g", noninf2),
    p.value = pval2
    ## estimate = 1
  ), class = "htest")

  list(
    critval.intersect = critval.intersect,
    test.intersect = test.int,
    test.1 = test.1,
    test.2 = test.2
  )

}


##' @description
##' Let \eqn{Y} denote the clinical outcome, \eqn{A} the binary treatment
##' variable, \eqn{X} baseline covariates, \(T\) the failure time,
##' and \(epsilon=1,2\) the cause of failure.
##' The following are our two target parameters
##' \deqn{E(Y|T>t, A=1)- E(Y|T>t, A=0)}
##' \deqn{P(T<t,\epsilon=1|A=1)- P(T<t,\epsilon=1|A=0)}
##' @title Estimation of mean clinical outcome truncated by event process
##' @param data data.frame
##' @param ymod model for clinical outcome given T>time
##' @param rmod model for missing data mechnaism for clinical outcome at T=time
##' @param amod treatment model (in RCT should just be 'a ~ 1')
##' @param eventmod Model for time-to-event process ('Event(time,status) ~ x')
##' @param time landmark time
##' @param cause primary event (in the 'status' variable of the 'Event'
##'   statement)
##' @param cens.code censoring code (0 default)
##' @param naive ff TRUE the unadjusted estimates ignoring baseline covariates
##'   is returned as the attribute 'naive'
##' @param ... additional arguments passed to lower level functions
##' @return estimate object
##' @author Klaus Kähler Holst
##' @examples
##' \dontrun{
##' mod1 <- predictor_glm(y ~ a * (x1 + x2))
##' mod2 <- predictor_glm(r ~ a * (x1 + x2), family = binomial)
##' a <- truncatedscore_estimate(
##'   data = dat,
##'   ymod = mod1,
##'   rmod = mod2,
##'   amod = a ~ 1,
##'   eventmod = mets::Event(time, status) ~ a * (x1+x2),
##'   time = 2
##' )
##'
##' s <- summary(a, noninf.t = -0.1)
##' print(s)
##' parameter(s)
##' }
##' @export
truncatedscore_estimate <- function(
                     data,
                     ymod,
                     rmod,
                     amod,
                     eventmod,
                     time,
                     cause = 1,
                     cens.code = 0,
                     naive = FALSE,
                     ...
                     ) {
  if (inherits(ymod, "formula")) {
    ymod <- predictor_glm(ymod)
  }
  if (inherits(rmod, "formula")) {
    rmod <- predictor_glm(rmod, family = binomial)
  }
  if (inherits(amod, "formula")) {
    amod <- predictor_glm(amod, family = binomial)
  }
  # Missing data model
  rmod$estimate(data)
  r <- rmod$response(data) == 1
  # Data with Y observed
  d1 <- data[which(r), , drop = TRUE]
  # Outcome model E(Y|R=1,X,A)
  ymod$estimate(d1)
  # Treatment model
  amod$estimate(data)
  a <- amod$response(data)
  a0 <- amod$response(data, eval = FALSE)
  treatment <- all.vars(update(amod$formula, ~1))
  alev <- c(a0[which(a==0)[1]], a0[which(a==1)[1]])
  y <- ymod$response(data, na.action=lava::na.pass0)
  tmpdata <- data
  est <- ic <-
    lab0 <- est.naive <-
      ic.naive <- c()
  for (aval in c(0, 1)) {
    Ia <- (a == aval)
    tmpdata[, treatment] <- alev[aval + 1]
    q1 <- ymod$predict(tmpdata)
    q2 <- rmod$predict(tmpdata)
    pa <- amod$predict(tmpdata)
    p <- pa[which(a == 1)[1]]
    if (aval == 0) pa <- 1 - pa
    idx1 <- which(Ia)
    pr <- mean(r[idx1]) # P(R=1|A=a)
    m1 <- mean(y[which(Ia & r)]) # E(Y|A=a,R=1)
    est.naive <- c(est.naive, m1)
    ic0 <- r * Ia / (pa * pr) * (y - m1)
    ic.naive <- cbind(ic.naive, ic0)
    ic1 <- ic0 -
      (a - p) * (aval - p) / (p * (1 - p) * pr) * (q1 - m1) * q2
    est1 <- m1 + mean(ic1)
    ic <- cbind(ic, ic1 - mean(ic1))
    est <- c(est, est1)
    lab0 <- c(lab0, sprintf("E(Y|T≥%d,A=%d)", time, aval))
  }
  lab0 <- c(lab0, "diff")
  res <- estimate(
    coef = c(est, est[2] - est[1]),
    IC = cbind(ic, ic[, 2] - ic[, 1]),
    labels = lab0,
    id = rownames(data)
  )
  data[, treatment] <- factor(data[, treatment])
  eventmod <- update(eventmod, as.formula(sprintf("~ %s+ .", treatment)))
  b <- mets::binregATE(eventmod,
    data = data, time = time,
    outcome = "cif", cause = cause, cens.code = cens.code,
    ...
    )
  lab <- paste0(gsub(
    "^treat", sprintf("Risk\\(T<%d|A=", time),
    names(b$riskDR)[1:2]), ")")
  lab <- c(lab, "riskdiff")
  best <- estimate(
    coef = c(b$riskDR, b$riskDR[2] - b$riskDR[1]),
    IC = NROW(b$riskDR.iid) * cbind(
      b$riskDR.iid,
      b$riskDR.iid[, 2] - b$riskDR.iid[, 1]
    ),
    id = rownames(data)
  )
  best <- estimate(best, labels=lab)
  res <- merge(res, best)
  if (naive) {
    eventmod0 <- update(eventmod, as.formula(sprintf("~ %s", treatment)))
    b0 <- mets::binregATE(eventmod0,
      data = data, time = time,
      outcome = "cif", cause = cause, cens.code = cens.code, ...
    )
    best0 <- estimate(
      coef = with(b0, c(riskDR, riskDR[2] - riskDR[1])),
      IC = with(b0, NROW(riskDR.iid) * cbind(
        riskDR.iid,
        riskDR.iid[, 2] - riskDR.iid[, 1]
      )),
      id = rownames(data), labels = lab
    )
    res0 <- estimate(
      coef = c(est.naive, est.naive[2] - est.naive[1]),
      IC = cbind(ic.naive, ic.naive[, 2] - ic.naive[, 1]),
      labels = lab0,
      id = rownames(data)
    )
    if (requireNamespace("cmprsk", quietly = TRUE)) {
      cc <- cmprsk::cuminc(data$time, data$status, data$a)
      F1.a0 <- cmprsk::timepoints(cc["0 1"], 2)
      F1.a1 <- cmprsk::timepoints(cc["1 1"], 2)
      lab <- paste0("cmprsk.", lab)
      ee <- estimate(
        coef = c(F1.a0$est, F1.a1$est),
        vcov = diag(c(F1.a0$var, F1.a1$var))
      ) |> estimate(cbind(c(1, 0, -1), c(0, 1, 1)), labels = lab)
    } else {
      ee <- NULL
    }
    res0 <- res0 + best0
    res <- structure(res, naive = list(res0, cmprsk = ee))
  }
  res$landmark.time <- time
  class(res) <- c("truncatedscore", class(res))
  return(res)
}

##' @export
summary.truncatedscore <- function(object,
                                   noninf.y = 0,
                                   noninf.t = 0,
                                   alpha = 0.05,
                                   parameter.sign = c(y = 1L, t = -1L),
                                   ...) {

  idx <- if (sign(parameter.sign[1]) < 0) c(1, 2) else c(2, 1)
  idx <- c(idx, if (sign(parameter.sign[2]) < 0) c(4, 5) else c(5, 4))
  nn <- names(coef(object))
  lab <- c(paste0(nn[idx[1]], " - ", nn[idx[2]]),
    paste0(nn[idx[3]], " - ", nn[idx[4]])
  )
  B <- matrix(0, 2, 6)
  B[1, idx[1:2]] <- c(1, -1)
  B[2, idx[3:4]] <- c(1, -1)
  est <- estimate(object, B)
  args <- list(
    thetahat1 = coef(est)[1],
    se1 = vcov(est)[1, 1]**.5,
    thetahat2 = coef(est)[2],
    se2 = vcov(est)[2, 2]**.5,
    corr = cov2cor(vcov(est))[1, 2],
    noninf1 = noninf.y,
    noninf2 = noninf.t,
    alpha = alpha
  )
  pval <- do.call(test_intersectsignedwald, args)
  res1 <- with(pval$test.1, cbind(estimate, statistic, p.value))
  rownames(res1) <- names(pval$test.1$estimate)
  res2 <- with(pval$test.2, cbind(estimate, statistic, p.value))
  rownames(res2) <- names(pval$test.2$estimate)
  res12 <- with(pval$test.intersect, cbind(NA, statistic, p.value))
  rownames(res12) <- "intersection"
  tests <- rbind(res1, res2, res12)
  res <- c(list(
    object = object,
    alpha = alpha,
    estimate = est,
    nonfinf.y = noninf.y,
    noninf.t = noninf.t,
    labels = lab,
    tests = tests
  ), pval)
  class(res) <- "summary.truncatedscore"
  return(res)
}

##' @export
print.summary.truncatedscore <- function(x, ...) {
  cat("\n")
  cli::cli_rule("Parameter estimates")
  print(x$object)
  cat("\n")
  cli::cli_rule("One-sided tests")
  cat("\nb1 = ", x$labels[1], "\n")
  print(x$test.1)
  cli::cli_h3("")
  cat("\nb2 = ", x$labels[2], "\n")
  print(x$test.2)
  cli::cli_rule("Intersection test")
  print(x$test.intersect)
  cli::cli_rule()
  cat("\n")
}

##' @export
parameter.summary.truncatedscore <- function(x, ...) {
  x$tests
}
##' @export
coef.summary.truncatedscore <- function(object, ...) {
  object$tests
}
