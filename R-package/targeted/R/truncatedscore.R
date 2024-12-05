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
testvals <- function(thetahat1,
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
                           prob.fct(SignWald.intersect, alpha, corr) + alpha, 1
  )
  pval1 <- ifelse(SignWald1 > 0, 0.5 * pchisq(SignWald1, 1, lower.tail = F), 1)
  pval2 <- ifelse(SignWald2 > 0, 0.5 * pchisq(SignWald2, 1, lower.tail = F), 1)
  cbind(
    SignWald.intersect,
    critval.intersect,
    pval.intersect,
    SignWald1,
    pval1,
    SignWald2,
    pval2
  )
}

##' @description
##' Let \eqn{Y} denote the clinical outcome, \eqn{A} the binary treatment variable,
##' \eqn{X} baseline covariates, \(T\) the failure time, and \(epsilon=1,2\)
##' the cause of failure. The following are our two target parameters
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
##' @export
truncatedscore_estimate <- function(
                     data,
                     ymod,
                     rmod,
                     amod,
                     eventmod,
                     time = 2,
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
    lab0 <- c(lab0, sprintf("E(Y|R=1,A=%g)", aval))
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
    "^treat", sprintf("Risk\\(T<%g|A=",time),
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
      laab <- paste0("cmprsk.", lab)
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
  class(res) <- c("truncatedscore", class(res))
  return(res)
}

##' @export
summary.truncatedscore <- function(object, ...) {
  return(object)
}
