#' Prediction of overall survival through Bayesian model averaging
#'
#' This function improves the accuracy of the OS forecast by combining joint
#' models developed based on each granular component of PFS through Bayesian
#' model averaging.
#'
#' @param data trial dataset with progression-free survival info (measurements
#'  of target lesions, time to non-target lesion, time to new lesion, time to
#'  death), covaraites, patient ID, randomization date, etc. Please see
#'  \code{\link{data}} for more detailed descriptions.
#' @param covariate a vector of column names of all the covariates.
#' @param longest_survival the longest possible survival of a patient in this
#'  trial.
#' @param n.MCMC.chain number of MCMC chains, default to 2.
#' @param no.adapt the number of iterations for adaptation.
#' @param burn.in number of burin-in iterations for MCMC.
#' @param MCMC.sample number of iterations for MCMC.
#' @param thin thinning interval for monitors.
#' @param method same as the `method` argument in \code{\link[runjags]{run.jags}},
#'  the method with which to call JAGS; probably a character vector specifying
#'  one of 'rjags', 'simple', 'interruptible', 'parallel', 'rjparallel',
#'  'background', 'bgparallel' or 'snow'. The 'rjags' and 'rjparallel' methods
#'  run JAGS using the rjags package, whereas other options do not require the
#'  rjags package and call JAGS as an external executable.
#' @param beta.prior.sd the standard deviation of the normal priors of all betas
#'  (coefficient of covariates)
#' @param alpha.prior the rate of the exponential priors of all alphas (the shape
#'  parameter)
#' @param eta.prior the rate of the exponential priors of all etas (correlations
#'  in the Clayton copula)
#' @param lambda.prior.sd the standard deviation of the normal priors of all
#'  lambdas
#' @param b.prior.sd the standard deviation of the normal priors of all random
#'  effects `b`
#' @param y.prior.sd the standard deviation of the normal prior of Y (the observed
#'  lesion measurement)
#' @param ... optional arguments that are passed to \code{runjags::run.jags()}
#'  function.
#'
#' @import ggplot2
#' @import runjags
#'
#' @return
#' \item{posterior_sample}{an \code{runjags} object generated through the \code{runjags::run.jags}
#'    function, which includes posterior samples of all the parameters and predicted
#'    survival times}
#' \item{rand_date}{the randomization dates of patients}
#' \item{covariate}{the covariates columns of the input dataset}
#' \item{OS}{the overall survival and its censored status provided in the input
#'    dataset}
#' \item{ID}{the patient ID}
#'
#'
#' @examples
#' \donttest{
#' result <- osbma(osbma::data, covariate = "trt")
#' }
#'
#' @references
#' (manuscript under review)
#' @seealso
#' \code{\link{predict.osbma}} \cr
#' \code{\link{plot.predict.osbma}} \cr
#' \code{\link[runjags]{run.jags}}
#'
#' @rdname osbma
#' @export

osbma <- function(data, covariate, longest_survival = 2000,
                  n.MCMC.chain = 2, no.adapt = 20000,
                  burn.in = 10000, MCMC.sample = 30000, thin = 3,
                  method = "rjparallel", beta.prior.sd = 1,
                  alpha.prior = 1, eta.prior = 1,
                  lambda.prior.sd = 1,
                  b.prior.sd = 1, y.prior.sd = 0.1, ...) {
  beta.prior.sd <- 1 / (beta.prior.sd)^2
  lambda.prior.sd <- 1 / (lambda.prior.sd)^2
  b.prior.sd <- 1 / (b.prior.sd)^2
  y.prior.sd <- 1 / (y.prior.sd)^2

  N <- nrow(data)
  delta_NL <- data$delta_NL
  delta_NT <- data$delta_NT
  delta_OS <- data$delta_OS
  delta_PFS <- as.numeric(data$delta_PFS)
  time_PFS <- data$time_PFS
  time_OS <- data$time_OS
  Y_os <- data$time_OS
  Y_2 <- data$time_NT
  Y_3 <- data$time_NL
  X <- data[covariate]
  Z_1 <- cbind(rep(1, N), X)
  Npar_0 <- ncol(X)
  Npar <- Npar_0 + 1
  Z_2 <- Z_1
  mu_colno <- stringr::str_detect(colnames(data), "mu.t")
  Y <- data[mu_colno]

  evl_time_colno <- stringr::str_detect(colnames(data), "evl_time.t")
  t <- data[evl_time_colno]

  evl_times <- ncol(Y)
  t.os <- numeric(N)
  t.cen <- numeric(N)
  T <- numeric(N)
  T1 <- numeric(N)
  p <- numeric(4)

  p[1] <- Npar_0 + 7 + Npar
  p[2] <- Npar + 2
  p[3] <- Npar + 2
  p[4] <- Npar + 1

  for (i in 1:N) {
    if (delta_OS[i] == 1) {
      t.os[i] <- Y_os[i]
      t.cen[i] <- 0
    } else {
      t.os[i] <- NA
      t.cen[i] <- Y_os[i]
    }
  }


  for (i in 1:N) {
    d <- Y[i, ]
    d <- d[!is.na(d)]
    T[i] <- length(d)
    Y[i, ] <- c(d, rep(0, evl_times - length(d)))
  }

  for (i in 1:N) {
    d <- t[i, ]
    d <- d[!is.na(d)]
    T1[i] <- length(d)
    t[i, ] <- c(d, rep(0, evl_times - length(d)))
  }

  Y <- cbind(Y, rep(0, N))
  t <- cbind(t, rep(0, N))

  t.os <- as.numeric(t.os)
  t.cen <- as.numeric(t.cen)
  t <- data.matrix(t)
  T <- as.numeric(T)
  Y <- data.matrix(Y)
  Y_os <- as.numeric(Y_os)
  delta_OS <- as.numeric(delta_OS)
  Y_2 <- as.numeric(Y_2)
  Y_3 <- as.numeric(Y_3)
  delta_NL <- as.numeric(delta_NL)
  delta_NT <- as.numeric(delta_NT)
  X <- data.matrix(X)
  Z_1 <- data.matrix(Z_1)
  Z_2 <- data.matrix(Z_2)
  b <- matrix(rep(rep(0, N), 2), ncol = 2)
  b_os.1 <- rep(10, N)
  os.2 <- rep(log(1), N)
  nt.2 <- rep(log(0.4), N)
  os.3 <- rep(log(1), N)
  nl.3 <- rep(log(0.4), N)

  ######### joint model ##########

  # bug files written to temporary directory on function call to satisfy CRAN
  # requirements of not accessing user's system files

  # "ospred_get_posterior_4JM.bug"
  ospred_get_posterior_4JM_file <- tempfile(fileext = ".bug")
  writeLines(ospred_get_posterior_4JM_text(), con = ospred_get_posterior_4JM_file)

  # "ospred_4JM.bug"
  ospred_4JM_file <- tempfile(fileext = ".bug")
  writeLines(ospred_4JM_text(), con = ospred_4JM_file)

  bugfile1 <- readLines(ospred_get_posterior_4JM_file)
  bugfile1 <- gsub(pattern = "beta.prior.sd", replacement = beta.prior.sd, x = bugfile1)
  bugfile1 <- gsub(pattern = "alpha.prior", replacement = alpha.prior, x = bugfile1)
  bugfile1 <- gsub(pattern = "eta.prior", replacement = eta.prior, x = bugfile1)
  bugfile1 <- gsub(pattern = "lambda.prior.sd", replacement = lambda.prior.sd, x = bugfile1)
  bugfile1 <- gsub(pattern = "y.prior.sd", replacement = y.prior.sd, x = bugfile1)
  bugfile2 <- gsub(pattern = "b.prior.sd", replacement = b.prior.sd, x = bugfile1)

  bugfile2_file <- tempfile(fileext = ".bug")
  writeLines(bugfile2, con = bugfile2_file)

  jag.model.name <- bugfile2_file
  t.cen.lim.wb.os <- round(cbind(ifelse(delta_OS == 1, 0, Y_os / 365), rep(longest_survival / 365, length(delta_OS))), 3) # round(ifelse(delta_OS == 1, longest_survival, Y_os)/365, 3)  #
  t.cen.lim <- round(cbind(ifelse(delta_OS == 1, 0, Y_os / 365), rep(longest_survival / 365, length(delta_OS))), 3) # round(ifelse(delta_OS == 1, longest_survival, Y_os)/365, 3) #
  t.cen.lim.nl <- round(cbind(ifelse(delta_NL == 1, 0, Y_3 / 365), rep(longest_survival / 365, length(delta_NL))), 3) # round(ifelse(delta_NL == 1, longest_survival, Y_3)/365, 3)  #
  y.lim_t <- y.lim_t.os <- rep(1, length(delta_OS)) # 1 - delta_OS #
  y.lim_wb <- y.lim_wb.os <- rep(1, length(delta_OS)) # 1 - delta_OS #
  y.lim_nl.os <- y.lim_nt.os <- rep(1, length(delta_OS)) # 1 - delta_OS #
  y.lim_nl <- rep(1, length(delta_NL)) # 1 - delta_NL  #
  y.lim_nt <- rep(1, length(delta_NT)) # 1 - delta_NT #
  Y.os.t_rep <- Y.os.nl_rep <- Y.os.nt_rep <- Y.os.wb_rep <- round(ifelse(t.cen == 0, NA, (t.cen + 1) / 365), 3)
  Y.nl_rep <- round(ifelse(delta_NL == 1, NA, (Y_3 + 1) / 365), 3)
  Y.nt_rep <- round(ifelse(delta_NT == 1, NA, (Y_2 + 1) / 365), 3)

  TL_ob_n <- which(delta_OS == 1)
  TL_c_n <- which(delta_OS == 0)

  posterior_sample <- NULL
  attempt <- 1
  while (is.null(posterior_sample) && attempt <= 5) {
    attempt <- attempt + 1

    try({
      posterior_sample <- runjags::run.jags(file.path(jag.model.name),
        data = list(
          N = N,
          longest_survival = round(longest_survival / 365, 3),
          TL_ob_n = TL_ob_n,
          TL_c_n = TL_c_n,
          Npar_0 = Npar_0,
          Npar = Npar,
          X = X,
          Z_1 = Z_1,
          Z_2 = Z_2,
          t.os = round(t.os / 365, 3),
          t.os.wb = round(t.os / 365, 3),
          t.cen.lim = t.cen.lim,
          t.cen.lim.wb.os = t.cen.lim.wb.os,
          t = round(t / 365, 3),
          T = T,
          Y = Y,
          Y_os = round(Y_os / 365, 3),
          Y_2 = round(Y_2 / 365, 3),
          Y_3 = round(Y_3 / 365, 3),
          delta_NT = delta_NT,
          delta_NL = delta_NL,
          delta_OS = delta_OS,
          evl_times = evl_times
        ),
        inits <- function() {
          list(
            Y.os.nl_rep = ifelse(is.na(Y.os.nl_rep), 1 / 365, Y.os.nl_rep),
            Y.os.nt_rep = ifelse(is.na(Y.os.nt_rep), 1 / 365, Y.os.nt_rep),
            t.os = Y.os.t_rep,
            t.os.wb = Y.os.wb_rep,
            Y.os.t_rep = ifelse(is.na(Y.os.t_rep), 1 / 365, Y.os.t_rep),
            Y.os.wb_rep = ifelse(is.na(Y.os.wb_rep), 1 / 365, Y.os.wb_rep)
          )
        },
        thin = thin,
        adapt = no.adapt,
        burnin = burn.in,
        n.chains = n.MCMC.chain,
        monitor = c(
          "beta_0", "beta_1", "lambda", "alpha_os.1", "eta_2",
          "alpha_2", "beta_2", "alpha_os.2", "beta_os.2", "eta_3",
          "alpha_3", "beta_3", "alpha_os.3", "beta_os.3",
          "alpha_os.4", "beta_os.4",
          "sigma_os.1", "sigma_os.2", "sigma_os.3",
          "sigma_b1", "sigma_b2", "rho"
        ),
        sample = MCMC.sample,
        method = method,
        ...
      )
    })
  }

  out_post <- summary(posterior_sample)

  # beta_0, beta_1, alpha_os.1, lambda, sigma_b1, sigma_b2, rho, sigma_os.1
  beta0_name <- paste0("beta_0[", seq(1, Npar, 1), "]")
  beta1_name <- paste0("beta_1[", seq(1, Npar, 1), "]")

  mg1 <- out_post[c(beta0_name, beta1_name, c(
    "alpha_os.1",
    "lambda", "sigma_b1", "sigma_b2", "rho", "sigma_os.1"
  )), "Mean"]
  seg1 <- out_post[c(beta0_name, beta1_name, c(
    "alpha_os.1",
    "lambda", "sigma_b1", "sigma_b2", "rho", "sigma_os.1"
  )), "SD"]

  beta_os_name <- paste0("beta_os.2[", seq(1, Npar, 1), "]")

  # beta_os.2, alpha_os.2, sigma_os.2
  mg2 <- out_post[c(beta_os_name, c(
    "alpha_os.2",
    "sigma_os.2"
  )), "Mean"]
  seg2 <- out_post[c(beta_os_name, c(
    "alpha_os.2",
    "sigma_os.2"
  )), "SD"]

  beta_os.3_name <- paste0("beta_os.3[", seq(1, Npar, 1), "]")

  # beta_os.3, alpha_os.3, sigma_os.3
  mg3 <- out_post[c(beta_os.3_name, c(
    "alpha_os.3",
    "sigma_os.3"
  )), "Mean"]
  seg3 <- out_post[c(beta_os.3_name, c(
    "alpha_os.3",
    "sigma_os.3"
  )), "SD"]

  beta_os.4_name <- paste0("beta_os.4[", seq(1, Npar, 1), "]")
  # beta_os.4, alpha_os.4
  mg4 <- out_post[c(beta_os.4_name, c("alpha_os.4")), "Mean"]
  seg4 <- out_post[c(beta_os.4_name, c("alpha_os.4")), "SD"]

  bugfile3 <- readLines(ospred_4JM_file)
  bugfile3 <- gsub(pattern = "beta.prior.sd", replacement = beta.prior.sd, x = bugfile3)
  bugfile3 <- gsub(pattern = "alpha.prior", replacement = alpha.prior, x = bugfile3)
  bugfile3 <- gsub(pattern = "eta.prior", replacement = eta.prior, x = bugfile3)
  bugfile3 <- gsub(pattern = "lambda.prior.sd", replacement = lambda.prior.sd, x = bugfile3)
  bugfile3 <- gsub(pattern = "y.prior.sd", replacement = y.prior.sd, x = bugfile3)
  bugfile4 <- gsub(pattern = "b.prior.sd", replacement = b.prior.sd, x = bugfile3)

  bugfile4_file <- tempfile(fileext = ".bug")
  writeLines(bugfile4, con = bugfile4_file)

  jag.model.name <- bugfile4_file
  t.cen.lim.wb.os <- round(cbind(ifelse(delta_OS == 1, 0, Y_os / 365), rep(longest_survival / 365, length(delta_OS))), 3) # round(ifelse(delta_OS == 1, longest_survival, Y_os)/365, 3) #
  t.cen.lim <- round(cbind(ifelse(delta_OS == 1, 0, Y_os / 365), rep(longest_survival / 365, length(delta_OS))), 3) # round(ifelse(delta_OS == 1, longest_survival, Y_os)/365, 3) #
  t.cen.lim.nl <- round(cbind(ifelse(delta_NL == 1, 0, Y_3 / 365), rep(longest_survival / 365, length(delta_NL))), 3) # round(ifelse(delta_NL == 1, longest_survival, Y_3)/365, 3) #
  y.lim_t <- y.lim_t.os <- rep(1, length(delta_OS)) # 1 - delta_OS #
  y.lim_wb <- y.lim_wb.os <- rep(1, length(delta_OS)) # 1 - delta_OS #
  y.lim_nl.os <- y.lim_nt.os <- rep(1, length(delta_OS)) # 1 - delta_OS #
  y.lim_nl <- rep(1, length(delta_NL)) # 1 - delta_NL  #
  y.lim_nt <- rep(1, length(delta_NT)) # 1 - delta_NT #
  Y.os.t_rep <- Y.os.nl_rep <- Y.os.nt_rep <- Y.os.wb_rep <- round(ifelse(t.cen == 0, NA, (t.cen + 1) / 365), 3)
  Y.nl_rep <- round(ifelse(delta_NL == 1, NA, (Y_3 + 1) / 365), 3)
  Y.nt_rep <- round(ifelse(delta_NT == 1, NA, (Y_2 + 1) / 365), 3)

  TL_ob_n <- which(delta_OS == 1)
  TL_c_n <- which(delta_OS == 0)

  posterior_sample <- NULL
  attempt <- 1
  while (is.null(posterior_sample) && attempt <= 3) {
    attempt <- attempt + 1

    try({
      posterior_sample <- runjags::run.jags(file.path(jag.model.name),
        data = list(
          N = N,
          longest_survival = round(longest_survival / 365, 3),
          TL_ob_n = TL_ob_n,
          TL_c_n = TL_c_n,
          Npar_0 = Npar_0,
          Npar = Npar,
          X = X,
          Z_1 = Z_1,
          Z_2 = Z_2,
          t.os = round(t.os / 365, 3),
          t.os.wb = round(t.os / 365, 3),
          t.cen.lim = t.cen.lim,
          t.cen.lim.wb.os = t.cen.lim.wb.os,
          t = round(t / 365, 3),
          T = T,
          Y = Y,
          Y_os = round(Y_os / 365, 3),
          Y_2 = round(Y_2 / 365, 3),
          Y_3 = round(Y_3 / 365, 3),
          delta_NT = delta_NT,
          delta_NL = delta_NL,
          delta_OS = delta_OS,
          evl_times = evl_times,
          mg1 = mg1,
          seg1 = seg1,
          mg2 = mg2,
          seg2 = seg2,
          mg3 = mg3,
          seg3 = seg3,
          mg4 = mg4,
          seg4 = seg4,
          p = p
        ),
        inits = function() {
          list(
            Y.os.nl_rep = ifelse(is.na(Y.os.nl_rep), 1 / 365, Y.os.nl_rep),
            Y.os.nt_rep = ifelse(is.na(Y.os.nt_rep), 1 / 365, Y.os.nt_rep),
            # Y.nl_rep = Y.nl_rep,
            # Y.nt_rep = Y.nt_rep,
            t.os = Y.os.t_rep,
            t.os.wb = Y.os.wb_rep,
            Y.os.t_rep = ifelse(is.na(Y.os.t_rep), 1 / 365, Y.os.t_rep),
            Y.os.wb_rep = ifelse(is.na(Y.os.wb_rep), 1 / 365, Y.os.wb_rep)
          )
        },
        thin = thin,
        adapt = no.adapt,
        burnin = burn.in,
        n.chains = n.MCMC.chain,
        monitor = c(
          "beta_0", "beta_1", "lambda", "alpha_os.1", "eta_2",
          "alpha_2", "beta_2", "alpha_os.2", "beta_os.2", "eta_3",
          "alpha_3", "beta_3", "alpha_os.3", "beta_os.3",
          "alpha_os.4", "beta_os.4",
          "Y.os.t_rep", "Y.os.nl_rep", "Y.os.nt_rep", "Y.os.wb_rep",
          "sigma_os.1", "sigma_os.2", "sigma_os.3",
          "sigma_b1", "sigma_b2", "rho",
          "Y.mix_rep", "w", "LL", "Pr", "g"
        ),
        sample = MCMC.sample,
        method = method, ...
      )
    })
  }

  JM_sample <- list(
    "posterior_sample" = posterior_sample$mcmc,
    "rand_date" = data$rand_date,
    "covariate" = X,
    "OS" = cbind(data$time_OS, data$delta_OS),
    "ID" = cbind(data$ID)
  )
  class(JM_sample) <- "osbma"
  return(JM_sample)
}

#' Summarizing osbma fits
#'
#' `summary` method for class "`osbma`"
#'
#' @param object an object of class "`osbma`", usually, a result of a call to
#'  \code{\link{osbma}}
#' @param ... further arguments included in the \link[runjags]{summary.runjags} function.
#'
#' @returns
#' The summary method returns a numeric matrix of summary statistics for each
#' variable. Please check \code{\link[runjags]{summary.runjags}} for more
#' details.
#'
#' @rdname osbma
#' @export
summary.osbma <- function(object, ...) {
  s_return <- summary(object$posterior_sample, ...)
  return(s_return)
}

#' @rdname osbma
#' @param x `osbma` object to print
#' @param ... further arguments. Not currently used.
#' @export

print.osbma <- function(x, ...) {
  print(x$posterior_sample)
}
#' overall survival prediction
#'
#' `predict` method for class "`osbma`"
#'
#' @param object an object of class "`osbma`", usually, a result of a call to
#'  \code{\link{osbma}}
#' @param quantile the quantile of the predicted overall survival
#' @param ... further arguments. Not currently used
#'
#' @returns
#' \item{prediction}{the predicted overall survival based on \code{osbma} fit}
#' \item{prediction_date}{the death date predicted by the \code{osbma} model}
#' \item{quantile}{the quantile specified by the user}
#' \item{covariate}{the columns of covariates provided in the input dataset}
#' \item{rand_date}{the randomization date of the patients in the trial}
#' \item{OS}{the overall survival and its censored status provided in the input
#'    dataset}
#' \item{ID}{the ID of the patients provided in the input dataset}
#' @rdname osbma
#' @export


predict.osbma <- function(object, quantile = 0.9, ...) {
  JM_median <- summary(object$posterior_sample)$quantiles[, 3]
  pred_col_JM <- stringr::str_detect(names(JM_median), "Y.mix_rep")
  pred_OS_JM <- JM_median[pred_col_JM]
  pred_OS_JM_date <- pred_OS_JM * 30 + object$rand_date
  JM_quantile <- coda::HPDinterval(object$posterior_sample,
    prob = quantile
  )[[1]][pred_col_JM, ]
  obj <- list(
    prediction = pred_OS_JM,
    prediction_date = pred_OS_JM_date,
    quantile = JM_quantile,
    covariate = object$covariate,
    rand_date = object$rand_date,
    OS = object$OS,
    ID = object$ID
  )
  class(obj) <- "predict.osbma"
  obj
}

#' plot based on prediction outcomes of an osbma object
#'
#' @param x an object of class "`predict.osbma`", usually, a result of a call to
#'  \code{\link{predict.osbma}}
#' @param trt.col.name the name of the assigned treatment column, default to `trt`
#' @param type type of plot the users want: "date" generates a plot with patients'
#'  survival timeline; "KM" generates a Kaplan-Meier curve of overall survival;
#'  a number n between 0 - N, where N is the total number of patients enrolled,
#'  generates a plot of the estimated time of the nth death in the trial.
#' @param ... further arguments. Not currently used.
#'
#' @returns
#' 3 types of plots
#' @rdname osbma
#' @export
plot.predict.osbma <- function(x, trt.col.name = "trt", type = "date", ...) {
  if (type == "KM") {
    # KM plot
    data <- as.data.frame(cbind(ifelse(x$OS[, 2] == 0, x$prediction, x$OS[, 1]), x$covariate[, trt.col.name]))
    colnames(data) <- c("os", "trt")
    trt1 <- subset(data, trt == 1)
    trt2 <- subset(data, trt == 2)

    # Time to death for OS
    fit1 <- survival::survfit(survival::Surv(os, rep(1, nrow(trt1))) ~ 1, data = trt1)
    fit2 <- survival::survfit(survival::Surv(os, rep(1, nrow(trt2))) ~ 1, data = trt2)

    fit <- list(trt_1 = fit1, trt_2 = fit2)

    survminer::ggsurvplot(
      fit,
      combine = TRUE,
      conf.int = TRUE,
      risk.table = "abs_pct",
      xlab = "Months",
      ylab = "Overall survival probability",
      legend.labs = c("trt1", "trt2"),
      title = "Kaplan-Meier Curve of Overall Survival"
    )
  } else if (type == "date") {
    # the nth death plot
    data <- data.frame(
      "date" = data.table::fifelse(x$OS[, 2] == 1, x$OS[, 1] * 30 + x$rand_date, as.Date(x$prediction_date)),
      "trt" = as.factor(x$covariate[, trt.col.name]),
      "rand_date" = x$rand_date,
      "ID" = x$ID
    )
    data_plot <- data[order(data$date), ]
    data_plot$death <- seq(1, length(x$prediction), 1)

    g <- ggplot(data_plot, aes(colour = trt, text = paste("ID:", ID, "\ndeath:", death, "\nrand.date:", rand_date, "\ndeath.date:", date))) +
      geom_linerange(data = data_plot, mapping = aes(x = death, ymin = rand_date, ymax = date)) +
      ylab("time") +
      ggtitle("Patients' Survival Timeline")

    plotly::ggplotly(g, tooltip = "text")
  } else {
    # predicted date range plot
    data <- data.frame(
      "date" = data.table::fifelse(x$OS[, 2] == 1, x$OS[, 1] * 30 + x$rand_date, as.Date(x$prediction_date)),
      "lower_bound" = x$quantile[, 1],
      "upper_bound" = x$quantile[, 2],
      "trt" = as.factor(x$covariate[, trt.col.name]),
      "rand_date" = x$rand_date,
      "ID" = x$ID
    )
    data_plot <- data[order(data$date), ]
    data_plot$death <- seq(1, length(x$prediction), 1)
    death.number <- as.integer(type)

    q <- ggplot(
      data_plot[death.number, ],
      aes(text = paste(
        "ID:", ID,
        "\ndeath:", death,
        "\nrand.date:", rand_date,
        "\ndeath.date:", date,
        "\nlower.bound:", lower_bound * 30 + rand_date,
        "\nupper.bound:", upper_bound * 30 + rand_date
      ))
    ) +
      geom_linerange(
        data = data_plot[death.number, ], mapping = aes(x = death.number, ymin = lower_bound * 30 + rand_date, ymax = upper_bound * 30 + rand_date),
        size = 1, color = "blue"
      ) +
      geom_point(
        data = data_plot[death.number, ], mapping = aes(x = death.number, y = date),
        size = 4, shape = 21, fill = "white"
      ) +
      ylab("time") +
      ggtitle(paste0("Predicted date of the ", death.number, "th death")) +
      xlab(death.number) +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    plotly::ggplotly(q, tooltip = "text")
  }
}
