# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

osbma <- function(data, covariate, longest.survival = 365 * 20 / 30,
                  n.MCMC.chain = 2, no.adapt = 80,
                  burn.in = 10, MCMC.sample = 100, thin = 10, method = "rjparallel", ...) {
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

  jag.model.name <- ospred_get_posterior_4JM_file
  # jag.model.name <- "ospred_get_posterior_4JM.bug"
  t.cen.lim <- cbind(ifelse(delta_OS == 1, 0, Y_os), rep(longest.survival, length(delta_OS)))
  t.cen.lim.wb.os <- ifelse(t.cen == 0, longest.survival, Y_os)
  t.cen.lim.nl <- cbind(ifelse(delta_NL == 1, 0, Y_3), rep(longest.survival, length(delta_NL))) # cbind(Y_3, rep(longest.survival, length(t.cen)))
  y.lim_t <- y.lim_t.os <- rep(1, length(delta_OS))
  y.lim_wb.os <- y.lim_wb <- 1 - delta_OS
  y.lim_nl.os <- y.lim_nt.os <- rep(1, length(delta_OS))
  y.lim_nl <- rep(1, length(delta_NL)) # 1 - delta_NL
  y.lim_nt <- rep(1, length(delta_NT))
  Y.os.t_rep <- Y.os.nl_rep <- Y.os.nt_rep <- Y.os.wb_rep <- ifelse(t.cen == 0, NA, t.cen + 1)
  Y.nl_rep <- ifelse(delta_NL == 1, NA, Y_3 + 1)
  Y.nt_rep <- ifelse(delta_NT == 1, NA, Y_2 + 1)


  tryCatch({
    posterior_sample <- runjags::run.jags(file.path(jag.model.name),
      data = list(
        N = N,
        Npar_0 = Npar_0,
        Npar = Npar,
        X = X,
        Z_1 = Z_1,
        Z_2 = Z_2,
        t.os = t.os,
        t.os.wb = t.os,
        t.cen.lim = t.cen.lim,
        t.cen.lim.wb.os = t.cen.lim.wb.os,
        y.lim_t = y.lim_t,
        y.lim_t.os = y.lim_t.os,
        y.lim_nt.os = y.lim_nt.os,
        # y.lim_nt = y.lim_nt,
        # y.lim_nl = y.lim_nl,
        y.lim_nl.os = y.lim_nl.os,
        y.lim_wb = y.lim_wb,
        y.lim_wb.os = y.lim_wb.os,
        t = t,
        T = T,
        Y = Y,
        Y_os = Y_os,
        Y_2 = Y_2,
        Y_3 = Y_3,
        delta_NT = delta_NT,
        delta_NL = delta_NL,
        delta_OS = delta_OS,
        evl_times = evl_times
      ),
      inits = function() {
        list(
          Y.os.t_rep = ifelse(is.na(Y.os.t_rep), 1, Y.os.t_rep),
          Y.os.nl_rep = ifelse(is.na(Y.os.nl_rep), 1, Y.os.nl_rep),
          Y.os.nt_rep = ifelse(is.na(Y.os.nt_rep), 1, Y.os.nt_rep),
          # Y.nl_rep = Y.nl_rep,
          # Y.nt_rep = Y.nt_rep,
          t.os = Y.os.t_rep,
          t.os.wb = Y.os.t_rep,
          Y.os.wb_rep = ifelse(is.na(Y.os.wb_rep), 1, Y.os.wb_rep)
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
        "sigma_os.3", "sigma_os.2", "sigma_os.1",
        "sigma_b1", "sigma_b2", "rho"
      ),
      sample = MCMC.sample,
      method = method, ...
    )
  })

  #     jag <- rjags::jags.model(file.path(jag.model.name),
  #       data = list(
  #         N = N,
  #         Npar_0 = Npar_0,
  #         Npar = Npar,
  #         X = X,
  #         Z_1 = Z_1,
  #         Z_2 = Z_2,
  #         t.os = t.os,
  #         t.os.wb = t.os,
  #         t.cen.lim = t.cen.lim,
  #         t.cen.lim.wb.os = t.cen.lim.wb.os,
  #         y.lim_t = y.lim_t,
  #         y.lim_t.os = y.lim_t.os,
  #         y.lim_nt.os = y.lim_nt.os,
  #         # y.lim_nt = y.lim_nt,
  #         # y.lim_nl = y.lim_nl,
  #         y.lim_nl.os = y.lim_nl.os,
  #         y.lim_wb = y.lim_wb,
  #         y.lim_wb.os = y.lim_wb.os,
  #         t = t,
  #         T = T,
  #         Y = Y,
  #         Y_os = Y_os,
  #         Y_2 = Y_2,
  #         Y_3 = Y_3,
  #         delta_NT = delta_NT,
  #         delta_NL = delta_NL,
  #         delta_OS = delta_OS,
  #         evl_times = evl_times
  #       ),
  #       inits <- function() {
  #         list(
  #           Y.os.t_rep = ifelse(is.na(Y.os.t_rep), 1, Y.os.t_rep),
  #           Y.os.nl_rep = ifelse(is.na(Y.os.nl_rep), 1, Y.os.nl_rep),
  #           Y.os.nt_rep = ifelse(is.na(Y.os.nt_rep), 1, Y.os.nt_rep),
  #           # Y.nl_rep = Y.nl_rep,
  #           # Y.nt_rep = Y.nt_rep,
  #           t.os = Y.os.t_rep,
  #           t.os.wb = Y.os.t_rep,
  #           Y.os.wb_rep = ifelse(is.na(Y.os.wb_rep), 1, Y.os.wb_rep)
  #         )
  #       },
  #       n.chains = n.MCMC.chain, n.adapt = no.adapt
  #     )
  #     update(jag, burn.in)
  #
  #     posterior_sample <- rjags::coda.samples(jag,
  #       c(
  #         "beta_0", "beta_1", "lambda", "alpha_os.1", "eta_2",
  #         "alpha_2", "beta_2", "alpha_os.2", "beta_os.2", "eta_3",
  #         "alpha_3", "beta_3", "alpha_os.3", "beta_os.3",
  #         "alpha_os.4", "beta_os.4",
  #         "Y.os.t_rep", "Y.os.nl_rep", "Y.os.nt_rep", "Y.os.wb_rep",
  #         "sigma_os.3", "sigma_os.2", "sigma_os.1",
  #         "sigma_b1", "sigma_b2", "rho"
  #       ),
  #       n.iter = MCMC.sample, thin = thin
  #     )
  #   },
  #   warning = function(war) {},
  #   error = function(err) {},
  #   finally = {}
  # )

  out_post <- summary(posterior_sample)

  # beta_0, beta_1, alpha_os.1, lambda, sigma_b1, sigma_b2, rho, sigma_os.1
  beta0_name <- paste0("beta_0[", seq(1, Npar, 1), "]")
  beta1_name <- paste0("beta_1[", seq(1, Npar, 1), "]")

  mg1 <- out_post[c(beta0_name, beta1_name, c(
    "alpha_os.1",
    "lambda", "sigma_b1", "sigma_b2", "rho", "sigma_os.1"
  )), "Median"]
  seg1 <- out_post[c(beta0_name, beta1_name, c(
    "alpha_os.1",
    "lambda", "sigma_b1", "sigma_b2", "rho", "sigma_os.1"
  )), "SD"]

  beta_os_name <- paste0("beta_os.2[", seq(1, Npar, 1), "]")

  # beta_os.2, alpha_os.2, sigma_os.2
  mg2 <- out_post[c(beta_os_name, c(
    "alpha_os.2",
    "sigma_os.2"
  )), "Median"]
  seg2 <- out_post[c(beta_os_name, c(
    "alpha_os.2",
    "sigma_os.2"
  )), "SD"]

  beta_os.3_name <- paste0("beta_os.3[", seq(1, Npar, 1), "]")

  # beta_os.3, alpha_os.3, sigma_os.3
  mg3 <- out_post[c(beta_os.3_name, c(
    "alpha_os.3",
    "sigma_os.3"
  )), "Median"]
  seg3 <- out_post[c(beta_os.3_name, c(
    "alpha_os.3",
    "sigma_os.3"
  )), "SD"]

  beta_os.4_name <- paste0("beta_os.4[", seq(1, Npar, 1), "]")
  # beta_os.4, alpha_os.4
  mg4 <- out_post[c(beta_os.4_name, c("alpha_os.4")), "Median"]
  seg4 <- out_post[c(beta_os.4_name, c("alpha_os.4")), "SD"]

  jag.model.name <- ospred_4JM_file
  # jag.model.name <- "ospred_4JM.bug"
  t.cen.lim <- cbind(ifelse(delta_OS == 1, 0, Y_os), rep(longest.survival, length(delta_OS)))
  t.cen.lim.wb.os <- ifelse(t.cen == 0, longest.survival, Y_os)
  t.cen.lim.nl <- cbind(ifelse(delta_NL == 1, 0, Y_3), rep(longest.survival, length(delta_NL))) # cbind(Y_3, rep(longest.survival, length(t.cen)))
  y.lim_t <- y.lim_t.os <- rep(1, length(delta_OS))
  y.lim_wb.os <- y.lim_wb <- 1 - delta_OS
  y.lim_nl.os <- y.lim_nt.os <- rep(1, length(delta_OS))
  y.lim_nl <- rep(1, length(delta_NL)) # 1 - delta_NL
  y.lim_nt <- rep(1, length(delta_NT))
  Y.os.t_rep <- Y.os.nl_rep <- Y.os.nt_rep <- Y.os.wb_rep <- ifelse(t.cen == 0, NA, t.cen + 1)
  Y.nl_rep <- ifelse(delta_NL == 1, NA, Y_3 + 1)
  Y.nt_rep <- ifelse(delta_NT == 1, NA, Y_2 + 1)
  # Y.os.t_rep = Y.os.nl_rep = Y.os.wb_rep = Y.nl_rep = t.cen + 1

  tryCatch({
    posterior_sample <- runjags::run.jags(file.path(jag.model.name),
      data = list(
        N = N,
        Npar_0 = Npar_0,
        Npar = Npar,
        X = X,
        Z_1 = Z_1,
        Z_2 = Z_2,
        t.os = t.os,
        t.os.wb = t.os,
        t.cen.lim = t.cen.lim,
        t.cen.lim.wb.os = t.cen.lim.wb.os,
        y.lim_t = y.lim_t,
        y.lim_t.os = y.lim_t.os,
        y.lim_nt.os = y.lim_nt.os,
        # y.lim_nt = y.lim_nt,
        # y.lim_nl = y.lim_nl,
        y.lim_nl.os = y.lim_nl.os,
        y.lim_wb = y.lim_wb,
        y.lim_wb.os = y.lim_wb.os,
        t = t,
        T = T,
        Y = Y,
        Y_os = Y_os,
        Y_2 = Y_2,
        Y_3 = Y_3,
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
          Y.os.t_rep = ifelse(is.na(Y.os.t_rep), 1, Y.os.t_rep),
          Y.os.nl_rep = ifelse(is.na(Y.os.nl_rep), 1, Y.os.nl_rep),
          Y.os.nt_rep = ifelse(is.na(Y.os.nt_rep), 1, Y.os.nt_rep),
          # Y.nl_rep = Y.nl_rep,
          # Y.nt_rep = Y.nt_rep,
          t.os = Y.os.t_rep,
          t.os.wb = Y.os.t_rep,
          Y.os.wb_rep = ifelse(is.na(Y.os.wb_rep), 1, Y.os.wb_rep)
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
        "sigma_os.3", "sigma_os.2", "sigma_os.1",
        "sigma_b1", "sigma_b2", "rho",
        "Y.mix_rep", "w"
      ),
      sample = MCMC.sample,
      method = method, ...
    )
  })
  #     jag <- rjags::jags.model(file.path(jag.model.name),
  #       data = list(
  #         N = N,
  #         Npar_0 = Npar_0,
  #         Npar = Npar,
  #         X = X,
  #         Z_1 = Z_1,
  #         Z_2 = Z_2,
  #         t.os = t.os,
  #         t.os.wb = t.os,
  #         t.cen.lim = t.cen.lim,
  #         t.cen.lim.wb.os = t.cen.lim.wb.os,
  #         y.lim_t = y.lim_t,
  #         y.lim_t.os = y.lim_t.os,
  #         y.lim_nt.os = y.lim_nt.os,
  #         # y.lim_nt = y.lim_nt,
  #         # y.lim_nl = y.lim_nl,
  #         y.lim_nl.os = y.lim_nl.os,
  #         y.lim_wb = y.lim_wb,
  #         y.lim_wb.os = y.lim_wb.os,
  #         t = t,
  #         T = T,
  #         Y = Y,
  #         Y_os = Y_os,
  #         Y_2 = Y_2,
  #         Y_3 = Y_3,
  #         delta_NT = delta_NT,
  #         delta_NL = delta_NL,
  #         delta_OS = delta_OS,
  #         evl_times = evl_times,
  #         mg1 = mg1,
  #         seg1 = seg1,
  #         mg2 = mg2,
  #         seg2 = seg2,
  #         mg3 = mg3,
  #         seg3 = seg3,
  #         mg4 = mg4,
  #         seg4 = seg4,
  #         p = p
  #       ),
  #       inits <- function() {
  #         list(
  #           Y.os.t_rep = ifelse(is.na(Y.os.t_rep), 1, Y.os.t_rep),
  #           Y.os.nl_rep = ifelse(is.na(Y.os.nl_rep), 1, Y.os.nl_rep),
  #           Y.os.nt_rep = ifelse(is.na(Y.os.nt_rep), 1, Y.os.nt_rep),
  #           # Y.nl_rep = Y.nl_rep,
  #           # Y.nt_rep = Y.nt_rep,
  #           t.os = Y.os.t_rep,
  #           t.os.wb = Y.os.t_rep,
  #           Y.os.wb_rep = ifelse(is.na(Y.os.wb_rep), 1, Y.os.wb_rep)
  #         )
  #       },
  #       n.chains = n.MCMC.chain, n.adapt = no.adapt
  #     )
  #
  #     update(jag, burn.in)
  #
  #     posterior_sample <- rjags::coda.samples(jag,
  #       c(
  #         "beta_0", "beta_1", "lambda", "alpha_os.1", "eta_2",
  #         "alpha_2", "beta_2", "alpha_os.2", "beta_os.2", "eta_3",
  #         "alpha_3", "beta_3", "alpha_os.3", "beta_os.3",
  #         "alpha_os.4", "beta_os.4",
  #         "Y.os.t_rep", "Y.os.nl_rep", "Y.os.nt_rep", "Y.os.wb_rep",
  #         "sigma_os.3", "sigma_os.2", "sigma_os.1",
  #         "sigma_b1", "sigma_b2", "rho",
  #         "Y.mix_rep", "w"
  #       ),
  #       n.iter = MCMC.sample, thin = thin
  #     )
  #   },
  #   warning = function(war) {},
  #   error = function(err) {},
  #   finally = {}
  # )

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



predict.osbma <- function(object, quantile, ...) {
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


plot.predict.osbma <- function(object, trt.col.name = "trt", type = "date", ...) {
  if (type == "KM") {
    # KM plot
    data <- as.data.frame(cbind(ifelse(object$OS[, 2] == 0, object$prediction, object$OS[, 1]), object$covariate[, trt.col.name]))
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
      xlab = "Months",
      ylab = "Overall survival probability",
      legend.labs = c("trt1", "trt2")
    )
  } else if (type == "date") {
    # the nth death plot
    data <- data.frame(
      "date" = data.table::fifelse(object$OS[, 2] == 1, object$OS[, 1] * 30 + object$rand_date, as.Date(object$prediction_date)),
      "trt" = as.factor(object$covariate[, trt.col.name]),
      "rand_date" = object$rand_date,
      "ID" = object$ID
    )
    data_plot <- data[order(data$date), ]
    data_plot$death <- seq(1, length(object$prediction), 1)
    # p <- ggplot(data_plot, aes(death, date, colour = trt, text = paste("rand.date:", rand_date, "\nID:", ID))) +
    #   geom_point()
    # plotly::ggplotly(p)

    g <- ggplot(data_plot, aes(colour = trt, text = paste("ID:", ID, "\ndeath:", death, "\nrand.date:", rand_date, "\ndeath.date:", date))) +
      geom_linerange(data = data_plot, mapping = aes(x = death, ymin = rand_date, ymax = date), size = 1) +
      ylab("time") +
      ggtitle("Patients' Survival Timeline")
    # geom_point(data=d, mapping=aes(x=drink, y=mean), size=4, shape=21, fill="white") +
    # opts(title="geom_linerange", plot.title=theme_text(size=40, vjust=1.5))
    plotly::ggplotly(g, tooltip = "text")
  } else {
    data <- data.frame(
      "date" = data.table::fifelse(object$OS[, 2] == 1, object$OS[, 1] * 30 + object$rand_date, as.Date(object$prediction_date)),
      "lower_bound" = object$quantile[, 1],
      "upper_bound" = object$quantile[, 2],
      "trt" = as.factor(object$covariate[, trt.col.name]),
      "rand_date" = object$rand_date,
      "ID" = object$ID
    )
    data_plot <- data[order(data$date), ]
    data_plot$death <- seq(1, length(object$prediction), 1)
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
