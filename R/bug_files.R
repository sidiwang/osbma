ospred_get_posterior_4JM_text <- function() {
  return(
    "################ Target Lesion #########################################
data
{
   for (i in 1:N)
   {
       Y_NLOS[i]<- 0
       Y_NTOS[i]<- 0
   }
}

model{

  C <- 1e+4

  tau <- pow(sigma_Y, -2)

  tau_os.1 <- pow(sigma_os.1, -2)

  tau_b1 <- pow(sigma_b1, -2)

  tau_b2 <- pow(sigma_b2, -2)

  tau_b2.b1 <- 1 / (1 - pow(rho, 2)) * tau_b2

  for (i in 1:N) {


    for(j in 1:T[i]) {


      Y[i, j] ~ dnorm(mu_1[i, j], tau)

      mu_1[i, j] <- inprod(X[i, 1:Npar_0], beta_0[1:Npar_0]) + t[i, j] * (beta_0[Npar_0 + 1] + b[i, 2]) + b[i, 1]

    }


    for (j in (T[i] + 1):(evl_times + 1)){

      mu_1[i, j] <- 0;

    }


    mu[i] <- sum(mu_1[i,]) / T[i]



    theta_os.1[i] <- exp(inprod(Z_1[i, 1:Npar], beta_1[1:Npar]) + lambda * mu[i] + b_os.1[i])

    Y.os.t_rep[i] ~ dweib(alpha_os.1, theta_os.1[i]) T(t.cen.lim[i, 1], longest_survival)


    b[i,1] ~ dnorm(0, tau_b1)

    b[i,2] ~ dnorm(mu_b2.b1[i], tau_b2.b1)

    mu_b2.b1[i] <- sigma_b2 * b[i,1] * rho / sigma_b1

    b_os.1[i] ~ dnorm(0, tau_os.1)

  }

  for (i in TL_ob_n){

    t.os[i] ~ dweib(alpha_os.1, theta_os.1[i]) T(, longest_survival)

  }

  for (i in TL_c_n){

    t.os[i] ~ dweib(alpha_os.1, theta_os.1[i]) T(t.cen.lim[i, 1], longest_survival)
  }




  ######################### Priors ###########################################


  for (i in 1:(Npar_0 + 1)) {
    beta_0[i] ~ dnorm(0.0, beta.prior.sd)
  }


  for (i in 1:Npar) {
    beta_1[i] ~ dnorm(0.0, beta.prior.sd)
  }


  alpha_os.1 ~ dexp(alpha.prior)

  lambda ~ dnorm(0, lambda.prior.sd)



  sigma_Y ~ dnorm(0, y.prior.sd) T(0,)

  sigma_b1 ~ dnorm(0, b.prior.sd) T(0,)

  sigma_b2 ~ dnorm(0, b.prior.sd) T(0,)

  sigma_os.1 ~ dnorm(0, b.prior.sd) T(0,)

  rho ~ dunif(-1, 1)


  ############################  Non target   #########################################

  tau_os.2 <- pow(sigma_os.2, -2)
  tau_nt.2 <- pow(sigma_nt.2, -2)


  for(i in 1:N) {



    NT.u[i] <- -pow(Y_2[i], alpha_NT.2b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_2[1:Npar])) * (-eta_2)
    NT.v[i] <- -pow(Y_os[i], alpha_2b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.2[1:Npar])) * (-eta_2)
    NT.LSE[i] <- log(exp(NT.u[i] - max(NT.u[i], NT.v[i], 0)) + exp(NT.v[i] - max(NT.u[i], NT.v[i], 0)) - exp(0 - max(NT.u[i], NT.v[i], 0))) + max(NT.u[i], NT.v[i], 0)

    L1.NT[i] <- L1.NT1[i] + L1.NT2[i] + L1.NT3[i] + L1.NT4[i]

    L1.NT1[i] <- log(eta_2 + 1) + (-pow(Y_2[i], (alpha_NT.2b[i])) * exp(inprod(Z_2[i, 1:Npar], beta_2[1:Npar])) - pow(Y_os[i], alpha_2b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.2[1:Npar]))) * (-(eta_2+1))
    L1.NT2[i] <- (-(2 * eta_2 + 1)/(eta_2)) *  NT.LSE[i]
    L1.NT3[i] <- log(alpha_NT.2b[i] * pow(Y_2[i], alpha_NT.2b[i] - 1)) + (inprod(Z_2[i, 1:Npar], beta_2[1:Npar]) - pow(Y_2[i], alpha_NT.2b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_2[1:Npar])))
    L1.NT4[i] <- log(alpha_2b[i] * pow(Y_os[i], alpha_2b[i] - 1)) + (inprod(Z_2[i, 1:Npar], beta_os.2[1:Npar]) - pow(Y_os[i], alpha_2b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.2[1:Npar])))


    L2.NT[i] <- L2.NT1[i] + L2.NT2[i] + L2.NT3[i]

    L2.NT1[i] <- (-pow(Y_2[i], alpha_NT.2b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_2[1:Npar]))) * (eta_2 + 1)
    L2.NT2[i] <- (-(eta_2 + 1)/(eta_2)) *  NT.LSE[i]
    L2.NT3[i] <- L1.NT3[i]


    L3.NT[i] <- L3.NT1[i] + L3.NT2[i] + L3.NT3[i]

    L3.NT1[i] <- (-pow(Y_os[i], alpha_2b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.2[1:Npar]))) * (eta_2 + 1)
    L3.NT2[i] <- (-(eta_2 + 1)/(eta_2)) *  NT.LSE[i]
    L3.NT3[i] <- L1.NT4[i]


    L4.NT[i] <- (-1/eta_2) *  NT.LSE[i]


    L.NT[i] <- L1.NT[i] * (delta_NT[i] * delta_OS[i]) + L2.NT[i] * (delta_NT[i] * (1 - delta_OS[i])) + L3.NT[i] * ((1 - delta_NT[i]) * delta_OS[i]) + L4.NT[i] * ((1 - delta_OS[i]) * (1 - delta_NT[i]))



    phi_NT[i] <- -L.NT[i] + C
    Y_NTOS[i] ~ dpois(phi_NT[i])


    os.2[i] ~ dnorm(0, tau_os.2)
    nt.2[i] ~ dnorm(0, tau_nt.2)

    b_os.2[i] <- exp(os.2[i])
    b_nt.2[i] <- exp(nt.2[i])


    alpha_2b[i] <- alpha_os.2 + b_os.2[i]
    alpha_NT.2b[i] <- alpha_2 + b_nt.2[i]


    ####### marginal model for OS ##########


    theta_os.2[i] <- exp(inprod(Z_2[i, 1:Npar], beta_os.2[1:Npar]))

    Y.os.nt_rep[i] ~ dweib(alpha_2b[i], theta_os.2[i]) T(t.cen.lim[i, 1], longest_survival)


  }


  ################## priors  ########################



  for (i in 1:Npar) {
    beta_2[i] ~ dnorm(0.0, beta.prior.sd)
  }

  for (i in 1:Npar) {
    beta_os.2[i] ~ dnorm(0.0, beta.prior.sd)
  }

  eta_2 ~ dexp(eta.prior)


  alpha_2 ~ dexp(alpha.prior)
  alpha_os.2 ~ dexp(alpha.prior)

  sigma_os.2 ~ dnorm(0, b.prior.sd) T(0,)
  sigma_nt.2 ~ dnorm(0, b.prior.sd) T(0,)


  ############################## New Lesion #######################################

  tau_os.3 <- pow(sigma_os.3, -2)
  tau_nl.3 <- pow(sigma_nl.3, -2)

  for(i in 1:N) {


    NL.u[i] <- -pow(Y_3[i], alpha_NL.3b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_3[1:Npar])) * (-eta_3)
    NL.v[i] <- -pow(Y_os[i], alpha_3b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.3[1:Npar])) * (-eta_3)
    NL.LSE[i] <- log(exp(NL.u[i] - max(NL.u[i], NL.v[i], 0)) + exp(NL.v[i] - max(NL.u[i], NL.v[i], 0)) - exp(0 - max(NL.u[i], NL.v[i], 0))) + max(NL.u[i], NL.v[i], 0)


    L1.NL[i] <- L1.NL1[i] + L1.NL2[i] + L1.NL3[i] + L1.NL4[i]

    L1.NL1[i] <- log(eta_3 + 1) + (-pow(Y_3[i], (alpha_NL.3b[i])) * exp(inprod(Z_2[i, 1:Npar], beta_3[1:Npar])) - pow(Y_os[i], alpha_3b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.3[1:Npar]))) * (-(eta_3+1))
    L1.NL2[i] <- (-(2 * eta_3 + 1)/(eta_3)) * NL.LSE[i]
    L1.NL3[i] <- log(alpha_NL.3b[i] * pow(Y_3[i], alpha_NL.3b[i] - 1)) + (inprod(Z_2[i, 1:Npar], beta_3[1:Npar]) - pow(Y_3[i], alpha_NL.3b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_3[1:Npar])))
    L1.NL4[i] <- log(alpha_3b[i] * pow(Y_os[i], alpha_3b[i] - 1)) + (inprod(Z_2[i, 1:Npar], beta_os.3[1:Npar]) - pow(Y_os[i], alpha_3b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.3[1:Npar])))


    L2.NL[i] <- L2.NL1[i] + L2.NL2[i] + L2.NL3[i]

    L2.NL1[i] <- (-pow(Y_3[i], alpha_NL.3b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_3[1:Npar]))) * (eta_3 + 1)
    L2.NL2[i] <- (-(eta_3 + 1)/(eta_3)) * NL.LSE[i]
    L2.NL3[i] <- L1.NL3[i]


    L3.NL[i] <- L3.NL1[i] + L3.NL2[i] + L3.NL3[i]

    L3.NL1[i] <- (-pow(Y_os[i], alpha_3b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.3[1:Npar]))) * (eta_3 + 1)
    L3.NL2[i] <- (-(eta_3 + 1)/(eta_3)) * NL.LSE[i]
    L3.NL3[i] <- L1.NL4[i]


    L4.NL[i] <- (-1/eta_3) * NL.LSE[i]


    L.NL[i] <- L1.NL[i] * (delta_NL[i] * delta_OS[i]) + L2.NL[i] * (delta_NL[i] * (1 - delta_OS[i])) + L3.NL[i] * ((1 - delta_NL[i]) * delta_OS[i]) + L4.NL[i] * ((1 - delta_OS[i]) * (1 - delta_NL[i]))



    phi_NL[i] <- -L.NL[i] + C
    Y_NLOS[i] ~ dpois(phi_NL[i])

    os.3[i] ~ dnorm(0, tau_os.3)
    nl.3[i] ~ dnorm(0, tau_nl.3)

    b_os.3[i] <- exp(os.3[i])
    b_nl.3[i] <- exp(nl.3[i])

    alpha_3b[i] <- alpha_os.3 + b_os.3[i]
    alpha_NL.3b[i] <- alpha_3 + b_nl.3[i]

    ######### marginal model for OS ##########


    theta_os.3[i] <- exp(inprod(Z_2[i, 1:Npar], beta_os.3[1:Npar]))

    Y.os.nl_rep[i] ~ dweib(alpha_3b[i], theta_os.3[i]) T(t.cen.lim[i, 1], longest_survival)


  }


  ################## priors  ########################



  for (i in 1:Npar)
  {
    beta_3[i] ~ dnorm(0.0, beta.prior.sd)
  }

  for (i in 1:Npar)
  {
    beta_os.3[i] ~ dnorm(0.0, beta.prior.sd)
  }

  eta_3 ~ dexp(eta.prior)


  alpha_3 ~ dexp(alpha.prior)
  alpha_os.3 ~ dexp(alpha.prior)


  sigma_os.3 ~ dnorm(0, b.prior.sd) T(0,)
  sigma_nl.3 ~ dnorm(0, b.prior.sd) T(0,)




  ################# weibull OS #####################

  for(i in 1:N) {

    theta_os.4[i] <- exp(inprod(Z_2[i, 1:Npar], beta_os.4[1:Npar]))

    Y.os.wb_rep[i] ~ dweib(alpha_os.4, theta_os.4[i]) T(t.cen.lim.wb.os[i, 1], longest_survival)

  }


  for (i in TL_ob_n){

   t.os.wb[i] ~ dweib(alpha_os.4, theta_os.4[i]) T(, longest_survival)

  }

  for (i in TL_c_n){

    t.os.wb[i] ~ dweib(alpha_os.4, theta_os.4[i]) T(t.cen.lim.wb.os[i, 1], longest_survival)

  }

  #### priors ####

  alpha_os.4 ~ dexp(alpha.prior)

  for (i in 1:Npar)
  {
    beta_os.4[i] ~ dnorm(0.0, beta.prior.sd)
  }


}
"
  )
}

ospred_4JM_text <- function() {
  return(
    "################ Target Lesion #########################################
data
{
   for (i in 1:N)
   {
       Y_NLOS[i]<- 0
       Y_NTOS[i]<- 0
   }
}

model{

  C <- 1e+4

  tau <- pow(sigma_Y, -2)

  tau_os.1 <- pow(sigma_os.1, -2)

  tau_b1 <- pow(sigma_b1, -2)

  tau_b2 <- pow(sigma_b2, -2)

  tau_b2.b1 <- 1 / (1 - pow(rho, 2)) * tau_b2

  for (i in 1:N) {


    for(j in 1:T[i]) {


      Y[i, j] ~ dnorm(mu_1[i, j], tau)

      mu_1[i, j] <- inprod(X[i, 1:Npar_0], beta_0[1:Npar_0]) + t[i, j] * (beta_0[Npar_0 + 1] + b[i, 2]) + b[i, 1]

    }


    for (j in (T[i] + 1):(evl_times + 1)){

      mu_1[i, j] <- 0;

    }


    mu[i] <- sum(mu_1[i,]) / T[i]



    theta_os.1[i] <- exp(inprod(Z_1[i, 1:Npar], beta_1[1:Npar]) + lambda * mu[i] + b_os.1[i])

    Y.os.t_rep[i] ~ dweib(alpha_os.1, theta_os.1[i]) T(t.cen.lim[i, 1], longest_survival)


    b[i,1] ~ dnorm(0, tau_b1)

    b[i,2] ~ dnorm(mu_b2.b1[i], tau_b2.b1)

    mu_b2.b1[i] <- sigma_b2 * b[i,1] * rho / sigma_b1

    b_os.1[i] ~ dnorm(0, tau_os.1)

  }

  for (i in TL_ob_n){

    t.os[i] ~ dweib(alpha_os.1, theta_os.1[i]) T(, longest_survival)

  }

  for (i in TL_c_n){

    t.os[i] ~ dweib(alpha_os.1, theta_os.1[i]) T(t.cen.lim[i, 1], longest_survival)
  }




  ######################### Priors ###########################################


  for (i in 1:(Npar_0 + 1)) {
    beta_0[i] ~ dnorm(0.0, beta.prior.sd)
  }


  for (i in 1:Npar) {
    beta_1[i] ~ dnorm(0.0, beta.prior.sd)
  }


  alpha_os.1 ~ dexp(alpha.prior)

  lambda ~ dnorm(0, lambda.prior.sd)



  sigma_Y ~ dnorm(0, y.prior.sd) T(0,)

  sigma_b1 ~ dnorm(0, b.prior.sd) T(0,)

  sigma_b2 ~ dnorm(0, b.prior.sd) T(0,)

  sigma_os.1 ~ dnorm(0, b.prior.sd) T(0,)

  rho ~ dunif(-1, 1)


  ############################  Non target   #########################################

  tau_os.2 <- pow(sigma_os.2, -2)
  tau_nt.2 <- pow(sigma_nt.2, -2)


  for(i in 1:N) {



    NT.u[i] <- -pow(Y_2[i], alpha_NT.2b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_2[1:Npar])) * (-eta_2)
    NT.v[i] <- -pow(Y_os[i], alpha_2b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.2[1:Npar])) * (-eta_2)
    NT.LSE[i] <- log(exp(NT.u[i] - max(NT.u[i], NT.v[i], 0)) + exp(NT.v[i] - max(NT.u[i], NT.v[i], 0)) - exp(0 - max(NT.u[i], NT.v[i], 0))) + max(NT.u[i], NT.v[i], 0)

    L1.NT[i] <- L1.NT1[i] + L1.NT2[i] + L1.NT3[i] + L1.NT4[i]

    L1.NT1[i] <- log(eta_2 + 1) + (-pow(Y_2[i], (alpha_NT.2b[i])) * exp(inprod(Z_2[i, 1:Npar], beta_2[1:Npar])) - pow(Y_os[i], alpha_2b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.2[1:Npar]))) * (-(eta_2+1))
    L1.NT2[i] <- (-(2 * eta_2 + 1)/(eta_2)) *  NT.LSE[i]
    L1.NT3[i] <- log(alpha_NT.2b[i] * pow(Y_2[i], alpha_NT.2b[i] - 1)) + (inprod(Z_2[i, 1:Npar], beta_2[1:Npar]) - pow(Y_2[i], alpha_NT.2b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_2[1:Npar])))
    L1.NT4[i] <- log(alpha_2b[i] * pow(Y_os[i], alpha_2b[i] - 1)) + (inprod(Z_2[i, 1:Npar], beta_os.2[1:Npar]) - pow(Y_os[i], alpha_2b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.2[1:Npar])))


    L2.NT[i] <- L2.NT1[i] + L2.NT2[i] + L2.NT3[i]

    L2.NT1[i] <- (-pow(Y_2[i], alpha_NT.2b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_2[1:Npar]))) * (eta_2 + 1)
    L2.NT2[i] <- (-(eta_2 + 1)/(eta_2)) *  NT.LSE[i]
    L2.NT3[i] <- L1.NT3[i]


    L3.NT[i] <- L3.NT1[i] + L3.NT2[i] + L3.NT3[i]

    L3.NT1[i] <- (-pow(Y_os[i], alpha_2b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.2[1:Npar]))) * (eta_2 + 1)
    L3.NT2[i] <- (-(eta_2 + 1)/(eta_2)) *  NT.LSE[i]
    L3.NT3[i] <- L1.NT4[i]


    L4.NT[i] <- (-1/eta_2) *  NT.LSE[i]


    L.NT[i] <- L1.NT[i] * (delta_NT[i] * delta_OS[i]) + L2.NT[i] * (delta_NT[i] * (1 - delta_OS[i])) + L3.NT[i] * ((1 - delta_NT[i]) * delta_OS[i]) + L4.NT[i] * ((1 - delta_OS[i]) * (1 - delta_NT[i]))



    phi_NT[i] <- -L.NT[i] + C
    Y_NTOS[i] ~ dpois(phi_NT[i])


    os.2[i] ~ dnorm(0, tau_os.2)
    nt.2[i] ~ dnorm(0, tau_nt.2)

    b_os.2[i] <- exp(os.2[i])
    b_nt.2[i] <- exp(nt.2[i])


    alpha_2b[i] <- alpha_os.2 + b_os.2[i]
    alpha_NT.2b[i] <- alpha_2 + b_nt.2[i]


    ####### marginal model for OS ##########


    theta_os.2[i] <- exp(inprod(Z_2[i, 1:Npar], beta_os.2[1:Npar]))

    Y.os.nt_rep[i] ~ dweib(alpha_2b[i], theta_os.2[i]) T(t.cen.lim[i, 1], longest_survival)


  }


  ################## priors  ########################



  for (i in 1:Npar) {
    beta_2[i] ~ dnorm(0.0, beta.prior.sd)
  }

  for (i in 1:Npar) {
    beta_os.2[i] ~ dnorm(0.0, beta.prior.sd)
  }

  eta_2 ~ dexp(eta.prior)


  alpha_2 ~ dexp(alpha.prior)
  alpha_os.2 ~ dexp(alpha.prior)

  sigma_os.2 ~ dnorm(0, b.prior.sd) T(0,)
  sigma_nt.2 ~ dnorm(0, b.prior.sd) T(0,)


  ############################## New Lesion #######################################

  tau_os.3 <- pow(sigma_os.3, -2)
  tau_nl.3 <- pow(sigma_nl.3, -2)

  for(i in 1:N) {


    NL.u[i] <- -pow(Y_3[i], alpha_NL.3b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_3[1:Npar])) * (-eta_3)
    NL.v[i] <- -pow(Y_os[i], alpha_3b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.3[1:Npar])) * (-eta_3)
    NL.LSE[i] <- log(exp(NL.u[i] - max(NL.u[i], NL.v[i], 0)) + exp(NL.v[i] - max(NL.u[i], NL.v[i], 0)) - exp(0 - max(NL.u[i], NL.v[i], 0))) + max(NL.u[i], NL.v[i], 0)


    L1.NL[i] <- L1.NL1[i] + L1.NL2[i] + L1.NL3[i] + L1.NL4[i]

    L1.NL1[i] <- log(eta_3 + 1) + (-pow(Y_3[i], (alpha_NL.3b[i])) * exp(inprod(Z_2[i, 1:Npar], beta_3[1:Npar])) - pow(Y_os[i], alpha_3b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.3[1:Npar]))) * (-(eta_3+1))
    L1.NL2[i] <- (-(2 * eta_3 + 1)/(eta_3)) * NL.LSE[i]
    L1.NL3[i] <- log(alpha_NL.3b[i] * pow(Y_3[i], alpha_NL.3b[i] - 1)) + (inprod(Z_2[i, 1:Npar], beta_3[1:Npar]) - pow(Y_3[i], alpha_NL.3b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_3[1:Npar])))
    L1.NL4[i] <- log(alpha_3b[i] * pow(Y_os[i], alpha_3b[i] - 1)) + (inprod(Z_2[i, 1:Npar], beta_os.3[1:Npar]) - pow(Y_os[i], alpha_3b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.3[1:Npar])))


    L2.NL[i] <- L2.NL1[i] + L2.NL2[i] + L2.NL3[i]

    L2.NL1[i] <- (-pow(Y_3[i], alpha_NL.3b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_3[1:Npar]))) * (eta_3 + 1)
    L2.NL2[i] <- (-(eta_3 + 1)/(eta_3)) * NL.LSE[i]
    L2.NL3[i] <- L1.NL3[i]


    L3.NL[i] <- L3.NL1[i] + L3.NL2[i] + L3.NL3[i]

    L3.NL1[i] <- (-pow(Y_os[i], alpha_3b[i]) * exp(inprod(Z_2[i, 1:Npar], beta_os.3[1:Npar]))) * (eta_3 + 1)
    L3.NL2[i] <- (-(eta_3 + 1)/(eta_3)) * NL.LSE[i]
    L3.NL3[i] <- L1.NL4[i]


    L4.NL[i] <- (-1/eta_3) * NL.LSE[i]


    L.NL[i] <- L1.NL[i] * (delta_NL[i] * delta_OS[i]) + L2.NL[i] * (delta_NL[i] * (1 - delta_OS[i])) + L3.NL[i] * ((1 - delta_NL[i]) * delta_OS[i]) + L4.NL[i] * ((1 - delta_OS[i]) * (1 - delta_NL[i]))



    phi_NL[i] <- -L.NL[i] + C
    Y_NLOS[i] ~ dpois(phi_NL[i])

    os.3[i] ~ dnorm(0, tau_os.3)
    nl.3[i] ~ dnorm(0, tau_nl.3)

    b_os.3[i] <- exp(os.3[i])
    b_nl.3[i] <- exp(nl.3[i])

    alpha_3b[i] <- alpha_os.3 + b_os.3[i]
    alpha_NL.3b[i] <- alpha_3 + b_nl.3[i]

    ######### marginal model for OS ##########


    theta_os.3[i] <- exp(inprod(Z_2[i, 1:Npar], beta_os.3[1:Npar]))

    Y.os.nl_rep[i] ~ dweib(alpha_3b[i], theta_os.3[i]) T(t.cen.lim[i, 1], longest_survival)


  }


  ################## priors  ########################



  for (i in 1:Npar)
  {
    beta_3[i] ~ dnorm(0.0, beta.prior.sd)
  }

  for (i in 1:Npar)
  {
    beta_os.3[i] ~ dnorm(0.0, beta.prior.sd)
  }

  eta_3 ~ dexp(eta.prior)


  alpha_3 ~ dexp(alpha.prior)
  alpha_os.3 ~ dexp(alpha.prior)


  sigma_os.3 ~ dnorm(0, b.prior.sd) T(0,)
  sigma_nl.3 ~ dnorm(0, b.prior.sd) T(0,)




  ################# weibull OS #####################

  for(i in 1:N) {

    theta_os.4[i] <- exp(inprod(Z_2[i, 1:Npar], beta_os.4[1:Npar]))

    Y.os.wb_rep[i] ~ dweib(alpha_os.4, theta_os.4[i]) T(t.cen.lim.wb.os[i, 1], longest_survival)

  }


  for (i in TL_ob_n){

   t.os.wb[i] ~ dweib(alpha_os.4, theta_os.4[i]) T(, longest_survival)

  }

  for (i in TL_c_n){

    t.os.wb[i] ~ dweib(alpha_os.4, theta_os.4[i]) T(t.cen.lim.wb.os[i, 1], longest_survival)

  }

  #### priors ####

  alpha_os.4 ~ dexp(alpha.prior)

  for (i in 1:Npar)
  {
    beta_os.4[i] ~ dnorm(0.0, beta.prior.sd)
  }


  ################ model averaging #############################

  for (i in 1:N) {

    Y.mix_rep[i] <- Y.os.t_rep[i] * w[1] + Y.os.nt_rep[i] * w[2] + Y.os.nl_rep[i] * w[3] + Y.os.wb_rep[i] * w[4]

    ### likelihoods of 4 joint models ####


    LL[i,1] <- delta_OS[i] * (log(alpha_os.1) + log(theta_os.1[i]) + (alpha_os.1 - 1) * log(Y_os[i])) - theta_os.1[i] * pow(Y_os[i], alpha_os.1)

    LL[i,2] <- delta_OS[i] * (log(alpha_2b[i]) + log(theta_os.2[i]) + (alpha_2b[i] - 1) * log(Y_os[i])) - theta_os.2[i] * pow(Y_os[i], alpha_2b[i])

    LL[i,3] <- delta_OS[i] * (log(alpha_3b[i]) + log(theta_os.3[i]) + (alpha_3b[i] - 1) * log(Y_os[i])) - theta_os.3[i] * pow(Y_os[i], alpha_3b[i])

    LL[i,4] <- delta_OS[i] * (log(alpha_os.4) + log(theta_os.4[i]) + (alpha_os.4 - 1) * log(Y_os[i])) - theta_os.4[i] * pow(Y_os[i], alpha_os.4)

  }

  #### probability of parameters given model j (j = 1:4) ####

  #### model 1 ####

  for (j in 1:(Npar_0 + 7 + Npar) ){

    pg1[j] <- 1 / pow(seg1[j], 2)

  }


  for (j in 1:(Npar_0 + 1)){

    Pr[j, 1] <- 0.5 * log(beta.prior.sd / 6.28) - 0.5 * beta.prior.sd * pow(beta_0[j], 2)

    g[j, 1] <- 0.5 * log(pg1[j] / 6.28) - 0.5 * pg1[j] * pow(beta_0[j] - mg1[j], 2)

  }

  for (j in (Npar_0 + 2):(Npar_0 + 1 + Npar)  ){

    Pr[j, 1] <- 0.5 * log(beta.prior.sd / 6.28) - 0.5 * beta.prior.sd * pow(beta_1[j - Npar_0 - 1], 2)

    g[j, 1] <- 0.5 * log(pg1[j] / 6.28) - 0.5 * pg1[j] * pow(beta_1[j - Npar_0 - 1] - mg1[j], 2)

  }



  Pr[Npar_0 + 2 + Npar, 1] <- -alpha_os.1 * alpha.prior + log(alpha.prior)

  g[Npar_0 + 2 + Npar, 1] <- -alpha_os.1 * mg1[Npar_0 + 2 + Npar] + log(mg1[Npar_0 + 2 + Npar])

  Pr[Npar_0 + 3 + Npar, 1] <- 0.5 * log(lambda.prior.sd / 6.28) - 0.5 * lambda.prior.sd * pow(lambda, 2)

  g[Npar_0 + 3 + Npar, 1] <- 0.5 * log(pg1[Npar_0 + 3 + Npar] / 6.28) - 0.5 * pg1[Npar_0 + 3 + Npar] * pow(lambda - mg1[Npar_0 + 3 + Npar], 2)

  Pr[Npar_0 + 4 + Npar, 1] <- 0.5 * log(b.prior.sd * 4 / 6.28) - 0.5 * b.prior.sd * pow(sigma_b1, 2)

  g[Npar_0 + 4 + Npar, 1] <- 0.5 * log(pg1[Npar_0 + 4 + Npar] * 4 / 6.28) - 0.5 * pg1[Npar_0 + 4 + Npar] * pow(sigma_b1 - mg1[Npar_0 + 4 + Npar], 2)

  Pr[Npar_0 + 5 + Npar, 1] <- 0.5 * log(b.prior.sd * 4 / 6.28) - 0.5 * b.prior.sd * pow(sigma_b2, 2)

  g[Npar_0 + 5 + Npar, 1] <- 0.5 * log(pg1[Npar_0 + 5 + Npar] * 4/ 6.28) - 0.5 * pg1[Npar_0 + 5 + Npar] * pow(sigma_b2 - mg1[Npar_0 + 5 + Npar], 2)

  Pr[Npar_0 + 6 + Npar, 1] <- log(0.5)

  g[Npar_0 + 6 + Npar, 1] <- log(1 / (seg1[Npar_0 + 6 + Npar] * sqrt(12)))

  Pr[Npar_0 + 7 + Npar, 1] <- 0.5 * log(b.prior.sd * 4/ 6.28) - 0.5 * b.prior.sd * pow(sigma_os.1, 2)

  g[Npar_0 + 7 + Npar, 1] <- 0.5 * log(pg1[Npar_0 + 7 + Npar] * 4/ 6.28) - 0.5 * pg1[Npar_0 + 7 + Npar] * pow(sigma_os.1 - mg1[Npar_0 + 7 + Npar], 2)


  #### model 2 ####

  for (j in 1:(Npar + 2)){

    pg2[j] <- 1 / pow(seg2[j], 2)

  }


  for (j in 1:Npar){

    Pr[j,2] <- 0.5 * log(beta.prior.sd / 6.28) - 0.5 * beta.prior.sd * pow(beta_os.2[j], 2)

    g[j,2] <- 0.5 * log(pg2[j] / 6.28) - 0.5 * pg2[j] * pow(beta_os.2[j] - mg2[j], 2)

  }



  Pr[Npar + 1, 2] <- -alpha_os.2 * alpha.prior + log(alpha.prior)

  g[Npar + 1, 2] <- -alpha_os.2 * mg2[Npar + 1] + log(mg2[Npar + 1])


  Pr[Npar + 2, 2] <- 0.5 * log(b.prior.sd / 6.28) - 0.5 * b.prior.sd * pow(sigma_os.2, 2)

  g[Npar + 2, 2] <- 0.5 * log(pg2[Npar + 2] / 6.28) - 0.5 * pg2[Npar + 2] * pow(sigma_os.2 - mg2[Npar + 2], 2)


  for (j in (Npar + 3):(Npar_0 + 7 + Npar)){

    Pr[j, 2] <- 0

    g[j, 2] <- 0

  }



  #### model 3 ####

  for (j in 1:(Npar + 2)){

    pg3[j] <- 1 / pow(seg3[j], 2)

  }

  for (j in 1:Npar){

    Pr[j, 3] <- 0.5 * log(beta.prior.sd / 6.28) - 0.5 * beta.prior.sd * pow(beta_os.3[j], 2)

    g[j, 3] <- 0.5 * log(pg3[j] / 6.28) - 0.5 * pg3[j] * pow(beta_os.3[j] - mg3[j], 2)

  }


  Pr[Npar + 1, 3] <- -alpha_os.3 * alpha.prior + log(alpha.prior)

  g[Npar + 1, 3] <- -alpha_os.3 * mg3[Npar + 1] + log(mg3[Npar + 1])


  Pr[Npar + 2, 3] <- 0.5 * log(b.prior.sd / 6.28) - 0.5 * b.prior.sd * pow(sigma_os.3, 2)

  g[Npar + 2, 3] <- 0.5 * log(pg3[Npar + 2] / 6.28) - 0.5 * pg3[Npar + 2] * pow(sigma_os.3 - mg3[Npar + 2], 2)


  for (j in (Npar + 3):(Npar_0 + 7 + Npar)){

    Pr[j, 3] <- 0

    g[j,3] <- 0

  }


    #### model 4 ####

  for (j in 1:(Npar + 1) ){

    pg4[j] <- 1 / pow(seg4[j], 2)

  }

  for (j in 1:Npar){

    Pr[j, 4] <- 0.5 * log(beta.prior.sd / 6.28) - 0.5 * beta.prior.sd * pow(beta_os.4[j], 2)

    g[j, 4] <- 0.5 * log(pg4[j] / 6.28) - 0.5 * pg4[j] * pow(beta_os.4[j] - mg4[j], 2)

  }


  Pr[Npar + 1, 4] <- -alpha_os.4 * alpha.prior + log(alpha.prior)

  g[Npar + 1, 4] <- -alpha_os.4 * mg4[Npar + 1] + log(mg4[Npar + 1])



  for (j in (Npar + 2):(Npar_0 + 7 + Npar)){

    Pr[j, 4] <- 0

    g[j,4] <- 0

  }




  #### total number of parameters in each model #####


  for (k in 1:4){

    PriorMod[k] <- 1/4

  }



  for (k in 1:4){

    SL[k] <- max(logR[k] - maxR, -500)

    log.Pstar[k] <- sum(LL[, k]) + sum(Pr[1:p[k], k]) + log(PriorMod[k])

    logR[k] <- log.Pstar[k] - sum(g[1:p[k], k])

    expSL[k] <- exp(SL[k])

    w[k] <- expSL[k] / sum(expSL[])


  }


  logR.sorted <- sort(logR[])
  maxR <- logR.sorted[4]



}
")
}

