library(boot)

#==============
# Basic values
#==============
nsite <- 60   # number of sites
nyear <- 20   # number of years
nreps <- 3    # number of replicates
ncovs <- 1    # number of environmental covariates
npcvs <- 2    # number of observational covariates
nspec <- 2    # number of species

# intercept and slopes for initial population size
beta0 <- matrix(c(5, 50, .5, .5), nspec, ncovs+1)
# standard deviation of stochasticity in initial population size
sigma0 <- rep(.1, nspec)
# intercept and slopes for local population growth rate
beta_phi <- matrix(c(0, 0, -.5, -.5, .5, -.5, 0, .5), nspec, nspec+ncovs+1)
# decay parameter for movement
kappa <- c(2, 2)
# standard deviation of stochasticity in population growth
sigma <- rep(.1, nspec)
# intercept and slopes for detection probability
beta_pobs <- matrix(c(logit(.67), logit(.67), .5, .5), nspec, npcvs+1)
# standard deviation of stochasticity in detection probability
sigma_pobs <- rep(.1, nspec)

#===============
# Simulate data
#===============
### Site locations
lat <- runif(nsite, 0, 5)
lon <- runif(nsite, 0, 5)

### Distance between sites
dist <- matrix(, nsite, nsite)
for (i in 1:nsite) {
  for (j in 1:nsite) {
    dist[i,j] <- sqrt((lat[i] - lat[j]) ^ 2 + (lon[i] - lon[j]) ^ 2)
  } # j
} # i

### Environmental covariates
x_mean <- matrix(rnorm(nsite * ncovs, 0, 1), nsite, ncovs)
x <- array(, dim=c(nsite, nyear, ncovs))
for (i in 1:nsite) {
  for (s in 1:ncovs) {
    x[i,,s] <- rnorm(nyear, x_mean[i,s], .2)
  } # s
} # i

### Metapopulation dynamics
lambda0 <- matrix(, nsite, nspec)
for (s in 1:nspec) {
  lambda0[,s] <- exp(cbind(1, x[,1,]) %*% beta0[s,] + rnorm(nsite, 0, sigma0[s]))
} # s

N <- array(, dim=c(nsite, nyear, nspec))
for (s in 1:nspec) {
  N[,1,s] <- rpois(nsite, lambda0[,s])
} # s

eta <- theta <- array(, dim=c(nsite, nsite, nspec))
for (s in 1:nspec) {
  eta[,,s] <- exp(-1 * kappa[s] * dist)
  theta[,,s] <- eta[,,s] / rowSums(eta[,,s])
} # s

phi <- lambda <- array(, dim=c(nsite, nyear-1, nspec))
for (t in 2:nyear) {
  for (s in 1:nspec) {
    phi[,t-1,s] <- exp(cbind(1, (N[,t-1,]- lambda0) / lambda0, x[,t,] - x[,1,]) %*% beta_phi[s,])
    lambda[,t-1,s] <- ((N[,t-1,s] * phi[,t-1,s]) %*% theta[,,s]) * exp(rnorm(nsite, 0, sigma[s]))
    N[,t,s] <- rpois(nsite, lambda[,t-1,s])
  } # s
} # t

### Count data
w <- array(rnorm(nsite * nyear * npcvs * nreps, 0, 1), dim=c(nsite, nyear, npcvs, nreps))
pobs <- y <- array(0, dim=c(nsite, nyear, nspec, nreps))
for (i in 1:nsite) {
  for (t in 1:nyear) {
    for (s in 1:nspec) {
      pobs[i,t,s,] <- inv.logit(cbind(1, t(w[i,t,,])) %*% beta_pobs[s,] + rnorm(nreps, 0, sigma_pobs[s]))
      y[i,t,s,] <- rbinom(nreps, N[i,t,s], pobs[i,t,s,])
    } # s
  } # t
} # i

#=======================
# Define MCMC algorithm
#=======================
# MCMC function
sdnm_mcmc <- function(y, x, w, dist, nmcmc) {

  # Setup variables
  nsite <- dim(y)[1]
  nyear <- dim(y)[2]
  nspec <- dim(y)[3]
  nreps <- dim(y)[4]
  ncovs <- dim(x)[3]
  npcvs <- dim(w)[3]
  ymax <- apply(y, 1:3, max)

  beta0_save <- array(0, dim=c(nspec, ncovs+1, nmcmc))
  sigma0_save <- matrix(0, nspec, nmcmc)
  beta_phi_save <- array(0, dim=c(nspec, nspec+ncovs+1, nmcmc))
  kappa_save <- matrix(0, nspec, nmcmc)
  sigma_save <- matrix(0, nspec, nmcmc)
  beta_pobs_save <- array(0, dim=c(nspec, npcvs+1, nmcmc))
  sigma_pobs_save <- matrix(0, nspec, nmcmc)
  N_save <- array(0, dim=c(nsite, nyear, nspec, nmcmc))

  # Priors
  beta0_mean <- matrix(0, nspec, ncovs+1)
  beta0_sd <- sqrt(1000)
  log_sigma0_mean <- 0
  log_sigma0_sd <- sqrt(1000)
  beta_phi_mean <- matrix(0, nspec, nspec+ncovs+1)
  beta_phi_sd <- sqrt(1000)
  log_kappa_mean <- 0
  log_kappa_sd <- sqrt(1000)
  log_sigma_mean <- 0
  log_sigma_sd <- sqrt(1000)
  beta_pobs_mean <- matrix(0, nspec, npcvs+1)
  beta_pobs_sd <- sqrt(1000)
  log_sigma_pobs_mean <- 0
  log_sigma_pobs_sd <- sqrt(1000)

  # Starting values
  N <- round((ymax + 1) / .5)
  beta0 <- matrix(0, nspec, ncovs+1)
  sigma0 <- rep(1, nspec)
  epsilon0 <- matrix(0, nsite, nspec)
  beta_phi <- matrix(0, nspec, nspec+ncovs+1)
  phi <- array(1, dim=c(nsite, nyear-1, nspec))
  kappa <- rep(1, nspec)
  eta <- theta <- array(, dim=c(nsite, nsite, nspec))
  for (s in 1:nspec) {
    eta[,,s] <- exp(-1 * kappa[s] * dist)
    theta[,,s] <- eta[,,s] / rowSums(eta[,,s])
  } # s
  sigma <- rep(1, nspec)
  epsilon <- array(0, dim=c(nsite, nyear-1, nspec))
  beta_pobs <- matrix(0, nspec, npcvs+1)
  sigma_pobs <- rep(1, nspec)
  epsilon_pobs <- array(0, dim=c(nsite, nyear, nspec, nreps))
  pobs <- array(.5, dim=c(nsite, nyear, nspec, nreps))

  # Tuning factors
  N_tune <- 1
  beta0_tune <- matrix(.035, nspec, ncovs+1)
  sigma0_tune <- .1
  epsilon0_tune <- .05
  beta_phi_tune <- matrix(c(.05, .1), nspec, nspec+ncovs+1)
  kappa_tune <- .1
  sigma_tune <- .1
  epsilon_tune <- .05
  beta_pobs_tune <- matrix(c(.02, .03), nspec, npcvs+1)
  sigma_pobs_tune <- .1
  epsilon_pobs_tune <- .05

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample N
    N_star <- array(rpois(nsite * nyear * nspec, N + N_tune), dim=c(nsite, nyear, nspec))
    lambda0 <- matrix(, nsite, nspec)
    mh1 <- mh2 <- array(, dim=c(nsite, nyear, nspec))
    for (s in 1:nspec) {
      lambda0[,s] <- exp(cbind(1, x[,1,]) %*% beta0[s,] + epsilon0[,s])
      mh1[,,s] <- apply(dbinom(y[,,s,], N_star[,,s], pobs[,,s,], log=T), 1:2, sum) + 
                  dpois(N_star[,,s], cbind(lambda0[,s], lambda[,,s]), log=T) + 
                  dpois(N[,,s], N_star[,,s] + N_tune, log=T)
      mh2[,,s] <- apply(dbinom(y[,,s,], N[,,s], pobs[,,s,], log=T), 1:2, sum) + 
                  dpois(N[,,s], cbind(lambda0[,s], lambda[,,s]), log=T) + 
                  dpois(N_star[,,s], N[,,s] + N_tune, log=T)
    } # s
    mh <- exp(mh1 - mh2)
    Nkeep <- ((mh > runif(nsite*nyear*nspec)) & (N_star >= ymax))
    N[Nkeep] <- N_star[Nkeep]

    ### Sample beta0
    beta0_star <- matrix(rnorm(nspec*(ncovs+1), beta0, beta0_tune), nspec, ncovs+1)
    lambda0 <- lambda0_star <- matrix(, nsite, nspec)
    for (s in 1:nspec) {
      lambda0[,s] <- exp(cbind(1, x[,1,]) %*% beta0[s,] + epsilon0[,s])
      lambda0_star[,s] <- exp(cbind(1, x[,1,]) %*% beta0_star[s,] + epsilon0[,s])
      mh1 <- sum(dpois(N[,1,s], lambda0_star[,s], log=T)) + 
             sum(dnorm(beta0_star[s,], beta0_mean[s,], beta0_sd, log=T))
      mh2 <- sum(dpois(N[,1,s], lambda0[,s], log=T)) + 
             sum(dnorm(beta0[s,], beta0_mean[s,], beta0_sd, log=T))
      mh <- exp(mh1 - mh2)
      if (mh > runif(1)) {
        beta0[s,] <- beta0_star[s,]
      }
    } # s

    # Sample sigma0
    sigma0_star <- exp(rnorm(nspec, log(sigma0), sigma0_tune))
    for (s in 1:nspec) {
      mh1 <- sum(dnorm(epsilon0[,s], 0, sigma0_star[s], log=TRUE)) + 
             dnorm(log(sigma0_star[s]), log_sigma0_mean, log_sigma0_sd, log=TRUE)
      mh2 <- sum(dnorm(epsilon0[,s], 0, sigma0[s], log=TRUE)) + 
             dnorm(log(sigma0[s]), log_sigma0_mean, log_sigma0_sd, log=TRUE)
      mh <- exp(mh1 - mh2)
      if (mh > runif(1)) {
        sigma0[s] <- sigma0_star[s]
      }
    } # s

    # Sample epsilon0
    epsilon0_star <- matrix(rnorm(nsite*nspec, epsilon0, epsilon0_tune), nsite, nspec)
    lambda0 <- lambda0_star <- matrix(, nsite, nspec)
    for (s in 1:nspec) {
      lambda0[,s] <- exp(cbind(1, x[,1,]) %*% beta0[s,] + epsilon0[,s])
      lambda0_star[,s] <- exp(cbind(1, x[,1,]) %*% beta0[s,] + epsilon0_star[,s])
    } # s
    mh1 <- dpois(N[,1,], lambda0_star, log=TRUE) + 
           dnorm(epsilon0_star, 0, sigma0, log=TRUE)
    mh2 <- dpois(N[,1,], lambda0, log=TRUE) + 
           dnorm(epsilon0, 0, sigma0, log=TRUE)
    mh <- exp(mh1 - mh2)
    tt <- matrix(runif(nsite*nspec), nsite, nspec)
    epsilon0[which(mh > tt)] <- epsilon0_star[which(mh > tt)]

    ### Sample beta_phi
    beta_phi_star <- matrix(rnorm(nspec*(nspec+ncovs+1), beta_phi, beta_phi_tune), nspec, nspec+ncovs+1)
    phi_star <- lambda <- lambda_star <- array(, dim=c(nsite, nyear-1, nspec))
    for (s in 1:nspec) {
      for (t in 2:nyear) {
        phi_star[,t-1,s] <- exp(cbind(1, (N[,t-1,]- lambda0) / lambda0, x[,t,] - x[,1,]) %*% beta_phi_star[s,])
      } # t
      lambda[,,s] <- t(t(N[,-nyear,s] * phi[,,s]) %*% theta[,,s]) * exp(epsilon[,t-1,s])
      lambda_star[,,s] <- t(t(N[,-nyear,s] * phi_star[,,s]) %*% theta[,,s]) * exp(epsilon[,t-1,s])
      mh1 <- sum(dpois(N[,-1,s], lambda_star[,,s], log=T)) + 
             sum(dnorm(beta_phi_star[s,], beta_phi_mean[s,], beta_phi_sd, log=T))
      mh2 <- sum(dpois(N[,-1,s], lambda[,,s], log=T)) + 
             sum(dnorm(beta_phi[s,], beta_phi_mean[s,], beta_phi_sd, log=T))
      mh <- exp(mh1 - mh2)
      if (mh > runif(1)) {
        beta_phi[s,] <- beta_phi_star[s,]
        phi[,,s] <- phi_star[,,s]
      }
    } # s

    ### Sample kappa
    kappa_star <- exp(rnorm(nspec, log(kappa), kappa_tune))
    eta_star <- theta_star <- array(, dim=c(nsite, nsite, nspec))
    lambda <- lambda_star <- array(, dim=c(nsite, nyear-1, nspec))
    for (s in 1:nspec) {
      eta_star[,,s] <- exp(-1 * kappa_star[s] * dist)
      theta_star[,,s] <- eta_star[,,s] / rowSums(eta_star[,,s])
      lambda[,,s] <- t(t(N[,-nyear,s] * phi[,,s]) %*% theta[,,s]) * exp(epsilon[,,s])
      lambda_star[,,s] <- t(t(N[,-nyear,s] * phi[,,s]) %*% theta_star[,,s]) * exp(epsilon[,,s])
    } # s
    mh1 <- sum(dpois(N[,-1,], lambda_star, log=T)) + 
           sum(dnorm(log(kappa_star), log_kappa_mean, log_kappa_sd, log=T))
    mh2 <- sum(dpois(N[,-1,], lambda, log=T)) + 
           sum(dnorm(log(kappa), log_kappa_mean, log_kappa_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      kappa <- kappa_star
      eta <- eta_star
      theta <- theta_star
    }

    # Sample sigma
    sigma_star <- exp(rnorm(nspec, log(sigma), sigma_tune))
    for (s in 1:nspec) {
      mh1 <- sum(dnorm(epsilon[,,s], 0, sigma_star[s], log=TRUE)) + 
             dnorm(log(sigma_star[s]), log_sigma_mean, log_sigma_sd, log=TRUE)
      mh2 <- sum(dnorm(epsilon[,,s], 0, sigma[s], log=TRUE)) + 
             dnorm(log(sigma[s]), log_sigma_mean, log_sigma_sd, log=TRUE)
      mh <- exp(mh1 - mh2)
      if (mh > runif(1)) {
        sigma[s] <- sigma_star[s]
      }
    } # s

    # Sample epsilon
    epsilon_star <- array(rnorm(nsite*(nyear-1)*nspec, epsilon, epsilon_tune), dim=c(nsite, nyear-1, nspec))
    lambda <- lambda_star <- array(, dim=c(nsite, nyear-1, nspec))
    for (s in 1:nspec) {
      lambda[,,s] <- t(t(N[,-nyear,s] * phi[,,s]) %*% theta[,,s]) * exp(epsilon[,,s])
      lambda_star[,,s] <- t(t(N[,-nyear,s] * phi[,,s]) %*% theta[,,s]) * exp(epsilon_star[,,s])
    } # s
    mh1 <- mh2 <- array(, dim=c(nsite, nyear-1, nspec))
    for (s in 1:nspec) {
      mh1[,,s] <- dpois(N[,-1,s], lambda_star[,,s], log=TRUE) + 
                  dnorm(epsilon_star[,,s], 0, sigma[s], log=TRUE)
      mh2[,,s] <- dpois(N[,-1,s], lambda[,,s], log=TRUE) + 
                  dnorm(epsilon[,,s], 0, sigma[s], log=TRUE)
    } # s
    mh <- exp(mh1 - mh2)
    tt <- array(runif(nsite*(nyear-1)*nspec), dim=c(nsite, nyear-1, nspec))
    epsilon[which(mh > tt)] <- epsilon_star[which(mh > tt)]

    ### Sample beta_pobs
    beta_pobs_star <- matrix(rnorm(nspec*(npcvs+1), beta_pobs, beta_pobs_tune), nspec, npcvs+1)
    pobs <- pobs_star <- array(0, dim=c(nsite,nyear,nspec,nreps))
    for (i in 1:nsite) {
      for (t in 1:nyear) {
        for (s in 1:nspec) {
          pobs[i,t,s,] <- inv.logit(cbind(1, t(w[i,t,,])) %*% beta_pobs[s,] + epsilon_pobs[i,t,s,])
          pobs_star[i,t,s,] <- inv.logit(cbind(1, t(w[i,t,,])) %*% beta_pobs_star[s,] + epsilon_pobs[i,t,s,])
        } # s
      } # t
    } # i
    for (s in 1:nspec) {
      mh1 <- sum(dbinom(y[,,s,], N[,,s], pobs_star[,,s,], log=T)) + 
             sum(dnorm(beta_pobs_star[s,], beta_pobs_mean[s,], beta0_sd, log=T))
      mh2 <- sum(dbinom(y[,,s,], N[,,s], pobs[,,s,], log=T)) + 
             sum(dnorm(beta_pobs[s,], beta_pobs_mean[s,], beta0_sd, log=T))
      mh <- exp(mh1 - mh2)
      if (mh > runif(1)) {
        beta_pobs[s,] <- beta_pobs_star[s,]
      }
    } # s

    # Sample sigma_pobs
    sigma_pobs_star <- exp(rnorm(nspec, log(sigma_pobs), sigma_pobs_tune))
    for (s in 1:nspec) {
      mh1 <- sum(dnorm(epsilon_pobs[,,s,], 0, sigma_pobs_star[s], log=TRUE)) + 
             dnorm(log(sigma_pobs_star[s]), log_sigma_pobs_mean, log_sigma_pobs_sd, log=TRUE)
      mh2 <- sum(dnorm(epsilon_pobs[,,s,], 0, sigma_pobs[s], log=TRUE)) + 
             dnorm(log(sigma_pobs[s]), log_sigma_pobs_mean, log_sigma_pobs_sd, log=TRUE)
      mh <- exp(mh1 - mh2)
      if (mh > runif(1)) {
        sigma_pobs[s] <- sigma_pobs_star[s]
      }
    } # s

    # Sample epsilon_pobs
    epsilon_pobs_star <- array(rnorm(nsite*nyear*nspec*nreps, epsilon_pobs, epsilon_pobs_tune), dim=c(nsite, nyear, nspec, nreps))
    pobs <- pobs_star <- array(0, dim=c(nsite,nyear,nspec,nreps))
    for (i in 1:nsite) {
      for (t in 1:nyear) {
        for (s in 1:nspec) {
          pobs[i,t,s,] <- inv.logit(cbind(1, t(w[i,t,,])) %*% beta_pobs[s,] + epsilon_pobs[i,t,s,])
          pobs_star[i,t,s,] <- inv.logit(cbind(1, t(w[i,t,,])) %*% beta_pobs[s,] + epsilon_pobs_star[i,t,s,])
        } # s
      } # t
    } # i
    mh1 <- mh2 <- array(, dim=c(nsite, nyear, nspec, nreps))
    for (s in 1:nspec) {
      mh1[,,s,] <- dbinom(y[,,s,], N[,,s], pobs_star[,,s,], log=TRUE) + 
                   dnorm(epsilon_pobs_star[,,s,], 0, sigma[s], log=TRUE)
      mh2[,,s,] <- dbinom(y[,,s,], N[,,s], pobs[,,s,], log=TRUE) + 
                   dnorm(epsilon_pobs[,,s,], 0, sigma[s], log=TRUE)
    } # s
    mh <- exp(mh1 - mh2)
    tt <- array(runif(nsite*nyear*nspec*nreps), dim=c(nsite, nyear, nspec, nreps))
    epsilon_pobs[which(mh > tt)] <- epsilon_pobs_star[which(mh > tt)]

    ### Save samples
    N_save[,,,k] <- N
    beta0_save[,,k] <- beta0
    sigma0_save[,k] <- sigma0
    beta_phi_save[,,k] <- beta_phi
    kappa_save[,k] <- kappa
    sigma_save[,k] <- sigma
    beta_pobs_save[,,k] <- beta_pobs
    sigma_pobs_save[,k] <- sigma_pobs

    setTxtProgressBar(txtProgressBar(min=0, max=nmcmc, style=3, width=50, char="+"), k)
  } # k

  # Write output
  list(N_save=N_save, 
       beta0_save=beta0_save, sigma0_save=sigma0_save, 
       beta_phi_save=beta_phi_save, kappa_save=kappa_save, sigma_save=sigma_save, 
       beta_pobs_save=beta_pobs_save, sigma_pobs_save=sigma_pobs_save)

} # sdnm_mcmc

#==========
# Run MCMC
#==========
nmcmc <- 100000
chain <- 5
out <- list()

for (i in 1:chain) {
  print(i)
  out[[i]] <- sdnm_mcmc(y=y, x=x, w=w, dist=dist, nmcmc=nmcmc)
  print(' ')
} # i

#==============
# Plot results
#==============
par(mfrow=c(6,5))
par(mar=c(1,2,1,1))
for (s in 1:nspec) {
  for (j in 1:(ncovs+1)) {
    plot(out[[1]]$beta0_save[s,j,], type='n')
    for (i in 1:chain) {
      lines(out[[i]]$beta0_save[s,j,], col=i+1)
    } # i
    abline(h=beta0[s,j], col=1)
  } # j

  for (j in 1:(nspec+ncovs+1)) {
    plot(out[[1]]$beta_phi_save[s,j,], type='n')
    for (i in 1:chain) {
      lines(out[[i]]$beta_phi_save[s,j,], col=i+1)
    } # i
    abline(h=beta_phi[s,j], col=1)
  } # k

  plot(out[[1]]$kappa_save[s,], type='n')
  for (i in 1:chain) {
    lines(out[[i]]$kappa_save[s,], col=i+1)
  } # i
  abline(h=kappa[s], col=1)

  for (j in 1:(npcvs+1)) {
    plot(out[[1]]$beta_pobs_save[s,j,], type='n')
    for (i in 1:chain) {
      lines(out[[i]]$beta_pobs_save[s,j,], col=i+1)
    } # i
    abline(h=beta_pobs[s,j], col=1)
  } # j

  plot(out[[1]]$sigma0_save[s,], type='n')
  for (i in 1:chain) {
    lines(out[[i]]$sigma0_save[s,], col=i+1)
  } # i
  abline(h=sigma0[s], col=1)

  plot(out[[1]]$sigma_save[s,], type='n')
  for (i in 1:chain) {
    lines(out[[i]]$sigma_save[s,], col=i+1)
  } # i
  abline(h=sigma[s], col=1)

  plot(out[[1]]$sigma_pobs_save[s,], type='n')
  for (i in 1:chain) {
    lines(out[[i]]$sigma_pobs_save[s,], col=i+1)
  } # i
  abline(h=sigma_pobs[s], col=1)
} # s


