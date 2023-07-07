rm(list= ls())
#'%ni%' <- function(x,y)!('%in%'(x,y))
library(parallel)
setwd("D:/Choennam Univ/FMD_korea/Code")

# 0. Set the seed
set.seed(1)

# 1. Setting-up the parallel computation
detectCores(logical= T)
n_core <- 10
cl <- makeCluster(n_core)
clusterExport(cl, c('n_core'))

# 2. Variables
n_mod <- 5
n_seq <- 12 #8
n_particle <- 200
max_iter <- n_particle / n_core
tuner <- 2.0
q_lim <- 0.5

rho_a <- c(31, 7, 13, 9) # Observed data
rho_b <- c(47, 9, 3, 7) # Observed data
rho_min <- 0.0
rho_max <- 1.0
beta_min <- 0.0
beta_max <- 1.5 # Calculated based on Orsel et al 2007
beta_pmu <- 0.6 # Calculated based on Dekker et al 2020 
beta_psd <- 0.1 # Calculated based on Dekker et al 2020
omega_min <- 0.0
omega_max <- 0.5 # Based on Kinsley et al 2016
psi_min <- 0.0
psi_max <- 1.0
alpha_min <- 0
alpha_max <- 1
gamma1_n <- 97 # Based on Kinsley et al 2016
gamma1_p <- 0.02 # Based on Kinsley et al 2016
gamma2_n <- 66 # Based on Kinsley et al 2016
gamma2_p <- 0.02 # Based on Kinsley et al 2016
gamma1_max <- 7
gamma3_min <- 3 # Based on Parida et al 2007
gamma3_max <- 7 # Based on Parida et al 2007
gamma_min <- 1 # Based on Kinsley et al 2016
gamma_max <- 9 # Based on Kinsley et al 2016

simperiod <- 20
n_farm <- 4
n_pigs <- c(825, 822, 854, 413)
tau_min <- 0
tau_max <- simperiod - 1
obs_ts <- list(c(0,0,0,0,0,0,0,0,0,0,18,17,17,22,5,10,16,1,9,2),
               c(0,0,0,0,0,0,0,0,0,16,47,36,35,24,23,15,3,6,7,3),
               c(0,0,0,0,0,0,0,0,0,0,9,3,4,3,21,5,4,10,20,11),
               c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,15,17,7))

m_l <- NULL

rho_l <- vector(mode= "list", length= n_mod)
tau_l <- vector(mode= "list", length= n_mod)
for (i in 1:n_mod) {
  rho_l[[i]] <- list(NULL,NULL,NULL,NULL)
  tau_l[[i]] <- list(NULL,NULL,NULL,NULL)
}
beta_l <- vector(mode= "list", length= n_mod)
omega_l <- vector(mode= "list", length= n_mod)
psi_l <- vector(mode= "list", length= n_mod)
alpha_l <- vector(mode= "list", length= n_mod)
gamma1_l <- vector(mode= "list", length= n_mod)
gamma2_l <- vector(mode= "list", length= n_mod)
gamma3_l <- vector(mode= "list", length= n_mod)

sim_ts <- rep(0, simperiod)

# 3. ABC-SMC for the initial sequence
t_init <- Sys.time()
diff_01 <- NULL; diff_02 <- NULL; diff_03 <- NULL; diff_04 <- NULL

for (j in 1:max_iter) {
  ms <- sample(c(1:n_mod), n_core, replace= T)
  rhos <- list(
    rbeta(n_core, rho_a[1], rho_b[1]),
    rbeta(n_core, rho_a[2], rho_b[2]),
    rbeta(n_core, rho_a[3], rho_b[3]),
    rbeta(n_core, rho_a[4], rho_b[4])
  )
  betas <- NULL
  while (length(betas) < n_core) {
    v <- rnorm(1, beta_pmu, beta_psd)
    if (v >= beta_min & v <= beta_max) {
      betas <- c(betas, v)
    }
  }
  omegas <- runif(n_core, omega_min, omega_max)
  psis <- runif(n_core, psi_min, psi_max)
  alphas <- runif(n_core, alpha_min, alpha_max)
  gamma1s <- gamma2s <- NULL
  while (length(gamma1s) < n_core) {
    gamma1 <- rbinom(1, gamma1_n, gamma1_p)
    gamma2 <- rbinom(1, gamma2_n, gamma2_p)
    if (gamma1 + gamma2 <= gamma_max & gamma1 <= gamma1_max & gamma1 >= gamma_min & gamma2 >= gamma_min) {
      gamma1s <- c(gamma1s, gamma1)
      gamma2s <- c(gamma2s, gamma2)
    }
  }
  gamma3s <- runif(n_core, gamma3_min, gamma3_max)
  taus <- list(runif(n_core, tau_min, tau_max), 
               runif(n_core, tau_min, tau_max),
               runif(n_core, tau_min, tau_max), 
               runif(n_core, tau_min, tau_max)
  )
  
  fid <- 1
  tmp <- clusterApplyLB(cl, 1:n_core, function(i, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas,
                                               gamma1s, gamma2s, gamma3s, taus, sim_ts) {
    m <- ms[i]
    n_pig <- n_pigs[fid]
    rho <- rhos[[fid]][i]
    beta <- betas[i]
    omega <- omegas[i]
    psi <- psis[i]
    alpha <- alphas[i]
    gamma1 <- round(gamma1s[i], 0)
    gamma2 <- round(gamma2s[i], 0)
    gamma3 <- round(gamma3s[i], 0)
    tau <- round(taus[[fid]][i], 0)
    seed <- ((j - 1) * n_core) + i
    
    dyn.load("FMD_sim.dll")
    v <- .C("_FMD_sim", 
            as.integer(i), as.integer(n_pig), as.integer(m), as.double(rho), as.double(beta), as.double(omega),
            as.double(psi), as.double(alpha), as.integer(gamma1), as.integer(gamma2), 
            as.integer(gamma3), as.integer(tau), as.integer(sim_ts))
    dyn.unload("FMD_sim.dll"); v
  }, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas, gamma1s, gamma2s, gamma3s, taus, sim_ts)
  diff_01 <- c(diff_01, unlist(lapply(lapply(tmp, unlist), function(x) {
    v <- tail(x, simperiod)
    sqrt(sum((obs_ts[[fid]][1:simperiod] - v)^2))
  })))
  fid <- 2
  tmp <- clusterApplyLB(cl, 1:n_core, function(i, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas,
                                               gamma1s, gamma2s, gamma3s, taus, sim_ts) {
    m <- ms[i]
    n_pig <- n_pigs[fid]
    rho <- rhos[[fid]][i]
    beta <- betas[i]
    omega <- omegas[i]
    psi <- psis[i]
    alpha <- alphas[i]
    gamma1 <- round(gamma1s[i], 0)
    gamma2 <- round(gamma2s[i], 0)
    gamma3 <- round(gamma3s[i], 0)
    tau <- round(taus[[fid]][i], 0)
    seed <- ((j - 1) * n_core) + i
    
    dyn.load("FMD_sim.dll")
    v <- .C("_FMD_sim", 
            as.integer(i), as.integer(n_pig), as.integer(m), as.double(rho), as.double(beta), as.double(omega),
            as.double(psi), as.double(alpha), as.integer(gamma1), as.integer(gamma2), 
            as.integer(gamma3), as.integer(tau), as.integer(sim_ts))
    dyn.unload("FMD_sim.dll"); v
  }, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas, gamma1s, gamma2s, gamma3s, taus, sim_ts)
  diff_02 <- c(diff_02, unlist(lapply(lapply(tmp, unlist), function(x) {
    v <- tail(x, simperiod)
    sqrt(sum((obs_ts[[fid]][1:simperiod] - v)^2))
  })))
  fid <- 3
  tmp <- clusterApplyLB(cl, 1:n_core, function(i, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas,
                                               gamma1s, gamma2s, gamma3s, taus, sim_ts) {
    m <- ms[i]
    n_pig <- n_pigs[fid]
    rho <- rhos[[fid]][i]
    beta <- betas[i]
    omega <- omegas[i]
    psi <- psis[i]
    alpha <- alphas[i]
    gamma1 <- round(gamma1s[i], 0)
    gamma2 <- round(gamma2s[i], 0)
    gamma3 <- round(gamma3s[i], 0)
    tau <- round(taus[[fid]][i], 0)
    seed <- ((j - 1) * n_core) + i
    
    dyn.load("FMD_sim.dll")
    v <- .C("_FMD_sim", 
            as.integer(i), as.integer(n_pig), as.integer(m), as.double(rho), as.double(beta), as.double(omega),
            as.double(psi), as.double(alpha), as.integer(gamma1), as.integer(gamma2), 
            as.integer(gamma3), as.integer(tau), as.integer(sim_ts))
    dyn.unload("FMD_sim.dll"); v
  }, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas, gamma1s, gamma2s, gamma3s, taus, sim_ts)
  diff_03 <- c(diff_03, unlist(lapply(lapply(tmp, unlist), function(x) {
    v <- tail(x, simperiod)
    sqrt(sum((obs_ts[[fid]][1:simperiod] - v)^2))
  })))
  fid <- 4
  tmp <- clusterApplyLB(cl, 1:n_core, function(i, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas,
                                               gamma1s, gamma2s, gamma3s, taus, sim_ts) {
    m <- ms[i]
    n_pig <- n_pigs[fid]
    rho <- rhos[[fid]][i]
    beta <- betas[i]
    omega <- omegas[i]
    psi <- psis[i]
    alpha <- alphas[i]
    gamma1 <- round(gamma1s[i], 0)
    gamma2 <- round(gamma2s[i], 0)
    gamma3 <- round(gamma3s[i], 0)
    tau <- round(taus[[fid]][i], 0)
    seed <- ((j - 1) * n_core) + i
    
    dyn.load("FMD_sim.dll")
    v <- .C("_FMD_sim", 
            as.integer(i), as.integer(n_pig), as.integer(m), as.double(rho), as.double(beta), as.double(omega),
            as.double(psi), as.double(alpha), as.integer(gamma1), as.integer(gamma2), 
            as.integer(gamma3), as.integer(tau), as.integer(sim_ts))
    dyn.unload("FMD_sim.dll"); v
  }, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas, gamma1s, gamma2s, gamma3s, taus, sim_ts)
  diff_04 <- c(diff_04, unlist(lapply(lapply(tmp, unlist), function(x) {
    v <- tail(x, simperiod)
    sqrt(sum((obs_ts[[fid]][1:simperiod] - v)^2))
  })))
  
  m_l <- c(m_l, ms)
  
  for (i in 1:n_mod) {
    v <- ms[i]
    beta_l[[v]] <- c(beta_l[[v]], betas[i])
    omega_l[[v]] <- c(omega_l[[v]], omegas[i])
    psi_l[[v]] <- c(psi_l[[v]], psis[i])
    alpha_l[[v]] <- c(alpha_l[[v]], alphas[i])
    gamma1_l[[v]] <- c(gamma1_l[[v]], gamma1s[i])
    gamma2_l[[v]] <- c(gamma2_l[[v]], gamma2s[i])
    gamma3_l[[v]] <- c(gamma3_l[[v]], gamma3s[i])
    
    for (k in 1:n_farm) {
      rho_l[[v]][[k]] <- c(rho_l[[v]][[k]], rhos[[k]][i])
      tau_l[[v]][[k]] <- c(tau_l[[v]][[k]], taus[[k]][i])
    }
  }
}
save.image(paste("fmd_modselect1.RData", sep= ""))
t_mid <- Sys.time(); print(1); print(t_mid - t_init); print(t_mid); Sys.sleep(0.01); gc()
diff_l01 <- diff_01; diff_l02 <- diff_02; diff_l03 <- diff_03; diff_l04 <- diff_04

# 4. ABC-SMC for the subsequent sequences
# 4.1. Information required for subsequent sequences
total_run <- n_particle

rho_w <- rho_l
rho_w <- rapply(rho_w, function(x) {ifelse(is.na(x), x, 1)}, how= "replace")
tau_w <- tau_l
tau_w <- rapply(tau_w, function(x) {ifelse(is.na(x), x, 1)}, how= "replace")

beta_w <- lapply(beta_l, function(x) {ifelse(is.na(x), x, 1)})
omega_w <- lapply(omega_l, function(x) {ifelse(is.na(x), x, 1)})
psi_w <- lapply(psi_l, function(x) {ifelse(is.na(x), x, 1)})
alpha_w <- lapply(alpha_l, function(x) {ifelse(is.na(x), x, 1)})
gamma1_w <- lapply(gamma1_l, function(x) {ifelse(is.na(x), x, 1)})
gamma2_w <- lapply(gamma2_l, function(x) {ifelse(is.na(x), x, 1)})
gamma3_w <- lapply(gamma3_l, function(x) {ifelse(is.na(x), x, 1)})

# 4.2. Running the simulation
for (seq in 2:n_seq) {
  m_tl <- NULL
  
  rho_tl <- vector(mode= "list", length= n_mod)
  tau_tl <- vector(mode= "list", length= n_mod)
  for (i in 1:n_mod) {
    rho_tl[[i]] <- list(NULL,NULL,NULL,NULL)
    tau_tl[[i]] <- list(NULL,NULL,NULL,NULL)
  }
  beta_tl <- vector(mode= "list", length= n_mod)
  omega_tl <- vector(mode= "list", length= n_mod)
  psi_tl <- vector(mode= "list", length= n_mod)
  alpha_tl <- vector(mode= "list", length= n_mod)
  gamma1_tl <- vector(mode= "list", length= n_mod)
  gamma2_tl <- vector(mode= "list", length= n_mod)
  gamma3_tl <- vector(mode= "list", length= n_mod)
  
  rho_tw <- vector(mode= "list", length= n_mod)
  tau_tw <- vector(mode= "list", length= n_mod)
  for (i in 1:n_mod) {
    rho_tw[[i]] <- list(NULL,NULL,NULL,NULL)
    tau_tw[[i]] <- list(NULL,NULL,NULL,NULL)
  }
  beta_tw <- vector(mode= "list", length= n_mod)
  omega_tw <- vector(mode= "list", length= n_mod)
  psi_tw <- vector(mode= "list", length= n_mod)
  alpha_tw <- vector(mode= "list", length= n_mod)
  gamma1_tw <- vector(mode= "list", length= n_mod)
  gamma2_tw <- vector(mode= "list", length= n_mod)
  gamma3_tw <- vector(mode= "list", length= n_mod)
  
  diff_tl01 <- NULL; diff_tl02 <- NULL; diff_tl03 <- NULL; diff_tl04 <- NULL
  
  # Measure the threshold for summary statistics
  cutoff_01 <- quantile(diff_l01, q_lim); cutoff_02 <- quantile(diff_l02, q_lim)
  cutoff_03 <- quantile(diff_l03, q_lim); cutoff_04 <- quantile(diff_l04, q_lim)
  
  # Calculate the standard deviation for perturbation kernel
  rho_sd <- vector(mode= "list", length= n_mod)
  tau_sd <- vector(mode= "list", length= n_mod)
  beta_sd <- NULL
  omega_sd <- NULL
  psi_sd <- NULL
  alpha_sd <- NULL
  gamma1_sd <- NULL
  gamma2_sd <- NULL
  gamma3_sd <- NULL
  
  df <- data.frame(table(m_l))
  df <- df[df$Freq > 1,]
  mods <- as.numeric(df$m_l)
  for (i in mods) {
    for (j in 1:n_farm) {
      rho_sd[[i]] <- c(rho_sd[[i]],
                       (tuner * sum(((rho_l[[i]][[j]] - mean(rho_l[[i]][[j]]))^2) * rho_w[[i]][[j]]) / sum(rho_w[[i]][[j]]))^0.5)
      tau_sd[[i]] <- c(tau_sd[[i]],
                       (tuner * sum(((tau_l[[i]][[j]] - mean(tau_l[[i]][[j]]))^2) * tau_w[[i]][[j]]) / sum(tau_w[[i]][[j]]))^0.5)
    }
    beta_sd <- c(beta_sd, 
                 (tuner * sum(((beta_l[[i]] - mean(beta_l[[i]]))^2) * beta_w[[i]]) / sum(beta_w[[i]]))^0.5)
    omega_sd <- c(omega_sd, 
                  (tuner * sum(((omega_l[[i]] - mean(omega_l[[i]]))^2) * omega_w[[i]]) / sum(omega_w[[i]]))^0.5)
    psi_sd <- c(psi_sd, 
                (tuner * sum(((psi_l[[i]] - mean(psi_l[[i]]))^2) * psi_w[[i]]) / sum(psi_w[[i]]))^0.5)
    alpha_sd <- c(alpha_sd, 
                  (tuner * sum(((alpha_l[[i]] - mean(alpha_l[[i]]))^2) * alpha_w[[i]]) / sum(alpha_w[[i]]))^0.5)
    gamma1_sd <- c(gamma1_sd, 
                   (tuner * sum(((gamma1_l[[i]] - mean(gamma1_l[[i]]))^2) * gamma1_w[[i]]) / sum(gamma1_w[[i]]))^0.5)
    gamma2_sd <- c(gamma2_sd, 
                   (tuner * sum(((gamma2_l[[i]] - mean(gamma2_l[[i]]))^2) * gamma2_w[[i]]) / sum(gamma2_w[[i]]))^0.5)
    gamma3_sd <- c(gamma3_sd, 
                   (tuner * sum(((gamma3_l[[i]] - mean(gamma3_l[[i]]))^2) * gamma3_w[[i]]) / sum(gamma3_w[[i]]))^0.5)
  }
  
  j <- 1
  while (length(m_tl) < n_particle) {
    ms <- sample(mods, n_core, replace= T)
    
    # Perturbate the general parameter values
    rhos <- list(NULL, NULL, NULL, NULL)
    taus <- list(NULL, NULL, NULL, NULL)

    betas <- NULL; k <- 1
    while (length(betas) < n_core) {
      beta_star <- sample(beta_l[[ms[k]]], size= 1, prob= beta_w[[ms[k]]])
      beta <- rnorm(1, beta_star, beta_sd[ms[k]]) 
      if (beta >= beta_min & beta <= beta_max) {
        betas <- c(betas, beta)
      }
    }
    omegas <- NULL; k <- 1
    while (length(omegas) < n_core) {
      omega_star <- sample(omega_l[[ms[k]]], size= 1, prob= omega_w[[ms[k]]])
      omega <- rnorm(1, omega_star, omega_sd[ms[k]]) 
      if (omega >= omega_min & omega <= omega_max) {
        omegas <- c(omegas, omega)
      }
    }
    psis <- NULL; k <- 1
    while (length(psis) < n_core) {
      psi_star <- sample(psi_l[[ms[k]]], size= 1, prob= psi_w[[ms[k]]])
      psi <- rnorm(1, psi_star, psi_sd[ms[k]]) 
      if (psi >= psi_min & psi <= psi_max) {
        psis <- c(psis, psi)
      }
    }
    alphas <- NULL; k <- 1
    while (length(alphas) < n_core) {
      alpha_star <- sample(alpha_l[[ms[k]]], size= 1, prob= alpha_w[[ms[k]]])
      alpha <- rnorm(1, alpha_star, alpha_sd[ms[k]]) 
      if (alpha >= alpha_min & alpha <= alpha_max) {
        alphas <- c(alphas, alpha)
      }
    }
    gamma1s <- gamma2s <- NULL; k <- 1
    while (length(gamma1s) < n_core) {
      gamma1_star <- sample(gamma1_l[[ms[k]]], size= 1, prob= gamma1_w[[ms[k]]])
      gamma1 <- rnorm(1, gamma1_star, gamma1_sd[ms[k]]) 
      gamma2_star <- sample(gamma2_l[[ms[k]]], size= 1, prob= gamma2_w[[ms[k]]])
      gamma2 <- rnorm(1, gamma2_star, gamma2_sd[ms[k]]) 
      if (gamma1 + gamma2 <= gamma_max & gamma1 <= gamma1_max & gamma1 >= gamma_min & gamma2 >= gamma_min) {
        gamma1s <- c(gamma1s, gamma1)
        gamma2s <- c(gamma2s, gamma2)
      }
    }
    gamma3s <- NULL; k <- 1
    while (length(gamma3s) < n_core) {
      gamma3_star <- sample(gamma3_l[[ms[k]]], size= 1, prob= gamma3_w[[ms[k]]])
      gamma3 <- rnorm(1, gamma3_star, gamma3_sd[ms[k]]) 
      if (gamma3 >= gamma3_min & gamma3 <= gamma3_max) {
        gamma3s <- c(gamma3s, gamma3)
      }
    }
    
    fid <- 1; k <- 1
    while (length(rhos[[fid]]) < n_core) {
      rho_star <- sample(rho_l[[ms[k]]][[fid]], size= 1, prob= rho_w[[ms[k]]][[fid]])
      rho <- rnorm(1, rho_star, rho_sd[[ms[k]]][fid]) 
      if (rho >= rho_min & rho <= rho_max) {
        rhos[[fid]] <- c(rhos[[fid]], rho)
      }
    }
    k <- 1
    while (length(taus[[fid]]) < n_core) {
      tau_star <- sample(tau_l[[ms[k]]][[fid]], size= 1, prob= tau_w[[ms[k]]][[fid]])
      tau <- rnorm(1, tau_star, tau_sd[[ms[k]]][fid]) 
      if (tau >= tau_min & tau <= tau_max) {
        taus[[fid]] <- c(taus[[fid]], tau)
      }
    }
    tmp <- clusterApplyLB(cl, 1:n_core, function(i, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas,
                                                 gamma1s, gamma2s, gamma3s, taus, sim_ts) {
      m <- ms[i]
      n_pig <- n_pigs[fid]
      rho <- rhos[[fid]][i]
      beta <- betas[i]
      omega <- omegas[i]
      psi <- psis[i]
      alpha <- alphas[i]
      gamma1 <- round(gamma1s[i], 0)
      gamma2 <- round(gamma2s[i], 0)
      gamma3 <- round(gamma3s[i], 0)
      tau <- round(taus[[fid]][i], 0)
      seed <- ((j - 1) * n_core) + i
      
      dyn.load("FMD_sim.dll")
      v <- .C("_FMD_sim", 
              as.integer(i), as.integer(n_pig), as.integer(m), as.double(rho), as.double(beta), as.double(omega),
              as.double(psi), as.double(alpha), as.integer(gamma1), as.integer(gamma2), 
              as.integer(gamma3), as.integer(tau), as.integer(sim_ts))
      dyn.unload("FMD_sim.dll"); v
    }, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas, gamma1s, gamma2s, gamma3s, taus, sim_ts)
    diff_01 <- NULL
    diff_01 <- c(diff_01, unlist(lapply(lapply(tmp, unlist), function(x) {
      v <- tail(x, simperiod)
      sqrt(sum((obs_ts[[fid]][1:simperiod] - v)^2))
    })))
    
    fid <- 2; k <- 1
    while (length(rhos[[fid]]) < n_core) {
      rho_star <- sample(rho_l[[ms[k]]][[fid]], size= 1, prob= rho_w[[ms[k]]][[fid]])
      rho <- rnorm(1, rho_star, rho_sd[[ms[k]]][fid]) 
      if (rho >= rho_min & rho <= rho_max) {
        rhos[[fid]] <- c(rhos[[fid]], rho)
      }
    }
    k <- 1
    while (length(taus[[fid]]) < n_core) {
      tau_star <- sample(tau_l[[ms[k]]][[fid]], size= 1, prob= tau_w[[ms[k]]][[fid]])
      tau <- rnorm(1, tau_star, tau_sd[[ms[k]]][fid]) 
      if (tau >= tau_min & tau <= tau_max) {
        taus[[fid]] <- c(taus[[fid]], tau)
      }
    }
    tmp <- clusterApplyLB(cl, 1:n_core, function(i, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas,
                                                 gamma1s, gamma2s, gamma3s, taus, sim_ts) {
      m <- ms[i]
      n_pig <- n_pigs[fid]
      rho <- rhos[[fid]][i]
      beta <- betas[i]
      omega <- omegas[i]
      psi <- psis[i]
      alpha <- alphas[i]
      gamma1 <- round(gamma1s[i], 0)
      gamma2 <- round(gamma2s[i], 0)
      gamma3 <- round(gamma3s[i], 0)
      tau <- round(taus[[fid]][i], 0)
      seed <- ((j - 1) * n_core) + i
      
      dyn.load("FMD_sim.dll")
      v <- .C("_FMD_sim", 
              as.integer(i), as.integer(n_pig), as.integer(m), as.double(rho), as.double(beta), as.double(omega),
              as.double(psi), as.double(alpha), as.integer(gamma1), as.integer(gamma2), 
              as.integer(gamma3), as.integer(tau), as.integer(sim_ts))
      dyn.unload("FMD_sim.dll"); v
    }, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas, gamma1s, gamma2s, gamma3s, taus, sim_ts)
    diff_02 <- NULL
    diff_02 <- c(diff_02, unlist(lapply(lapply(tmp, unlist), function(x) {
      v <- tail(x, simperiod)
      sqrt(sum((obs_ts[[fid]][1:simperiod] - v)^2))
    })))
    
    fid <- 3; k <- 1
    while (length(rhos[[fid]]) < n_core) {
      rho_star <- sample(rho_l[[ms[k]]][[fid]], size= 1, prob= rho_w[[ms[k]]][[fid]])
      rho <- rnorm(1, rho_star, rho_sd[[ms[k]]][fid]) 
      if (rho >= rho_min & rho <= rho_max) {
        rhos[[fid]] <- c(rhos[[fid]], rho)
      }
    }
    k <- 1
    while (length(taus[[fid]]) < n_core) {
      tau_star <- sample(tau_l[[ms[k]]][[fid]], size= 1, prob= tau_w[[ms[k]]][[fid]])
      tau <- rnorm(1, tau_star, tau_sd[[ms[k]]][fid]) 
      if (tau >= tau_min & tau <= tau_max) {
        taus[[fid]] <- c(taus[[fid]], tau)
      }
    }
    tmp <- clusterApplyLB(cl, 1:n_core, function(i, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas,
                                                 gamma1s, gamma2s, gamma3s, taus, sim_ts) {
      m <- ms[i]
      n_pig <- n_pigs[fid]
      rho <- rhos[[fid]][i]
      beta <- betas[i]
      omega <- omegas[i]
      psi <- psis[i]
      alpha <- alphas[i]
      gamma1 <- round(gamma1s[i], 0)
      gamma2 <- round(gamma2s[i], 0)
      gamma3 <- round(gamma3s[i], 0)
      tau <- round(taus[[fid]][i], 0)
      seed <- ((j - 1) * n_core) + i
      
      dyn.load("FMD_sim.dll")
      v <- .C("_FMD_sim", 
              as.integer(i), as.integer(n_pig), as.integer(m), as.double(rho), as.double(beta), as.double(omega),
              as.double(psi), as.double(alpha), as.integer(gamma1), as.integer(gamma2), 
              as.integer(gamma3), as.integer(tau), as.integer(sim_ts))
      dyn.unload("FMD_sim.dll"); v
    }, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas, gamma1s, gamma2s, gamma3s, taus, sim_ts)
    diff_03 <- NULL
    diff_03 <- c(diff_03, unlist(lapply(lapply(tmp, unlist), function(x) {
      v <- tail(x, simperiod)
      sqrt(sum((obs_ts[[fid]][1:simperiod] - v)^2))
    })))
    
    fid <- 4; k <- 1
    while (length(rhos[[fid]]) < n_core) {
      rho_star <- sample(rho_l[[ms[k]]][[fid]], size= 1, prob= rho_w[[ms[k]]][[fid]])
      rho <- rnorm(1, rho_star, rho_sd[[ms[k]]][fid]) 
      if (rho >= rho_min & rho <= rho_max) {
        rhos[[fid]] <- c(rhos[[fid]], rho)
      }
    }
    k <- 1
    while (length(taus[[fid]]) < n_core) {
      tau_star <- sample(tau_l[[ms[k]]][[fid]], size= 1, prob= tau_w[[ms[k]]][[fid]])
      tau <- rnorm(1, tau_star, tau_sd[[ms[k]]][fid]) 
      if (tau >= tau_min & tau <= tau_max) {
        taus[[fid]] <- c(taus[[fid]], tau)
      }
    }
    tmp <- clusterApplyLB(cl, 1:n_core, function(i, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas,
                                                 gamma1s, gamma2s, gamma3s, taus, sim_ts) {
      m <- ms[i]
      n_pig <- n_pigs[fid]
      rho <- rhos[[fid]][i]
      beta <- betas[i]
      omega <- omegas[i]
      psi <- psis[i]
      alpha <- alphas[i]
      gamma1 <- round(gamma1s[i], 0)
      gamma2 <- round(gamma2s[i], 0)
      gamma3 <- round(gamma3s[i], 0)
      tau <- round(taus[[fid]][i], 0)
      seed <- ((j - 1) * n_core) + i
      
      dyn.load("FMD_sim.dll")
      v <- .C("_FMD_sim", 
              as.integer(i), as.integer(n_pig), as.integer(m), as.double(rho), as.double(beta), as.double(omega),
              as.double(psi), as.double(alpha), as.integer(gamma1), as.integer(gamma2), 
              as.integer(gamma3), as.integer(tau), as.integer(sim_ts))
      dyn.unload("FMD_sim.dll"); v
    }, j, ms, fid, n_pigs, rhos, betas, omegas, psis, alphas, gamma1s, gamma2s, gamma3s, taus, sim_ts)
    diff_04 <- NULL
    diff_04 <- c(diff_04, unlist(lapply(lapply(tmp, unlist), function(x) {
      v <- tail(x, simperiod)
      sqrt(sum((obs_ts[[fid]][1:simperiod] - v)^2))
    })))
    
    idx <- diff_01 <= cutoff_01 & diff_02 <= cutoff_02 & diff_03 <= cutoff_03 & diff_04 <= cutoff_04
    if (length(m_tl) + sum(idx) > n_particle) {
      v <- length(m_tl) + sum(idx) - n_particle
      ids <- which(idx == T)[1:v]
    } else {
      ids <- which(idx)
    }
    
    if (sum(idx) > 0) {
      m_tl <- c(m_tl, ms[ids])
      tmp <- ms[ids]
      
      for (i in c(1:length(tmp))) {
        v <- tmp[i]
        beta_tl[[v]] <- c(beta_tl[[v]], betas[ids[i]])
        beta_tw[[v]] <- c(beta_tw[[v]], dnorm(betas[ids[i]], beta_pmu, beta_psd) / 
                            sum(beta_w[[v]] * dnorm(betas[ids[i]], beta_l[[v]], beta_sd[v])))
        omega_tl[[v]] <- c(omega_tl[[v]], omegas[ids[i]])
        omega_tw[[v]] <- c(omega_tw[[v]], 1 / 
                             sum(omega_w[[v]] * dnorm(omegas[ids[i]], omega_l[[v]], omega_sd[v])))
        psi_tl[[v]] <- c(psi_tl[[v]], psis[ids[i]])
        psi_tw[[v]] <- c(psi_tw[[v]], 1 / 
                           sum(psi_w[[v]] * dnorm(psis[ids[i]], psi_l[[v]], psi_sd[v])))
        alpha_tl[[v]] <- c(alpha_tl[[v]], alphas[ids[i]])
        alpha_tw[[v]] <- c(alpha_tw[[v]], 1 / 
                             sum(alpha_w[[v]] * dnorm(alphas[ids[i]], alpha_l[[v]], alpha_sd[v])))
        gamma1_tl[[v]] <- c(gamma1_tl[[v]], gamma1s[ids[i]])
        gamma1_tw[[v]] <- c(gamma1_tw[[v]], dbinom(round(gamma1s[ids[i]], 0), gamma1_n, gamma1_p) / 
                         sum(gamma1_w[[v]] * dnorm(gamma1s[ids[i]], gamma1_l[[v]], gamma1_sd[v])))
        gamma2_tl[[v]] <- c(gamma2_tl[[v]], gamma2s[ids[i]])
        gamma2_tw[[v]] <- c(gamma2_tw[[v]], dbinom(round(gamma2s[ids[i]], 0), gamma2_n, gamma2_p) / 
                         sum(gamma2_w[[v]] * dnorm(gamma2s[ids[i]], gamma2_l[[v]], gamma2_sd[v])))
        gamma3_tl[[v]] <- c(gamma3_tl[[v]], gamma3s[ids[i]])
        gamma3_tw[[v]] <- c(gamma3_tw[[v]], 1 / 
                              sum(gamma3_w[[v]] * dnorm(gamma3s[ids[i]], gamma3_l[[v]], gamma3_sd[v])))
        
        for (k in 1:n_farm) {
          rho_tl[[v]][[k]] <- c(rho_tl[[v]][[k]], rhos[[k]][ids[i]])
          rho_tw[[v]][[k]] <- c(rho_tw[[v]][[k]], dbeta(rhos[[k]][ids[i]], rho_a[k], rho_b[k]) / 
                                  sum(rho_w[[v]][[k]] * dnorm(rhos[[k]][ids[i]], rho_l[[v]][[k]], rho_sd[[v]][k])))
          tau_tl[[v]][[k]] <- c(tau_tl[[v]][[k]], taus[[k]][ids[i]])
          tau_tw[[v]][[k]] <- c(tau_tw[[v]][[k]], 1 / 
                                  sum(tau_w[[v]][[k]] * dnorm(taus[[k]][ids[i]], tau_l[[v]][[k]], tau_sd[[v]][k])))
        }
      }
      
      diff_tl01 <- c(diff_tl01, diff_01[idx]); diff_tl02 <- c(diff_tl02, diff_02[idx])
      diff_tl03 <- c(diff_tl03, diff_03[idx]); diff_tl04 <- c(diff_tl04, diff_04[idx])
      
      if (length(m_tl) %% 100 == 0) {
        print(length(m_tl)); Sys.sleep(0.001)
      }
    }
    
    total_run <- total_run + n_core
    j <- j + 1
  }
  m_l <- m_tl
  beta_l <- beta_tl
  omega_l <- omega_tl
  psi_l <- psi_tl
  alpha_l <- alpha_tl
  gamma1_l <- gamma1_tl
  gamma2_l <- gamma2_tl
  gamma3_l <- gamma3_tl
  beta_w <- beta_tw
  omega_w <- omega_tw
  psi_w <- psi_tw
  alpha_w <- alpha_tw
  gamma1_w <- gamma1_tw
  gamma2_w <- gamma2_tw
  gamma3_w <- gamma3_tw
  rho_l <- rho_tl
  tau_l <- tau_tl
  rho_w <- rho_tw
  tau_w <- tau_tw

  diff_l01 <- diff_tl01; diff_l02 <- diff_tl02
  diff_l03 <- diff_tl03; diff_l04 <- diff_tl04
  
  save.image(paste("fmd_modselect", seq, ".RData", sep= ""))
  t_mid <- Sys.time(); print(seq); print(t_mid - t_init); print(t_mid); Sys.sleep(0.001); gc()
}
stopCluster(cl)

# Analysis
rm(list= ls())
load("D:/Choennam Univ/FMD_Korea/Code/fmd_modselect1.RData")
table(m_l)
dat <- data.frame("seq"= 1, table(m_l))
load("D:/Choennam Univ/FMD_Korea/Code/fmd_modselect2.RData"); c(cutoff_01,cutoff_02,cutoff_03,cutoff_04)
table(m_l)
dat <- rbind(dat, data.frame("seq"= 2, table(m_l)))
load("D:/Choennam Univ/FMD_Korea/Code/fmd_modselect3.RData"); c(cutoff_01,cutoff_02,cutoff_03,cutoff_04)
table(m_l)
dat <- rbind(dat, data.frame("seq"= 3, table(m_l)))
load("D:/Choennam Univ/FMD_Korea/Code/fmd_modselect4.RData"); c(cutoff_01,cutoff_02,cutoff_03,cutoff_04)
table(m_l)
dat <- rbind(dat, data.frame("seq"= 4, table(m_l)))
load("D:/Choennam Univ/FMD_Korea/Code/fmd_modselect5.RData"); c(cutoff_01,cutoff_02,cutoff_03,cutoff_04)
table(m_l)
dat <- rbind(dat, data.frame("seq"= 5, table(m_l)))
load("D:/Choennam Univ/FMD_Korea/Code/fmd_modselect6.RData"); c(cutoff_01,cutoff_02,cutoff_03,cutoff_04)
table(m_l)
dat <- rbind(dat, data.frame("seq"= 6, table(m_l)))
load("D:/Choennam Univ/FMD_Korea/Code/fmd_modselect7.RData"); c(cutoff_01,cutoff_02,cutoff_03,cutoff_04)
table(m_l)
dat <- rbind(dat, data.frame("seq"= 7, table(m_l)))
load("D:/Choennam Univ/FMD_Korea/Code/fmd_modselect8.RData"); c(cutoff_01,cutoff_02,cutoff_03,cutoff_04)
table(m_l)
dat <- rbind(dat, data.frame("seq"= 8, table(m_l)))
load("D:/Choennam Univ/FMD_Korea/Code/fmd_modselect9.RData"); c(cutoff_01,cutoff_02,cutoff_03,cutoff_04)
table(m_l)
dat <- rbind(dat, data.frame("seq"= 9, table(m_l)))
load("D:/Choennam Univ/FMD_Korea/Code/fmd_modselect10.RData"); c(cutoff_01,cutoff_02,cutoff_03,cutoff_04)
table(m_l)
dat <- rbind(dat, data.frame("seq"= 10, table(m_l)))
load("D:/Choennam Univ/FMD_Korea/Code/fmd_modselect11.RData"); c(cutoff_01,cutoff_02,cutoff_03,cutoff_04)
table(m_l)
dat <- rbind(dat, data.frame("seq"= 11, table(m_l)))
load("D:/Choennam Univ/FMD_Korea/Code/fmd_modselect12.RData"); c(cutoff_01,cutoff_02,cutoff_03,cutoff_04)
table(m_l)
dat <- rbind(dat, data.frame("seq"= 12, table(m_l)))

library(ggplot2)
ggplot(data= dat, aes(x= m_l, y= Freq)) +
  geom_bar(stat= "identity") +
  facet_wrap(~seq, nrow= 2)
