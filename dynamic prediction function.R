
#### dynamic prediction function ####

unitlevel_DP_JM <- function(object, 
                            object_summary,
                            surv_data, 
                            long_data, 
                            type = c("SurvProb", "Density"), 
                            idVar = "id", 
                            pred_ij, 
                            simulate = TRUE, 
                            survTimes, 
                            pred.time, 
                            LeftTrunc_var = NULL, 
                            MCMCnum = 200L, 
                            CI.levels = c(0.025, 0.975), 
                            log = FALSE, 
                            scale = 1.6, 
                            init.b = NULL, 
                            seed = 1L) {
  
  
  TT <- long_data$followup_year 

  #### prepare data ####
  timeVar <- "followup_year"
  hasScale <- T # biomarker is continuous, has measurement error
  gammas <- object_summary[object_summary$variable %in% 
                                     c("gamma[1]","gamma[2]","gamma[3]","gamma[4]","gamma[5]","gamma[6]"),2][[1]]
  betas1 <- object_summary[object_summary$variable %in% 
                                     c("beta1[1]","beta1[2]","beta1[3]","beta1[4]","beta1[5]","beta1[6]","beta1[7]"),2][[1]]
  betas2 <- object_summary[object_summary$variable %in% 
                                     c("beta2[1]","beta2[2]","beta2[3]","beta2[4]","beta2[5]","beta2[6]","beta2[7]"),2][[1]]
  xi <- 1
  alpha1 <- object_summary[object_summary$variable %in% c("alpha1"),2][[1]]
  alpha2 <- object_summary[object_summary$variable %in% c("alpha2"),2][[1]]
  sigma_e <- object_summary[object_summary$variable %in% c("sigma_e"),2][[1]]
  sigma_bi <- object_summary[object_summary$variable %in% c("sigma_bi[1]","sigma_bi[2]"),2][[1]]
  sigma_bij <- object_summary[object_summary$variable %in% c("sigma_bij[1]", "sigma_bij[2]","sigma_bij[3]","sigma_bij[4]"),2][[1]]
  sigma_ri <- object_summary[object_summary$variable=="sigma_ri",2][[1]]
  R_bi <- object_summary[object_summary$variable=="R_bi[1,2]",2][[1]]
  R_bij <- object_summary[object_summary$variable %in% 
                                    c("R_bij[1,2]", "R_bij[1,3]", "R_bij[1,4]","R_bij[2,3]","R_bij[2,4]","R_bij[3,4]"),2][[1]]
  
  var_ri = sigma_ri^2
  var_e = sigma_e^2
  
  length_bi <- length(sigma_bi)
  Cor_bi <- diag(1, length_bi)
  Cor_bi[lower.tri(Cor_bi)] <- R_bi
  Cor_bi <- Cor_bi + t(Cor_bi) - diag(diag(Cor_bi))  # Make symmetric
  var_bi <- diag(sigma_bi) %*% Cor_bi %*%  diag(sigma_bi)
  
  length_bij <- length(sigma_bij)
  Cor_bij <- diag(1, length_bij)
  Cor_bij[lower.tri(Cor_bij)] <- R_bij
  Cor_bij <- Cor_bij + t(Cor_bij) - diag(diag(Cor_bij))  # Make symmetric
  var_bij <- diag(sigma_bij) %*% Cor_bij %*%  diag(sigma_bij)
  
  max.time <- max(TT)
  
  
  
  
  ## filter ij's who are at risk at pred.time
  pred_data_surv <- surv_data %>% filter(eventtime > pred.time)
  pred_data_long <- long_data %>% 
    filter((id_tooth %in% pred_data_surv$id_tooth) & (followup_year <= pred.time)) %>% group_by(id_tooth) 
  
  ## filter ij's who are at risk and being censored at pred.time
  pred_data_surv_censored <- surv_data %>% filter(status == 0)
  pred_data_long_censored <- long_data %>% 
    filter((id_tooth %in% pred_data_surv_censored$id_tooth)) %>% group_by(id_tooth) 
  
  
  
  ## the specific ij of interest to make dynamic prediction
  pred_data_long_ij <- pred_data_long %>% filter(id_tooth == pred_ij)
  
  id <- pred_data_long_ij$id_tooth
  X_long <- c(pred_data_long_ij$sex, pred_data_long_ij$std_age, pred_data_long_ij$smoke, 
              pred_data_long_ij$molar, pred_data_long_ij$lower)
  X_surv <- c(pred_data_long_ij$sex[1], pred_data_long_ij$std_age[1], pred_data_long_ij$smoke[1], 
              pred_data_long_ij$molar[1], pred_data_long_ij$lower[1])
  #  Z_long <- data.frame(3)
  yA <- log(pred_data_long_ij$max_ppd)
  yB <- pred_data_long_ij$mobility
  previous.time <- pred_data_long_ij$followup_year  # observed times before prediction 
  

  n <- length(TT) # number of total observation time
  n.tp <- length(pred.time)
  ncz <- 7
  
  #### useful functions ####
  
  rmvt <- function (n, mu, Sigma, df) {
    p <- length(mu)
    if (is.list(Sigma)) {
      ev <- Sigma$values
      evec <- Sigma$vectors
    } else {
      ed <- eigen(Sigma, symmetric = TRUE)
      ev <- ed$values
      evec <- ed$vectors
    }
    X <- drop(mu) + tcrossprod(evec * rep(sqrt(pmax(ev, 0)), each = p), 
                               matrix(rnorm(n * p), n)) / rep(sqrt(rchisq(n, df)/df), each = p)
    if (n == 1L) drop(X) else t.default(X)
  }
  
  dmvt <- function (x, mu, Sigma = NULL, invSigma = NULL, df, log = FALSE, prop = TRUE) {
    if (!is.numeric(x)) 
      stop("'x' must be a numeric matrix or vector")
    if (!is.matrix(x)) 
      x <- rbind(x)
    p <- length(mu)
    if (is.null(Sigma) && is.null(invSigma))
      stop("'Sigma' or 'invSigma' must be given.")
    if (!is.null(Sigma)) {
      if (is.list(Sigma)) {
        ev <- Sigma$values
        evec <- Sigma$vectors
      } else {
        ed <- eigen(Sigma, symmetric = TRUE)
        ev <- ed$values
        evec <- ed$vectors
      }
      invSigma <- evec %*% (t(evec)/ev)
      if (!prop)
        logdetSigma <- sum(log(ev))
    } else {
      if (!prop)
        logdetSigma <- c(-determinant(invSigma)$modulus)
    }
    ss <- x - rep(mu, each = nrow(x))
    quad <- rowSums((ss %*% invSigma) * ss)/df
    if (!prop)
      fact <- lgamma((df + p)/2) - lgamma(df/2) - 0.5 * (p * (log(pi) + 
                                                                log(df)) + logdetSigma)
    if (log) {
      if (!prop) as.vector(fact - 0.5 * (df + p) * log(1 + quad)) else as.vector(- 0.5 * (df + p) * log(1 + quad))
    } else {
      if (!prop) as.vector(exp(fact) * ((1 + quad)^(-(df + p)/2))) else as.vector(((1 + quad)^(-(df + p)/2)))
    }
  }
  
  densLong_norm <- function (y, eta.y, scale, log = FALSE, data) {
    dnorm(x = y, mean = eta.y, sd = scale, log = log)
  }
  
  dmvnorm <- function (x, mu, Sigma = NULL, invSigma = NULL, log = FALSE, prop = TRUE) {
    if (!is.matrix(x))
      x <- rbind(x)
    if (is.matrix(mu)) {
      if (nrow(mu) != nrow(x))
        stop("incorrect dimensions for 'mu'.")
      p <- ncol(mu)
    } else {
      p <- length(mu)
      mu <- rep(mu, each = nrow(x))
    }
    if (is.null(Sigma) && is.null(invSigma))
      stop("'Sigma' or 'invSigma' must be given.")
    if (!is.null(Sigma)) {
      if (is.list(Sigma)) {
        ev <- Sigma$values
        evec <- Sigma$vectors
      } else {
        ed <- eigen(Sigma, symmetric = TRUE)
        ev <- ed$values
        evec <- ed$vectors            
      }
      invSigma <- evec %*% (t(evec) / ev)
      if (!prop)
        logdetSigma <- sum(log(ev))
    } else {
      if (!prop)
        logdetSigma <- - determinant(as.matrix(invSigma))$modulus
    }
    ss <- x - mu
    quad <- 0.5 * rowSums((ss %*% invSigma) * ss)
    if (!prop)
      fact <- - 0.5 * (p * log(2 * pi) + logdetSigma)
    if (log) {
      if (!prop) as.vector(fact - quad) else as.vector(- quad)
    } else {
      if (!prop) as.vector(exp(fact - quad)) else as.vector(exp(- quad))
    }
  }
  
  densLong_binomial <- function (y, eta.y, scale, log = FALSE, data) {
    y*eta.y - log(1+exp(eta.y))
  }
  
  
  log.posterior.b <- function (b, yA, yB, X_long, X_surv, previous.time, 
                               betas1, betas2, gammas, xi, alpha1, alpha2, 
                               var_e, var_bi, var_bij, var_ri) {
    # b=c(bi1, bi2, bij10, bij11, bij20, bij21, ri)
    mu.yA <- matrix(c(rep(1, length(previous.time)), 
                      X_long, previous.time), ncol = 7, byrow = F) %*% betas1 + 
      matrix(c(rep(1, length(previous.time)), rep(1, length(previous.time)), 
               previous.time), ncol = 3, byrow = F) %*% c(b[1], b[3:4])
    
    mu.yB <- matrix(c(rep(1, length(previous.time)), 
                      X_long, previous.time), ncol = 7, byrow = F) %*% betas2 + 
      matrix(c(rep(1, length(previous.time)), rep(1, length(previous.time)), 
               previous.time), ncol = 3, byrow = F) %*% c(b[2], b[5:6])
    
    
    logYA <- dnorm(yA, mu.yA, sd=sqrt(var_e), log =T)
    log.p.yA <- sum(logYA)
    
    logYB <- yB*mu.yB - log(1+exp(mu.yB))
    log.p.yB <- sum(logYB)
    
    
    sk <- c(-0.9914554, -0.9491079, -0.8648644, -0.7415312, -0.5860872,
            -0.4058452, -0.2077850,  0.0000000,  0.2077850,  0.4058452,
            0.5860872,  0.7415312,  0.864864,  0.9491079,  0.9914554)
    wk <- c(0.02293532, 0.06309209, 0.10479001, 0.14065326, 0.16900473,
            0.19035058, 0.20443294, 0.20948214, 0.20443294, 0.19035058,
            0.16900473, 0.14065326, 0.10479001, 0.06309209, 0.02293532)
    length_GK <- length(sk)
    
    
    log.survival <- 0
    for (j in 1:15){
      sk_scaled_j = (sk[j]+1)*pred.time/2
      wk_scaled_j = wk[j]*pred.time/2
      log.survival <- log.survival - wk_scaled_j*xi*sk_scaled_j^(xi-1)*exp((
        alpha1*betas1[4]+alpha1*b[4]+alpha2*betas2[4]+alpha2*b[6])*sk_scaled_j)*
        exp(b[7] + c(1,X_surv) %*% gammas + alpha1*(c(1, X_surv) %*% betas1[1:6]+b[1]+b[3])
            + alpha2*(c(1, X_surv) %*% betas2[1:6]+b[2]+b[5]))
    }
    
    
    log.p.bi <- dmvnorm(b[1:2], mu=rep(0, 2), Sigma=var_bi, log=T, prop=F)
    log.p.bij <- dmvnorm(b[3:6], mu=rep(0, 4), Sigma=var_bij, log=T, prop=F)
    log.p.ri <- dnorm(b[7], mean=0, sd=sigma_ri, log=T)
    
    log.p.yA + log.p.yB + log.survival + log.p.bi + log.p.bij + log.p.ri
  }

  
  S_t <- function(t, X_surv, betas1, betas2, gammas, 
                  xi, alpha1, alpha2,
                  b_Ai,
                  b_Bi,
                  b_Aij0,
                  b_Aij1,
                  b_Bij0,
                  b_Bij1,
                  b_ri){
  
    sk <- c(-0.9914554, -0.9491079, -0.8648644, -0.7415312, -0.5860872,
            -0.4058452, -0.2077850,  0.0000000,  0.2077850,  0.4058452,
            0.5860872,  0.7415312,  0.864864,  0.9491079,  0.9914554)
    wk <- c(0.02293532, 0.06309209, 0.10479001, 0.14065326, 0.16900473,
            0.19035058, 0.20443294, 0.20948214, 0.20443294, 0.19035058,
            0.16900473, 0.14065326, 0.10479001, 0.06309209, 0.02293532)
    length_GK <- length(sk)
    
    log_S_t <- 0
    for (j in 1:15) {
      sk_scaled_j = (sk[j]+1)*t/2
      wk_scaled_j = wk[j]*t/2
      
      log_S_t <- log_S_t - wk_scaled_j*xi*sk_scaled_j^(xi-1)*exp((
        alpha1*betas1[7]+alpha1*b_Aij1+alpha2*betas2[7]+alpha2* b_Bij1)*sk_scaled_j)*
        exp(b_ri + c(1,X_surv) %*% gammas + alpha1*(c(1, X_surv) %*% betas1[1:6]+b_Ai+ b_Aij0)
            + alpha2*(c(1, X_surv) %*% betas2[1:6]+b_Bi+b_Bij0))
    }
    exp(log_S_t)
    
  }
  
  
  
  
  #### calculate the Empirical Bayes estimates and the variance ####
  gammas.new <- gammas
  betas1.new <- betas1
  betas2.new <- betas1
  xi.new <- xi
  alpha1.new <- alpha1
  alpha2.new <- alpha2
  var_e.new <- var_e
  var_bi.new <- var_bi
  var_bij.new <- var_bij
  var_ri.new <- var_ri
  
  ff <- function (b) -log.posterior.b(b, yA=yA, yB=yB, X_long=X_long, 
                                      X_surv=X_surv, previous.time=previous.time, 
                                      betas1.new, betas2.new, gammas.new, 
                                      xi.new, 
                                      alpha1.new,
                                      alpha2.new,
                                      var_e.new, 
                                      var_bi.new,
                                      var_bij.new,
                                      var_ri.new)
  #  start <- if (is.null(init.b)) rep(0, 7) else init.b[i, ]
  start <- rep(0,7)
  
  modes.b <- matrix(0, n.tp, ncz)
  invVars.b <- Vars.b <- vector("list", n.tp)
  
  opt <- try(optim(start, ff, method = "BFGS", hessian = TRUE), silent = TRUE)
  modes.b[1,] <- opt$par
  invVars.b[[1]] <- opt$hessian/scale
  Vars.b[[1]] <- scale * solve(opt$hessian)
  
  
  
  #### dynamic prediction (CI) #### 
  set.seed(123)
  out <- matrix(data=NA, nrow = MCMCnum, ncol = length(survTimes))
  success.rate <- matrix(FALSE, MCMCnum, n.tp)
  b.old <- b.new <- modes.b
  b.matrix <-  matrix(data=NA, nrow = MCMCnum, ncol = 7)
  
  if (n.tp == 1)
    dim(b.old) <- dim(b.new) <- c(1L, ncz)
  
  mcmc <- object
  mcmc_colname <- names(mcmc)
  mcmc_colname_no_b <- mcmc_colname[!grepl("^bi", mcmc_colname)]
  mcmc_no_b <- data.frame(mcmc[mcmc_colname_no_b])
  samples <- sample(length(mcmc$`beta1[1]`), MCMCnum)
  mcmc <- mcmc_no_b[samples,]
  
  proposed.b <- mapply(rmvt, n = MCMCnum, mu = split(modes.b, row(modes.b)), Sigma = Vars.b, df = 4,
                       SIMPLIFY = FALSE)
  proposed.b[] <- lapply(proposed.b, function (x) if (is.matrix(x)) x else rbind(x))
  dmvt.proposed <- mapply(dmvt, x = proposed.b, mu = modes.b,
                          Sigma = Vars.b, MoreArgs = list(df = 4, log = TRUE), 
                          SIMPLIFY = FALSE)
  
  for (m in 1:MCMCnum) { 
    
    
    # Step 1: extract parameter values
    betas1.new <- c(mcmc$beta1.1.[m],mcmc$beta1.2.[m],mcmc$beta1.3.[m],mcmc$beta1.4.[m],mcmc$beta1.5.[m],mcmc$beta1.6.[m],mcmc$beta1.7.[m])
    betas2.new <- c(mcmc$beta2.1.[m],mcmc$beta2.2.[m],mcmc$beta2.3.[m],mcmc$beta2.4.[m],mcmc$beta2.5.[m],mcmc$beta2.6.[m],mcmc$beta2.7.[m])
    gammas.new <- c(mcmc$gamma.1.[m], mcmc$gamma.2.[m], mcmc$gamma.3.[m],mcmc$gamma.4.[m], mcmc$gamma.5.[m], mcmc$gamma.6.[m])
    xi.new <- 1
    alpha1.new <- mcmc$alpha1[m]
    alpha2.new <- mcmc$alpha2[m]
    
    sigma_e <- mcmc$sigma_e[m]
    sigma_bi <- c(mcmc$sigma_bi.1.[m], mcmc$sigma_bi.2.[m])
    sigma_bij <- c(mcmc$sigma_bij.1.[m], mcmc$sigma_bij.2.[m], mcmc$sigma_bij.3.[m], mcmc$sigma_bij.4.[m])
    sigma_ri <- mcmc$sigma_ri[m]  
    R_bi <- mcmc$R_bi.1.2.[m]
    R_bij <- c(mcmc$R_bij.1.2.[m], mcmc$R_bij.1.3.[m], mcmc$R_bij.1.4.[m], 
               mcmc$R_bij.2.3.[m], mcmc$R_bij.2.4.[m], mcmc$R_bij.3.4.[m])
    
    var_ri.new = sigma_ri^2
    var_e.new = sigma_e^2
    
    length_bi <- length(sigma_bi)
    Cor_bi <- diag(1, length_bi)
    Cor_bi[lower.tri(Cor_bi)] <- R_bi
    Cor_bi <- Cor_bi + t(Cor_bi) - diag(diag(Cor_bi))  # Make symmetric
    var_bi.new <- diag(sigma_bi) %*% Cor_bi %*%  diag(sigma_bi)
    
    length_bij <- length(sigma_bij)
    Cor_bij <- diag(1, length_bij)
    Cor_bij[lower.tri(Cor_bij)] <- R_bij
    Cor_bij <- Cor_bij + t(Cor_bij) - diag(diag(Cor_bij))  # Make symmetric
    var_bij.new <- diag(sigma_bij) %*% Cor_bij %*%  diag(sigma_bij)
    
    
    SS <- vector("list", n.tp)
    CIF <- vector("list", n.tp)
    for (i in seq_len(n.tp)) {
      # Step 2: simulate new random effects values
      p.b <- proposed.b[[i]][m, ]
      
      dmvt.old <- dmvt(x = b.old[i,], mu = modes.b[i, ], invSigma = invVars.b[[i]], 
                       df = 4, log = TRUE)
      
      dmvt.prop <- dmvt.proposed[[i]][m]
      a <- min(exp(log.posterior.b(p.b, yA=yA, yB=yB, X_long=X_long, X_surv=X_surv, 
                                   previous.time=previous.time, 
                                   betas1.new, betas2.new, gammas.new, xi.new, 
                                   alpha1.new, alpha2.new, 
                                   var_e.new, var_bi.new, var_bij.new, var_ri.new) 
                   + dmvt.old - 
                     log.posterior.b(b.old[i,], yA=yA, yB=yB, X_long=X_long, X_surv=X_surv, 
                                     previous.time=previous.time, 
                                     betas1.new, betas2.new, gammas.new, xi.new, 
                                     alpha1.new, alpha2.new, 
                                     var_e.new, var_bi.new, var_bij.new, var_ri.new)  
                   - dmvt.prop), 1)
      
      ind <- runif(1) <= a
      success.rate[m, i] <- ind
      if (!is.na(ind) && ind)
        b.new[i,] <- p.b
      b_Ai.new <- b.new[1]
      b_Bi.new <- b.new[2]
      b_Aij0.new <- b.new[3]
      b_Aij1.new <- b.new[4]
      b_Bij0.new <- b.new[5]
      b_Bij1.new <- b.new[6]
      b_ri.new <- b.new[7]
      
      
      
      # Step 3: compute Pr(T > t_k | T > t_{k - 1}; theta.new, b.new)
      S.last <- S_t(pred.time, X_surv, betas1.new, betas2.new, gammas.new, xi.new, 
                    alpha1.new, alpha2.new, 
                    b_Ai.new, b_Bi.new, b_Aij0.new, b_Aij1.new, b_Bij0.new, b_Bij1.new, b_ri.new)
      
      S.last <- rep(S.last, length(survTimes))
      S.pred <- numeric(length(survTimes))
      for (l in seq_along(S.pred))
        S.pred[l] <- S_t(survTimes[l], X_surv, betas1.new, betas2.new, gammas.new, xi.new, 
                         alpha1.new, alpha2.new, 
                         b_Ai.new, b_Bi.new, b_Aij0.new, b_Aij1.new, b_Bij0.new, b_Bij1.new, b_ri.new)
      
      SS[[i]] <- S.pred/S.last
      CIF[[i]] <- 1 - SS[[i]]
    }
    b.old <- b.new
    b.matrix[m,] <- b.new
    out[m,] <- SS[[i]]
  }
  
  pred.cil<- c()
  pred.ciu <- c()
  pred.median <- c()
  pred.mean <- c()
  for (i in 1:length(survTimes)) {
    pred.cil[i] <- quantile(out[,i], probs = CI.levels[1])
    pred.ciu[i] <- quantile(out[,i], probs = CI.levels[2])
    pred.median[i] <- median(out[,i])
    pred.mean[i] <- mean(out[,i])
    
  }
  
  

  #### return output ####
  result <- cbind(
    times = survTimes,
    "Mean" = colMeans(out, na.rm = TRUE),
    "Median" = pred.median,
    "Lower" = pred.cil,
    "Upper" = pred.ciu
  )
  
  return(result)
}


