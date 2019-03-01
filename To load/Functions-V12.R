
# add the factor half in front of the smoothness penalty
# see whether the covariance matrix coincides with the one from the GAM output

library(mgcv)
library(Matrix)

# Function modified for any dimention of covariate matrix Z  - March 1, 2018

# Simulation function
BSMethGammSim <- function(n, posit, theta.0, beta, random.eff = F, mu.e=0, 
                          sigma.ee=1,p0 = 0.003, p1 = 0.9, X, Z,binom.link="logit"){
  if( !is.matrix((Z)) ) message ("covariate Z is not a matrix")
#  if( !is.matrix(beta) ) message ("the covariate effect parameter beta is not a matrix")
  
  if( !(nrow(X)==nrow(Z) & nrow(X) == n) ) message("Both X and Z should have n rows")
  if( !(ncol(X)==nrow(theta.0) & ncol(X) ==nrow (beta) & ncol(X)==length(posit) )) message ("The columns of X should be the same as length of beta theta.0 and posit; They all equals to the number of CpGs")
  if( ncol(beta)!= ncol(Z)) message("beta and Z should have the same dimentions")
  
  # the random effect term
  if(random.eff == T){
    my.e <- rnorm(n, mean=mu.e, sd = sqrt(sigma.ee))
  }else{
    my.e <- rep(mu.e, n)
  }
  
  my.theta <- t(sapply(1:n, function(i){
    theta.0 + rowSums(sapply(1:ncol(Z), function(j){Z[i,j] * beta[,j]})) + my.e[i]
  }))
  
  
  # Transform my.theta to my.pi for each (i, j)
  
  my.pi <- t(sapply(1:nrow(my.theta), function(i){
    #exp(my.theta[i,])/(1+exp(my.theta[i,]))
    binomial(link=binom.link)$linkinv(my.theta[i,])
  }))
  #  Generate S-ij based on the my.pi and my.Z
  #---------------------------------------------------#
  my.S <- my.pi
  for ( i in 1:nrow(my.S)){
    for ( j in 1:ncol(my.S)){
      my.S[i,j] <- rbinom (1, size = X[i,j], prob= my.pi[i,j])
    }
  }
  #---------------------------------------------------#
  # Generate Y-ij based on the S-ij and the error rate (1-p1) and p0
  #---------------------------------------------------#
  my.Y <- my.S
  for ( i in 1:nrow(my.Y)){
    for ( j in 1:ncol(my.Y)){
      my.Y[i,j] <- sum(rbinom(my.S[i,j], size =1, prob=p1)) +
        sum(rbinom(X[i,j]-my.S[i,j], size = 1, prob=p0))
    }
  }
  out = list(S = my.S, Y = my.Y, theta = my.theta, pi = my.pi)
}
# EM estimate function
BSMethEM <- function (data, p0 = 0.003, p1 = 0.9, epsilon = 10^(-6),
                      epsilon.lambda = 10^(-3), maxStep = 200, n.k , detail=T, binom.link = "logit",
                      method="REML"){
  if( (ncol(data)-4) ==1){
    Z <- matrix(data[,-(1:4)], nrow = nrow(data))
  }else{
    Z <- data[,-(1:4)]
  }
  colnames(Z) <- colnames(data)[-c(1:4)]
  # The gam formula corresponding to the Z
  formula.z.part <- sapply(1:(ncol(data)-4), function(i){ 
    paste0("s(Posit, k = n.k[",i+1,"], fx=F, bs=\"cr\", by = Z[,", i, "])"  ) } )
  my.covar.fm <- paste(c("s(Posit, k=n.k[1], fx=F, bs=\"cr\")", formula.z.part), collapse="+")
  
  # Fit gam for the initial value
  gam.int <- gam(as.formula( paste0("Y/X ~", my.covar.fm)), family =binomial(link = binom.link),weights=X, data = data, method = method)
  # Estimates
  old.pi.ij <- gam.int$fitted.values;old.par <-gam.int$coefficients;lambda <- gam.int$sp 
  # Update
  out <-  BSMethEMUpdate (data, old.pi.ij, p0 = p0, p1 = p1, n.k=n.k, binom.link = binom.link, method = method)
  new.par<-out$par; new.lambda <- out$lambda;new.pi.ij <- out$pi.ij
  i <- 1; Est.points <- rbind(c(old.par, lambda), c(new.par, new.lambda))
# Do the iteration
  # The stopping criterion is that estimator a and b are close enough
  # I exclude the criterio that lambda-a and lambda-b are close
  while ( sum((old.par - new.par)^2) > epsilon & i < maxStep & sum((lambda-new.lambda)^2) > epsilon.lambda){
    i <- i +1; 
    old.par <- new.par
    old.pi.ij <- new.pi.ij
    
    out <-  BSMethEMUpdate (data, old.pi.ij, p0 = p0, p1 = p1, n.k = n.k, binom.link = binom.link, method = method)
    new.par<-out$par; new.lambda <- out$lambda;new.pi.ij <- out$pi.ij
    
    Est.points <- rbind(Est.points, c(new.par,new.lambda))
    if(detail){
      print(paste0("iteration", i))
    }
  }
  FinalGamObj <- out$gam.Obj
 
  #--------------------------------------------
  # Calculate SE of alpha
  #-------------------------------------------
  # the part to calculate the SE of alpha
  my.design.matrix <-  model.matrix.gam(FinalGamObj);Y <- data$Y ; X <- data$X
  pred.pi <- predict.gam(FinalGamObj, type ="response" )
  # Q1: the second partial derivative w.r.t alpha^2 
  # Q2: the second derivative w.r.t alpha & alpha_star
  
 
  #Q1 <- matrix(NA, nrow =length(new.par), ncol = length(new.par))
  #for ( l in 1: nrow(Q1)){
    #for (m in 1:ncol(Q1)){
  #  Q1[l, m] <-  sum(-X * pred.pi * (1-pred.pi) * my.design.matrix[, m]*my.design.matrix[,l] ) 
  #  }
  #}
  res <- outer( 1:length(new.par), 1:length(new.par), Vectorize(function(l,m)  sum(-X * pred.pi * (1-pred.pi) * my.design.matrix[, m]*my.design.matrix[,l] )) )
  smoth.mat <-sapply(1:(ncol(Z)+1), function(i){FinalGamObj$smooth[[i]]$S[[1]] * new.lambda[i]})
  smoth.mat[[length(smoth.mat) + 1]] <- 0  # assume the lambda for the constant of the interaction is 0 -- no penalization
  span.penal.matrix <- as.matrix(bdiag( smoth.mat[c(length(smoth.mat), (1:(length(smoth.mat)-1)))] ))
    Q1_with_lambda <- res - span.penal.matrix  # delete the factor 2
    Q1_no_lambda <- res
  
  Q2 <- outer(1:length(new.par), 1:length(new.par), Vectorize(function(l,m){
    term1 <- Y*p1*p0/(p1*pred.pi + p0 * (1-pred.pi))^2 + (X-Y)*(1-p1)*(1-p0)/((1-p1)*pred.pi + (1-p0) * (1-pred.pi))^2
    sum(term1 * pred.pi * (1-pred.pi) * my.design.matrix[,m] * my.design.matrix[,l])
  }))
  #-------------------------------------------------------------------------
  # from the variance-covariance of alpha, to report
  # 1. var(beta_p(t))
  # 2. Report the chi-square statistic for each covariate
  # 3. Report the p-value for each covariate
  # Added in June 13-18, 2018. 
  #-------------------------------------------------------------------------
  
  # ------ 3: estimate of beta(t) --- # 
  uni.pos <- unique(data$Posit); uni.id <- match(uni.pos, data$Posit)
  BZ <- my.design.matrix[uni.id, 1:n.k[1]]
  BZ.beta = lapply(1:ncol(Z), function(i){mgcv::smooth.construct(s(Posit, k = n.k[i+1], fx = F, bs = "cr") , data = data[uni.id,], knots = NULL)$X })
  
  cum_s <-cumsum(n.k) 
  alpha.sep <- lapply(1:ncol(Z), function(i){new.par[ (cum_s[i]+1):cum_s[i+1]]}); alpha.0 <- new.par[1:n.k[1]]
  
  Beta.out<-  cbind( BZ %*% alpha.0, sapply(1:ncol(Z), function(i){BZ.beta[[i]] %*% alpha.sep[[i]]}))
  
  colnames(Beta.out) <- c("Intercept", colnames(Z))
  
  
  # Simply check to see 
  
  
  
  
  tryCatch({
    var.cov.alpha <-solve(-(Q1_with_lambda + Q2))
 
  #-------1: Estimate SE of beta(t) --- move to another function plot_BSMethEM ?
  #------    keep this part of the code here, for now.
  
  SE.out <- vector("list", (ncol(Z)+1)); names(SE.out) <- c("Intercept", colnames(Z))
 
 
  # SE of beta.0
 
  var.alpha.0 <- var.cov.alpha[1:n.k[1], 1:n.k[1]]
  var.beta.0 <-  BZ %*% var.alpha.0 %*% t(BZ)
  #SE.out[[1]] <- sqrt(diag(var.beta.0))
 
  # SE of the effect of Zs [beta.1(t), beta.2(t), beta.3(t) ...]
 
  var.alpha.sep <- lapply(1:ncol(Z), function(i){var.cov.alpha[ (cum_s[i]+1):cum_s[i+1],(cum_s[i]+1):cum_s[i+1]]})
  var.beta <- lapply(1:ncol(Z), function(i){BZ.beta[[i]] %*% var.alpha.sep[[i]] %*% t(BZ.beta[[i]])})
  
  SE.out <- cbind(sqrt(diag(var.beta.0)), sapply(1:ncol(Z), function(i){sqrt(diag(var.beta[[i]]))}))
  rownames(SE.out) <- uni.pos; colnames(SE.out) <-  c("Intercept", colnames(Z))
  SE.pos <- uni.pos
   
  
  # -------2: The chi-square statistics and P values for each covariate
   
   chi.sq <- c( t(alpha.0) %*% solve(var.alpha.0) %*% alpha.0, sapply(1:ncol(Z), function(i){t(alpha.sep[[i]]) %*% solve(var.alpha.sep[[i]]) %*% alpha.sep[[i]]})) 
   #pvalue <- 1- pchisq(chi.sq, df = n.k)
   
   edf <- FinalGamObj$edf1
   
   df = sapply(1:(ncol(Z)+1), function(i){sum(edf[(cum_s[i]-n.k[i]+1):cum_s[i]])})
   pvalue <- 1- pchisq(chi.sq, df = df)
   #pvalue <- 1- pchisq(chi.sq, df = n.k)
   pvalue.log <-  pchisq(chi.sq, df = df, log.p = T, lower.tail = F) 
   
  names(pvalue.log) <- names(pvalue) <- names(chi.sq) <- c("Intercept", colnames(Z))
   # ----- try other covariance matrix and other edf
   #vcov.gam(FinalGamObj)  # results coincide with Vp
   #all(var.cov.alpha - FinalGamObj$Vp < 0.00000000000001)
   #all(var.cov.alpha - FinalGamObj$Ve < 0.00000000000001)
   #all(var.cov.alpha - FinalGamObj$Vc < 0.00000000000001)
   
   #res.ve <- chi.sq.function(FinalGamObj$Ve,alpha.0 = alpha.0, alpha.sep = alpha.sep, n.k = n.k, edf = edf);
   #res.vp <- chi.sq.function(FinalGamObj$Vp,alpha.0 = alpha.0, alpha.sep = alpha.sep, n.k = n.k, edf = edf);
   #res.vc <- chi.sq.function(FinalGamObj$Vc,alpha.0 = alpha.0, alpha.sep = alpha.sep, n.k = n.k, edf = edf)
   
   #chi.sq.wald <- cbind(res.ve[, "chi.sq"], res.vp[, "chi.sq"], res.vc[, "chi.sq"])
   #p.value.wald <- cbind(res.ve[, "pvalue"], res.vp[, "pvalue"], res.vc[, "pvalue"])
   
   #colnames(chi.sq.wald) <- colnames(p.value.wald) <-  c("Ve", "Vp", "Vc")
   
   
  
   out = list(est = new.par, lambda = new.lambda, est.pi = new.pi.ij, ite.points = Est.points, FinalGamObj=FinalGamObj,
              cov1 = var.cov.alpha, chi.sq = chi.sq, pvalue = pvalue, pvalue.log = pvalue.log, SE.out = SE.out, SE.pos = SE.pos, Beta.out = Beta.out )
   return(out)
  },
  error = function(err){
    return (out = list(est = new.par, lambda = new.lambda, est.pi = new.pi.ij, ite.points = Est.points, FinalGamObj=FinalGamObj,
                       cov1 = diag(1, length(new.par)), Beta.out = Beta.out ) )
  }
  )
  # New stuff stopped here
  #---------------------------------------------------------------------------
  
 # tryCatch({var.cov.lambda <-solve(-(Q1_with_lambda + Q2)) 
 #  out = list(est = new.par, lambda = new.lambda, est.pi = new.pi.ij, ite.points = Est.points, FinalGamObj=FinalGamObj,
 #            cov1 = var.cov.lambda, chi.sq = chi.sq, pvalue = pvalue, SE.out = SE.out, SE.pos = SE.pos)
 #  return(out)
 # },
 #error = function(err){return (out = list(est = new.par, lambda = new.lambda, est.pi = new.pi.ij, ite.points = Est.points, FinalGamObj=FinalGamObj,
 #                                          cov1 = diag(1, length(new.par)))) }
#)
  
}
# EM estimate update function; the update step used in the "BSMethEM" function
BSMethEMUpdate <- function (data, pi.ij, p0 , p1 , n.k, binom.link, method){
  if( !(nrow(data)==length(pi.ij)) ) message("The row of data should be compatible with the length of initial value pi.ij")
  # The E-step
  # Calculate the "posterior" probability
  eta.1 <-  p1 *pi.ij /(p1*pi.ij + p0*(1-pi.ij)) # posterior probability given an observed methylated rates, what is the probability that the reads are truely methylated
  eta.0 <-  (1-p1) * pi.ij /((1-p1) * pi.ij + (1-p0)*(1-pi.ij))
  
  Y <- data$Y ; X <- data$X
  E.S <- Y * eta.1 + (X-Y) * eta.0
  
  # The M step
  if(ncol(Z)==1){
    Z <- matrix(data[,-(1:4)], nrow = nrow(data))
  }else{
    Z <- data[,-(1:4)]
  }
  # The gam formula corresponding to the Z
  formula.z.part <- sapply(1:(ncol(data)-4), function(i){ 
    paste0("s(Posit, k = n.k[",i+1,"], fx=F, bs=\"cr\", by = Z[,", i, "])"  ) } )
  my.covar.fm <- paste(c("s(Posit, k=n.k[1], fx=F, bs=\"cr\")", formula.z.part), collapse="+")
  
  gam.int.see <- mgcv::gam(as.formula( paste0("E.S/X ~", my.covar.fm)), family = binomial(link=binom.link),weights=X,
                           data = data, method=method)
  
  out <- list( pi.ij = gam.int.see$fitted.values, par = gam.int.see$coefficients,
               lambda = gam.int.see$sp, gam.Obj = gam.int.see)
  return(out)
}

BSMethGAM <- function (data,n.k , detail=T, binom.link = "logit", method ="REML"){
  if(ncol(Z)==1){
    Z <- matrix(my.span.dat[,-(1:4)], nrow = nrow(my.span.dat))
  }else{
    Z <- my.span.dat[,-(1:4)]
  }
  
  formula.z.part <- sapply(1:(ncol(data)-4), function(i){ 
    paste0("s(Posit, k = n.k[",i+1,"], fx=F, bs=\"cr\", by = Z[,", i, "])"  ) } )
  my.covar.fm <- paste(c("s(Posit, k=n.k[1], fx=F, bs=\"cr\")", formula.z.part), collapse="+")
  
  # Fit gam for the initial value
  gam.int <- gam( as.formula( paste0("Y/X ~", my.covar.fm)), family =binomial(link = binom.link),weights=X, data = data)
  
  return(gam.int)
}


chi.sq.function <- function(var.cov,alpha.0, alpha.sep, n.k, edf ){
  cum_s <-cumsum(n.k)
  var.alpha.0 <- var.cov[1:n.k[1], 1:n.k[1]]
  var.alpha.sep <- lapply(1:(length(n.k)-1), function(i){var.cov[ (cum_s[i]+1):cum_s[i+1],(cum_s[i]+1):cum_s[i+1]]})
  chi.sq <- c( t(alpha.0) %*% solve(var.alpha.0) %*% alpha.0, sapply(1:(length(n.k)-1), function(i){t(alpha.sep[[i]]) %*% solve(var.alpha.sep[[i]]) %*% alpha.sep[[i]]})) 
  
  pvalue <- 1- pchisq(chi.sq, df = sapply(1:((length(n.k)-1)+1), function(i){sum(edf[(cum_s[i]-n.k[i]+1):cum_s[i]])}))
  cbind(chi.sq, pvalue)
}
