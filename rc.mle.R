rc.mle <- function(formula,data,mu.e,B){
  default_opts <- callr::r(function(){options()})
  options(warn=-1,scipen=0,digits = 3)
  #
  Y <- model.response(data = model.frame(formula,data))
  X <- model.matrix(formula, data)

  p <- ncol(X)

  W <- data[,p]

  op.rlike <- optim(fn=rlike,gr=NULL,par=c(mean(W),var(W)),
                    method="BFGS",hessian=F,w=W)

  mu.star <- (mu.e-op.rlike$par[1])/(mu.e*(1-op.rlike$par[1]))
  xhat <- 1-(1-data[,p])*mu.star
  data.RC <- cbind(data[,1:(p-1)],xhat)


  func <- betareg(formula=y~z+xhat|1, data=data.RC)

  #print(func$fitted.values)

  ep.phi <- sqrt(diag(func$vcov))[p+1]
  ep.beta <- sqrt(diag(func$vcov))[1:p]
  ep <- c(ep.phi,ep.beta)

  np <- length(func$coefficients$mean)+1

  theta.b <- matrix(NA,nrow=B,ncol=np)

  n <- length(Y)

  for(b in 1:B){
    b.i <- sample(1:n,size=n,replace=TRUE)
    data.boot <-data[b.i,]

    W.boot <- data.boot[,p]
    op.rlike.boot <- optim(fn=rlike,gr=NULL,par=c(mean(W.boot),var(W.boot)),
                      method="BFGS",hessian=F,w=W.boot)

    mu.star.boot <- (mu.e-op.rlike.boot$par[1])/(mu.e*(1-op.rlike.boot$par[1]))
    xhat.boot <- 1-(1-data.boot[,p])*mu.star.boot
    data.RC.boot <- cbind(data.boot[,1:(p-1)],xhat.boot)

    mod.boot <- betareg(formula=y~z+xhat.boot|1,data=data.RC.boot)

    theta.b[b,] <- c(mod.boot$coefficients$precision,mod.boot$coefficients$mean)

  }

  temp <- rep(apply(theta.b,2,mean),each=B)

  theta.est <- matrix(temp,nrow=B,ncol=np,byrow=FALSE)

  argu <- (theta.b-theta.est)^2

  epRC <- sqrt((1/(B-1))*apply(argu,2,sum))

  est <- c(func$coefficients$precision,func$coefficients$mean)

  zstat.RC <- est/epRC

  pvalue.RC <- 2*(1-pnorm(abs(zstat.RC),mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))

  final.rc <- data.frame(est,ep,epRC,zstat.RC,pvalue.RC)

  rownames(final.rc) <- c("(phi)","(Intercept)",names(data.RC[,2:ncol(data.RC)]))
  colnames(final.rc) <- c("Estimate", "Std. Error","Std. Error(boot)", "z value", "Pr(>|z|)")


  if(!is.na(func$pseudo.r.squared))
    cat("\nPseudo R-squared:", func$pseudo.r.squared,"\n")

  print(final.rc,class = TRUE)
  
  options(default_opts)

}
