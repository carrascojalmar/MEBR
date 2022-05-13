rlike <- function(nTheta,w){
  mu.w <- nTheta[1]
  phi.w <- nTheta[2]

  shape1 <- mu.w*phi.w
  shape2 <- (1-mu.w)*phi.w

  func <- dbeta(w,shape1,shape2,log=TRUE)

  op <- sum(func)
  return(-op)
}
integrand <- function(x,formula,theta,theta0,mu.e,dataSet){

  q <- length(dataSet)

  Y <- dataSet[1]
  X <- cbind(1,dataSet[2:(q-1)])
  W <- dataSet[q]

  phi <- theta[1]
  beta <- theta[2:(q + 1)]

  mu.w <- theta0[1]
  phi.w <- theta0[2]

  ll <- sapply(x,function(xi)
  {
    eta <- X%*%beta[1:(q-1)]+beta[q]*xi
    rho <- phi

    mu <- exp(eta)/(1+exp(eta))
    phi <- exp(rho)

    shape.A <- mu*phi
    shape.B <- (1-mu)*phi

    temp1 <- dbeta(Y,shape.A,shape.B)

    #.........................................................

    mu.star <- (mu.e-mu.w)/(mu.e*(1-mu.w))
    phi.star <- (1-mu.w)*phi.w

    temp2 <- ((1-xi)^(mu.star*phi.star-1))*((xi-W)^((1-mu.star)*phi.star-1))

    func<-temp1*temp2

    return(func)
  })

  return(ll)
}
l.pl <- function(theta,formula,mu.e,dataSet){

  DataCompleto <- dataSet

  Y <- model.response(data = model.frame(formula,dataSet))
  X <- model.matrix(formula, dataSet)

  p <- dim(X)[2]
  W <- X[,p]

  op.rlike <- optim(fn=rlike,gr=NULL,par=c(mean(W),var(W)),
                    method="BFGS",hessian=T,w=W)

  n <- length(Y)

  temp <- matrix(NA,nrow=n)

  for(i in 1:n){
    dataSet<-as.numeric(dataSet[i,])
    temp0 <- integrate(integrand,lower=dataSet[p],upper=1,
                       formula=formula,theta=theta,
                       theta0=op.rlike$par,dataSet=dataSet,
                       mu.e=mu.e,stop.on.error = FALSE)$value
    temp[i] <- -op.rlike$value+log(temp0)
    dataSet <- DataCompleto
  }

  ll<-sum(temp)
  return(-ll)
}
lvero <- function(theta,formula,mu.e,dataSet){

  Y <- model.response(data = model.frame(formula,dataSet))
  X <- model.matrix(formula, dataSet)

  nP <- dim(X)[2]+1
  nN <- 2

  theta1 <- theta[1:nP]

  a <- nP+1
  b <- nP+nN

  theta2 <- theta[a:b]

  alldata <- dataSet

  temp <- matrix(NA,nrow=length(dataSet[,1]))

  for(i in 1:length(dataSet[,1]))
  {
    dataSet<-as.numeric(dataSet[i,])
    temp0 <- integrate(integrand,lower=dataSet[ncol(X)],
                       upper=1,formula=formula,theta=theta1,
                       theta0=theta2,mu.e=mu.e, dataSet=dataSet,
                       stop.on.error = FALSE)$value

    mu.w.s <- theta2[1]
    phi.w.s <- theta2[2]

    shape1 <- mu.w.s*phi.w.s
    shape2 <- (1-mu.w.s)*phi.w.s

    temp[i] <- dbeta(Y,shape1,shape2,log=T)+log(temp0)
    dataSet <- alldata
  }

  ll<-sum(temp)
  return(-ll)
}
npl.mle <- function(formula,mu.e,data){
  default_opts <- callr::r(function(){options()})
  options(warn=-1,scipen=0,digits = 3)
  #
  dataSet <- data

  #
  pp <- ncol(dataSet)
  w.temp <- dataSet[,pp]
  op.temp <- optim(fn=rlike,par=c(mean(w.temp),var(w.temp)),method="BFGS",
                   hessian=T,w=w.temp)
  #

  op.aux <- betareg(formula=y~z+w|1,data=dataSet)
  par.ini <- c(op.aux$coefficients$precision,op.aux$coefficients$mean)
  func <- optim(fn=l.pl,formula=formula,
      par=par.ini,method="BFGS",hessian=T,
      mu.e=mu.e,dataSet=dataSet)
  #

  nP <- length(func$par)
  nN <- 2
  li <- nP+1
  ls <- nP+nN

  I.11 <- func$hessian
  S.11 <- op.temp$hessian

  I.12 <- hessian(lvero,x=c(func$par,op.temp$par), formula=formula,
                  mu.e=mu.e,dataSet=dataSet)[1:nP,li:ls]
  I.11.inv <- solve(I.11)

  V <- I.11.inv+I.11.inv%*%(I.12%*%solve(S.11)%*%t(I.12))%*%I.11.inv
  ep <- sqrt(diag(V))
  zstat <- func$par/ep
  pvalue <- 2*(1-pnorm(abs(zstat),mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))

  final.pseudo <- data.frame(func$par,ep,zstat,pvalue)

  lab.temp <- names(dataSet)[-c(1,ncol(dataSet))]

  rownames(final.pseudo) <- c("(phi)","(Intercept)", lab.temp ,"x")
  colnames(final.pseudo) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  mu.star <- (mu.e-op.temp$par[1])/(mu.e*(1-op.temp$par[1]))

  xhat <- 1-(1-w.temp)*mu.star

  X.temp <- cbind(model.matrix(formula, dataSet)[,1:(pp-1)],xhat)

  eta.NPL <- X.temp%*%func$par[2:(pp+1)]
  g.y <- log(dataSet[,1]/(1-dataSet[,1]))
  R.squared.PL <- (cor(eta.NPL,g.y))^2

  if(!is.na(R.squared.PL))
    cat("\nPseudo R-squared:", R.squared.PL,"\n")

  print(final.pseudo,class = TRUE)

  options(default_opts)
}


