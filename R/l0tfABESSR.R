
subMatrix <- function(X, A) {
  return(as.matrix(X[, A]))
}


#' @title Create a Difference Matrix using the ABESS approach
#' @description This function creates a difference matrix of a specified size and order based on ABESS idea.
#' @param n The number of rows and columns of the matrix.
#' @param q The order of the difference (0 or 1).
#' @return A matrix of order \code{n} representing the difference of order \code{q}.
#' @examples
#' \dontrun{
#' DiffMatrix(4, 0)
#' DiffMatrix(10, 1)
#' }
#' @export
DiffMatrix <- function(n, q) {
  if(q == 0) {
    phi <- matrix(1, n, n)
    for(i in 1:n) {
      for(j in 1:n) {
        if(i < j) phi[i, j] <- 0
      }
    }
  } else if(q == 1) {
    phi <- matrix(1, n, n)
    for(i in 1:n) {
      for(j in 2:n) {
        if(i == j) {
          phi[i, j] <- 1
        } else {
          phi[i, j] <- max(i - j + 1, 0)
        }
      }
    }
  }
  return(phi)
}


SplicingsSparseR <- function(A=A,I=I,k=k,y=y,q=q,tau=0){
  n <- length(y)
  All <- 1:n
  S <- 1:(q+1)
  X <- DiffMatrix(n,q)
  kmax <- k
  Al <- union(S,A)
  beta <- as.numeric(n)
  beta[Al] <- solve(t(subMatrix(X,Al))%*%subMatrix(X,Al))%*%t(subMatrix(X,Al))%*%y
  beta[I] <- 0
  eps1 <- as.numeric(n)
  eps2 <- as.numeric(n)
  for(i in 1:n){
    h <- t(subMatrix(X,i))%*%subMatrix(X,i)/n
    d <- t(subMatrix(X,i))%*%(y-X%*%beta)/n
    eps1[i] <- h*(beta[i])^2/2
    eps2[i] <- h/2*(d/h)^2
  }
  L0 <- t(y-X%*%beta)%*%(y-X%*%beta)/(2*n)
  L <- L0
  if(length(A)>length(I)){
    kmax <- length(I)
  }
  for(j in 1:kmax){
    A0 <- as.vector(A[order(eps1[A],decreasing = F)[1:j]])
    I0 <- as.vector(I[order(eps2[I],decreasing = T)[1:j]])
    Anew <- c(setdiff(A,A0),I0)
    Inew <- c(setdiff(I,I0),A0)
    Alnew <- union(Anew,S)
    betanew <- as.numeric(n)
    betanew[Alnew] <- solve(t(subMatrix(X,Alnew))%*%subMatrix(X,Alnew))%*%t(subMatrix(X,Alnew))%*%y
    betanew[Inew] <- 0
    Lnew <- t(y-X%*%betanew)%*%(y-X%*%betanew)/(2*n)
    L0 <- L
    if(L-Lnew>0){
      A <- Anew
      I <- Inew
      beta <- betanew
      L <- Lnew
      for(i in 1:n){
        h <- t(subMatrix(X,i))%*%subMatrix(X,i)/n
        d <- t(subMatrix(X,i))%*%(y-X%*%beta)/n
        eps1[i] <- h*(beta[i])^2/2
        eps2[i] <- h/2*(d/h)^2
      }
      if(L0-L<tau){
        return(list(A=A,I=I,beta=beta))
      }
    }
  }
  return(list(A=A,I=I,beta=beta))
}


L0tfkSparseR <- function(A0=A0,y=y,q=q,max.step=50,k=k){
  n <- length(y)
  All <- 1:n
  S <- 1:(q+1)
  I0 <- setdiff(setdiff(All,S),A0)
  X <- DiffMatrix(n,q)
  for(i in 1:max.step){
    fit <- SplicingsSparseR(A=A0,I=I0,k=k,y=y,q=q)
    A <- fit$A
    I <- fit$I
    if(identical(A,A0) & identical(I,I0)){
      break
    }else{
      A0 <- A
      I0 <- I
    }
  }
  Al <- union(A,S)
  beta <- as.numeric(n)
  beta[Al] <- solve(t(subMatrix(X,Al))%*%subMatrix(X,Al))%*%t(subMatrix(X,Al))%*%y
  beta[I] <- 0
  yhat <- X%*%beta
  return(list(A=A,I=I,yhat=yhat,beta=beta))
}


#' @title Sequential L0 Trend Filtering based on ABESS with R
#' @description This function implements a sequential L0 trend filtering algorithm for sparse signal recovery. It iteratively selects the most significant variables based on a given criterion, such as the Bayesian Information Criterion (BIC).
#' @param y A numeric vector of response variables.
#' @param kmax The maximum number of variables to include in the model.
#' @param q The order of the difference matrix used in the trend filtering.
#' @return A list containing the following elements:
#' \itemize{
#'   \item{kopt}{The optimal number of variables chosen by the algorithm.}
#'   \item{beta}{Coefficients of the model for each iteration up to \code{kmax}.}
#'   \item{A.all}{A list of indices of variables included in the model for each iteration up to \code{kmax}.}
#'   \item{A}{The final set of indices of variables included in the model.}
#'   \item{yhat}{Predicted values using the final model.}
#' }
#' @examples
#' \dontrun{
#' y <- rnorm(100)
#' result <- l0tfABESSR(y = y, kmax = 10, q = 1)
#' plot(y)
#' lines(result$yhat, col = "blue")
#' }
#' @export
l0tfABESSR <- function(y=y,kmax=kmax,q=q){
  n <- length(y)
  All <- 1:n
  S <- 1:(q+1)
  A0 <- NULL
  I0 <- setdiff(setdiff(All,S),A0)
  k <- kmax
  X <- DiffMatrix(n,q)
  betaAll <- NULL
  mse <- as.numeric(k)
  bic <- as.numeric(k)
  A.all <- list()
  for(j in 1:k){
    Al <- union(S,A0)
    beta <- as.numeric(n)
    beta[Al] <- solve(t(subMatrix(X,Al))%*%subMatrix(X,Al))%*%t(subMatrix(X,Al))%*%y
    beta[I0] <- 0
    eps2 <- as.numeric(n)
    for(i in 1:n){
      h <- t(subMatrix(X,i))%*%subMatrix(X,i)/n
      d <- t(subMatrix(X,i))%*%(y-X%*%beta)/n
      eps2[i] <- h/2*(d/h)^2
    }
    A0 <- union(A0,I0[which.max(eps2[I0])])
    result <- L0tfkSparseR(A0=A0,y=y,q=q,k=j)
    A.all[[j]] <- result$A
    betaAll <- cbind(betaAll,result$beta)
    A0 <- result$A
    I0 <- setdiff(setdiff(All,S),A0)
  }
  for(i in 1:k){
    mse[i] <- mean((y-X%*%betaAll[,i])^2)
    bic[i] <- n*log(mse[i]) + 2*log(n)*(i+q+1)
  }
  kopt <- which.min(bic)
  beta <- betaAll[,kopt]
  A <- sort(as.vector(A.all[[kopt]]))
  yhat <- X%*%beta
  return(list(kopt=kopt,beta=beta,A.all=A.all,A=A,yhat=yhat))
}



#' @title Simulate 'Blocks' Data
#' @description This function simulates data based on a 'blocks' signal, commonly used in benchmarking signal processing methods.
#' @param n Number of data points to generate.
#' @param sigma Standard deviation of the noise added to the signal.
#' @param seed An optional seed for random number generation to make results reproducible.
#' @return A list containing the simulated data (y), the true signal without noise (y0), the x-coordinates (x), change points (tau), and indices of change points (SetA).
#' @examples
#' \dontrun{
#' simData <- SimuBlocks(n = 100, sigma = 0.1)
#' plot(simData$x, simData$y)
#' lines(simData$x, simData$y0, col = "red")
#' }
#' @import stats 
#' @export
SimuBlocks <- function (n, sigma = 0.1, seed=NA) 
{
  if (!is.na(seed)) set.seed(seed)
  x = seq(1/n, 1,length.out = n)
  tau <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 
           0.76, 0.78, 0.81)
  h <- c(0, 4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)/5
  A = sapply(tau, function(z) which(x>=z)[1])
  tau1=c(0,tau,1)
  y0 = 0*x
  for (j in 1:(length(A)+1)) y0[x>tau1[j] & x<=tau1[j+1]] = h[j]
  y <- y0 + sigma*rnorm(n)
  return(list(y = y, x = x, y0 = y0, tau = tau, SetA = A))  
}


#' @title Simulate 'Wave' Data
#' @description Generates data based on a 'wave' pattern. Useful for testing and benchmarking in signal processing and statistical analysis.
#' @param n Number of data points.
#' @param sigma Noise level.
#' @param seed Random seed for reproducibility.
#' @return A list of the generated y values, the true signal (y0), x-coordinates (x), change points (tau), and indices of change points (SetA).
#' @examples
#' \dontrun{
#' simWave <- SimuWave(n = 100, sigma = 0.1)
#' plot(simWave$x, simWave$y)
#' lines(simWave$x, simWave$y0, col = "blue")
#' }
#' @import stats 
#' @export
SimuWave <- function (n, sigma = 0.1, seed=NA) 
{
  if (!is.na(seed)) set.seed(seed)
  x = seq(1/n, 1,length.out = n)
  tau = c(256, 512, 768, 1024, 1152, 1280, 1344, 1408)/1472
  nknot = length(tau)
  A = sapply(tau, function(z) which(x>=z)[1])
  tau=x[A]
  tau1 = c(0, tau, 1)
  h =  cumsum(c(1, (-1)^(1:nknot)*(1:nknot)))
  phi = rep(1, n)
  for (j in 1:(nknot+1)) phi = cbind(phi,pmin(pmax(x-tau1[j], 0), tau1[j+1]-tau1[j]))
  y0 = as.vector(phi%*%c(0, h))
  y <- y0 + sigma*rnorm(n)
  return(list(y = y, x = x, y0 = y0, tau = tau, SetA = A))  
}


#' @title L0 Trend Filtering is calculated using the ABESS approach and R functions
#' @description This function applies L0 trend filtering to simulated data. It supports different types of simulated data, such as 'Blocks' and 'Wave'.
#' @param sigma The standard deviation of the noise added to the simulated data.
#' @param dgm The type of data generation model to use. Possible values are "Blocks" or "Wave".
#' @param n The number of data points to be generated in the simulation.
#' @param seed An optional seed for random number generation to make the simulation reproducible.
#' @return The result of the L0 trend filtering algorithm, typically including the optimal number of variables chosen, the coefficients of the model, and other relevant information.
#' @examples
#' \dontrun{
#' # Simulating 'Blocks' data
#' blocksData <- SimuBlocks(n = 100, sigma = 0.1, seed=123)
#' resultBlocks <- Simul0tfABESSR(sigma = 0.1, dgm = "Blocks", n = 100, seed=123)
#' plot(blocksData$y0)
#' lines(resultBlocks$yhat,col="red")
#' 
#' # Simulating 'Wave' data
#' waveData <- SimuWave(n = 100, sigma = 0.1, seed=123)
#' resultWave <- Simul0tfABESSR(sigma = 0.1, dgm = "Wave", n = 100, seed=123)
#' plot(waveData$y0)
#' lines(resultWave$yhat,col="red")
#' }
#' @export
Simul0tfABESSR <- function(sigma=0.1,dgm="Blocks",n=n,seed=NA){
  if (dgm=="Blocks"){
    data = SimuBlocks(n=n, sigma=sigma, seed=seed)
    q=0; kmax=20;
  }
  if (dgm=="Wave"){
    data = SimuWave(n=n, sigma=sigma, seed=seed)
    q=1; kmax=20;
  }
  resL0 = l0tfABESSR(as.numeric(data$y), kmax=kmax, q=q)
  return(resL0)
}


#' @title Benchmark l0tfABESSR and l0tfABESSC functions
#' @name BenchmarksABESS
#' @description Use R package \code{microbenchmark} to compare the performance of R functions (\code{l0tfABESSR}) and RCpp functions (\code{l0tfABESSC}).
#' @examples
#' \dontrun{
#' tm1 <- microbenchmark::microbenchmark(
#'   rnR = Simul0tfABESSR(sigma=0.1,dgm="Blocks",n=200),
#'   rnC = Simul0tfABESSC(sigma=0.1,dgm="Blocks",n=200),
#'   times = 100
#' )
#' print(summary(tm1)[,c(1,3,5,6)])
#' tm2 <- microbenchmark::microbenchmark(
#'   rnR = Simul0tfABESSR(sigma=0.1,dgm="Wave",n=200),
#'   rnC = Simul0tfABESSC(sigma=0.1,dgm="Wave",n=200),
#'   times = 100
#' )
#' print(summary(tm2)[,c(1,3,5,6)])
#' }
#' @import ggplot2
#' @import boot 
#' @import bootstrap
#' @import DAAG
#' @import coda
#' @import knitr
#' @import microbenchmark
#' @import Rcpp 
#' @import RcppArmadillo
#' @useDynLib SA23204173
NULL

