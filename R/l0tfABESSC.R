
SplicingsSparseC <- function(A=A,I=I,k=k,y=y,q=q,tau=0){
  n <- length(y)
  All <- 1:n
  S <- 1:(q+1)
  X <- DiffMatrix(n,q)
  kmax <- k
  Al <- union(S,A)
  beta <- as.numeric(n)
  Alc <- Al-1
  beta[Al] <- solveEquation(X,Alc,y)
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
    Alnewc <- Alnew-1
    betanew[Alnew] <- solveEquation(X,Alnewc,y)
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

L0tfkSparseC <- function(A0=A0,y=y,q=q,max.step=50,k=k){
  n <- length(y)
  All <- 1:n
  S <- 1:(q+1)
  I0 <- setdiff(setdiff(All,S),A0)
  X <- DiffMatrix(n,q)
  for(i in 1:max.step){
    fit <- SplicingsSparseC(A=A0,I=I0,k=k,y=y,q=q)
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
  Alc <- Al-1
  beta[Al] <- solveEquation(X,Alc,y)
  beta[I] <- 0
  yhat <- X%*%beta
  return(list(A=A,I=I,yhat=yhat,beta=beta))
}


#' @title Sequential L0 Trend Filtering based on ABESS with Rcpp
#' @description Given the response data, the maximum number of detected knots, and the order, we obtain the solution to the L0 trend filtering problem based on the idea of ABESS.
#' @param y  Response data.
#' @param kmax Maximum number of detected knots for sequential trend filtering.
#' @param q Order and integer used in creating the differential matrix, influencing the filtering process.
#' @return A list containing the optimal number of detected knots (kopt), the optimized beta coefficients, the full history of active sets (A.all), the final active set (A), and the predicted values (yhat) after the trend filtering process.
#' @examples
#' \dontrun{
#' y <- rnorm(100)
#' result <- l0tfABESSC(y = y, kmax = 10, q = 1)
#' plot(y)
#' lines(result$yhat, col = "red")
#' }
#' @export
l0tfABESSC <- function(y=y,kmax=kmax,q=q){
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
    Alc <- Al-1
    beta[Al] <- solveEquation(X,Alc,y)
    beta[I0] <- 0
    eps2 <- as.numeric(n)
    for(i in 1:n){
      h <- t(subMatrix(X,i))%*%subMatrix(X,i)/n
      d <- t(subMatrix(X,i))%*%(y-X%*%beta)/n
      eps2[i] <- h/2*(d/h)^2
    }
    A0 <- union(A0,I0[which.max(eps2[I0])])
    result <- L0tfkSparseC(A0=A0,y=y,q=q,k=j)
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

#' @title L0 trend filtering is calculated using the ABESS approach and Rcpp functions
#' @description Based on the concept of ABESS, the L0 trend filtering solution obtained using Rcpp functions estimates piecewise constant and piecewise linear functions.
#' @param sigma The standard deviation of the noise added to the simulated data.
#' @param dgm The type of data generation model to use. Possible values are "Blocks" or "Wave".
#' @param n The number of data points to be generated in the simulation.
#' @param seed An optional seed for random number generation to make the simulation reproducible.
#' @return The result of the L0 trend filtering algorithm, typically including the optimal number of variables chosen, the coefficients of the model, and other relevant information.
#' @examples
#' \dontrun{
#' # Simulating 'Blocks' data
#' blocksData <- SimuBlocks(n = 100, sigma = 0.1, seed=123)
#' resultBlocks <- Simul0tfABESSC(sigma = 0.1, dgm = "Blocks", n = 100, seed=123)
#' plot(blocksData$y0)
#' lines(resultBlocks$yhat,col="red")
#' 
#' # Simulating 'Wave' data
#' waveData <- SimuWave(n = 100, sigma = 0.1, seed=123)
#' resultWave <- Simul0tfABESSC(sigma = 0.1, dgm = "Wave", n = 100, seed=123)
#' plot(waveData$y0)
#' lines(resultWave$yhat,col="red")
#' }
#' @import RcppArmadillo
#' @import Rcpp 
#' @export
Simul0tfABESSC <- function(sigma=0.1,dgm="Blocks",n=n,seed=NA){
  if (dgm=="Blocks"){
    data = SimuBlocks(n=n, sigma=sigma, seed=seed)
    q=0; kmax=20;
  }
  if (dgm=="Wave"){
    data = SimuWave(n=n, sigma=sigma, seed=seed)
    q=1; kmax=20;
  }
  resL0 = l0tfABESSC(as.numeric(data$y), kmax=kmax, q=q)
  return(resL0)
}


