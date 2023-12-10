#' @title Generate a Difference Matrix based on AMIAS method
#' @description This function generates a matrix used for calculating differences of a certain order. It is particularly useful in numerical methods or for creating specific matrix patterns. The matrix created by this function can be used in various applications where a difference operation of order 'q' is needed over 'n' data points.
#' @param n The number of data points.
#' @param q The order of the difference to be computed.
#' @return A matrix of size (n-q-1) by n, with elements arranged to facilitate the calculation of a 'q'-th order difference.
#' @examples
#' \dontrun{
#' diffMatrix <- DiffMat(n = 10, q = 2)
#' print(diffMatrix)
#' }
#' @export
DiffMat <- function(n, q) {
  X <- matrix(0, (n-q-1), n)
  for (i in 1:(n-q-1)) {
    for (j in 1:n) {
      if (j >= i && j <= i + q + 1) {
        X[i, j] <- (-1)^(j - i + q - 1) * choose(q + 1, j - i)
      }
    }
  }
  return(X)
}

L0kAMIASR <- function(A0=A0,y=y,max.step=50,k=k,q=q){
  n <- length(y)
  rho <- n^(q+1)
  D <- DiffMat(n,q)
  m <- nrow(D);p <- ncol(D)
  All <- 1:m
  A0 <- A0
  I0 <- setdiff(All,A0)
  u <- as.numeric(m);v <- as.numeric(m);epi <- as.numeric(m)
  for(i in 1:max.step){
    u[I0] <- solve(D[I0,]%*%t(D[I0,]))%*%D[I0,]%*%y
    u[A0] <- 0
    v[I0] <- 0
    v[A0] <- D[A0,]%*%y-D[A0,]%*%t(D[I0,])%*%u[I0]
    epi <- abs(v+u/rho)
    A <- order(epi,decreasing = T)[1:k]
    I <- setdiff(All,A)
    if(identical(A,A0) & identical(I,I0)){
      break
    }else{
      A0 <- A
      I0 <- I
    }
  }
  beta <- as.numeric(p)
  beta <- y-t(D[I,])%*%u[I]
  return(list(A=A,I=I,beta=beta))
}

#' @title L0 trend filtering using AMIAS method with R functions
#' @description This function implements an AMIAS algorithm using L0-norm minimization for sparse representations. The main objective of the function is to identify the most relevant variables in a dataset by minimizing the L0-norm. The function iteratively selects the most significant variables and computes the corresponding model parameters.
#' @param y A numeric vector representing the data from which the model is to be identified.
#' @param kmax The maximum number of variables to be included in the model. The algorithm will stop when this number is reached.
#' @param q The order of the difference to be calculated, used within the 'DiffMat' function call inside.
#' @return A list containing several components: 'kopt' - the optimal number of knots based on BIC, 'beta' - the model parameters at the optimal iteration, 'A.all' - a list containing the selected variables at each iteration, 'A' - the final set of selected variables at the optimal iteration.
#' @examples
#' \dontrun{
#' result <- l0tfAMIASR(y = rnorm(100), kmax = 10, q = 2)
#' print(result$kopt)
#' print(result$beta)
#' 
#' }
#' @export
l0tfAMIASR <- function(y=y,kmax=kmax,q=q){
  n <- length(y)
  D <- DiffMat(n,q)
  m <- nrow(D); p <- ncol(D)
  All <- 1:m
  A0 <- NULL
  I0 <- setdiff(All,A0)
  betaAll <- NULL
  mse <- as.numeric(kmax)
  bic <- as.numeric(kmax)
  A.all <- list()
  for(j in 1:kmax){
    u <- as.numeric(m)
    u[I0] <- solve(D[I0,]%*%t(D[I0,]))%*%D[I0,]%*%y
    u[A0] <- 0
    A0 <- union(A0,I0[which.max(abs(u[I0]))])
    result <- L0kAMIASR(A0=A0,y=y,k=j,q=q)
    A.all[[j]] <- sort(as.vector(result$A))
    betaAll <- cbind(betaAll,result$beta)
    A0 <- result$A
    I0 <- setdiff(All,A0)
  }
  for(i in 1:kmax){
    mse[i] <- mean((y-betaAll[,i])^2)
    bic[i] <- n*log(mse[i]) + 2*log(n)*i
  }
  kopt <- which.min(bic)
  beta <- betaAll[,kopt]
  A <- A.all[[kopt]]
  return(list(kopt=kopt,beta=beta,A.all=A.all,A=A))
}



#' @title L0 trend filtering is calculated using the AMIAS approach and R functions
#' @description This function applies L0 trend filtering to simulated data using an AMIAS approach and sparse representation. The function supports different types of simulated data generation models, such as 'Blocks' and 'Wave'. 
#' @param sigma The standard deviation of the noise added to the simulated data.
#' @param dgm The type of data generation model to use. Possible values are "Blocks" or "Wave". These models determine the characteristics of the simulated data.
#' @param n The number of data points to be generated in the simulation.
#' @param seed An optional seed for random number generation to make the simulation reproducible.
#' @return The result of applying the L0 trend filtering algorithm using AMIAS on the simulated data, typically including the optimal number of variables chosen, the coefficients of the model, and other relevant information.
#' @examples
#' \dontrun{
#' # Simulating 'Blocks' data
#' blocksData <- SimuBlocks(n = 100, sigma = 0.1, seed=123)
#' resultBlocks <- Simul0tfAMIASR(sigma = 0.1, dgm = "Blocks", n = 100, seed=123)
#' plot(blocksData$y0)
#' lines(resultBlocks$yhat,col="red")
#' 
#' # Simulating 'Wave' data
#' waveData <- SimuWave(n = 100, sigma = 0.1, seed=123)
#' resultWave <- Simul0tfAMIASR(sigma = 0.1, dgm = "Wave", n = 100, seed=123)
#' plot(waveData$y0)
#' lines(resultWave$yhat,col="red")
#' }
#' @export
Simul0tfAMIASR <- function(sigma = 0.1, dgm = "Blocks", n = n, seed = NA) {
  if (dgm=="Blocks"){
    data = SimuBlocks(n=n, sigma=sigma, seed=seed)
    q=0; kmax=20;
  }
  if (dgm=="Wave"){
    data = SimuWave(n=n, sigma=sigma, seed=seed)
    q=1; kmax=20;
  }
  resL0 = l0tfAMIASR(as.numeric(data$y), kmax=kmax, q=q)
  return(resL0)
}






