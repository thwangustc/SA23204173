## -----------------------------------------------------------------------------
  library(ggplot2)
  library(stats)
  library(boot) 
  library(bootstrap)
  library(DAAG)
  library(coda)
  library(microbenchmark)
  library(Rcpp)

## ----eval=TRUE----------------------------------------------------------------
SimuDoppler <- function(n, sigma = 0.1, seed=NA){
  if (!is.na(seed)) set.seed(seed)
  x <- seq(1/n, 1,length.out = n)
  y0 <- sqrt(x*(1-x))*sin(2*pi*(1+0.05)/(x+0.05))
  y <- y0 + sigma*rnorm(n)
  #y <- y0 + sigma*rt(n, 4)
  return(list(y = y, y0 = y0, x=x))
}

## ----eval=TRUE----------------------------------------------------------------
  n <- 1000
  Doppler <- SimuDoppler(n,sigma=0.1)
  data <- cbind(Doppler$x,Doppler$y)
  colnames(data) <- c("X","Y")
  rownames(data) <- 1:n
  data[1:10,]
  knitr::kable(data[1:10,])

## ----eval=TRUE----------------------------------------------------------------
  data[1:10,]
  knitr::kable(data[1:10,])

## -----------------------------------------------------------------------------
  plot(Doppler$x,Doppler$y,type="l",xlab="X",ylab="Y")

## -----------------------------------------------------------------------------
  my_sample <- function(x,size,prob=NULL){
    if(length(x)==1) x <- 1:x
    n <- length(x)
    if(is.null(prob)) prob <- rep(1/n,n)
    cp <- cumsum(prob)
    m <- size
    U <- runif(m)
    r <- x[findInterval(U,cp)+1]
    return(r)
  }

## ----fig.height=4-------------------------------------------------------------
  x <- 10
  size <- 1e5
  X <- my_sample(x,size)
  Y <- sample(x,size,replace=TRUE)
  counts <- table(c(rep(0,size),rep(1,size)),c(X, Y))
  barplot(counts, main="Sampling distribution generation eg1",
  xlab="X",ylab='Count', col=c("darkblue","red"),
 	legend = c('Inverse','Sample'), beside=TRUE)
  
  x <- 81:85
  size <- 1e5
  X <- my_sample(x,size)
  Y <- sample(x,size,replace=TRUE)
  counts <- table(c(rep(0,size),rep(1,size)),c(X, Y))
  barplot(counts, main="Sampling distribution generation eg2",
  xlab="X",ylab='Count', col=c("darkblue","red"),
 	legend = c('Inverse','Sample'), beside=TRUE)
  
  x <- c(89,98,100)
  size <- 1e5
  X <- my_sample(x,size)
  Y <- sample(x,size,replace=TRUE)
  counts <- table(c(rep(0,size),rep(1,size)),c(X, Y))
  barplot(counts, main="Sampling distribution generation eg3",
  xlab="X",ylab='Count', col=c("darkblue","red"),
 	legend = c('Inverse','Sample'), beside=TRUE)
  
  x <- c("man","woman")
  prob <- c(0.4,0.6)
  size <- 1e5
  X <- my_sample(x,size,prob=prob)
  Y <- sample(x,size,replace=TRUE,prob=prob)
  counts <- table(c(rep(0,size),rep(1,size)),c(X, Y))
  barplot(counts, main="Sampling distribution generation eg4",
  xlab="X",ylab='Count', col=c("darkblue","red"),
 	legend = c('Inverse','Sample'), beside=TRUE)
  
  x <- c(89,98,100,101)
  prob <- c(0.2,0.3,0.1,0.4)
  size <- 1e5
  X <- my_sample(x,size,prob=prob)
  Y <- sample(x,size,replace=TRUE,prob=prob)
  counts <- table(c(rep(0,size),rep(1,size)),c(X, Y))
  barplot(counts, main="Sampling distribution generation eg5",
  xlab="X",ylab='Count', col=c("darkblue","red"),
 	legend = c('Inverse','Sample'), beside=TRUE)

## -----------------------------------------------------------------------------
  my_Laplace <- function(n){
    x <- runif(n)
    r <- ifelse(x>=1/2,-log(2-2*x),log(2*x))
    return(r)
  }
  
  library(ggplot2)
  x <- my_Laplace(1000)
  mydata <- data.frame(x)
  ggplot(mydata, aes(x = x))+ geom_histogram(aes(y=after_stat(density)))+geom_density()

## -----------------------------------------------------------------------------
  n <- 10000
  x <- seq(-10,10,length=n)
  y0 <- function(x){
    ifelse(x>0,0.5*exp(-x),0.5*exp(x))
  }
  y1 <- my_Laplace(n)
  data <- data.frame(y1)
  ggplot(data,aes(x=y1))+ geom_histogram(aes(y=after_stat(density)),alpha = 0.4, bins = 30)+geom_density(color=2)+stat_function(fun=y0)
  
  n <- 100000
  x <- seq(-10,10,length=n)
  y1 <- my_Laplace(n)
  data <- data.frame(y1)
  ggplot(data,aes(x=y1))+ geom_histogram(aes(y=after_stat(density)),alpha = 0.4, bins = 30)+geom_density(color=2)+stat_function(fun=y0)

## -----------------------------------------------------------------------------
  my_beta32 <- function(n){
    y <- numeric(n)#生成的Beta分布记录
    j <- 0
    i <- 0#Beta分布的指标集合
    while (i < n) {
    u <- runif(1)
    j <- j + 1
    x <- runif(1) 
    if ( 27/4*x^2*(1-x) > u) {
      i <- i + 1
      y[i] <- x
      }
    }
    return(y)
  }

## -----------------------------------------------------------------------------
  library(ggplot2)
  n <- 100
  X <- my_beta32(n)
  Y <- rbeta(n,3,2)
  counts <- table(c(rep(0,n),rep(1,n)),c(X, Y))
  data <- data.frame(values = c(X,Y),
                   group = c(rep("Acceptance-rejection", 100),
                             rep("Beta", 100)))
  ggplot(data, aes(x = values, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 50)

## -----------------------------------------------------------------------------
  n <- 1000
  X <- my_beta32(n)
  Y <- rbeta(n,3,2)
  counts <- table(c(rep(0,n),rep(1,n)),c(X, Y))
  data <- data.frame(values = c(X,Y),
                   group = c(rep("Acceptance-rejection", 100),
                             rep("Beta", 100)))
  ggplot(data, aes(x = values, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 50)
  
  n <- 10000
  X <- my_beta32(n)
  Y <- rbeta(n,3,2)
  counts <- table(c(rep(0,n),rep(1,n)),c(X, Y))
  data <- data.frame(values = c(X,Y),
                   group = c(rep("Acceptance-rejection", 100),
                             rep("Beta", 100)))
  ggplot(data, aes(x = values, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 50)
  
  n <- 100000
  X <- my_beta32(n)
  Y <- rbeta(n,3,2)
  counts <- table(c(rep(0,n),rep(1,n)),c(X, Y))
  data <- data.frame(values = c(X,Y),
                   group = c(rep("Acceptance-rejection", 100),
                             rep("Beta", 100)))
  ggplot(data, aes(x = values, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 50)
  
  n <- 1000000
  X <- my_beta32(n)
  Y <- rbeta(n,3,2)
  counts <- table(c(rep(0,n),rep(1,n)),c(X, Y))
  data <- data.frame(values = c(X,Y),
                   group = c(rep("Acceptance-rejection", 100),
                             rep("Beta", 100)))
  ggplot(data, aes(x = values, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 50)

## -----------------------------------------------------------------------------
  my_knernel <- function(n){
    i <- 0
    y <- as.numeric(n)
    while(i < n){
      i <- i+1
      U1 <- runif(1,min=-1,max=1)
      U2 <- runif(1,min=-1,max=1)
      U3 <- runif(1,min=-1,max=1)
      if(abs(U3) >= max(abs(U1),abs(U2))){
      y[i] <- U2
    }else{
      y[i] <- U3
      }
    }
    return(y)
  }

## -----------------------------------------------------------------------------
  n <- 10000
  x <- seq(-1,1,length=n)
  y0 <- function(x){
    3/4*(1-x^2)
  }
  y1 <- my_knernel(n)
  data <- data.frame(y1)
  ggplot(data,aes(x=y1))+ geom_density()+stat_function(fun=y0,color=2)+geom_histogram(aes(y=after_stat(density)),alpha = 0.4, bins = 50)
  
  n <- 100000
  x <- seq(-1,1,length=n)
  y1 <- my_knernel(n)
  data <- data.frame(y1)
  ggplot(data,aes(x=y1))+ geom_density()+stat_function(fun=y0,color=2)+geom_histogram(aes(y=after_stat(density)),alpha = 0.4, bins = 50)
  
  n <- 1000000
  x <- seq(-1,1,length=n)
  y1 <- my_knernel(n)
  data <- data.frame(y1)
  ggplot(data,aes(x=y1))+ geom_density()+stat_function(fun=y0,color=2)+geom_histogram(aes(y=after_stat(density)),alpha = 0.4, bins = 50)

## -----------------------------------------------------------------------------
  set.seed(100)
  K <- 100
  n <- 1000000
  result <- as.numeric(K)
  rho <- c(0.5,0.8,1) 
  rhomin <- 0
  pi1 <- as.numeric(K)
  pi2 <- as.numeric(K)
  pi3 <- as.numeric(K)
  for(i in 1:K){
    sample1 <- rbinom(1, n, prob=2*rho[1]/pi)#rho=0.5
    sample2 <- rbinom(1, n, prob=2*rho[2]/pi)#rho=0.8
    sample3 <- rbinom(1, n, prob=2*rho[3]/pi)#rho=1
    pi1[i] <- 2*rho[1]*n/sample1
    pi2[i] <- 2*rho[2]*n/sample2
    pi3[i] <- 2*rho[3]*n/sample3
  }
  var(pi1)
  var(pi2)
  var(pi3)

## -----------------------------------------------------------------------------
  set.seed(100)
  B <- 1000000
  x <- runif(B)
  V0 <- var(exp(x))/B
  V0

## -----------------------------------------------------------------------------
  set.seed(100)
  B <- 1000000
  x <- runif(B/2)
  V <- var(exp(x)+exp(1-x))/(2*B)
  V

## -----------------------------------------------------------------------------
  V0/V

## -----------------------------------------------------------------------------
  M <- 10000; k <- 1000 # what if k is larger?
  r <- M/k #replicates per stratum
  N <- 50 #number of times to repeat the estimation
  T2 <- numeric(k)
  est <- matrix(0, N, 2)
  g<-function(x)exp(-x)/(1+x*x)*(x>0)*(x<1)
  for (i in 1:N) {
  est[i, 1] <- mean(g(runif(M)))
  for(j in 1:k)T2[j]<-mean(g(runif(M/k,(j-1)/k,j/k)))
  est[i, 2] <- mean(T2)
  }
  round(apply(est,2,mean),4)
  round(apply(est,2,sd),5)

## -----------------------------------------------------------------------------
  f_1 <- function(x){
    return(4/pi/(1+x^2))
  }

  f_2 <- function(x){
    return(exp(1/2)*exp(-x/2)/2)
  }
  
  g <- function(x){
    return(x^2/sqrt(2*pi)*exp(-x^2/2))
  }
  
  g_1 <- function(x){
    return((g(x))^2/f_1(x))
  }
  
  g_2 <- function(x){
    return((g(x))^2/f_2(x))
  }
  
  integrate(g_1,1,Inf)$value
  integrate(g_2,1,Inf)$value

## -----------------------------------------------------------------------------
  B <- 1e6
  U <- runif(B)
  X <- 1-2*log(1-U)
  mean(g(X)/f_2(X))
  
  integrate(g,1,Inf)$value

## -----------------------------------------------------------------------------
  a <- c(0,0.2,0.4,0.6,0.8,1)
  hattheta <- as.numeric(5)
  vartheta <- as.numeric(5)
  M <- 1e4
  k <- 5
  m <- M/k
  X <- as.numeric(m)
  
  f <- function(x){
    return(exp(-x)/(1-exp(-1)))
  }
  
  g <- function(x){
    return(exp(-x)/(1+x^2))
  }
  
  for(i in 1:5){
    U <- runif(m)
    f0 <- function(x){
      fi0 <- integrate(f,a[i],a[i+1])$value
      return(f(x)/fi0)
    }
    f0x <- function(x){
      return(integrate(f0,a[i],x)$value)
    }
    fi <- function(u,x){
      return(f0x(x)-u)
    }
    for(j in 1:m){
      X[j] <- uniroot(fi,c(a[i],a[i+1]),u=U[j])$root
      X[j] <- g(X[j])/f0(X[j])
    }
    hattheta[i] <- mean(X)
    vartheta[i] <- var(X)
  }
  sum(hattheta)
  sqrt(sum(vartheta))

## -----------------------------------------------------------------------------
  set.seed(123)
  B <- 1e6
  n <- 20
  alpha <- 0.05
  s <- 0
  
  for(i in 1:B){
    x <- rchisq(n,2)
    UCL <- (n-1) * var(x) / qchisq(alpha, df=n-1)
    if(UCL>4){
      s <- s+1
    }
  }
  s/B

## -----------------------------------------------------------------------------
  B <- 1e6
  n <- 20
  alpha <- 0.05
  s <- 0
  
  for(i in 1:B){
    x <- rnorm(n,0,2)
    UCL <- mean(x)+qt(0.95,n-1)*sd(x)/sqrt(n)
    if(UCL>0){
      s <- s+1
    }
  }
  s/B

## -----------------------------------------------------------------------------
  set.seed(123)
  B <- 1e6
  n <- 20
  alpha <- 0.05
  s <- 0
  
  for(i in 1:B){
    x <- rchisq(n,2)
    UCL <- mean(x)+qt(0.95,n-1)*sd(x)/sqrt(n)
    if(UCL>2){
      s <- s+1
    }
  }
  s/B

## -----------------------------------------------------------------------------
  ##卡方分布
  set.seed(123)
  B <- 1e6
  n <- 20
  alpha <- 0.05
  s <- 0
  
  for(i in 1:B){
    x <- rchisq(n,1)
    UCL <- mean(x)+qt(0.975,n-1)*sd(x)/sqrt(n)
    UCR <- mean(x)-qt(0.975,n-1)*sd(x)/sqrt(n)
    if(UCL<1 | UCR>1){
      s <- s+1
    }
  }
  s/B

## -----------------------------------------------------------------------------
  ##均匀分布
  set.seed(123)
  B <- 1e6
  n <- 20
  alpha <- 0.05
  s <- 0
  
  for(i in 1:B){
    x <- runif(n,0,2)
    UCL <- mean(x)+qt(0.975,n-1)*sd(x)/sqrt(n)
    UCR <- mean(x)-qt(0.975,n-1)*sd(x)/sqrt(n)
    if(UCL<1 | UCR>1){
      s <- s+1
    }
  }
  s/B

## -----------------------------------------------------------------------------
  ##指数分布
  set.seed(123)
  B <- 1e6
  n <- 20
  alpha <- 0.05
  s <- 0
  
  for(i in 1:B){
    x <- rexp(n,1)
    UCL <- mean(x)+qt(0.975,n-1)*sd(x)/sqrt(n)
    UCR <- mean(x)-qt(0.975,n-1)*sd(x)/sqrt(n)
    if(UCL<1 | UCR>1){
      s <- s+1
    }
  }
  s/B

## -----------------------------------------------------------------------------
  m <- 1000
  m1 <- 950
  m2 <- 50
  p <- as.numeric(m)
  set.seed(123)
  
  p[1:950] <- runif(m1)
  p[951:1000] <- rbeta(m2,0.1,1)
  pbon <- p.adjust(p,method="bonferroni")
  pbh <- p.adjust(p,method="BH")
  sum(pbon<0.1)
  sum(pbh<0.1)

## -----------------------------------------------------------------------------
  m <- 1000
  m1 <- 950
  m2 <- 50
  M <- 1000
  
  p <- as.numeric(m)
  FWER <- as.numeric(M)
  FWER0 <- as.numeric(M)
  FDR <- as.numeric(M)
  FDR0 <- as.numeric(M)
  TPR <- as.numeric(M)
  TPR0 <- as.numeric(M)
  
  for(i in 1:M){
    set.seed(i)
    p[1:950] <- runif(m1)
    p[951:1000] <- rbeta(m2,0.1,1)
    pbon <- p.adjust(p,method="bonferroni")
    pbh <- p.adjust(p,method="BH")
    
    FWER[i] <- ifelse(sum(pbon[1:950]<0.1)==0,0,1)
    FWER0[i] <- ifelse(sum(pbh[1:950]<0.1)==0,0,1)
    
    FDR[i] <- length(intersect(which(pbon<0.1),1:950))/sum(pbon<0.1)
    FDR0[i] <- length(intersect(which(pbh<0.1),1:950))/sum(pbh<0.1)
    
    TPR[i] <- length(intersect(which(pbon<0.1),951:1000))/50
    TPR0[i] <- length(intersect(which(pbh<0.1),951:1000))/50
  }
  mean(FWER);mean(FWER0)
  mean(FDR);mean(FDR0)
  mean(TPR);mean(TPR0)

## -----------------------------------------------------------------------------
  bootstrapexp <- function(n,lambda,B,M){
    set.seed(1)
    x <- rexp(n,lambda)
    lambdaresult <- as.numeric(M)
    seresult <- as.numeric(M)
    
    for(i in 1:M){
      lambdahat <- as.numeric(B)
      sehat <- as.numeric(B)
      for(b in 1:B){
        xstar <- sample(x,n,replace=TRUE)
        lambdahat[b] <- 1/mean(xstar)
      }
      seresult[i] <- sd(lambdahat)
      lambdaresult[i] <- mean(lambdahat)-1/mean(x)
    }
    bias <- mean(lambdaresult)
    se <- mean(seresult)
    return(list(bias=bias,se=se))
  }
  
  lambda <- 2
  B <- 1000
  M <- 1000
  n1 <- 5
  n2 <- 10
  n3 <- 20
  
  bootstrapexp(n1,lambda,B,M)
  bootstrapexp(n2,lambda,B,M)
  bootstrapexp(n3,lambda,B,M)

## -----------------------------------------------------------------------------
  x <- c(576,635,558,578,666,580,555,661,651,605,653,575,545,572,594)
  y <- c(339,330,281,303,344,307,300,343,336,313,312,274,276,288,296)
  rhat <- cor(x,y)
  
  bootse <- function(LSAT,GPA,R){
    m <- length(LSAT)
    th <- as.numeric(R)
    for(i in 1:R){
      A <- sample(1:m, size = m, replace = TRUE)
      #B <- sample(1:m, size = m, replace = TRUE)
      th[i] <- cor(LSAT[A],GPA[A])
      #th[i] <- cor(LSAT[A],GPA[B])
    }
    return(sd(th))
  }
  
  set.seed(100)
  B <- 1000
  n <- length(x)
  R <- t <- as.numeric(B)
  Rsd <- as.numeric(B)
  LSAT <- GPA <- NULL
  for (b in 1:B) {
    i <- sample(1:n, size = n, replace = TRUE)
    LSAT <- x[i]
    GPA <- y[i]
    R[b] <- cor(LSAT, GPA)
    Rsd[b] <- bootse(LSAT,GPA,R = 100)
  }
  t <- (R-rhat)/Rsd
  rhat-quantile(t,0.025)*sd(R)
  rhat-quantile(t,0.975)*sd(R)
  

## ----warning=FALSE------------------------------------------------------------
  data <- c(3,5,7,18,43,85,91,98,100,130,230,487)
  mean(data)
  lambda.boot <- function(data, ind) {
    mean(data[ind])
  }
  boot.obj <- boot(data, statistic = lambda.boot, R = 1000)
  boot.ci(boot.obj, type = c("basic", "norm", "perc","bca"))

## -----------------------------------------------------------------------------
  data(scor, package = "bootstrap")
  n <- nrow(scor)
  thetaeigenvalue <- function(dat,ind){
    data <- dat[ind,]
    sigma <- eigen(cov(data))
    return(max(sigma$val)/sum(sigma$val))
  }
  theta.hat <- thetaeigenvalue(scor,1:n)
  theta.jack <- numeric(n)
  for(i in 1:n){
    theta.jack[i] <- thetaeigenvalue(scor,setdiff(1:n,i))
  }
  bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
  se.jack <- sqrt((n-1)*mean((theta.jack - mean(theta.jack))^2))
  round(c(original=theta.hat,bias.jack=bias.jack,se.jack=se.jack),4)

## ----warning=FALSE------------------------------------------------------------
  attach(ironslag)
  n <- length(magnetic) 
  e1 <- e2 <- e3 <- e4 <- numeric(n)
  c <- NULL
  k <- 1
  for(i in 1:n){
    for(j in 1:n){
      if(j <= i){
        c[[k]] <- c(i,j)
        k <- k+1
      }
    }
  }
  m <- choose(n,2)
  for (k in 1:m) {
    y <- magnetic[-c[[k]]]
    x <- chemical[-c[[k]]]
    J1 <- lm(y ~ x)
    yhat1_1 <- J1$coef[1] + J1$coef[2] * chemical[c[[k]][1]]
    yhat1_2 <- J1$coef[1] + J1$coef[2] * chemical[c[[k]][2]]
    e1[k] <- (magnetic[c[[k]][1]]-yhat1_1)^2 + (magnetic[c[[k]][2]]-yhat1_2)^2
    
    J2 <- lm(y ~ x + I(x^2))
    yhat2_1 <- J2$coef[1] + J2$coef[2] * chemical[c[[k]][1]] +
    J2$coef[3] * chemical[c[[k]][1]]^2
    yhat2_2 <- J2$coef[1] + J2$coef[2] * chemical[c[[k]][2]] +
    J2$coef[3] * chemical[c[[k]][2]]^2
    e2[k] <- (magnetic[c[[k]][1]] - yhat2_1)^2 + (magnetic[c[[k]][2]] - yhat2_2)^2
    
    J3 <- lm(log(y) ~ x)
    logyhat3_1 <- J3$coef[1] + J3$coef[2] * chemical[c[[k]][1]]
    yhat3_1 <- exp(logyhat3_1)
    logyhat3_2 <- J3$coef[1] + J3$coef[2] * chemical[c[[k]][2]]
    yhat3_2 <- exp(logyhat3_2)
    e3[k] <- (magnetic[c[[k]][1]] - yhat3_1)^2+(magnetic[c[[k]][2]] - yhat3_2)^2
    
    J4 <- lm(log(y) ~ log(x))
    logyhat4_1 <- J4$coef[1] + J4$coef[2] * log(chemical[c[[k]][1]])
    yhat4_1 <- exp(logyhat4_1)
    logyhat4_2 <- J4$coef[1] + J4$coef[2] * log(chemical[c[[k]][2]])
    yhat4_2 <- exp(logyhat4_2)
    e4[k] <- (magnetic[c[[k]][1]] - yhat4_1)^2+(magnetic[c[[k]][2]] - yhat4_2)^2
  }
  
  c(mean(e1), mean(e2), mean(e3), mean(e4))

## -----------------------------------------------------------------------------
  set.seed(123)
  attach(chickwts)
  x <- sort(as.vector(weight[feed == "soybean"]))
  y <- sort(as.vector(weight[feed == "linseed"]))
  detach(chickwts)
  n <- length(x)
  m <- length(y)
  z <- c(x, y) 
  k <- n+m
  R <- 1000
  D <- numeric(R) 
  options(warn = -1)
  
  ecdf <- function(vec,x){
    n <- length(x)
    R <- as.numeric(n)
    for(i in 1:n){
      R[i] <- mean(vec<=x[i])
    }
    return(R)
  }
  
  W <- function(vec1,vec2,x,y){
    f1 <- ecdf(vec1,x)
    f2 <- ecdf(vec1,y)
    g1 <- ecdf(vec2,x)
    g2 <- ecdf(vec2,y)
    return(sum(f1-g1)^2+sum(f2-g2)^2)
  }
  
  D0 <- m*n/(m+n)^2*W(x,y,x,y)
  for (i in 1:R) {
    K <- sample(k, size = n, replace = FALSE)
    x1 <- z[K]
    y1 <- z[-K] 
    D[i] <- m*n/(m+n)^2*W(x1,y1,x1,y1)
  }
  p <- mean(c(D0, D) >= D0)
  options(warn = 0)
  p

## -----------------------------------------------------------------------------
  set.seed(123)
  attach(chickwts)
  x <- sort(as.vector(weight[feed == "sunflower"]))
  y <- sort(as.vector(weight[feed == "linseed"]))
  detach(chickwts)
  n <- length(x)
  m <- length(y)
  z <- c(x, y) 
  k <- n+m
  R <- 1000
  D <- numeric(R) 
  options(warn = -1)
  
  ecdf <- function(vec,x){
    n <- length(x)
    R <- as.numeric(n)
    for(i in 1:n){
      R[i] <- mean(vec<=x[i])
    }
    return(R)
  }
  
  W <- function(vec1,vec2,x,y){
    f1 <- ecdf(vec1,x)
    f2 <- ecdf(vec1,y)
    g1 <- ecdf(vec2,x)
    g2 <- ecdf(vec2,y)
    return(sum(f1-g1)^2+sum(f2-g2)^2)
  }
  
  D0 <- m*n/(m+n)^2*W(x,y,x,y)
  for (i in 1:R) {
    K <- sample(k, size = n, replace = FALSE)
    x1 <- z[K]
    y1 <- z[-K] 
    D[i] <- m*n/(m+n)^2*W(x1,y1,x1,y1)
  }
  p <- mean(c(D0, D) >= D0)
  options(warn = 0)
  p

## -----------------------------------------------------------------------------
  set.seed(123)
  count5 <- function(x, y) {
    X <- x - mean(x)
    Y <- y - mean(y)
    outx <- sum(X > max(Y)) + sum(X < min(Y))
    outy <- sum(Y > max(X)) + sum(Y < min(X))
    return((max(c(outx, outy))))
  }

  n1 <- 20
  n2 <- 30
  n <- n1+n2
  mu1 <- mu2 <- 0
  sigma1 <- sigma2 <- 1
  m <- 1000
  B <- 100
  p <- as.numeric(m)
  
  for(i in 1:m){
    x <- rnorm(n1, mu1, sigma1)
    y <- rnorm(n2, mu2, sigma2)
    z <- c(x,y)
    x <- x - mean(x) 
    y <- y - mean(y)
    D0 <- count5(x,y)
    D <- as.numeric(m)
  
    for(j in 1:B){
      K <- sample(n, size = n1, replace = FALSE)
      x1 <- z[K]
      y1 <- z[-K] 
      D[j] <- count5(x1,y1)
    }
    p[i] <- mean(c(D0, D) >= D0)
    options(warn = 0)
  }
  mean(p<=0.0625)

## -----------------------------------------------------------------------------
  set.seed(123)
  n1 <- 20
  n2 <- 50
  n <- n1+n2
  mu1 <- mu2 <- 0
  sigma1 <- sigma2 <- 1
  m <- 1000
  B <- 100
  p <- as.numeric(m)
  
  for(i in 1:m){
    x <- rnorm(n1, mu1, sigma1)
    y <- rnorm(n2, mu2, sigma2)
    z <- c(x,y)
    x <- x - mean(x) 
    y <- y - mean(y)
    D0 <- count5(x,y)
    D <- as.numeric(m)
  
    for(j in 1:B){
      K <- sample(n, size = n1, replace = FALSE)
      x1 <- z[K]
      y1 <- z[-K] 
      D[j] <- count5(x1,y1)
    }
    p[i] <- mean(c(D0, D) >= D0)
    options(warn = 0)
  }
  mean(p<=0.0625)

## -----------------------------------------------------------------------------
  logisticroot <- function(N,b1,b2,b3,f0){
    X1 <- rpois(N,1)
    X2 <- rexp(N,1)
    X3 <- rbinom(N,1,0.5)
    z <- function(a,X1,X2,X3,b1,b2,b3,f0){
      p <- as.numeric(N)
      p <- 1/(1+exp(-a-b1*X1-b2*X2-b3*X3))
      return(mean(p)-f0)
    }
    return(uniroot(z,c(-1000,1000),X1=X1,X2=X2,X3=X3,b1=b1,b2=b2,b3=b3,f0=f0)$root)
  }
  ##示例
  logisticroot(10000,1,2,3,0.5)

## -----------------------------------------------------------------------------
  N <- 1e4;b1 <- 0;b2 <- 1;b3 <- -1
  set.seed(123)
  logisticroot(N,b1,b2,b3,0.1)
  logisticroot(N,b1,b2,b3,0.01)
  logisticroot(N,b1,b2,b3,0.001)
  logisticroot(N,b1,b2,b3,0.0001)

## -----------------------------------------------------------------------------
  N <- 1e4;b1 <- 0;b2 <- 1;b3 <- -1
  B <- 100
  set.seed(123)
  m <- seq(-10,-0.0001,length=100)
  y <- as.numeric(B)
  for(i in 1:B){
    y[i] <- logisticroot(N,b1,b2,b3,exp(m[i]))
  }
  y

## -----------------------------------------------------------------------------
  dat <- data.frame(X = m, Y = y)
  ggplot(data = dat, mapping = aes(x = X, y = Y)) + geom_line() + ylab("a")+xlab("log(f0)")

## -----------------------------------------------------------------------------
  c <- 0.001
  x <- seq(-10,-0.0001,length=1000)
  y <- log(c/(exp(-x)-1))
  plot(x,y,type="l")

## -----------------------------------------------------------------------------
  B <- 1e3
  sigma <- 0.5
  X <- as.numeric(B)
  X[1] <- 0
  k <- 0 ##k代表拒绝的次数
  Z <- rnorm(B,0,sigma)
  U <- runif(B)
  for(i in 2:B){
    Y <- X[i-1]+Z[i]
    if(U[i] <= exp(-abs(Y)+abs(X[i-1]))){
      X[i] <- Y
    }else{
      X[i] <- X[i-1]
      k <- k+1
    }
  }
  par(fig = c(0.1, 0.95, 0.1, 0.9))
  plot(X,xlab=bquote(sigma==0.5),type="l",ylab="X")
  print(k/B)

## -----------------------------------------------------------------------------
  B <- 1e3
  sigma <- 1
  X <- as.numeric(B)
  X[1] <- 0
  k <- 0 ##k代表拒绝的次数
  Z <- rnorm(B,0,sigma)
  U <- runif(B)
  for(i in 2:B){
    Y <- X[i-1]+Z[i]
    if(U[i] <= exp(-abs(Y)+abs(X[i-1]))){
      X[i] <- Y
    }else{
      X[i] <- X[i-1]
      k <- k+1
    }
  }
  par(fig = c(0.1, 0.95, 0.1, 0.9))
  plot(X,xlab=bquote(sigma==1),type="l",ylab="X")
  print(k/B)

## -----------------------------------------------------------------------------
  B <- 1e3
  sigma <- 2.5
  X <- as.numeric(B)
  X[1] <- 0
  k <- 0 ##k代表拒绝的次数
  Z <- rnorm(B,0,sigma)
  U <- runif(B)
  for(i in 2:B){
    Y <- X[i-1]+Z[i]
    if(U[i] <= exp(-abs(Y)+abs(X[i-1]))){
      X[i] <- Y
    }else{
      X[i] <- X[i-1]
      k <- k+1
    }
  }
  par(fig = c(0.1, 0.95, 0.1, 0.9))
  plot(X,xlab=bquote(sigma==2.5),type="l",ylab="X")
  print(k/B)

## -----------------------------------------------------------------------------
  B <- 1e3
  sigma <- 16
  X <- as.numeric(B)
  X[1] <- 0
  k <- 0 ##k代表拒绝的次数
  Z <- rnorm(B,0,sigma)
  U <- runif(B)
  for(i in 2:B){
    Y <- X[i-1]+Z[i]
    if(U[i] <= exp(-abs(Y)+abs(X[i-1]))){
      X[i] <- Y
    }else{
      X[i] <- X[i-1]
      k <- k+1
    }
  }
  par(fig = c(0.1, 0.95, 0.1, 0.9))
  plot(X,xlab=bquote(sigma==16),type="l",ylab="X")
  print(k/B)

## -----------------------------------------------------------------------------
  N <- 5000 
  burn <- 1000 #burn-in length
  X <- matrix(0, N, 2) 
  rho <- 0.9 
  mu1 <- 0
  mu2 <- 0
  sigma1 <- 1
  sigma2 <- 1
  s1 <- sqrt(1-rho^2)*sigma1
  s2 <- sqrt(1-rho^2)*sigma2
  X[1, ] <- c(mu1, mu2) 
  for (i in 2:N) {
    x2 <- X[i-1, 2]
    m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
    X[i, 1] <- rnorm(1, m1, s1)
    x1 <- X[i, 1]
    m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
    X[i, 2] <- rnorm(1, m2, s2)
  }
  b <- burn + 1
  x <- X[b:N, ]

## -----------------------------------------------------------------------------
  vecx <- x[,1]
  vecy <- x[,2]
  fit <- lm(vecy~vecx)
  res <- fit$residuals
  qqnorm(res)
  qqline(res)

## -----------------------------------------------------------------------------
  value <- c(vecy,vecx)
  group <- c(rep(1,4000),rep(0,4000))
  bartlett.test(value~group)

## ----warning=FALSE------------------------------------------------------------
  set.seed(100)
  ## genchain为生成链的函数
  genchain <- function(sigma, N, X1) {
    f <- function(x, sigma) {
      if (any(x < 0)) return (0)
      stopifnot(sigma > 0)
      return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
    }
    x <- numeric(N)
    x[1] <- X1
    u <- runif(N)
    for (i in 2:N) {
      xt <- x[i-1]
      y <- rchisq(1, df = xt)
      num <- f(y, sigma) * dchisq(xt, df = y)
      den <- f(xt, sigma) * dchisq(y, df = xt)
      if (u[i] <= num/den){
        x[i] <- y 
      }else {
        x[i] <- xt
      }
    }
    return(x)
  }

  sigma <- 4 
  k <- 5 #number of chains to generate
  n <- 1000
  b <- 100 #burn-in length
  x0 <- rchisq(5,df=1)
  X <- matrix(0, nrow=k, ncol=n)
  for (i in 1:k){
    X[i, ] <- genchain(sigma, n, x0[i]) 
  }
  mc1 <- as.mcmc(X[1,])
  mc2 <- as.mcmc(X[2,])
  mc3 <- as.mcmc(X[3,])
  mc4 <- as.mcmc(X[4,])
  mc5 <- as.mcmc(X[5,])
  gelman.diag(mcmc.list(mc1,mc2,mc3,mc4,mc5))

## -----------------------------------------------------------------------------
  Gelman.Rubin <- function(psi) {
    psi <- as.matrix(psi)
    n <- ncol(psi)
    k <- nrow(psi)
    psi.means <- rowMeans(psi) 
    B <- n * var(psi.means) 
    psi.w <- apply(psi, 1, "var") 
    W <- mean(psi.w) 
    v.hat <- W*(n-1)/n + (B/n) 
    r.hat <- v.hat / W 
    return(r.hat)
  }

  psi <- t(apply(X, 1, cumsum))
  for (i in 1:nrow(psi)){
    psi[i,] <- psi[i,] / (1:ncol(psi))
  }
  print(Gelman.Rubin(psi))

## -----------------------------------------------------------------------------
  direct_MLE <- function(u, v) {
    f <- function(lambda, u, v) {
        n <- length(u)
        r <- as.numeric(n)
        for(i in 1:n){
          r[i] <- (-u[i]*exp(-lambda*u[i])+v[i]*exp(-lambda*v[i]))/(exp(-lambda*u[i])-exp(-lambda*v[i]))
        }
        return(sum(r))
    }
    opt <- uniroot(f,c(0,5),u=u,v=v)$root
    return(opt)
  }

  u <- c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
  v <- c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3)
  direct_MLE(u,v)

## -----------------------------------------------------------------------------
  EM_MLE <- function(u, v, tol = 1e-6, max_iter = 1000) {
    n <- length(u)
    lambda <- 1  
    for (i in 1:max_iter) {
        E_X <- (exp(-lambda * u) * u - exp(-lambda * v) * v + 1/lambda * (exp(-lambda * u) - exp(-lambda * v))) / (exp(-lambda * u) - exp(-lambda * v))
        lambda_new <- n / sum(E_X)
        if (abs(lambda_new - lambda) < tol) {
            break
        }
        lambda <- lambda_new
    }
    return(lambda)
  }
  EM_MLE(u, v)

## -----------------------------------------------------------------------------
  solve.game <- function(A) {
    #solve the two player zero-sum game by simplex method
    #optimize for player 1, then player 2
    #maximize v subject to ...
    #let x strategies 1:m, and put v as extra variable
    #A1, the <= constraints
    min.A <- min(A)
    A <- A - min.A #so that v >= 0
    max.A <- max(A)
    A <- A / max(A)
    m <- nrow(A)
    n <- ncol(A)
    it <- n^3
    a <- c(rep(0, m), 1) #objective function
    A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
    b1 <- rep(0, n)
    A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
    b3 <- 1
    sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
    maxi=TRUE, n.iter=it)
    #the ’solution’ is [x1,x2,...,xm | value of game]
    #minimize v subject to ...
    #let y strategies 1:n, with v as extra variable
    a <- c(rep(0, n), 1) #objective function
    A1 <- cbind(A, rep(-1, m)) #constraints <=
    b1 <- rep(0, m)
    A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
    b3 <- 1
    sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
    maxi=FALSE, n.iter=it)
    soln <- list("A" = A * max.A + min.A,
                 "x" = sx$soln[1:m],
                 "y" = sy$soln[1:n],
                 "v" = sx$soln[m+1] * max.A + min.A)
    soln
  }

## -----------------------------------------------------------------------------
  library(boot)
  A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
        2,0,0,0,-3,-3,4,0,0,
        2,0,0,3,0,0,0,-4,-4,
        -3,0,-3,0,4,0,0,5,0,
        0,3,0,-4,0,-4,0,5,0,
        0,3,0,0,4,0,-5,0,-5,
        -4,-4,0,0,0,5,0,0,6,
        0,0,4,-5,-5,0,0,0,6,
        0,0,4,0,0,5,-6,-6,0), 9, 9)
  m <- solve.game(A)
  print(round(cbind(m$x, m$y), 6))
  t(m$x)%*%A%*%(m$y)

## -----------------------------------------------------------------------------
  B <- A+2
  m <- solve.game(B)
  print(round(cbind(m$x, m$y), 6))
  t(m$x)%*%B%*%(m$y)
  
  C <- 2*A
  n <- solve.game(C)
  print(round(cbind(n$x, n$y), 6))
  t(n$x)%*%C%*%(n$y)

## -----------------------------------------------------------------------------
  ## 举例
  A <- NULL
  A[[1]] <- c("stat","2023")
  A[[2]] <- 24
  A[[3]] <- list(3,4)
  
  print(A)
  unlist(A)
  as.vector(A)

## -----------------------------------------------------------------------------
  y <- c(1,5,6,9,12,9)
  dim(y)

## -----------------------------------------------------------------------------
  X <- matrix(1:6,nrow=3)
  is.matrix(X)
  is.array(X)

## -----------------------------------------------------------------------------
  x <- c("add","baa","caa")
  y <- 1:3
  data <- data.frame(x,y)
  data
  as.matrix(data)

## -----------------------------------------------------------------------------
  # 创建一个有0行但有2列的数据框
  dat1 <- data.frame(列1 = numeric(0), 列2 = character(0))
  dim(dat1)
  
  # 创建一个有0列的数据框
  dat2 <- data.frame()
  dim(dat2)
  # 这将创建一个既没有行也没有列的空数据框。

## -----------------------------------------------------------------------------
  scale01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    (x - rng[1]) / (rng[2] - rng[1])
  }

## -----------------------------------------------------------------------------
  ## 假设data是我们的数据框，我们将scale01函数应用于数据框的每一列
  data <- data.frame(X=c(1,2,5,10),Y=c(3,5,6,9))
  sapply(data, scale01)
  
  ## 假设data1是我们的数据框，我们仅将scale01函数应用于数据框的数值列
  data2 <- data.frame(X=c(1,2,5,10),Y=c(3,5,6,9),Z=c("a","l","n","aa"))
  numcols <- sapply(data2, is.numeric)
  sapply(data2[,numcols], scale01)

## -----------------------------------------------------------------------------
  ## 计算数值型数据框的每列的标准差
  data <- data.frame(X=c(1,2,5,10),Y=c(3,5,6,9))
  vapply(data, sd, numeric(1))
  
  ## 计算混合型数据框的每个数值列的标准差
  data2 <- data.frame(X=c(1,2,5,10),Y=c(3,5,6,9),Z=c("a","l","n","aa"))
  numcols <- vapply(data2, is.numeric, logical(1))
  vapply(data2[numcols], sd, numeric(1))

## -----------------------------------------------------------------------------
  gibbssampler <- function(n, a, b, N) {
    x <- numeric(N)
    y <- numeric(N)
    y[1] <- runif(1)  ## 初始值
    ## 交替迭代
    for (i in 2:N) {
      x[i] <- rbinom(1, n, y[i-1])
      y[i] <- rbeta(1, x[i]+a, n-x[i]+b)
    }
    return(data.frame(x, y))
  }

## ----warning=FALSE------------------------------------------------------------
  cppFunction('
  DataFrame gibbssampler_Rcpp(int n, double a, double b, int N) {
    NumericVector x(N);
    NumericVector y(N);
    y[0] = R::runif(0, 1);  // 初始值
    // 交替迭代
    for (int i = 1; i < N; ++i) {
      x[i] = R::rbinom(n, y[i-1]);
      y[i] = R::rbeta(x[i]+a, n-x[i]+b);
    }
    return DataFrame::create(Named("x")=x, Named("y")=y);
  }
  ')

## -----------------------------------------------------------------------------
  ## 使用microbenchmark包比较这两个函数的性能。
  n <- 20
  a <- 2
  b <- 1
  N <- 1000

  benchmark <- microbenchmark(
    R = gibbssampler(n, a, b, N),
    Rcpp = gibbssampler_Rcpp(n, a, b, N),
    times = 30
  )
  print(benchmark)

