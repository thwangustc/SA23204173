## -----------------------------------------------------------------------------
  library(SA23204173)
  library(microbenchmark)

## -----------------------------------------------------------------------------
  data1 <- SimuBlocks(200,sigma=0.1,seed=123)
  plot(data1$x,data1$y,xlab="X",ylab="Y")
  data2 <- SimuWave(200,sigma=0.1,seed=123)
  plot(data2$x,data2$y,col="red",xlab="X",ylab="Y")

## -----------------------------------------------------------------------------
  AB1 <- l0tfABESSR(y=data1$y,kmax=20,q=0)
  AM1 <- l0tfAMIASR(y=data1$y,kmax=20,q=0)
  plot(data1$x,data1$y,xlab="X",ylab="Y")
  lines(data1$x,AB1$yhat,col="red")
  lines(data1$x,AM1$beta,col="blue")
  
  AB2 <- l0tfABESSR(y=data2$y,kmax=20,q=1)
  AM2 <- l0tfAMIASR(y=data2$y,kmax=20,q=1)
  plot(data2$x,data2$y,xlab="X",ylab="Y")
  lines(data2$x,AB2$yhat,col="red")
  lines(data2$x,AM2$beta,col="blue")

## -----------------------------------------------------------------------------
 set.seed(123)
 tm1 <- microbenchmark::microbenchmark(
   ABR = Simul0tfABESSR(sigma=0.1,dgm="Blocks",n=60),
   ABC = Simul0tfABESSC(sigma=0.1,dgm="Blocks",n=60),
   AMR = Simul0tfAMIASR(sigma=0.1,dgm="Blocks",n=60),
   times = 30
  )
 print(summary(tm1)[,c(1,3,5,6)])
 tm2 <- microbenchmark::microbenchmark(
   ABR = Simul0tfABESSR(sigma=0.1,dgm="Wave",n=60),
   ABC = Simul0tfABESSC(sigma=0.1,dgm="Wave",n=60),
   AMR = Simul0tfAMIASR(sigma=0.1,dgm="Wave",n=60),
   times = 30
  )
 print(summary(tm2)[,c(1,3,5,6)])

