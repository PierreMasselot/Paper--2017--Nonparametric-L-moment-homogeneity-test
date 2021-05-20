######################################################################################
#
#            Nonparametric L-moment Homogeneity test
#
#     Authors : Pierre Masselot, Fateh Chebana, Taha B.M.J. Ouarda
#
#                     Date : July 2016
#
#	         Dependencies : packages 'lmomco', 'lmom'
#
######################################################################################

Lmom.homog.test <- function(X, Nsim = 500, type = "B", typenorm = "2", center = type %in% c("B","Yr","Ys"), printInfo = T)
{
  if (!type %in% c("HW","B","M","Yr","Ys")) stop("Type does not exist. Use one of: 'HW','B','M','Y'")
  N <- length(X) # number of sites
  # Verification of the coherence between sites
  Pvec <- vector("numeric",N)
  X <- lapply(X,as.matrix)
  P <- sapply(X,ncol)
  if (length(unique(P)) > 1) stop("All sites have not the same number of variables.") else P <- unique(P)
  if (printInfo){
     nomtest <- switch(type,HW = "Hosking-Wallis",M = "Permutations",B = "Bootstrap",Ys = "Polya",Yr = "Polya")
     print(sprintf("Performing the %s test with %s variables.",nomtest,P))
  }
  nsite <- sapply(X,function(x) nrow(as.matrix(x))) # record length of each site
  nsitecum <- c(0,cumsum(nsite)) # cumulated record length of each site
  Vobs <- computeVstat(X,typenorm=typenorm) # Test statistic on observed data
  # Centering the sites if requested
  if (center){
     mutot <- apply(matrix(sapply(X,apply,2,mean),nrow=P,ncol=N),1,mean)
     Xcen <- lapply(X,scale,T,F)
     Xcen <- lapply(Xcen,sweep,2,mutot,"+")
  } else {
     Xcen <- X
  }
  # Pooling of all sites
  Xtot <- do.call(rbind,Xcen)
  # Simulations
  Xsim <- vector("list",Nsim)
  Vsim <- vector("numeric",Nsim)
  # HW test
  if (type == "HW"){
     if (P > 2) stop("HW test works for 1 and 2 variables")
     tautot <- vector("list",P)
     taubar <- matrix(NA,P,4)
     paramKappa <- matrix(NA,P,4)
     for (j in 1:P){
         tautot[[j]] <- matrix(NA,N,4)
         for (i in 1:N){
              tautot[[j]][i,] <- samlmu(X[[i]][,j])
         }
         taubar[j,] <- apply(tautot[[j]],2,function(x) (nsite%*%x)/sum(nsite))
         paramKappa[j,] <- pelkap(lmom = taubar[j,])
     }
     if (P == 2){
        # Parameter of the logistic Gumbel copula
        tauKen <- sapply(X,function(x)cor(x[,1],x[,2],method="kendall"))
        tauMean <- (nsite%*%tauKen)/sum(nsite)
        mparam <- tauMean/(1-tauMean) + 1
     } 
     # Kappa simulation
     for (b in 1:Nsim){
         if (printInfo && b%%10 == 0) print(sprintf("Simulation %s / %s",b,Nsim))
         if (P == 2){
            # Simulation using the copula
            Usim <- simulateLogisticCopula(sum(nsite),mparam)
         } else {
            Usim <- as.matrix(runif(sum(nsite)))
         }
         Xsimtot <- matrix(NA,sum(nsite),P)
         for (j in 1:P) Xsimtot[,j] <- quakap(Usim[,j], paramKappa[j,])
         Xsim[[b]] <- split(as.data.frame(Xsimtot),cut(1:sum(nsite),nsitecum,labels=1:N))
         Vsim[b] <- computeVstat(Xsim[[b]],typenorm=typenorm) # Test statistic
     }     
  }
  # Permutations test
  if (type == "M"){
     # Sampling
     for (b in 1:Nsim){
         if (printInfo && b%%10 == 0) print(sprintf("Simulation %s / %s",b,Nsim))
         Xsimtot <- Xtot[sample(1:sum(nsite),sum(nsite)),]
         Xsim[[b]] <- split(as.data.frame(Xsimtot),cut(1:sum(nsite),nsitecum,labels=1:N))
         Vsim[b] <- computeVstat(Xsim[[b]],typenorm=typenorm) # Test statistic
     }   
  }
  # Bootstrap test
  if (type == "B"){
     # Sampling
     for (b in 1:Nsim){
         if (printInfo && b%%10 == 0) print(sprintf("Simulation %s / %s",b,Nsim))
         Xsimtot <- Xtot[sample(1:sum(nsite),sum(nsite),replace=T),]
         Xsim[[b]] <- split(as.data.frame(Xsimtot),cut(1:sum(nsite),nsitecum,labels=1:N))
         Vsim[b] <- computeVstat(Xsim[[b]],typenorm=typenorm) # Test statistic
     } 
  }
  # Polya test
  if (type == "Ys"){
     # Sampling
     for (b in 1:Nsim){
         if (printInfo && b%%10 == 0) print(sprintf("Simulation %s / %s",b,Nsim))
         Xsim[[b]] <- vector("list",N)
         for (j in 1:N) Xsim[[b]][[j]] <- Xtot[PolyaResampling(1:sum(nsite),nsite[j]),,drop=F]
         Vsim[b] <- computeVstat(Xsim[[b]],typenorm=typenorm) # Test statistic
     } 
  }
  if (type == "Yr"){
     # Preparing the Polya sample
     polyasamp <- Xtot[PolyaResampling(1:sum(nsite),sum(nsite)),,drop=F]
     # Sample using the booststrap on the polya sample
     for (b in 1:Nsim){
         if (printInfo && b%%10 == 0) print(sprintf("Simulation %s / %s",b,Nsim))
         Xsimtot <- polyasamp[sample(1:sum(nsite),sum(nsite),replace=T),,drop=F]
         Xsim[[b]] <- split(as.data.frame(Xsimtot),cut(1:sum(nsite),nsitecum,labels=1:N))
         Vsim[b] <- computeVstat(Xsim[[b]],typenorm=typenorm) # Test statistic
     }
  }
  # Final decision
  Hstat <- (Vobs - mean(Vsim)) / sd(Vsim)
  p.value <- sum(Vsim>Vobs)/Nsim
  return(list(Vobs = Vobs, Vsim = Vsim, H = Hstat, p.value = p.value))
}

#############################################################################
#############################################################################
#  Internally used functions

#---- Multivariate Vstat -----
computeVstat <- function(X, typenorm = "2")
# X : list. The sites, as in the function Lmom.homog.test ;
# typenorm : type of norm to use (see help for the function 'norm' in the base package).
{
  N <- length(X) # Number of sites
  nsite <- sapply(X,function(x) nrow(as.matrix(x))) # record length pf each site
  Tbar <- 0 # Mean of Lcomoment matrices
  Tlist <- vector("list",N)
  for (i in 1:N){
      Tlist[[i]] <- Lcomoment.matrix(as.data.frame(X[[i]]),k=2)$matrix / apply(X[[i]],2,mean) # Lcomoment ratio matrix
      Tbar <- Tbar + nsite[i]*Tlist[[i]] # their mean
  }
  Tbar <- Tbar / sum(nsite)
  Vstat <- 0
  for (i in 1:N){
      Vstat <- Vstat + nsite[i]*norm(Tlist[[i]] - Tbar, type=typenorm)^2 # compute the Vstat from Lcomoment ratio matrices
  }
  Vstat <- sqrt(Vstat/sum(nsite))
  return(Vstat)
}

#----- Simulate from the algorithm of Ghoudi et al (1998) ----
simulateLogisticCopula <- function(n,m)
# n : numeric. Number of observations to generate ;
# m : numeric. Parameter of the copula.
{
  Z <- Qntl_GL(runif(n),m)
  AZ <- A_GL(Z,m)
  pZ <- 1 - 1/m
  B <- rbinom(n,1,pZ)
  u1 <- runif(n)
  u2 <- runif(n)
  W <- vector("numeric",n)
  for (i in 1:n) W[i] <- ifelse(B[i]==1,u1[i],u1[i]*u2[i])
  U1 <- W^(Z / AZ)
  U2 <- W^((1-Z) / AZ)
  return(cbind(U1,U2))
}


Qntl_GL <- function(u,m){
  return((u^(1/m))/(u^(1/m)+(1-u)^(1/m)))
}

A_GL <- function(u,m){
  return((u^m + (1-u)^m)^(1/m))
}

#------  Polya sampling -------
PolyaResampling <- function(x,n)
# x : sample from which the resampling is done
# n : the number of records to sample
{
  xsamp <- vector(class(x[1]),n)
  for (i in 1:n){
      xsamp[i] <- sample(x,1)
      x <- c(x,xsamp[i])
  }
  return(xsamp) 
}
