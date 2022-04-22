#Utils

#' LBMbar
#'
#' @param X numeric vector or matrix
#'
#' @return 1-X
#'
#' @examples
#' a <- c(0.1,0.2,0.7)
#' LBMbar(a)

LBMbar <- function(X) {
  X.bar <- 1 - X
  return(X.bar)
}



#'Softmax
#'
#' @param x numeric vector or matrix
#'
#' @return softmax function applied to the vector or matrix
#'
#' @examples
#' a<- c(10,20,30)
#' .softmax(a)
.softmax <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}




#' Check boundaries
#'
#' @param x numeric vector or matrix
#' @param zero smallest tolerance allowed
#'
#' @return round element too close of zero or one to a value in ]0,1[
#' @examples
#' a<- c(0.2,0.5,0)
#' log(a)
#' log(check_boundaries(a))
check_boundaries <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x > 1 - zero] <- 1 - zero
  x[x <     zero] <-     zero
  return(x)
}


#' Quadratic form
#'
#' @param A  square matrix
#' @param x  vector
#'
#' @return quadratic form associated to the matrix A applied in x
#'
quad_form <- function(A,x) {t(x) %*% A %*% x}


#' Logistic function
#'
#' @param x numeric vector
#'
#' @return logistic function applied to a vector

#'
logistic <- function(x) {1/(1 + exp(-x))}


#' Logit
#'
#' @param x numeric vector
#'
#' @return logit of each element of the vector
logit    <- function(x) {log(x/(1 - x))}


#' Clustering indicator
#'
#' @param clustering vector of labels
#' @param k integer of maximum number of labels
#'
#' @return clustering indicator with k column
#' @examples a <- c(1,1,1,2,2,2,3,3)
#' clustering_indicator(a,3)
#' clustering_indicator(a,4)
clustering_indicator <- function(clustering,k) {
  nBlocks <- k
  nNodes  <- length(clustering)
  Z <- matrix(0,nNodes, nBlocks)
  Z[cbind(seq.int(nNodes), clustering)] <- 1
  return(Z)
}




#' Clustering initialization
#'
#' @param connectivity Binary matrix
#' @param Q1 number of clusters for rows
#' @param Q2 number of clusters for columns
#' @param type type of initialization : "hierarchical_clust", "spectral_clust" or "kmeans_clust"
#'
#' @return list containing an initial clustering for rows and columns
#' @import stats
#' @examples a<- matrix(0,10,10)
#' library(CoOPLBM)
#' a[1:5,1:5] = runif(25)<0.9
#' a[6:10,1:5] = runif(25)<0.5
#' a[1:5,6:10] = runif(25)<0.3
#' a[6:10,6:10] = runif(25)<0.1
#' print("AAAAAAAAAAAAAAAAAAAAAA")
#' print(a)
#' print(clustinit_LBM(a,2,2))

clustinit_LBM<-function(connectivity,Q1,Q2,type="hierarchical_clust"){
  if (type=="hierarchical_clust"){
    if (Q1 > 1) {
      D <- as.matrix(dist(connectivity, method = "manhattan"))
      D <- as.dist(ape::additive(D))
      cl01 <- cutree(hclust(D, method = "ward.D2"), Q1)
    } else {
      cl01 <- rep(1L,nrow(connectivity))
    }
    if (Q2 > 1) {
      D <- as.matrix(dist(t(connectivity), method = "manhattan"))
      D <- as.dist(ape::additive(D))
      cl02 <- cutree(hclust(D, method = "ward.D2"), Q2)
    } else {
      cl02 <- rep(1L,nrow(t(connectivity)))
    }
    cl0=list(cl01,cl02)
    return(cl0)
  }
  else if (type=="spectral_clust"){
    n1=dim(connectivity)[1]
    n2=dim(connectivity)[2]
    if (Q1 > 1) {
      degrees = apply(connectivity, 1,sum)
      degrees = (1/degrees)%x%matrix(1,1,n2)
      L=degrees*connectivity
      dec=svd(L)
      U = dec$u[,1:Q1]
      km=kmeans(U,Q1)
      tau=rdist(km$centers,U)
      cl01=apply(tau,2,which.max)
    }
    else{
      cl01 <- rep(1L,nrow(connectivity))
    }
    if (Q2 > 1) {
      degrees = apply(connectivity, 2,sum)
      degrees = (1/degrees)%x%matrix(1,1,n1)
      L=degrees*t(connectivity)
      dec=svd(L)
      U = dec$u[,1:Q2]
      km=kmeans(U,Q2)
      tau=rdist(km$centers,U)
      cl02=apply(tau,2,which.max)
    }
    else{
      cl02 <- rep(1L,ncol(connectivity))
    }
    cl0=list(cl01,cl02)
    return(cl0)
  }
  else if (type=="kmeans_clust"){
    if (Q1 > 1) {
      D  <- as.matrix(dist(connectivity, method = "euclidean"))
      cl01 <- as.integer(kmeans(ape::additive(D), Q1, nstart = 50, iter.max = 100)$cl)
    } else {
      cl01 <- rep(1L, nrow(connectivity))
    }
    if (Q2 > 1) {
      D  <- as.matrix(dist(t(connectivity), method = "euclidean"))
      cl02 <- as.integer(kmeans(ape::additive(D), Q2, nstart = 50, iter.max = 100)$cl)
    } else {
      cl02 <- rep(1L, nrow(t(connectivity)))
    }
    cl0=list(cl01,cl02)
    return(cl0)
  }
}

#' Membership to clustering
#'
#' @param membership a clustering indicator matrix
#'
#' @return the clustering associated to the matrix

membertoclust<-function(membership){
  return(apply(membership,1,which.max))
}

#' LBM_plot
#'
#' @param models
#'
#' @return plot of the ICL of the models contained in models
#' @importFrom graphics plot
LBM_plot<-function(models){
  modnames=names(models)
  Q=c()
  ICL=c()
  for (na in modnames){
    q = sum(as.numeric(strsplit(na,split="-")[[1]]))
    Q=c(Q,q)
    ICL<-c(ICL,models[[na]]$ICL)

  }

  plot(Q,ICL,xlab="Q1+Q2",ylab="ICL")
}


#' Memberships
#'
#' @param tau
#'
#' @return hard clustering indicator of a soft clustering
#'
#' @examples a<- matrix(0,5,2)
#' a[,1] = runif(5)
#' a[,2] = 1-a[,1]
#' print(a)
#' print(memberships(a))
memberships = function(tau) {
  dimtau=dim(tau)
  mem=matrix(0,dimtau[1],dimtau[2])
  for (i in 1:dimtau[1]){
    a=which.max(tau[i,])
    mem[i,a]=1
  }
  return(mem)
}



#' Update tau for VEM
#'
#' @param connectivity binary matrix of connectivity
#' @param alpha1 proportion of groups by row
#' @param alpha2 proportion of groups by column
#' @param pi matrix of probabilities of connections between groups
#' @param tau1 soft clustering of rows
#' @param tau2 soft clustering of columns
#'
#' @return new estimation of tau

LBM_update_tau<-function(connectivity,alpha1,alpha2,pi,tau1,tau2){
  barconnectivity=LBMbar(connectivity)
  Q1=length(alpha1)
  Q2=length(alpha2)
  if ((Q1>1)&(Q2>1)){
    newtau1<-connectivity %*% tau2 %*% t(log(pi)) + barconnectivity %*% tau2 %*% t(log(1 - pi))
    newtau1 <- t(apply(sweep(newtau1, 2, log(alpha1), "+"), 1, .softmax))
    newtau2<-t(connectivity) %*% tau1 %*% log(pi) + t(barconnectivity) %*% tau1 %*% log(1 - pi)
    newtau2 <- t(apply(sweep(newtau2, 2, log(alpha2), "+"), 1, .softmax))
  }
  else if ((Q1>1)&(Q2==1)){
    newtau1<-connectivity %*% tau2 %*% t(log(pi)) + barconnectivity %*% tau2 %*% t(log(1 - pi))
    newtau1 <- t(apply(sweep(newtau1, 2, log(alpha1), "+"), 1, .softmax))
    newtau2=tau2
  }
  else if ((Q1==1)&(Q2>1)){
    newtau2<-t(connectivity) %*% tau1 %*% log(pi) + t(barconnectivity) %*% tau1 %*% log(1 - pi)
    newtau2 <- t(apply(sweep(newtau2, 2, log(alpha2), "+"), 1, .softmax))
    newtau1=tau1
  }
  else if ((Q1==1)&(Q2==1)){
    newtau1=tau1
    newtau2=tau2
  }
  return(list(newtau1,newtau2))
}


#' Update pi for VEM
#'
#' @param connectivity binary matrix of connectivity
#' @param tau1 soft clustering of rows
#' @param tau2 soft clustering of columns
#'
#' @return new estimation of pi
LBM_update_pi = function(connectivity,tau1,tau2) {
  N1<-dim(connectivity)[1]
  N2<-dim(connectivity)[2]
  pi    <- check_boundaries(t(tau1)%*%connectivity%*%tau2 / t(tau1)%*%matrix(1,N1,N2)%*%tau2)
  return(pi)
}


#' Update alpha for VEM
#'
#' @param tau soft clustering (row or column)
#'
#' @return new estimation of alpha (row or column)
LBM_update_alpha = function(tau) {
  alpha <- check_boundaries(colMeans(tau))
  return(alpha)
}


#' Update alpha for Gibbs
#'
#' @param Z clustering (row or column)
#' @param Q number of group (row or column)
#'
#' @return new estimation of alpha or beta

LBM_update_alpha3 = function(Z,Q) {
  a = factor(Z, levels= 1:Q)
  return(check_boundaries(prop.table(table(a))))
}


#' Update pi for Gibbs
#'
#' @param connectivity binary connectivity matrix
#' @param Z1 row clustering
#' @param Z2 column clustering
#' @param Q1 number of groups for rows
#' @param Q2 number of groups for columns
#'
#' @return new value of pi given the parameters
#' @export
LBM_update_pi3 = function(connectivity,Z1,Z2,Q1,Q2){
  a=factor(Z1, levels=1:Q1)
  b=factor(Z2, levels=1:Q2)
  pi = check_boundaries(matrix(table(expand.grid(a,b)[connectivity==1,])/table(expand.grid(a,b)),Q1,Q2))
  return(pi)
}


#' Update lambda and mu
#'
#' @param rowSumsR sum of observation by row
#' @param colSumsR sum of observation by column
#' @param connectivity binary connectivity matrix
#' @param fixPointIter number of iteration for the fixed point algorithm
#'
#' @return list containing the estimation of lambda_i, mu_j x G, and the matrix lambda_i x mu_j x G
#' @export
LBM_update_lambda_mu = function(rowSumsR,colSumsR,connectivity,fixPointIter=3){
  lambda_i = rep(1,dim(connectivity)[1])
  for (k in 1:fixPointIter){
    lambda_i[rowSumsR>0] = c(rowSumsR[rowSumsR>0]/(connectivity[rowSumsR>0,colSumsR>0])%*%t(colSumsR[colSumsR>0]/ (lambda_i[rowSumsR>0]%*%(connectivity[rowSumsR>0,colSumsR>0]))))
    lambda_i[rowSumsR==0]=0
  }
  lambda_i = lambda_i/max(lambda_i)
  mu_j=c(colSumsR/t(lambda_i)%*%(connectivity))
  mu_j[colSumsR==0]=0
  res = list(lambda_i= lambda_i, mu_j = mu_j, lambda_mu = lambda_i%*%t(mu_j))
  res
}


#' Update connectivity
#'
#' @param V binary connectivity matrix of the observed matrix
#' @param Z1 clustering by row
#' @param Z2  clustering by column
#' @param pi probability of connection between groups
#' @param lambda_mu matrix containing the product lambda_i x mu_j x G
#'
#' @return simulation of a connectivity matrix with missing links added
#' @importFrom Rlab rbern


LBM_update_connectivity3 = function(V,Z1,Z2,pi,lambda_mu){
  pi_ij = pi[Z1,Z2][V==0]
  exp_lambda_mu = exp(-lambda_mu[V==0])
  P_ij = pi_ij*exp_lambda_mu/((1-pi_ij)+pi_ij*exp_lambda_mu)
  simulation = rbern(length(P_ij),p=P_ij)
  connectivity = V
  connectivity[V==0] = simulation
  connectivity

}


#' Update clustering
#'
#' @param R matrix of observation
#' @param alpha1 proportion of clustering by row
#' @param alpha2 proportion of clustering by column
#' @param pi matrix of probabilities of connections between groups
#' @param lambda_mu matrix containing the product lambda_i x mu_j x G
#' @param Z1 clustering by row
#' @param Z2 clustering by column
#' @param lfactorialR matrix of log factorial of the R matrix
#'
#' @return simulate new values for Z1 and Z2 given the parameters
#' @export
LBM_update_Z<-function(R,alpha1,alpha2,pi,lambda_mu,Z1,Z2,lfactorialR){
  Q1=length(alpha1)
  Q2=length(alpha2)
  n1 = dim(R)[1]
  n2 = dim(R)[2]


  if ((Q1>1)&(Q2>1)){
    logpi = log(pi)
    R0 = R==0
    R1 = !R0
    log_alpha1 = log(alpha1)
    log_alpha2 = log(alpha2)
    log_lambda_mu = log(lambda_mu$lambda_mu)
    exp_lambda_mu = exp(-lambda_mu$lambda_mu)
    A=(R1*(-lambda_mu$lambda_mu-lfactorialR+R*log_lambda_mu))

    newZ1 = matrix(rowSums(A,na.rm=T),n1,Q1)
    for (k in 1:Q1){
      newZ1[,k]=newZ1[,k] +log_alpha1[k] ++R1%*%logpi[k,Z2]+ rowSums(R0*log(1-(sweep((1- exp_lambda_mu),2,pi[k,Z2],"*"))))
    }
    newZ1=apply(newZ1,1,function(x){sample(1:Q1,size=1,prob=.softmax(x))})

    newZ2 = matrix(colSums(A,na.rm=T),n2,Q2)
    for (l in 1:Q2){
      newZ2[,l]=newZ2[,l] +log_alpha2[l] +logpi[newZ1,l]%*%R1+colSums(R0*log(1-(sweep((1- exp_lambda_mu),1,pi[newZ1,l],"*"))))
    }
    newZ2 = apply(newZ2,1,function(x){sample(1:Q2,size=1,prob=.softmax(x))})
  }
  else if ((Q1>1)&(Q2==1)){

    logpi = log(pi)
    R0 = R==0
    R1 = !R0
    log_alpha1 = log(alpha1)
    log_alpha2 = log(alpha2)
    log_lambda_mu = log(lambda_mu$lambda_mu)
    exp_lambda_mu = exp(-lambda_mu$lambda_mu)
    A=(R1*(-lambda_mu$lambda_mu-lfactorialR+R*log_lambda_mu))
    newZ1 = matrix(rowSums(A,na.rm=T),n1,Q1)
    for (k in 1:Q1){
      newZ1[,k]=newZ1[,k] +log_alpha1[k] ++R1%*%logpi[k,Z2]+ rowSums(R0*log(1-(sweep((1- exp_lambda_mu),2,pi[k,Z2],"*"))))
    }
    newZ1=apply(newZ1,1,function(x){sample(1:Q1,size=1,prob=.softmax(x))})
    newZ2= Z2

  }
  else if ((Q1==1)&(Q2>1)){
    logpi = log(pi)
    R0 = R==0
    R1 = !R0
    log_alpha1 = log(alpha1)
    log_alpha2 = log(alpha2)
    log_lambda_mu = log(lambda_mu$lambda_mu)
    exp_lambda_mu = exp(-lambda_mu$lambda_mu)
    A=(R1*(-lambda_mu$lambda_mu-lfactorialR+R*log_lambda_mu))

    newZ2 = matrix(colSums(A,na.rm=T),n2,Q2)
    for (l in 1:Q2){
      newZ2[,l]=newZ2[,l] +log_alpha2[l] +logpi[Z1,l]%*%R1+colSums(R0*log(1-(sweep((1- exp_lambda_mu),1,pi[Z1,l],"*"))))
    }

    newZ2 = apply(newZ2,1,function(x){sample(1:Q2,size=1,prob=.softmax(x))})
    newZ1 = Z1

  }
  else if ((Q1==1)&(Q2==1)){
    newZ1=Z1
    newZ2=Z2
  }
  return(list(newZ1,newZ2))
}


#' Update connectivity
#'
#' @param V binary connectivity matrix
#' @param Z1 clustering of rows
#' @param Z2 clustering of columns
#' @param pi matrix of probabilities of connections between groups
#' @param lambda_mu matrix of all the coefficient lambda x mu x G
#'
#' @return probability of M = 1 given that V = 0. Only return the values where V = 0
#' @export
LBM_connectivity_prob3= function(V,Z1,Z2,pi,lambda_mu){
  pi_ij = pi[Z1,Z2][V==0]
  exp_lambda_mu = exp(-lambda_mu[V==0])
  P_ij = pi_ij*exp_lambda_mu/((1-pi_ij)+pi_ij*exp_lambda_mu)
  P_ij}


#' ICL for LBM
#'
#' @param members1 clustering indicator of rows
#' @param members2 clustering indicator of column
#' @param alpha1 proportion of group by row
#' @param alpha2 proportion of group by column
#' @param pi matrix of probabilities of connections between groups
#' @param connectivity binary connectivity matrix
#'
#' @return estimation of ICL for classic LBM
#' @export
LBM_ICL<-function(members1,members2,alpha1,alpha2,pi, connectivity){
  Q1=length(alpha1)
  Q2=length(alpha2)
  N1=dim(connectivity)[1]
  N2=dim(connectivity)[2]
  barconnectivity=LBMbar(connectivity)
  logL=sum(members1%*%log(alpha1))+sum(members2%*%log(alpha2))+sum((t(members1)%*%connectivity%*%members2)*log(pi))+sum((t(members1)%*%barconnectivity%*%members2)*log(1-pi))
  pena = log(N1)*(Q1-1)/2 + log(N2)*(Q2-1) + (Q1*Q2)/2*log(N1*N2)
  return(logL-pena)
}


#' ICL for CoOP-LBM (binary)
#'
#' @param members1 clustering indicator of rows
#' @param members2 clustering indicator of column
#' @param alpha1 proportion of group by row
#' @param alpha2 proportion of group by column
#' @param pi matrix of probabilities of connections between groups
#' @param R observed connectivity matrix
#' @param lambdamu matrix of lambda_i x lambda_j
#'
#' @return estimation of ICL for EDD LBM using the likelihood not taking into account the Poisson distribution
#' @export
LBM_ICL_2<-function(members1,members2,alpha1,alpha2,pi, R,lambdamu){
  V = data.matrix(1*(R>0))
  Q1=length(alpha1)
  Q2=length(alpha2)
  N1=dim(V)[1]
  N2=dim(V)[2]
  pi_ij = pi[membertoclust(members1),membertoclust(members2)]*(1-exp(-lambdamu))
  L1 = log(pi_ij)
  L2 = log(1-pi_ij)
  L3 = ifelse(V,L1,L2)

  logL=sum(members1%*%log(alpha1))+sum(members2%*%log(alpha2))+sum(L3)
  pena=log(N1)*(Q1-1)/2 + log(N2)*(Q2-1) + (Q1*Q2+N1+N2-1)/2*log(N1*N2)
  return(logL-pena)
}


#' ICL for CoOP-LBM (Poisson)
#'
#' @param members1 clustering indicator of rows
#' @param members2 clustering indicator of column
#' @param alpha1 proportion of group by row
#' @param alpha2 proportion of group by column
#' @param pi matrix of probabilities of connections between groups
#' @param R observed connectivity matrix
#' @param lambdamu matrix of lambda_i x lambda_j
#'
#' @return estimation of ICL for EDD LBM using the likelihood taking into account the Poisson distribution
#' @export
LBM_ICL_3<-function(members1,members2,alpha1,alpha2,pi, R,lambdamu){
  V = 1*(R>0)
  Q1=length(alpha1)
  Q2=length(alpha2)
  N1=dim(V)[1]
  N2=dim(V)[2]
  pi_ij = pi[membertoclust(members1),membertoclust(members2)]
  L1 = sum((log(pi_ij)-lambdamu)[V>0])
  L2 = sum(R*log(lambdamu) - lfactorial(R))
  L3 = sum(1-pi_ij[V==0]*(1-exp(-lambdamu[V==0])))

  logL=sum(members1%*%log(alpha1))+sum(members2%*%log(alpha2))+L1 +L2 + L3
  pena=log(N1)*(Q1-1)/2 + log(N2)*(Q2-1) + (Q1*Q2+N1+N2-1)/2*log(N1*N2)
  return(logL-pena)
}






#' Variational ICL for EDD-LBM
#'
#' @param members1 soft clustering of rows
#' @param members2 soft clustering of column
#' @param alpha1 proportion of group by row
#' @param alpha2 proportion of group by column
#' @param pi matrix of probabilities of connections between groups
#' @param R observed connectivity matrix
#' @param lambdamu matrix of lambda_i x lambda_j
#'
#' @return Variational estimation of ICL for EDD LBM using the likelihood  not taking into account the Poisson distribution
#' @export
LBM_ICL_2.2<-function(members1,members2,alpha1,alpha2,pi, R,lambdamu){
  V = data.matrix(1*(R>0))
  Q1=length(alpha1)
  Q2=length(alpha2)
  N1=dim(V)[1]
  N2=dim(V)[2]
  L = 0

  for ( i in 1:Q1){
    for (j in 1:Q2){

      pi_ij = pi[i,j] *(1 -exp(-lambdamu))
      pQ1Q2 = (members1[,i]%*%t(members2[,j]))
      L1 = log(pi_ij)
      L2 = log(1-pi_ij)
      L3 = ifelse(V,L1,L2)
      L=L+sum(L3*pQ1Q2)

    }
  }
  logL=L+sum(members1%*%log(alpha1))+sum(members2%*%log(alpha2))
  pena=log(N1)*(Q1-1)/2 + log(N2)*(Q2-1) + (Q1*Q2+N1+N2-1)/2*log(N1*N2)
  return(logL-pena)
}



NODF = function(Matrix){
  M = Matrix[rowSums(Matrix)>0,colSums(Matrix)>0]
  n1 =dim(M)[1]
  n2 = dim(M)[2]
  Mbis = 1*(M>0)
  M2 = Mbis[,order(colSums(Mbis),decreasing=T)]
  M2 = M2[order(rowSums(Mbis),decreasing=T),]

  CS=colSums(M2)
  Ncol=0
  for (i in 1:(n2-1)){
    for (j in (i+1):n2){
      if (CS[i]>CS[j]){
        Ncol = Ncol + M2[,i]%*%(M2[,j]/sum(M2[,j]))
      }
    }
  }

  RS=rowSums(M2)
  Nrow=0
  for (i in 1:(n1-1)){
    for (j in (i+1):n1){
      if (RS[i]>RS[j]){
        Nrow = Nrow + M2[i,]%*%(M2[j,]/sum(M2[j,]))
      }
    }
  }

  list(row =Nrow/(n1*(n1-1)/2), col =Ncol/(n2*(n2-1)/2), matrix = (Ncol +Nrow)/((n1*(n1-1)/2)+(n2*(n2-1)/2))  )
}













