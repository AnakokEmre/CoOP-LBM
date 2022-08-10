
#LBM
#' Variational Expectation Maximization for LBM
#'
#' @param connectivity Binary matrix of connectivity
#' @param Q1 Number of row clusters
#' @param Q2 Number of column clusters
#' @param Z1 Initial clustering of rows
#' @param Z2 Initial clustering of columns
#' @param param List of parameters
#'
#' @return Estimated LBM parameters and clustering
#' @export
#'
#' @import parallel

fit_supervised_LBM<-function(connectivity,Q1,Q2,Z1=c(),Z2=c(),estimOptions=list()) {
  param <- list(
    maxIter       = 50,
    fixPointIter  = 3,
    threshold     = 1e-3,
    ICL_function  = LBM_ICL,
    initMethod = "hierarchical_clust"
  )
  param[names(estimOptions)] <- estimOptions

  if ((length(Z1)==0)|(length(Z2)==0)){
    cl0=clustinit_LBM(connectivity,Q1,Q2,type = param$initMethod)
    if(length(Z1)==0){
      Z1 = cl0[[1]]
    }
    if (length(Z2)==0){
      Z2 = cl0[[2]]
    }
  }

  tau1 = clustering_indicator(Z1,Q1)
  tau2 = clustering_indicator(Z2,Q2)

  delta     <- vector("numeric", param$maxIter)
  pi=matrix(0.5,Q1,Q2)
  i <- 0;
  cond <- FALSE

  while (!cond) {

    i <- i + 1

    ## ______________________________________________________
    ## M-step
    #
    # update the parameters of the SBM (a.k.a alpha and pi)
    alpha1=LBM_update_alpha(tau1)
    alpha2=LBM_update_alpha(tau2)

    new_pi=LBM_update_pi(connectivity,tau1,tau2)
    ## ______________________________________________________
    ## Variational E-Step
    #
    for (k in seq.int(param$fixPointIter)) {
      tau=LBM_update_tau(connectivity,alpha1,alpha2,new_pi,tau1,tau2)
      tau1=tau[[1]]
      tau2=tau[[2]]
    }
    ## Check convergence
    delta[i] <- sqrt(sum((new_pi - pi)^2)) / sqrt(sum((pi)^2))
    cond     <- (i > param$maxIter) |  (delta[i] < param$threshold)
    pi=new_pi
  }

  members1=memberships(tau1)
  members2=memberships(tau2)
  icl=LBM_ICL(members1,members2,alpha1,alpha2,pi,connectivity)
  res=list()
  res$tau1=tau1
  res$tau2=tau2
  res$pi=pi
  res$alpha1=alpha1
  res$alpha2=alpha2
  res$ICL=icl
  res$membership1=members1
  res$membership2=members2
  res$cluster1 = membertoclust(members1)
  res$cluster2 = membertoclust(members2)
  return(res)
}



#' Corrected Observation Process for Latent Block Model
#'
#' @param R Counting data connectivity matrix
#' @param Q1 Number of row clusters
#' @param Q2 Number of columns clusters
#' @param Z1 Initial clustering of rows
#' @param Z2 Initial clustering of columns
#' @param estimOptions List of parameters
#'
#' @return Estimated LBM parameters, clustering, lambda, mu and G for a given number of groups.
#' @export
#'
fit_supervised_CoOP_LBM<-function(R,Q1,Q2,Z1=c(),Z2=c(),estimOptions=list()) {

  param <- list(
    maxIter       = 50,
    maxHeat       = 50,
    fixPointIter  = 3,
    cores         = 1,
    ICL_function  = LBM_ICL_3,
    initMethod="hierarchical_clust"
  )
  param[names(estimOptions)] <- estimOptions

  connectivity = 1*(R>0)
  if ((length(Z1)==0)|(length(Z2)==0)){
    cl0=clustinit_LBM(connectivity,Q1,Q2,type = param$initMethod)
    if(length(Z1)==0){
      Z1 = cl0[[1]]
    }
    if (length(Z2)==0){
      Z2 = cl0[[2]]
    }
  }


  n1 = dim(R)[1]
  n2 = dim(R)[2]

  V = connectivity
  sumV = sum(V)
  rowSumsR = rowSums(R)
  colSumsR = colSums(R)
  lfactorialR = lfactorial(R)

  i <- 0;
  cond <- TRUE


  pi=LBM_update_pi3(connectivity ,Z1,Z2,Q1,Q2)


  while (cond) {

    i <- i + 1
    ## ______________________________________________________
    ## M-step
    alpha1=LBM_update_alpha3(Z1,Q1)
    alpha2=LBM_update_alpha3(Z2,Q2)
    lambda_mu=LBM_update_lambda_mu(rowSumsR ,colSumsR ,connectivity ,fixPointIter = param$fixPointIter)

    ## ______________________________________________________
    ##simulation step

    connectivity = LBM_update_connectivity3(V,Z1,Z2,pi,lambda_mu$lambda_mu)

    ## ______________________________________________________
    ## M-step

    pi=LBM_update_pi3(connectivity,Z1,Z2,Q1,Q2)

    ## ______________________________________________________
    ##simulation step

    new_Z = LBM_update_Z(R,alpha1,alpha2,pi,lambda_mu,Z1,Z2,lfactorialR)

    Z1 = new_Z[[1]]
    Z2 = new_Z[[2]]

    ## ______________________________________________________
    ## condition

    cond     <- (i < (param$maxHeat+param$maxIter))


    if (i == param$maxHeat){
      res_alpha1 = alpha1
      res_alpha2 = alpha2
      res_mem1 = clustering_indicator(Z1,Q1)
      res_mem2 = clustering_indicator(Z2,Q2)
      res_pi = pi
      res_lambda = lambda_mu$lambda_i
      res_mu = lambda_mu$mu_j
      j=1
    }
    if (i> param$maxHeat){
      j= j+1
      res_alpha1 = (1-1/(j+1))*res_alpha1 + alpha1/(j+1)
      res_alpha2 = (1-1/(j+1))*res_alpha2 + alpha2/(j+1)
      res_mem1 = (1-1/(j+1))*res_mem1 + clustering_indicator(Z1,Q1)/(j+1)
      res_mem2 = (1-1/(j+1))*res_mem2 + clustering_indicator(Z2,Q2)/(j+1)
      res_pi = (1-1/(j+1))*res_pi + pi/(j+1)
      res_lambda = (1-1/(j+1))*res_lambda + lambda_mu$lambda_i/(j+1)
      res_mu = (1-1/(j+1))*res_mu+ lambda_mu$mu_j/(j+1)
    }



  }


  icl=param$ICL_function(res_mem1,res_mem2,res_alpha1,res_alpha2,res_pi,R,(res_lambda)%*%t(res_mu))
  res=list()
  res$pi=res_pi
  res$alpha1=as.numeric(res_alpha1)
  res$alpha2=as.numeric(res_alpha2)
  res$ICL=icl
  res$membership1=res_mem1
  res$membership2=res_mem2
  res$cluster1 = membertoclust(res_mem1)
  res$cluster2 = membertoclust(res_mem2)
  res$lambda = res_lambda/max(res_lambda)
  res$mu = res_mu/max(res_mu)
  res$lambda_mu_G = res$lambda%*%t(res_mu)
  res$G = max(res_mu)
  res$connectivity_prob=V
  res$connectivity_prob[V==0] = LBM_connectivity_prob3(V,res$cluster1,res$cluster2,res$pi,res$lambda_mu_G)
  res$observation_prob = V
  res$observation_prob[V==0] = (res$pi[res$cluster1,res$cluster2][V==0] * (1-exp(-res$lambda_mu_G[V==0])))
  res$row_coverage = rowSums(V)/rowSums(res$connectivity_prob)
  res$col_coverage = colSums(V)/colSums(res$connectivity_prob)

  return(res)
}







#' forward exploration
#' @noRd
#' @param models list of models
#' @param k1 number of row clusters
#' @param k2 number of column clusters
#' @param connectivity connectivity matrix
#' @param f method (VEM or EDD-VEM)
#' @param param list of parameters
#' @export
#' @return Forward exploration of a given list of model
forward_explo<-function(models,k1,k2,connectivity,f,param){
  cl01<-membertoclust(models[[paste(as.character(k1),as.character(k2),sep="-")]]$membership1)
  cl02<-membertoclust(models[[paste(as.character(k1),as.character(k2),sep="-")]]$membership2)
  n1=length(unique(cl01))
  n2=length(unique(cl02))
  #print(n1)
  #print(n2)
  best_one=list()
  if (n1==k1&n2==k2){
    candidates <- mclapply(1:n1, function(j) {
      cl1 <- cl01
      J  <- which(cl1 == j)
      if (length(J) > 1) {
        J1 <- base::sample(J, floor(length(J)/2))
        J2 <- setdiff(J, J1)
        cl1[J1] <- j;
        cl1[J2] <- n1 + 1
        model=f(connectivity,n1+1,k2,cl1,cl02,param)
      }
      else {
        model=models[[paste(as.character(k1+1),as.character(k2),sep="-")]]
      }
      return(model)
    },mc.cores = param$cores)
    best_one[[1]] <- candidates[[which.max(sapply(candidates, function(candidate) candidate$ICL))]]

    candidates <- mclapply(1:n2, function(j) {
      cl2 <- cl02
      J  <- which(cl2 == j)
      if (length(J) > 1) {
        J1 <- base::sample(J, floor(length(J)/2))
        J2 <- setdiff(J, J1)
        cl2[J1] <- j;
        cl2[J2] <- n2 + 1
        model=f(connectivity,k1,k2+1,cl01,cl2,param)
      }
      else {
        model=models[[paste(as.character(k1),as.character(k2+1),sep="-")]]
      }
      return(model)
    },mc.cores = param$cores)
    best_one[[2]] <- candidates[[which.max(sapply(candidates, function(candidate) candidate$ICL))]]
  }
  else{
    best_one[[1]]=models[[paste(as.character(k1+1),as.character(k2),sep="-")]]
    best_one[[2]]=models[[paste(as.character(k1),as.character(k2+1),sep="-")]]
  }
  return(best_one)
}

#' backward exploration
#' @noRd
#' @param models list of models
#' @param k1 number of row clusters
#' @param k2 number of column clusters
#' @param connectivity connectivity matrix
#' @param f method (VEM or EDD-VEM)
#' @param param list of parameters
#' @export
#'
#' @return Forward exploration of a given list of model
backward_explo<-function(models,k1,k2,connectivity,f,param){
  cl01<-membertoclust(models[[paste(as.character(k1),as.character(k2),sep="-")]]$membership1)
  cl02<-membertoclust(models[[paste(as.character(k1),as.character(k2),sep="-")]]$membership2)
  cl1 <- factor(cl01)
  cl2 <- factor(cl02)
  n1=nlevels(cl1)
  n2=nlevels(cl2)
  #print(n1)
  #print(n2)
  best_one=list()
  if (n1==k1&n2==k2){
    candidates<-mclapply(combn(n1, 2, simplify = FALSE),function(couple){
      cl_fusion1 <- cl1
      levels(cl_fusion1)[which(levels(cl_fusion1) == paste(couple[1]))] <- paste(couple[2])
      levels(cl_fusion1) <- as.character(1:(n1-1))
      cl_fusion1<-as.numeric(cl_fusion1)
      model=f(connectivity,n1-1,k2,cl_fusion1,cl02,param)
      return(model)
    },mc.cores = param$cores)
    best_one[[1]] <- candidates[[which.max(sapply(candidates, function(candidate) candidate$ICL))]]

    candidates<-mclapply(combn(n2, 2, simplify = FALSE),function(couple){
      cl_fusion2 <- cl2
      levels(cl_fusion2)[which(levels(cl_fusion2) == paste(couple[1]))] <- paste(couple[2])
      levels(cl_fusion2) <- as.character(1:(n2-1))
      cl_fusion2<-as.numeric(cl_fusion2)
      model=f(connectivity,k1,n2-1,cl01,cl_fusion2,param)
      return(model)
    },mc.cores = param$cores)
    best_one[[2]] <- candidates[[which.max(sapply(candidates, function(candidate) candidate$ICL))]]
  }
  else{
    best_one[[1]]=models[[paste(as.character(k1-1),as.character(k2),sep="-")]]
    best_one[[2]]=models[[paste(as.character(k1),as.character(k2-1),sep="-")]]
  }
  return(best_one)
}


#' Main algorithm for LBM estimation using binary data
#'
#' @param connectivity Binary connectivity matrix
#' @param estimOptions List of options about estimation
#' @param exploOptions List of options about exploration
#'
#' @return List of models fitted on different number of cluster
#' @import robber
#' @export
#'
#' @examples
#' n1 = 100
#' n2 = 100
#' Q1=3
#' Q2=3
#' alpha1=c(.25,0.25, .5)
#' alpha2=c(0.10,0.4,0.5)
#' P <- matrix(c(0.9,0.6,0.4,0.7,0.5,0.3,0.5,0.3,0.1), Q1, Q2)
#'simulation1=robber::simulate_lbm(P,alpha1,alpha2,n1,n2)
#'M = simulation1$A
#'Z1 =simulation1$Z
#'Z2 =simulation1$W
#'G= 300
#'lambda_i =rbeta(n1,0.3,1.5)
#'mu_j = rbeta(n2,0.3,1.5)
#'lambda_i = lambda_i/max(lambda_i)
#'mu_j = mu_j/max(mu_j)
#'N0=lambda_i%*%t(mu_j)
#'N0 = N0*G
#'N=matrix(rpois(n1*n2,N0),nrow=n1)
#'R = M*N
#'obsrow = rowSums(R)>0
#'obscol = colSums(R)>0
#'R_obs = R[obsrow,obscol]
#'M_obs = M[obsrow,obscol]
#' R = M*N
#' models=fit_unsupervised_LBM(R,exploOptions=list(plot=F))
fit_unsupervised_LBM<-function(connectivity,estimOptions=list(),exploOptions=list()){

  current_estimOptions <- list(
    maxIter       = 50,
    fixPointIter  = 3,
    threshold     = 1e-3,
    cores         = 1
  )
  current_estimOptions[names(estimOptions)] <- estimOptions


  current_exploOptions <- list(
    initMethod="hierarchical_clust",
    verbosity     = 1,
    plot          = TRUE,
    maxExplo      = 1.5,
    maxGroups     = 10,
    reinitialize  = FALSE

  )
  current_exploOptions[names(exploOptions)] <- exploOptions


  k<-c(1,1)
  models<-list()
  max=-1e20
  whmax=c(0,0)
  gr=1
  cond=TRUE
  V = 1*(connectivity>0)

  while(cond){
    if (gr==1){
      name=unlist(lapply(1:k[2],function(k2){paste(as.character(k[1]),as.character(k2),sep="-")}))
      mods=mclapply(1:k[2],function(k2){

        if (current_exploOptions$verbosity){print(paste('k={',k[1],',',k2,'}'))}
        cl0=clustinit_LBM(V,k[1],k2,current_exploOptions$initMethod)
        Z1=cl0[[1]]
        Z2=cl0[[2]]
        model<-fit_supervised_LBM(V,k[1],k2,Z1,Z2,current_estimOptions)
        return(model)
      },mc.cores = current_estimOptions$cores)
      k2max=which.max(sapply(mods, function(mod) mod$ICL))

      names(mods)<-name
      models<-c(models,mods)
      if (models[[paste(as.character(k[1]),as.character(k2max),sep="-")]]$ICL>max){
        whmax=c(k[1],k2max)
        max=models[[paste(as.character(k[1]),as.character(k2max),sep="-")]]$ICL
      }
      if (current_exploOptions$plot){LBM_plot(models)}
    }
    else if (gr==2){
      name=unlist(lapply(1:k[1],function(k1){paste(as.character(k1),as.character(k[2]),sep="-")}))
      mods=mclapply(1:k[1],function(k1){
        if (current_exploOptions$verbosity){print(paste('k={',k1,',',k[2],'}'))}

        cl0=clustinit_LBM(V,k1,k[2],current_exploOptions$initMethod)
        Z1=cl0[[1]]
        Z2=cl0[[2]]
        model<-fit_supervised_LBM(V,k1,k[2],Z1,Z2,current_estimOptions)
        return(model)
      },mc.cores = current_estimOptions$cores)
      k1max=which.max(sapply(mods, function(mod) mod$ICL))
      names(mods)<-name
      models<-c(models,mods)
      if (models[[paste(as.character(k1max),as.character(k[2]),sep="-")]]$ICL>max){
        whmax=c(k1max,k[2])
        max=models[[paste(as.character(k1max),as.character(k[2]),sep="-")]]$ICL
      }
      if (current_exploOptions$plot){LBM_plot(models)}
    }

    cond=(
      (k[1]<4)|
        (k[1]<round((current_exploOptions$maxExplo*whmax[1])+0.1))|
        (k[2]<4)|
        (k[2]<round((current_exploOptions$maxExplo*whmax[2])+0.1)))&
      ((k[1]<current_exploOptions$maxGroups)&
         (k[2]<current_exploOptions$maxGroups))

    if ((k[1]<max(4,round((current_exploOptions$maxExplo*whmax[1])+0.1)))&
        (k[2]<max(4,round((current_exploOptions$maxExplo*whmax[2])+0.1)))){
      gr=which.min(k)
      k[which.min(k)]<-k[which.min(k)]+1
    }
    else if((k[1]>=max(4,round((current_exploOptions$maxExplo*whmax[1])+0.1)))&
            (k[2]<max(4,round((current_exploOptions$maxExplo*whmax[2])+0.1)))){
      k[2]<-k[2]+1
      gr=2
    }
    else if((k[1]<max(4,round((current_exploOptions$maxExplo*whmax[1])+0.1)))&
            (k[2]>=max(4,round((current_exploOptions$maxExplo*whmax[2])+0.1)))){
      k[1]<-k[1]+1
      gr=1
    }
  }
  if (current_exploOptions$reinitialize==TRUE){
    max2=max
    it<-0
    cond<-TRUE
    if (current_exploOptions$verbosity){print(k)}
    while(cond&(it<3)){
      it<-it+1
      for (k1 in 1:(k[1]-1)){
        for (k2 in 1:(k[2]-1)){
          if (current_exploOptions$verbosity){print(paste('forward','k={',k1,',',k2,'}'))}
          model<-forward_explo(models,k1,k2,V,fit_supervised_LBM,current_estimOptions)

          if (model[[1]]$ICL>models[[paste(as.character(k1+1),as.character(k2),sep="-")]]$ICL){
            models[[paste(as.character(k1+1),as.character(k2),sep="-")]]=model[[1]]

            if (models[[paste(as.character(k1+1),as.character(k2),sep="-")]]$ICL>max2){
              max2=models[[paste(as.character(k1+1),as.character(k2),sep="-")]]$ICL
            }
          }

          if (model[[2]]$ICL>models[[paste(as.character(k1),as.character(k2+1),sep="-")]]$ICL){
            models[[paste(as.character(k1),as.character(k2+1),sep="-")]]=model[[2]]

            if (models[[paste(as.character(k1),as.character(k2+1),sep="-")]]$ICL>max2){
              max2=models[[paste(as.character(k1),as.character(k2+1),sep="-")]]$ICL
            }
          }
          if (current_exploOptions$plot){LBM_plot(models)}
        }
      }
      for (k1 in c(k[1]:3)){
        for (k2 in c(k[2]:3)){
          if (current_exploOptions$verbosity){print(paste('backward','k={',k1,',',k2,'}'))}
          model<-backward_explo(models,k1,k2,V,fit_supervised_LBM,current_estimOptions)

          if (model[[1]]$ICL>models[[paste(as.character(k1-1),as.character(k2),sep="-")]]$ICL){
            models[[paste(as.character(k1-1),as.character(k2),sep="-")]]=model[[1]]

            if (models[[paste(as.character(k1-1),as.character(k2),sep="-")]]$ICL>max2){
              max2=models[[paste(as.character(k1-1),as.character(k2),sep="-")]]$ICL
            }
          }
          if (model[[2]]$ICL>models[[paste(as.character(k1),as.character(k2-1),sep="-")]]$ICL){
            models[[paste(as.character(k1),as.character(k2-1),sep="-")]]=model[[2]]

            if (models[[paste(as.character(k1),as.character(k2-1),sep="-")]]$ICL>max2){
              max2=models[[paste(as.character(k1),as.character(k2-1),sep="-")]]$ICL
            }
          }
          if (current_exploOptions$plot){LBM_plot(models)}
        }
      }
      if (max2>max){
        max=max2
      }
      else{
        cond=FALSE
      }
    }
  }
  print(paste0("LBM : Best model has been fitted with " ,names(best_ICL(models))," groups."))
  return(models)
}


#' Main algorithm for CoOP-LBM estimation using counting data
#'
#' @param connectivity Counting data connectivity matrix
#' @param estimOptions List of options about estimation
#' @param exploOptions List of options about exploration
#'
#' @return List of model with different number of cluster
#' @import robber
#' @export
#'
#' @examples
#' n1 = 100
#' n2 = 100
#' Q1=3
#' Q2=3
#' alpha1=c(.25,0.25, .5)
#' alpha2=c(0.10,0.4,0.5)
#' P <- matrix(c(0.9,0.6,0.4,0.7,0.5,0.3,0.5,0.3,0.1), Q1, Q2)
#' simulation1=robber::simulate_lbm(P,alpha1,alpha2,n1,n2)
#' M = simulation1$A
#' Z1 =simulation1$Z
#' Z2 =simulation1$W
#' G= 300
#' lambda_i =rbeta(n1,0.3,1.5)
#' mu_j = rbeta(n2,0.3,1.5)
#'lambda_i = lambda_i/max(lambda_i)
#'mu_j = mu_j/max(mu_j)
#'N0=lambda_i%*%t(mu_j)
#'N0 = N0*G
#'N=matrix(rpois(n1*n2,N0),nrow=n1)
#'R = M*N
#'obsrow = rowSums(R)>0
#'obscol = colSums(R)>0
#'R_obs = R[obsrow,obscol]
#'M_obs = M[obsrow,obscol]
#' R = M*N
#' models=fit_unsupervised_CoOP_LBM(R,exploOptions=list(plot=F))
#'

fit_unsupervised_CoOP_LBM<-function(connectivity,estimOptions=list(),exploOptions=list()){

  current_estimOptions <- list(
    maxIter       = 50,
    maxHeat       = 50,
    fixPointIter  = 3,
    cores         = 1,

    ICL_function  = LBM_ICL_3
  )
  current_estimOptions[names(estimOptions)] <- estimOptions


  current_exploOptions <- list(
    initMethod="hierarchical_clust",
    verbosity     = 1,
    plot          = TRUE,
    maxExplo      = 1.5,
    maxGroups     = 10,
    reinitialize  = FALSE

  )
  current_exploOptions[names(exploOptions)] <- exploOptions


  k<-c(1,1)
  models<-list()
  max=-1e20
  whmax=c(0,0)
  gr=1
  cond=TRUE
  V = 1*(connectivity>0)

  while(cond){
    if (gr==1){
      name=unlist(lapply(1:k[2],function(k2){paste(as.character(k[1]),as.character(k2),sep="-")}))
      mods=mclapply(1:k[2],function(k2){

        if (current_exploOptions$verbosity){print(paste('k={',k[1],',',k2,'}'))}

        cl0=clustinit_LBM(V,k[1],k2,current_exploOptions$initMethod)
        tau1=cl0[[1]]
        tau2=cl0[[2]]
        model<-fit_supervised_CoOP_LBM(connectivity,k[1],k2,tau1,tau2,current_estimOptions)
        return(model)
      },mc.cores = current_estimOptions$cores)

      k2max=which.max(sapply(mods, function(mod) mod$ICL))


      names(mods)<-name
      models<-c(models,mods)
      if (models[[paste(as.character(k[1]),as.character(k2max),sep="-")]]$ICL>max){
        whmax=c(k[1],k2max)
        max=models[[paste(as.character(k[1]),as.character(k2max),sep="-")]]$ICL
      }

      if (current_exploOptions$plot){LBM_plot(models)}

    }
    else if (gr==2){
      name=unlist(lapply(1:k[1],function(k1){paste(as.character(k1),as.character(k[2]),sep="-")}))
      mods=mclapply(1:k[1],function(k1){

        if (current_exploOptions$verbosity){print(paste('k={',k1,',',k[2],'}'))}

        cl0=clustinit_LBM(V,k1,k[2],current_exploOptions$initMethod)
        tau1=cl0[[1]]
        tau2=cl0[[2]]
        model<-fit_supervised_CoOP_LBM(connectivity,k1,k[2],tau1,tau2,current_estimOptions)
        return(model)
      },mc.cores = current_estimOptions$cores)

      k1max=which.max(sapply(mods, function(mod) mod$ICL))

      names(mods)<-name
      models<-c(models,mods)

      if (models[[paste(as.character(k1max),as.character(k[2]),sep="-")]]$ICL>max){
        whmax=c(k1max,k[2])
        max=models[[paste(as.character(k1max),as.character(k[2]),sep="-")]]$ICL
      }

      if (current_exploOptions$plot){LBM_plot(models)}
    }

    cond=(
      (k[1]<4)|
        (k[1]<round((current_exploOptions$maxExplo*whmax[1])+0.1))|
        (k[2]<4)|
        (k[2]<round((current_exploOptions$maxExplo*whmax[2])+0.1)))&
      ((k[1]<current_exploOptions$maxGroups)&
         (k[2]<current_exploOptions$maxGroups))

    if ((k[1]<max(4,round((current_exploOptions$maxExplo*whmax[1])+0.1)))&
        (k[2]<max(4,round((current_exploOptions$maxExplo*whmax[2])+0.1)))){
      gr=which.min(k)
      k[which.min(k)]<-k[which.min(k)]+1
    }
    else if((k[1]>=max(4,round((current_exploOptions$maxExplo*whmax[1])+0.1)))&
            (k[2]<max(4,round((current_exploOptions$maxExplo*whmax[2])+0.1)))){
      k[2]<-k[2]+1
      gr=2
    }
    else if((k[1]<max(4,round((current_exploOptions$maxExplo*whmax[1])+0.1)))&
            (k[2]>=max(4,round((current_exploOptions$maxExplo*whmax[2])+0.1)))){
      k[1]<-k[1]+1
      gr=1
    }
  }
  if (current_exploOptions$reinitialize==TRUE){
    max2=max
    it<-0
    cond<-TRUE
    if (current_exploOptions$verbosity){print(k)}
    while(cond&(it<3)){
      it<-it+1
      for (k1 in 1:(k[1]-1)){
        for (k2 in 1:(k[2]-1)){
          if (current_exploOptions$verbosity){print(paste('forward','k={',k1,',',k2,'}'))}
          model<-forward_explo(models,k1,k2,connectivity,fit_supervised_CoOP_LBM,current_estimOptions)

          if (model[[1]]$ICL>models[[paste(as.character(k1+1),as.character(k2),sep="-")]]$ICL){
            models[[paste(as.character(k1+1),as.character(k2),sep="-")]]=model[[1]]

            if (models[[paste(as.character(k1+1),as.character(k2),sep="-")]]$ICL>max2){
              max2=models[[paste(as.character(k1+1),as.character(k2),sep="-")]]$ICL
            }
          }

          if (model[[2]]$ICL>models[[paste(as.character(k1),as.character(k2+1),sep="-")]]$ICL){
            models[[paste(as.character(k1),as.character(k2+1),sep="-")]]=model[[2]]

            if (models[[paste(as.character(k1),as.character(k2+1),sep="-")]]$ICL>max2){
              max2=models[[paste(as.character(k1),as.character(k2+1),sep="-")]]$ICL
            }
          }
          if (current_exploOptions$plot){LBM_plot(models)}
        }
      }
      for (k1 in c(k[1]:3)){
        for (k2 in c(k[2]:3)){
          if (current_exploOptions$verbosity){print(paste('backward','k={',k1,',',k2,'}'))}
          model<-backward_explo(models,k1,k2,connectivity,fit_supervised_CoOP_LBM,current_estimOptions)

          if (model[[1]]$ICL>models[[paste(as.character(k1-1),as.character(k2),sep="-")]]$ICL){
            models[[paste(as.character(k1-1),as.character(k2),sep="-")]]=model[[1]]

            if (models[[paste(as.character(k1-1),as.character(k2),sep="-")]]$ICL>max2){
              max2=models[[paste(as.character(k1-1),as.character(k2),sep="-")]]$ICL
            }
          }
          if (model[[2]]$ICL>models[[paste(as.character(k1),as.character(k2-1),sep="-")]]$ICL){
            models[[paste(as.character(k1),as.character(k2-1),sep="-")]]=model[[2]]
            if (models[[paste(as.character(k1),as.character(k2-1),sep="-")]]$ICL>max2){
              max2=models[[paste(as.character(k1),as.character(k2-1),sep="-")]]$ICL
            }
          }
          if (current_exploOptions$plot){LBM_plot(models)}
        }
      }
      if (max2>max){
        max=max2
      }
      else{
        cond=FALSE
      }
    }
  }
  print(paste0("CoOP-LBM : Best model has been fitted with " ,names(best_ICL(models))," groups."))
  return(models)
}




#' Best ICL
#'
#' @param models a list of model with their ICL
#'
#' @return index of the model which has the best ICL
#' @export
#'
#' @examples
#'
#' n1 = 100
#' n2 = 100
#' Q1=3
#' Q2=3
#' alpha1=c(.25,0.25, .5)
#' alpha2=c(0.10,0.4,0.5)
#' P <- matrix(c(0.9,0.6,0.4,0.7,0.5,0.3,0.5,0.3,0.1), Q1, Q2)
#'simulation1=robber::simulate_lbm(P,alpha1,alpha2,n1,n2)
#'M = simulation1$A
#'Z1 =simulation1$Z
#'Z2 =simulation1$W
#'G= 300
#'lambda_i =rbeta(n1,0.3,1.5)
#'mu_j = rbeta(n2,0.3,1.5)
#'lambda_i = lambda_i/max(lambda_i)
#'mu_j = mu_j/max(mu_j)
#'N0=lambda_i%*%t(mu_j)
#'N0 = N0*G
#'N=matrix(rpois(n1*n2,N0),nrow=n1)
#'R = M*N
#'obsrow = rowSums(R)>0
#'obscol = colSums(R)>0
#'R_obs = R[obsrow,obscol]
#'M_obs = M[obsrow,obscol]
#' R = M*N
#' models=fit_unsupervised_CoOP_LBM(R,exploOptions=list(plot=F))
#' print(best_ICL(models))
best_ICL = function(models){
  ICL=unlist(lapply(models, '[[','ICL'))
  k=which.max(ICL)
  return(k)
}







