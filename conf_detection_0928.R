library(RDRToolbox)
library(GPfit)
library(pracma)
library(scales)
library(dHSIC)
conf_detection<-function(X,Y,N_initial_steps=5,N_projection_steps=10,N_fminsearch_iter=100,level=.05,N_neighbours=10,rescale_T=TRUE){
  # Inital Dimensionality Reduction
  # Initial guess: Isomap to get initial T.
  data=cbind(X,Y)
  rescale_T_vec=rep(rescale_T,dim(data)[1])
  T_est=Isomap(data, dims = 1, k=N_neighbours,verbose=FALSE)[[1]] # T hat, from two dimensional to one dimension
  for (k in 1:N_initial_steps){
    T_est=ifelse(rescale_T_vec,rescale(T_est),T_est)
    gp1=GP_fit(T_est, X)
    gp2=GP_fit(T_est, Y)
    fminsearch_output=fminsearch(f=function(t,gp1,gp2){
      gp1_prediction=predict.GP(gp1,t)
      gp2_prediction=predict.GP(gp2,t)
      return(sum((X-gp1_prediction$Y_hat)^2+(Y-gp2_prediction$Y_hat)^2))
    }, maxiter = N_fminsearch_iter, x0=T_est,gp1=gp1,gp2=gp2)
    T_est=fminsearch_output$xval
  }
  
  ##########################################################################
  # Dependence Criterion and its Minimization
  k=0
  Conti=TRUE
  while((k<N_projection_steps) & Conti){
    k=k+1
    T_est=ifelse(rescale_T_vec,rescale(T_est),T_est)
    gp1=GP_fit(T_est, X)
    gp2=GP_fit(T_est, Y)
    fminsearch_output=fminsearch(f=function(t,gp1,gp2){
      gp1_prediction=predict.GP(gp1,t)
      gp2_prediction=predict.GP(gp2,t)
      Nx=X-gp1_prediction$Y_hat
      Ny=Y-gp2_prediction$Y_hat
      return(dhsic.test(Nx,Ny,method="gamma")$statistic+dhsic.test(Nx,t,method="gamma")$statistic+
               dhsic.test(Ny,t,method="gamma")$statistic)
    }, maxiter = N_fminsearch_iter, x0=T_est,gp1=gp1,gp2=gp2)
    T_est=fminsearch_output$xval
    gp1_prediction=predict.GP(gp1,T_est)
    gp2_prediction=predict.GP(gp2,T_est)
    Nx=X-gp1_prediction$Y_hat
    Ny=Y-gp2_prediction$Y_hat
    pHSIC_Nx_Ny=dhsic.test(Nx,Ny,method="gamma")$p.value
    pHSIC_Nx_T_est=dhsic.test(Nx,T_est,method="gamma")$p.value
    pHSIC_Ny_T_est=dhsic.test(Ny,T_est,method="gamma")$p.value
    Conti=!((pHSIC_Nx_Ny>level) & (pHSIC_Nx_T_est>level) & (pHSIC_Ny_T_est>level))
  }
  Var_Nx=var(Nx)
  Var_Ny=var(Ny)
  list(confounder_detected=!Conti,T_est=T_est,Nx=Nx,Ny=Ny,projection_k=k,pHSIC_Nx_Ny=pHSIC_Nx_Ny, 
       pHSIC_Nx_T_est= pHSIC_Nx_T_est,pHSIC_Ny_T_est=pHSIC_Ny_T_est,ratio_Nx_Ny=Var_Nx/Var_Ny,
       u=gp1,v=gp2)
}
