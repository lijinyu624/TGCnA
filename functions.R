library(clusterProfiler)
library(fda)
library(ggplot2)
library(ggpubr)
library(Matrix)
library(pdfCluster)
library(WGCNA)
library(matrixcalc)
library(MASS)
library(parallel)
library(doParallel)
library(robustbase)
library(reshape2)

TGCN = function(tvec, numreps, obs)
{
    nt = length(numreps)
    
    ## centralize the original data
    res.mtx_FullSample = Normalize.res.mtx(obs,numreps)
    
    ## choose optimal K (number of component)
    chooseK_out  = chooseK(res.mtx_FullSample, numreps, nt)
    mdl = chooseK_out$logMDL
    aic = chooseK_out$logAIC
    K 	= which.min(mdl)
    print(paste("The estimated number of components K:",K))
    
    ## get eigen vector by K
    GetFactorRes = GetFactor(res.mtx_FullSample, numreps, K)
    U 	= GetFactorRes$u
    if(is.vector(U)){U=matrix(U)}
    
    ## obtain smoothed D square
    D_2 = B_Spline_for_D2(res.mtx_FullSample, numreps, tvec, U, nbasis=4)
    
    corrlist = list()
    for(tt in 1:nt)
    {
      
      TimePoint = tt
      print(paste("Time Point :",TimePoint))
      ## slice time t D square
      Dt_2 = D_2[,TimePoint]
      if(length(Dt_2)>1){Dt_2 = diag(D_2[,TimePoint])}	  
      
      ## slice time t X matrix
      Xt = res.mtx_FullSample[,(sum(numreps[1:TimePoint-1])+1):sum(numreps[1:TimePoint])]
      
      ## obtain time t total covariance
      sigt2 = Xt%*%t(Xt)/(ncol(Xt)-1)
      
      ## obtain time t total covariance
      if(length(Dt_2)==1){St = sigt2 - (U*Dt_2)%*%t(U)}
      if(length(Dt_2)>1){St = sigt2 - U%*%Dt_2%*%t(U)}
      
      ## low rank matrix
      if(length(Dt_2)==1){LowRank = (U*Dt_2)%*%t(U)}
      if(length(Dt_2)>1){LowRank  = U%*%Dt_2%*%t(U)}
      
      ## obtain time t residual covariance
      St = sigt2 - LowRank
      
      ## initialize threshold matrix
      taut = matrix(0,nrow(St),ncol(St))
      
      ## get sparse component
      St_shrink = getSparseComponent(sigt2, St, taut)
      
      ## Return the percentage of negative eigenvalues of reconstructed covariance
      sigt_2 = U%*%Dt_2%*%t(U) + St_shrink
      
      ## obtain correlation matrix
      corrlist[[TimePoint]] = cov2cor(sigt_2)
      
    }

    return(corrlist)
}


#### Sparsification
getSparseComponent = function(sigt2, St, taut)
{  
    c   = 0.5
    eta = 4
    St_shrink = matrix(0,nrow(St),ncol(St))
    ## generate sparse residual covariance
    for(i in 1:nrow(St))
    {
      for(j in 1:ncol(St))
      {
        if(i!=j)
        {
          # calculate thrshold matrix
          taut[i,j] = c*(sigt2[i,i]*sigt2[j,j])^0.5 
          # if a non-diagonal element is greater than the threshold
          if(abs(St[i,j])>taut[i,j]) 
          {
            # use shrinkage function for non-diagondal elements greater than threshold
            St_shrink[i,j] = St[i,j]*(1-taut[i,j]/abs(St[i,j]))^eta 
          }
        }
        # if a diagonal element is greater than 0
        else if(St[i,j] >= 0)
        {
          # if a St diagonal element is greater than 0, reserve it.
          St_shrink[i,j] = St[i,j] 
        }
        else
        {
          # if a St diagonal element is smaller than 0, set it to 0.
          St_shrink[i,j] = 0 
        }
      }
    }
    
    return(St_shrink)
}


#### TGCN Clustering
TGCN_Cluster = function(corrlist, newdata)
{
  nt = length(corrlist)
  geneID = c(1:nrow(newdata))
  ClustRes = lapply(1:nt, function(tt) Cluster(corrlist[[tt]], 
                                                geneID,
                                                softPower     = 6, 
                                                type 			    = 'unsigned', 
                                                method 		    = c('average'), 
                                                minModuleSize = 30))
  
  mergeRes = lapply(1:nt, function(tt) mergemods(corrlist[[tt]], 
                                                nt, 
                                                ClustRes[[tt]], 
                                                MEDissThres = 0.25,
                                                newdata))
  
  modres   = lapply(1:nt, function(tt) mergeRes[[tt]]$moduleLabels)
    
  return(modres)
}


#### select the best number of components
chooseK = function(res.mtx, numreps,nt)
{
  n = ncol(res.mtx)
  p = nrow(res.mtx)
  
  svd = svd(res.mtx/sqrt(sum(numreps)-length(numreps)))

  cl  = makeCluster(detectCores()-1)
  registerDoParallel(cl)
  output = foreach(K = 1:(n-nt-1))%dopar%
  {
    print(K)
    U = as.matrix(svd$u[,1:K])
    D = svd$d[1:K]
    if(length(D)>1){D = diag(D)}	  
    
    OmegaS = res.mtx%*%t(res.mtx)/(sum(numreps)-length(numreps))
    sigma2 = mean(diag(OmegaS - U%*%(D^2)%*%t(U)))
    d      = sum(numreps)-nt
    l      = svd$d^2  + sigma2
    rho_K = prod(l[(K+1):d]^(1/(d-K)))/(1/(d-K)*sum(l[(K+1):d]))	
    list(log(-n*(d - K)*log(rho_K) + K/2*(2*d - K)*log(n)), log(-2*n*(d - K)*log(rho_K)+2*K*(2*d-K)))
  }	
  
  stopCluster(cl)
  output = unlist(output)
  logMDL = output[seq(1,length(output),2)]
  logAIC = output[seq(2,length(output),2)]
  return(list('logMDL' = logMDL,'logAIC' = logAIC))
}



#### Smooth the eigenvalue
B_Spline_for_D2 = function(Sample, numreps, tvec, U, nbasis)
{
  nt = length(numreps)
  K  = ncol(U)
  if(is.null(K)){K=1}
  logtvec = log(tvec)
  varmtx  = matrix(0, K, nt)
  
  # get eigenvalue
  for(tt in 1:length(numreps))
  {
    eta.tt = Sample[,(1+sum(numreps[0:(tt-1)])):sum(numreps[1:tt])]
    cov.tt = t(eta.tt)%*%U
    varmtx[,tt] = diag(t(cov.tt)%*%cov.tt)/(numreps[tt]-1)
  }
  
  # use b-spline method to smooth varmtx
  bs.basis = create.bspline.basis(rangeval = range(logtvec), nbasis = nbasis)
  vart     = t(sapply(1:nrow(varmtx), function(ii) exp(eval.fd(logtvec, smooth.basis(argvals = logtvec, y = log(varmtx[ii,]),fdParobj = bs.basis)$fd))))
  
  return(vart)
}



### get factors using pc method
GetFactor = function(res.mtx, numreps, K){
  eig.w = svd(res.mtx/sqrt(sum(numreps)-length(numreps)))
  u = eig.w$u[,1:K] 
  d = eig.w$d[ 1:K]
  v = eig.w$v[,1:K]
  return(list('u' = u, 'd' = d, 'v' = v))
}



#### Cluster for close module 
Cluster = function(cormtx, geneID, softPower, method = c('average'), type, minModuleSize)
{
  adjacency 	  = adjacency.fromSimilarity(cormtx*0.999999, type = type, power = softPower);
  dissTOM 		  = 1-TOMsimilarity(adjacency);
  geneTree 		  = hclust(as.dist(dissTOM), method = method);
  minModuleSize = minModuleSize;
  dynamicMods 	= cutreeDynamic(dendro   = geneTree, 
                               distM     = dissTOM,
                               deepSplit = 3, 
                               pamRespectsDendro = FALSE,
                               minClusterSize    = minModuleSize);
  genelist = list()
  for (i in 1:(max(dynamicMods)))
  {
    genelist[[i]] = geneID[c(which(as.numeric(dynamicMods)==i))]
  }
  return(list('dissTOM' = dissTOM, 'dynamicMods' = dynamicMods, 'geneTree' = geneTree, 'genelist' = genelist))
}



#### merge close modules
mergemods = function(cormtx, nt, ClustRes, MEDissThres, newdata)
{
  merge = mergeCloseModules(cormtx, labels2colors(ClustRes$dynamicMods), cutHeight = MEDissThres, verbose = 3)
  mergecolor = merge$color
  colorOrder = c("grey", standardColors(100));
  moduleLabels = match(mergecolor, colorOrder)-1
  genelist.merge = list()
  for (i in 1:(max(moduleLabels)))
  {
    genelist.merge[[i]]= row.names(newdata)[c(which(as.numeric(moduleLabels)==i))]
  }
  genelist.merge = genelist.merge[lapply(genelist.merge,length)!=0] 
  return(list('merge' = merge, 'mergecolor' = mergecolor, 'moduleLabels' = moduleLabels, 'genelist.merge' = genelist.merge))
}



#### normalize residual matrix by dividing the first singular value d1 at each stage 
Normalize.res.mtx = function(res.mtx,numreps)
{
  p  = nrow(res.mtx)
  nt = length(numreps)
  
  d1 = sapply(1:nt, function(tt) 
    (svd(res.mtx[,(1+sum(numreps[0:(tt-1)])):sum(numreps[1:tt])])$d)[1])
  
  nor.res.mtx = lapply(1:nt, function(tt)
    res.mtx[,(1+sum(numreps[0:(tt-1)])):sum(numreps[1:tt])]/d1[tt]*sqrt(p))
  
  nor.res.mtx = do.call(cbind, nor.res.mtx)
  
  return(nor.res.mtx)
}



