# TGCnA
Temporal Gene Coexpression Network Analysis (TGCnA) framework for the
transcriptomic time-course data. doi: https://doi.org/10.1101/359612

The mathematical nature of TGCnA is the joint modeling of multiple
covariance matrices across time points using a “low-rank plus sparse" framework, in which the network
similarity across time points is explicitly modeled in the low-rank component. Tests on both simulation
data and a real transcriptomic data set showed that TGCnA improved the covariance estimation loss and
identiﬁed more robust and interpretable gene modules.

- input: data matrix, tvec, numreps
- output: estimated correlation matrix at each time point, clustering results at each time by using WGCNA

#### Example
#### The example data set is a 2000 by 75 matrix of 2000 genes, 5 time points and 15 replications for each time point with 5 clusters 
#load("Data_t5_c5_r15.rdata")
### set the time point vector and number of replicates at each time point
tvec    = c( 1, 2, 3, 4, 5) # set the time-lapses from the beginning
numreps = c(15,15,15,15,15) # set the number of replication for each time point


#### Return estimated correlation list results for each time point
correlationMatrix = NULL
correlationMatrix = TGCN(tvec, numreps, Data_t5_c5_r15)


#### (Opotional) Return clustering results list for each time point
newdata = as.data.frame(Data_t5_c5_r15)
modres  = TGCN_Cluster(correlationMatrix, newdata)
