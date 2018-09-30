source("C:/Users/Yutong Lai/Dropbox/TempCoexNet/draft_doNOTdelete/code/functionsNew.R")

#### Example
#### The example data set is a 2000 by 75 matrix of 2000 genes, 5 time points and 15 replications for each time point with 5 clusters 
load("C:/Users/Yutong Lai/Dropbox/TempCoexNet/draft_doNOTdelete/code/Data_t5_c5_r15.rdata")
tvec    = c( 1, 2, 3, 4, 5) # set the time-lapses from the beginning
numreps = c(15,15,15,15,15) # set the number of replication for each time point


#### Returen correlation list for each time point
correlationMatrix = NULL
correlationMatrix = TGCN(tvec, numreps, Data_t5_c5_r15)


#### (Opotional) Return cluster list for each time point
newdata = as.data.frame(Data_t5_c5_r15)
modres  = TGCN_Cluster(correlationMatrix, newdata)
