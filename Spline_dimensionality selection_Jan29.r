set.seed(123)
####
#setwd("~/Research/wb_application/topics")
setwd("C:/Users/CZhao/Dropbox/Research/PairCorrelationFunction/data application/wb_application/topics")
ptm<-proc.time()
#library(profvis)
library(dplyr)
library(MASS)
library(doSNOW)
library(foreach)
library(orthogonalsplinebasis)
library(doParallel)
library(data.table)
workpath = getwd()
# specify the parameters won't change during this computation
cpu_num =4
TT=30
R=0.12
# small value tolerance level
small_value = 10^-5
#
# n.h = 5
# h.min = 0.008
h.max = 0.016
n.L = 5
L.min = 4
L.max = 20
L.list = seq(L.min,L.max,by= 4)
R.max = R+h.max
#h.list = seq(h.min,h.max,by= 0.002)
# -----------------------------------------------------------------------------------------------#
# load data 
load(paste0(workpath,"/data/users_process.RData")) #?????????
# random select 10 % users
sample_10 = sample(1:length(process),length(process)*0.10)
processN = process
process = process[sample_10]
#load  data into r 
N  =length(process)
#ptm<-proc.time()
process = sapply(process, function(x) unique(x))
# calculate distance #######
#ptm<-proc.time()
event_num_N = sapply(process,function(x) length(unlist(x)))
event_process  = unlist(process)
l_event_process = length(event_process)
# set up parallel backends
cl = makeCluster(cpu_num, type = "SOCK")
registerDoSNOW(cl)
## generate a series of spline basis  
for (l in 1:n.L){
  L=L.list[l]
  knots<-expand.knots(seq(0,R,length.out = L-2)) # L basis function
  assign(paste0("basi_",L),OrthogonalSplineBasis(knots, order=4, keep.duplicates=FALSE))
  save(list = c(paste0("basi_",L)),file = paste0("basi_",L,".RData")) # make sure list equal to a string vector 
}


##
loc_all = 1:l_event_process
sum_N_i_vec= sapply(1:N, function(x) sum(event_num_N[1:x]))
start_vec = sum_N_i_vec-event_num_N+1
# calculate distance for different process
#ptm<-proc.time()
start_vec = c(start_vec,l_event_process+1) # add one more value to cheat for last process point js
#################################################
ptm<-proc.time()
distVDR_Cova = foreach(j = 1:l_event_process,.packages = c("orthogonalsplinebasis")) %dopar% {
  N_i =  which.max(j <start_vec)-1
  sum_N_i = sum_N_i_vec[N_i]
  start   = start_vec[N_i]
  loc_N_i =start:sum_N_i
  start_point = event_process[j]
  # different process distance
  event_ij =loc_all[-(start:sum_N_i)]
  events_ij_dp = abs(start_point-event_process[-(start:sum_N_i)])
  retain_logi_dp= events_ij_dp<=R
  dist_retain_dp=events_ij_dp[retain_logi_dp] 
  start_pro = rep(j,length(dist_retain_dp))
  end_pro = event_ij[retain_logi_dp]
  A_l=NULL
  dev  = 1/(TT-dist_retain_dp)
  for(l in 1:n.L){
    L=L.list[l]
    knots<-expand.knots(seq(0,R,length.out = L-2)) # L basis function
    # orthogonalize the basi functions 
    basi<-OrthogonalSplineBasis(knots, order=4, keep.duplicates=FALSE)
    psi.mat = evaluate(basi,dist_retain_dp)
    A_l[[l]] =t(psi.mat) %*% diag(dev) %*% (psi.mat)
  }
  #print(paste0("calculating with j = ",j))
  list(dist_retain_dp,start_pro,end_pro,A_l)
}
proc.time()-ptm
##
stopCluster(cl)
distVD = unlist(sapply(distVDR_Cova, function(x)x[[1]]))
covaVD_1 =unlist( sapply(distVDR_Cova, function(x)x[[2]]))
covaVD_2 = unlist(sapply(distVDR_Cova, function(x)x[[3]]))
matA_list = lapply(distVDR_Cova, function(x)x[[4]])
matA_l_list = NULL
for(l in 1:n.L){
  L=L.list[l]
  matA= matrix(rep(0,L*L),nrow = L)
  for(j in 1:l_event_process){
    matA = matA_list[[j]][[l]]+matA
  }
  matA_l_list[[l]]=matA/2
} # return matA_l_list: A matrix list in different dimension L
## remove duplicated cases
distVD = distVD[unique(duplicated(distVD))]
covaVD_1 = covaVD_1[unique(duplicated(distVD))]
covaVD_2 = covaVD_2[unique(duplicated(distVD))]
dataVD = data.table(distVD,covaVD_1,covaVD_2)
#matA = Reduce('+',matA_list)/2
# reorder the whole table
dataVD = dataVD[order(distVD)] # order distvd 
distVD = dataVD$distVD
covaVD_1 = dataVD$covaVD_1 # start point for dp distance  
covaVD_2 = dataVD$covaVD_2 # end point for dp distance
#n_event = unlist(lapply(1:N, function(x) rep(x,event_num_N[x])))
dist_dp_matrix = data.table(covaVD_1,covaVD_2,distVD)
colnames(dist_dp_matrix) = c("start","end","distance")
dp_len =length(distVD)

rm(dataVD,distVDR_Cova,covaVD_1,covaVD_2,distVD)
gc()
### distance and position from same process
cl = makeCluster(cpu_num, type = "SOCK")
registerDoSNOW(cl)
distVS_cal=foreach(i = 1:N, .combine='rbind',.packages = "data.table") %dopar% {
  obser=process[[i]]
  distSP_i =dist(obser, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  distM=as.matrix(distSP_i)
  distM[!lower.tri(distM)] =0
  position_i =which(as.matrix(distM)!=0,arr.ind = T)
  colnames(position_i)=c("start","end")
  dataVS = data.table(as.vector(distSP_i),position_i)
  dataVS = dataVS[V1<=R,]
  dataVS=dataVS[order(V1)]
  dataVS
} # return a data frame with dim l_sp*3
colnames(distVS_cal)= c("distance","start","end")
dist_sp_matrix = distVS_cal[,c(2,3,1)]
sp_len =dim(dist_sp_matrix)[1]
stopCluster(cl)
### calculate b vector in different dimension 
b_list=NULL
for(l in 1:n.L){
  L=L.list[l]
  knots<-expand.knots(seq(0,R,length.out = L-2)) # L basis function
  basi<-OrthogonalSplineBasis(knots, order=4, keep.duplicates=FALSE)
  # orthogonalize the basi functions 
  vecb_mat =  evaluate(basi,dist_sp_matrix$distance)
  for(k in 1:L){
    vecb_mat[,k]=vecb_mat[,k]/((TT-dist_sp_matrix$distance))
  }
  b_list[[l]]=colSums(vecb_mat)
}

rm(distVS_cal)
gc()
####################################################################################################
cl= makeCluster(cpu_num)
registerDoParallel(cl)
###-----------------------lscv ------------------------------------------###
lscv.vec = rep(0,n.L)
#----------------------- M_1(L) --------------------------------------------#
ptm<-proc.time()
M_1_hat = foreach(j = 1:dp_len,.combine = '+',.packages=c('data.table','MASS','dplyr','orthogonalsplinebasis')) %dopar% {
#profvis({
for(j in 1:500){
  position_uv= as.numeric(dist_dp_matrix[j,1:2])
  diff_uv =as.numeric(dist_dp_matrix[j,3])
  #leave pairs within same process
  LeavePair_1= dist_sp_matrix[start==position_uv[1] | end== position_uv[2]][,distance]
  #leave pairs between different process
  LeavePair_2= dist_dp_matrix[start==position_uv[1] | end== position_uv[2]][,distance]
  M_1_uv  = rep(0,n.L)
  # calculate the b vector which calculated with pairs being left-out \
  dev1 = TT-LeavePair_1
  dev2  = TT-LeavePair_2
  for (l in 1:n.L){
    L=L.list[l]
    knots<-expand.knots(seq(0,R,length.out = L-2)) # L basis function
    basi<-OrthogonalSplineBasis(knots, order=4, keep.duplicates=FALSE)
    #
    vecA = evaluate(basi,LeavePair_2)
    if(length(LeavePair_1)>0){
      ## basi function
      vecb_mat =  evaluate(basi,LeavePair_1)
      for(ll in 1:L){
        vecA[,ll]=vecA[,ll]/dev2
        vecb_mat[,ll]=vecb_mat[,ll]/dev1
      }
      A_out = t(vecA) %*%(vecA)/(N-1)
      b_out=colSums(vecb_mat)
    }else {
        for(ll in 1:L){
          vecA[,ll]=vecA[,ll]/dev2
        }
        A_out = t(vecA) %*%(vecA)/(N-1) 
        b_out = rep(0,L)}
    #
    A=matA_l_list[[l]]/(N-1) -A_out
    b=b_list[[l]] -b_out
    #
    theta.est = ginv(A, tol = sqrt(.Machine$double.eps)) %*% as.numeric(b)
    diff_uv=max(diff_uv,small_value)
    pcf_uv_est = evaluate(basi,diff_uv)%*% theta.est
    # weight function 
    weight = 1/(diff_uv*(TT-diff_uv)) ##????????
    # calculate M_1_hat function
    M_1_uv[l]= weight*(pcf_uv_est)^2  
    }
  M_1_uv 
}
proc.time() -ptm

#})
#------------------------ M_2(h) --------------------------------------------#
# i = 1:sp_len
ptm  <-proc.time() 
M_2_hat = foreach(i = 1:sp_len,.combine = '+',.packages=c('data.table','MASS','dplyr','orthogonalsplinebasis')) %dopar% {
  position_uv= as.numeric(dist_sp_matrix[i,1:2])
  diff_uv =as.numeric(dist_sp_matrix[i,3])
  #leave pairs within same process
  LeavePair_1= dist_sp_matrix[start==position_uv[1] | end== position_uv[2]][,distance]
  #leave pairs between different process
  LeavePair_2= dist_dp_matrix[start==position_uv[1] | end== position_uv[2]][,distance]
  M_2_uv  = rep(0,n.L)
  # calculate the b vector which calculated with pairs being left-out \
  dev1 = TT-LeavePair_1
  dev2  = TT-LeavePair_2
  for (l in 1:n.L){
    L=L.list[l]
    knots<-expand.knots(seq(0,R,length.out = L-2)) # L basis function
    basi<-OrthogonalSplineBasis(knots, order=4, keep.duplicates=FALSE)
    # calculate the leave one pair out matrix A and vector b
    vecb_mat =  evaluate(basi,LeavePair_1)
    if(length(LeavePair_2)>0){
      ## basi function
      vecA = evaluate(basi,LeavePair_2)
      for(ll in 1:L){
        vecA[,ll]=vecA[,ll]/dev2
        vecb_mat[,ll]=vecb_mat[,ll]/dev1
      }
      A_out = t(vecA) %*%(vecA)/(N-1)
      b_out=colSums(vecb_mat)
    }else {
      for(ll in 1:L){
        vecb_mat[,ll]=vecb_mat[,ll]/dev1
      }
      A_out = matrix(rep(0,L*L),nrow=L)
      b_out=colSums(vecb_mat)}
    # use original A matrix and B matrix minus the leave one pair out A and b 
    A=matA_l_list[[l]]/(N-1) -A_out
    b=b_list[[l]] -b_out
    # get the estimate of theta and the pcf at ditance which is the leave one pair distance u-v
    theta.est = ginv(A, tol = sqrt(.Machine$double.eps)) %*% as.numeric(b)
    diff_uv=max(diff_uv,small_value)
    pcf_uv_est = evaluate(basi,diff_uv)%*% theta.est
    # weight function 
    weight = 1/(diff_uv*(TT-diff_uv))
    # calculate M_1_hat function
    M_2_uv[l]= weight*(pcf_uv_est)
    }
  M_2_uv
}
proc.time()-ptm
#################
stopCluster(cl)
# calculate LSCV in total 
M_L = M_1_hat/(N*(N-1))-(2*M_2_hat)/N
# return M_h
lscv.vec=M_L
  
rm(M_1_hat,M_2_hat,RetainPair.1,RetainPair.2,RetainPair.1.temp,RetainPair.2.temp,RetainPair.1.logi,RetainPair.2.logi,position_sp,position_dp,dist_matrix)
gc()

###########################################################################
###########  Minimize the LSCV to get the best dimension L ################
L.choice = L.list[which.min(lscv.vec)]
out.sim = list(L.choice,L.list,lscv.vec)
save(out.sim,file=paste0("spline dimension select for data application.RData"))
#proc.time()-ptm

proc.time()-ptm