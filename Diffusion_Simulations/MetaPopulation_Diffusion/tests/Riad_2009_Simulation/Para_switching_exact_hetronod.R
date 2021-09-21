Para_switching_exact_hetronod<-function(delta_vec,gamma_1_vec,gamma_2_vec,beta1,beta2,p0){
  M=4;#number of compartments
  

  N=length(gamma_1_vec);
  
  A_d=array(0,c(M,M,N));# node base transition matrix
  A_d[1,2,]<-gamma_1_vec; A_d[3,4,]<-gamma_1_vec;
  A_d[2,1,]<-gamma_2_vec; A_d[4,3,]<-gamma_2_vec;
  A_d[3,1,]<-delta_vec; A_d[4,2,]<-delta_vec;
  

  
 
  list(A_d,beta1,beta2,p0);
}