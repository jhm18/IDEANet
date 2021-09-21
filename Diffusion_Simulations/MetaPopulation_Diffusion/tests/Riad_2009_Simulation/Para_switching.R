Para_switching<-function(delta,gamma_1,gamma_2,bet,betprim){
  M=4;#number of compartments
  l=3;#number of layers
  q=matrix(0,1,3);# matrix of influencer compartment for each layer which has the dimension of 1 by number of layers
  q[1]=3;q[2]=4;q[3]=4;
  
  
  A_d=matrix(0,M,M);# node base transition matrix
  A_d[1,2]<-gamma_1; A_d[3,4]<-gamma_1;
  A_d[2,1]<-gamma_2; A_d[4,3]<-gamma_2;
  A_d[3,1]<-delta; A_d[4,2]<-delta;
  
  A_b=array(0,c(M,M,l));#edgebase transittion array for different layers
  
  A_b[1,3,1]<-bet;
  A_b[2,4,1]<-bet;
  
  A_b[1,3,2]<-bet;
  A_b[2,4,2]<-bet;
  
 
  A_b[2,4,3]<-betprim;
  
 
  list(M,q,l,A_d,A_b);
}