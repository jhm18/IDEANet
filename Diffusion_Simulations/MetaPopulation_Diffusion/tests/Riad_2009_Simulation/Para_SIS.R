Para_SIS<-function(delta,beta){
  M=2;#number of compartments
  q=matrix(2,1,1);# matrix of influencer compartment for each layer which has the dimension of 1 by number of layers
  l=length(q);
  A_d=matrix(0,M,M);# node base transition matrix
  A_d[2,1]<-delta;
  A_b=array(0,c(M,M,l));#edgebase transittion array for different layers
  A_b[1,2,1]<-beta;
  list(M,q,l,A_d,A_b);
}