Para_SAIS_2layer<-function(delta,beta,beta_a,kappa,mu){
  M=3;#number of compartments
  q=matrix(2,1,2);# matrix of influencer compartment for each layer which has the dimension of 1 by number of layers
  l=length(q);
  A_d=matrix(0,M,M);# node base transition matrix
  A_d[2,1]<-delta;
  A_b=array(0,c(M,M,l));#edgebase transittion array for different layers
  A_b[1,2,1]<-beta;
  A_b[1,3,1]<-kappa;
  A_b[3,2,1]<-beta_a;
  A_b[1,3,2]<-mu;
  list(M,q,l,A_d,A_b);
}