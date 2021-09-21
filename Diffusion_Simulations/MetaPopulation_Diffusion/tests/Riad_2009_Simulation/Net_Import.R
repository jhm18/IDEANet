Net_Import<-function(File,N){
  L=scan(File);
  lL=length(L);
  i=seq(3,lL,3)
  L1=L[i-2];
  L2=L[i-1];
  L3=L[i];
 ft=NeighborhoodDataWD(N,L1,L2,L3);
 I1=matrix(ft[[2]],1,N);
 I2=matrix(ft[[3]],1,N);
 list(ft[1],I1,I2); 
}

