NetGen_ER<-function(N,p){
s=N*N;
  l=0;
L1<-vector("integer",s);
L2<-vector("integer",s);

  for (i in 1:(N-1)){
    for(j in (i+1):N){ 
          if(runif(1)<=p){
            l=l+1;
            L1[l]=i; L2[l]=j;
            l=l+1;
            L1[l]=j; L2[l]=i;
            };
    };
  };
L1<-L1[-(l+1):-s];
L2<-L2[-(l+1):-s]; 
L3<-rep(1,l);
ft=NeighborhoodDataWD(N,L1,L2,L3);
I1=matrix(ft[[2]],1,N);
I2=matrix(ft[[3]],1,N);
list(ft[1],I1,I2); 
}