NetCmbn<-function(NetSet,N){
  l=length(NetSet);
  I1=matrix(0,l,N);
  I2=matrix(0,l,N);
  Neigh=list();
  for(i in 1:l){
    Neigh[[i]]<-NetSet[[i]][[1]][[1]];
        I1[i,]<-NetSet[[i]][[2]][1,];
        I2[i,]<-NetSet[[i]][[3]][1,];
         };
list(Neigh,I1,I2);
  }

  