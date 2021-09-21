Post_NodesState<-function(x0,M,N,ts,n_index,j_index,lasteventnumber,timstp,Runtime){
  T<-c(0,cumsum(ts[1:lasteventnumber]));
  timnu<-floor(Runtime/timstp);
  Tr<-0:timnu*timstp;
  NodStt<-matrix(as.integer(0),timnu+1,N);
  dum<-as.integer(x0);
  mj=1;lif=0;
  for(i in 1:lasteventnumber){
    for(j in mj:(timnu+1)){
      if(T[i]<=Tr[j]&&Tr[j]<T[i+1]){NodStt[j,]<-dum;lif<-j;};
      if(T[i+1]<=Tr[j]){mj=j;break;}
    }
    dum[n_index[i]]<-j_index[i];
  };
  if(lif<(timnu+1)){for(j in (lif+1):(timnu+1)){NodStt[j,]<-dum;}};
  list(Tr,NodStt); 
}



  
