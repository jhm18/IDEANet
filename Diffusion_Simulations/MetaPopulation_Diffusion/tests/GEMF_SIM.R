GEMF_SIM<-function(Para,Net,x0,maxNumevent,Runtime,N){
  M=Para[[1]];q=Para[[2]];L=Para[[3]];A_d=Para[[4]];
  A_b=Para[[5]]; Neigh=Net[[1]];I1=Net[[2]];I2=Net[[3]];
  
  
  bil<-matrix(0,L,M);
  for(i in 1:L){bil[i,]=rowSums(A_b[,,i]);};
  
  bi<-list();
  temp<-matrix(0,M,L);
  for (i in 1:M){
    for(j in 1:L){ temp[,j]<-A_b[i,,j]};
    bi[[i]]<-temp;
  };
  di<-matrix(0,M,1);
  di[,1]<-rowSums(A_d);
  
  X<-vector("integer",N);X<-x0;
  
  Rn<-vector("numeric",N);
  Nq=matrix(0,L,N);
  
  for(n in 1:N){
    for(l in 1:L){
      if(I1[l,n]!=0&X[n]==q[1,l]){     
         
        for(c in I1[l,n]:I2[l,n]){ 
          m<-Neigh[[l]][1,c];
          w<-Neigh[[l]][2,c];
          Nq[l,m]<-Nq[l,m]+w;  
          };
        
         };
    };
    };
  
  for(n in 1:N){
  Rn[n]<-sum(bil[,X[n]]*Nq[,n])+di[X[n]];};
  R=sum(Rn);
  
  
  
  
  ts<-vector("numeric",maxNumevent);
  n_index<-vector("numeric",maxNumevent);
  i_index<-vector("numeric",maxNumevent);
  j_index<-vector("numeric",maxNumevent);
  
  s=0;Tf=0;
  dum<-matrix(0,L,1);
  
  while((s<maxNumevent)&&(R>=1e-6)&&(Tf<Runtime)){
    s<-s+1; 
    ts[s]<--log(runif(1))/R;
    
    ns<-sample(1:N,1,prob=Rn);
    is<-X[ns];
    dum[,1]<-Nq[,ns];
    pf<-A_d[is,]+as.vector(bi[[is]]%*%dum);
    
    js<-sample(1:M,1,prob=pf);
    
    n_index[s]<-ns;
    j_index[s]<-js;
    i_index[s]<-is;
    
    X[ns]<-js;
    
    R<-R-Rn[ns];
    Rn[ns]<-sum(bil[,js]*Nq[,ns])+di[js];
    R<-R+Rn[ns]
    
    for(l in 1:L){
      
      
        if(q[1,l]==js&I1[l,ns]!=0){
        for(c in I1[l,ns]:I2[l,ns]){ 
         n<-Neigh[[l]][1,c];
         w<-Neigh[[l]][2,c];
                                                 
          Nq[l,n]<-Nq[l,n]+w;
          Rn[n]<-Rn[n]+bil[l,X[n]]*w;
          R<-R+bil[l,X[n]]*w;
        }
      };
      
      
      if(q[1,l]==is&I1[l,ns]!=0){
        for(c in I1[l,ns]:I2[l,ns]){ 
          n<-Neigh[[l]][1,c];
          w<-Neigh[[l]][2,c];
                    
          Nq[l,n]<-max(0,Nq[l,n]-w);
           
          Rn[n]<-max(0,Rn[n]-bil[l,X[n]]*w);
          
          R<-R-bil[l,X[n]]*w;
        
        }
      };
     
    };
    
    Tf<-Tf+ts[s];
  };
  lasteventnumber<-s;
  if((Tf<Runtime)&& (s==maxNumevent)){print("increase maxNumevent inorder to match Runtime");};
  list(ts,n_index,i_index,j_index,Tf,lasteventnumber);
}