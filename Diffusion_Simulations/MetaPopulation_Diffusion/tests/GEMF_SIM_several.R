GEMF_SIM_several<-function(Para,Net,x0,maxNumevent,Runtime,N,numrun){
  M=Para[[1]];q=Para[[2]];L=Para[[3]];A_d=Para[[4]];
  A_b=Para[[5]]; Neigh=Net[[1]];I1=Net[[2]];I2=Net[[3]];
  x0=as.integer(x0);
  
  

  
  
  bil<-matrix(0,L,M);
  for(i in 1:L){bil[i,]=rowSums(A_b[,,i]);};
  
  
  bi<-list();
  temp<-matrix(0,M,L);
  for (i in 1:M){
    for(j in 1:L){ temp[,j]<-A_b[i,,j]};
    bi[[i]]<-temp;
  };
  di<-matrix(0,M,1);
  di[,1]<-rowSums(A_d)
  
  
  
  Nq=matrix(0,L,N);
  
  Rn<-vector("numeric",N);
  
  
  ts<-matrix(0,maxNumevent,numrun);
  n_index<-matrix(as.integer(0),maxNumevent,numrun);
  i_index<-matrix(as.integer(0),maxNumevent,numrun);
  j_index<-matrix(as.integer(0),maxNumevent,numrun);
  Tfvec<-vector("numeric",numrun);
  lasteventnumbervec<-vector("numeric",numrun);
  

  
  
 
  
  
  for(n in 1:N){
    for(l in 1:L){
      if(I1[l,n]!=0&x0[n]==q[1,l]){  #if node n is in influenecr state of layer l and it has some neigbor in layer l   
        
        for(c in I1[l,n]:I2[l,n]){ 
          m<-Neigh[[l]][1,c];
          w<-Neigh[[l]][2,c];
          Nq[l,m]<-Nq[l,m]+w;  # nq(l,m) is the number of neigbor of node m in layer l that are in influencer state of layer l
        };
        
      };
    };
  };
  
  for(n in 1:N){
    Rn[n]<-sum(bil[,x0[n]]*Nq[,n])+di[x0[n]];};
  
  
  X<-vector("integer",N);
  
  Nqvar=matrix(0,L,N);
  
  Rnvar<-vector("numeric",N);
  
  pf<-vector("numeric",M);
  
  
  
  ######################################function f1  does the simulation
  
  f1<-function(trin){
   
    
   X<<-x0;
   Nqvar<<-Nq;
   Rnvar<<-Rn;
  
    
  
    R<-sum(Rn);
    s<-0;
    dum<-matrix(as.integer(0),L,1);
    
    Tf<-0;
    while((s<maxNumevent)&&(R>=1e-6)&&(Tf<Runtime)){
      s<-s+1; 
      
      tper<--log(runif(1))/R;
      
      ts[s,trin]<<-tper;
     
      ns<-sample(1:N,1,prob=Rnvar);
      is<-X[ns];
      dum[,1]<-Nqvar[,ns];
      pf<<-A_d[is,]+as.vector(bi[[is]]%*%dum);
      
      js<-sample(1:M,1,prob=pf);
      
      n_index[s,trin]<<-ns;
      j_index[s,trin]<<-js;
      i_index[s,trin]<<-is;
      
      X[ns]<<-js;
     
     R<-R-Rnvar[ns];
      Rnvar[ns]<<-sum(bil[,js]*Nqvar[,ns])+di[js];
      R<-R+Rnvar[ns]
      
    
      for(l in 1:L){
        
        
        if(q[1,l]==js&I1[l,ns]!=0){
          for(c in I1[l,ns]:I2[l,ns]){
            n<-Neigh[[l]][1,c];
            w<-Neigh[[l]][2,c];
            Nqvar[l,n]<<-Nqvar[l,n]+w;
            Rnvar[n]<<-Rnvar[n]+bil[l,X[n]]*w;
            R<-R+bil[l,X[n]]*w;
          }
          
        };
        
        if(q[1,l]==is&I1[l,ns]!=0){
          for(c in I1[l,ns]:I2[l,ns]){
            n<-Neigh[[l]][1,c];
            w<-Neigh[[l]][2,c];
                    
               Nqvar[l,n]<<-max(Nqvar[l,n]-w,0);
            
              Rnvar[n]<<-max(Rnvar[n]-bil[l,X[n]]*w,0);
           
              R<-R-bil[l,X[n]]*w;
                
          }
          
      
        };
        
      };
      
      Tf<-Tf+tper;
    };
    lasteventnumbervec[trin]<<-s;
    Tfvec[trin]<<-Tf;
    if((Tf<Runtime)&& (s==maxNumevent)){print("increase maxNumevent inorder to match Runtime");};
  };
  
  
  
  
  
  #######running the simulation
  
 
  for(i in 1:(numrun)){
    
   
    f1(i);
    
    print(i);};
  
  
 ############ 
  
  
  
  
  list(ts,n_index,i_index,j_index,Tfvec,lasteventnumbervec);}