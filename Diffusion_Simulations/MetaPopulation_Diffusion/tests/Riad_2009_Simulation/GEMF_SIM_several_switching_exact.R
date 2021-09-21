GEMF_SIM_several_switching_exact<-function(Para,IN1,IN2,Nei1,Nei2,x0,maxNumevent,Runtime,N,numrun){
 
  
  M=5;
  x0=as.integer(x0);
  

  
A_d=Para[[1]];
beta1=Para[[2]];
beta2=Para[[3]];
p0=Para[[4]]  
Nei2[4,]<-p0;


dum=dim(Nei2)[2];
actvec<-vector("integer",dum)+as.integer(10);

for(n in 1:N){
  if(IN2[1,n]!=0){
 for(ind in IN2[1,n]:IN2[2,n]){
   m=Nei2[1,ind];
   p00=Nei2[4,ind];
   corind=Nei2[3,ind];
   act=((x0[n]==2||x0[n]==4)&&(x0[m]==2||x0[m]==4))&&(runif(1)<=p00);
 if (actvec[ind]==10){
   actvec[ind]<-as.integer(act);
   actvec[corind]<-as.integer(act)
   }
   
   }
  }
  
}

 
  
  
 di<-matrix(0,M,1);
  di[,1]<-rowSums(A_d)
  
  
  

  
  Rn<-vector("numeric",N);
  h1<-vector("numeric",N);
  h2<-vector("numeric",N);
  
  
  ts<-matrix(0,maxNumevent,numrun);
  n_index<-matrix(as.integer(0),maxNumevent,numrun);
  i_index<-matrix(as.integer(0),maxNumevent,numrun);
  j_index<-matrix(as.integer(0),maxNumevent,numrun);
  Tfvec<-vector("numeric",numrun);
  lasteventnumbervec<-vector("numeric",numrun);
  

  
  
 
  
  

  
  for(n in 1:N){
    
      if((x0[n]==3||x0[n]==4) & IN1[1,n]!=0){  #if node n is in state 3 or 4 and it has neigbor in layer 1   
        
        for(c in IN1[1,n]:IN1[2,n]){ 
          m<-Nei1[1,c];
          w<-Nei1[2,c];
          h1[m]<-h1[m]+w;  # h1(m) is the number of neigbor of node m in layer 1 that are in state 3 or 4
        };
        
      };
   
  };
  
  
  for(n in 1:N){
    
    if((x0[n]==4)&IN2[1,n]!=0){  #if node n is in state 4 and it has potential neigbor in layer 2  
      
      for(c in IN2[1,n]:IN2[2,n]){ 
        m<-Nei2[1,c];
        w<-Nei2[2,c];
        act<-actvec[c];
        h2[m]<-h2[m]+w*act;  # h2(m) is the number of neigbor of node m in layer 2 that are in state  4 and link is active
      };
      
    };
    
  };
  
 
  
  
  
    for(n in 1:N){
    Rn[n]<-di[x0[n]];
    
    if (x0[n]==1){Rn[n]<-Rn[n]+beta1*h1[n]}
    if(x0[n]==2){Rn[n]<-Rn[n]+beta1*h1[n]+beta2*h2[n]}
};
  
  
  
  
  
  
  
  
  X<-vector("integer",N);
  
  h1var=vector("numeric",N);
  h2var=vector("numeric",N);
  Rnvar<-vector("numeric",N);
  actvecvar<-actvec;
  pf<-vector("numeric",5);
  
  
  
  ###########function f1  does the simulation
  
  f1<-function(trin){
   
    
   X<<-x0;
   h1var<<-h1;
   h2var<<-h2;
   Rnvar<<-Rn;
   actvecvar<<-actvec;
    
  
    R<-sum(Rn);
    s<-0;
    Tf<-0;
    while((s<maxNumevent)&&(R>=1e-6)&&(Tf<Runtime)){
      s<-s+1; 
      
      tper<--log(runif(1))/R;
      
      ts[s,trin]<<- tper;
     
      ns<-sample(1:N,1,prob=Rnvar);
      is<-X[ns];
      
      pf<<-A_d[is,]
      
      if (is==1){pf<<-pf+c(0,0,beta1*h1var[ns],0,0)}
      if(is==2){pf<<-pf+c(0,0,0,beta1*h1var[ns]+beta2*h2var[ns],0)}
      
      js<-sample(1:M,1,prob=pf);
      
      n_index[s,trin]<<-ns;
      j_index[s,trin]<<-js;
      i_index[s,trin]<<-is;
      
      X[ns]<<-js;
     
      
      #update_____________________________________
      
      if (is==2 && js==1){#s2 to s1
        
        R<-max(R-Rnvar[ns],0);
        Rnvar[ns]<<-di[js]+beta1*h1var[ns];
        R<-R+Rnvar[ns];
        h2var[ns]<<-0;
        
        if (IN2[1,ns]!=0){
        for(c in IN2[1,ns]:IN2[2,ns]){
          actvecvar[c]<<-as.integer(0);
          actvecvar[Nei2[3,c]]<<-as.integer(0);
          
        }
          }
        }
        
    #-----------------------------------------    
        if (is==1 && js==2){#s1 to s2
          
          if (IN2[1,ns]!=0){
          for(c in IN2[1,ns]:IN2[2,ns]){
           
             act<-as.integer(((X[Nei2[1,c]]==2||X[Nei2[1,c]]==4)&&runif(1)<=Nei2[4,c]))
            
            actvecvar[c]<<-act;
            actvecvar[Nei2[3,c]]<<-act;
            
            #h2var[ns]<<-h2var[ns]+(X[Nei2[1,c]]==4)*act*Nei2[2,Nei2[3,c]];
            
            if(identical(numeric(0), (X[Nei2[1,c]]==4)*act*Nei2[2,Nei2[3,c]])){
              
              h2var[ns]<<-h2var[ns]+0
            }
            
            if(!identical(numeric(0), (X[Nei2[1,c]]==4)*act*Nei2[2,Nei2[3,c]])){
              
              h2var[ns]<<-h2var[ns]+ (X[Nei2[1,c]]==4)*act*Nei2[2,Nei2[3,c]];
            }
            
            
          }
          }
          
          R<-max(R-Rnvar[ns],0);
          Rnvar[ns]<<-di[js]+beta1*h1var[ns]+beta2*h2var[ns];
          R<-R+Rnvar[ns];
          
        
        }
     
    
 #-----------------------------------------    
      if (is==3 && js==4){#I1 to I2
        
        if (IN2[1,ns]!=0){
          for(c in IN2[1,ns]:IN2[2,ns]){
            
            act<-as.integer(((X[Nei2[1,c]]==2||X[Nei2[1,c]]==4)&&runif(1)<=Nei2[4,c]))
            
            actvecvar[c]<<-act;
            actvecvar[Nei2[3,c]]<<-act;
            
            #h2var[ns]<<-h2var[ns]+(X[Nei2[1,c]]==4)*act*Nei2[2,Nei2[3,c]];
            
            if(identical(numeric(0), (X[Nei2[1,c]]==4)*act*Nei2[2,Nei2[3,c]])){
              
              h2var[ns]<<-h2var[ns]+0
            }
            
            if(!identical(numeric(0), (X[Nei2[1,c]]==4)*act*Nei2[2,Nei2[3,c]])){
              
              h2var[ns]<<-h2var[ns]+ (X[Nei2[1,c]]==4)*act*Nei2[2,Nei2[3,c]];
            }
            
            
            h2var[Nei2[1,c]]<<-h2var[Nei2[1,c]]+act*Nei2[2,c];
            
            if(X[Nei2[1,c]]==2){Rnvar[Nei2[1,c]]<<-Rnvar[Nei2[1,c]]+beta2*act*Nei2[2,c];
            R<-R+beta2*act*Nei2[2,c];
            }
            
            }
        }
        
        R<-max(R-Rnvar[ns],0);
        Rnvar[ns]<<-di[js];
        R<-R+Rnvar[ns];
        
        
      }
 #-----------------------------------------    
      if (is==4 && js==3){#I2 to I1
        
        if (IN2[1,ns]!=0){
          for(c in IN2[1,ns]:IN2[2,ns]){
            
            
            h2var[Nei2[1,c]]<<-max(h2var[Nei2[1,c]]-actvecvar[c]*Nei2[2,c],0);
            
            if(X[Nei2[1,c]]==2){Rnvar[Nei2[1,c]]<<-max(Rnvar[Nei2[1,c]]-beta2*actvecvar[c]*Nei2[2,c],0);
            R<-max(R-beta2*actvecvar[c]*Nei2[2,c],0);
            };
            
            
            actvecvar[c]<<-as.integer(0);
            actvecvar[Nei2[3,c]]<<-as.integer(0);
            h2var[ns]<<-0;
            

          }
        }
        
        R<-max(R-Rnvar[ns],0);
        Rnvar[ns]<<-di[js];
        R<-R+Rnvar[ns];
        
        
      }    
      
      
      
      #-----------------------------------------    
      if (is==1 && js==3){#s1 to I1
        
        if (IN1[1,ns]!=0){
          for(c in IN1[1,ns]:IN1[2,ns]){
            
            
            h1var[Nei1[1,c]]<<-h1var[Nei1[1,c]]+Nei1[2,c];
            
            
            if(X[Nei1[1,c]]==1 || X[Nei1[1,c]]==2){Rnvar[Nei1[1,c]]<<-Rnvar[Nei1[1,c]]+beta1*Nei1[2,c];
            R<-R+beta1*Nei1[2,c];
            };
            
            
          }
        }
        
        R<-max(R-Rnvar[ns],0);
        Rnvar[ns]<<-di[js];
        R<-R+Rnvar[ns];
        
        
      }   
      
      
    #  -----------------------------------------    
      if (is==3 && js==5){#I1 to R
        
        if (IN1[1,ns]!=0){
          for(c in IN1[1,ns]:IN1[2,ns]){
            
            
            h1var[Nei1[1,c]]<<-max(h1var[Nei1[1,c]]-Nei1[2,c],0);
            
            
            if(X[Nei1[1,c]]==1 || X[Nei1[1,c]]==2){Rnvar[Nei1[1,c]]<<-max(Rnvar[Nei1[1,c]]-beta1*Nei1[2,c],0);
            R<-max(R-beta1*Nei1[2,c],0);
            }
            
            
          }
        }
        
        R<-max(R-Rnvar[ns],0);
        Rnvar[ns]<<-di[js];
        R<-R+Rnvar[ns];
        
        
      }   
      
      
      
      
      
      
      
      #-----------------------------------------    
      if (is==2 && js==4){#s2 to I2
        
        if (IN1[1,ns]!=0){
          for(c in IN1[1,ns]:IN1[2,ns]){
            
            
            h1var[Nei1[1,c]]<<-h1var[Nei1[1,c]]+Nei1[2,c];
            
            
            if(X[Nei1[1,c]]==1||X[Nei1[1,c]]==2){Rnvar[Nei1[1,c]]<<-Rnvar[Nei1[1,c]]+beta1*Nei1[2,c];
            R<-R+beta1*Nei1[2,c];
            }
            
            
          }
        }
        
        
        if (IN2[1,ns]!=0){
          for(c in IN2[1,ns]:IN2[2,ns]){
            
            
            h2var[Nei2[1,c]]<<-h2var[Nei2[1,c]]+actvecvar[c]*Nei2[2,c];
            
            
            if(X[Nei2[1,c]]==2 && actvecvar[c]==1){Rnvar[Nei2[1,c]]<<-Rnvar[Nei2[1,c]]+beta2*Nei2[2,c];
            R<-R+beta2*Nei2[2,c];
            }
            
            
          }
        }
        
        
        
        
        R<-max(R-Rnvar[ns],0);
        Rnvar[ns]<<-di[js];
        R<-R+Rnvar[ns];
        
        
      }   
      
      #-----------------------------------------    
      if (is==4 && js==5){#I2 to R
        
        if (IN1[1,ns]!=0){
          for(c in IN1[1,ns]:IN1[2,ns]){
            
            
            h1var[Nei1[1,c]]<<-max(h1var[Nei1[1,c]]-Nei1[2,c],0);
            
            
            if(X[Nei1[1,c]]==1||X[Nei1[1,c]]==2){Rnvar[Nei1[1,c]]<<-max(Rnvar[Nei1[1,c]]-beta1*Nei1[2,c],0);
            R<-max(R-beta1*Nei1[2,c],0);
            };
            
            
          }
        }
        
        
        if (IN2[1,ns]!=0){
          for(c in IN2[1,ns]:IN2[2,ns]){
            
            
            h2var[Nei2[1,c]]<<-max(h2var[Nei2[1,c]]-actvecvar[c]*Nei2[2,c],0);
            
            
            if(X[Nei2[1,c]]==2 && actvecvar[c]==1){Rnvar[Nei2[1,c]]<<-max(Rnvar[Nei2[1,c]]-beta2*Nei2[2,c],0);
            R<-max(R-beta2*Nei2[2,c],0);
            }
            actvecvar[c]<<-0;
            
          };
        }
        
        
        
        
        R<-max(R-Rnvar[ns],0);
        Rnvar[ns]<<-di[js];
        R<-R+Rnvar[ns];
        
        
      }    
      
      
     
     
   
      
      Tf<-Tf+tper;
    };
    lasteventnumbervec[trin]<<-s;
    Tfvec[trin]<<-Tf;
    if((Tf<Runtime)&& (s==maxNumevent)){print("increase maxNumevent inorder to match Runtime");};
  print(R);
  print(sum(Rnvar));
    };
  
  
  
  
  
  #######running the simulation
  
 
  for(i in 1:(numrun)){
    
   
    f1(i);
    
    print(i);};
  
  
 ############ 
  
  
  
  
  list(ts,n_index,i_index,j_index,Tfvec,lasteventnumbervec);}


