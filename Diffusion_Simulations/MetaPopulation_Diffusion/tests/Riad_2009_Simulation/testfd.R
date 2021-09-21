source("NeighborhoodDataWD.R");
source("Net_Import.R");
source("Post_Population.R");
source("Post_NodesState.R");
source("GEMF_SIM_several_switching_exact.R");
source("NetGen_Geo.R");
source("NetGen_ER.R");
source("Para_switching_exact.R");
source("samp.R");

#nets--------------------------------------------------------------

N=500;
Net1=Net_Import("a1wdedges.txt",N);

# Net2=Net_Import("a2wdedges.txt",N);


I1=Net1[[2]]; I2=Net1[[3]];
 Nei1=Net1[[1]][[1]];






# Nei2=rbind(Nei2a,Nei2b)


#initial conditions__________________________________________________________________________
M=4;

x0=matrix(as.integer(3),1,N); 


  #simulation terminator 
  maxNumevent=1000000;Runtime=40;
  
  Para=Para_switching_exact(delta=4,gamma_1=4,gamma_2=4,beta1=1,beta2=0,p0=0)
  

numrun=1;
trin=1;



A_d=Para[[1]];
beta1=Para[[2]];

di<-matrix(0,M,1);
di[,1]<-rowSums(A_d)
bi=rbind(c(0,0,beta1,0),c(0,0,0,beta1),c(0,0,0,0),c(0,0,0,0))
bil=c(beta1,beta1,0,0);






ts<-matrix(0,maxNumevent,numrun);
n_index<-matrix(as.integer(0),maxNumevent,numrun);
i_index<-matrix(as.integer(0),maxNumevent,numrun);
j_index<-matrix(as.integer(0),maxNumevent,numrun);
Tfvec<-vector("numeric",numrun);
lasteventnumbervec<-vector("numeric",numrun);








h1<-vector("numeric",N);

for(n in 1:N){
  
  if(I1[n]!=0 && (x0[n]==3 || x0[n]==4) ){  #if node n is in state 3 or 4 and it has neigbor in layer 1   
    
    for(c in I1[n]:I2[n]){ 
      m<-Nei1[1,c];
     
      h1[m]<-h1[m]+1;  # h1(m) is the number of neigbor of node m in layer 1 that are in state 3 or 4
    };
    
  };
  
};





Rn<-vector("numeric",N);


for(n in 1:N){
  Rn[n]<-di[x0[n]]+bil[x0[n]]*h1[n];
  
};








xvar<-vector("integer",N);

h1var<-vector("numeric",N);

Rnvar<-vector("numeric",N);





###########function f1  does the simulation



  
  xvar<-x0;
  h1var<-h1;
 
  Rnvar<-Rn;

  
  
  R<-sum(Rn);
  s<-0;
  Tf<-0;
  while((s<maxNumevent)&&(R>=1e-6)&&(Tf<Runtime)){
    s<-s+1; 
    
    tper<--log(runif(1))/R;
    
    ts[s,trin]<- tper;
   
    
    ns<-samp(Rnvar);
   
     is<-xvar[ns];
    
    pf<-A_d[is,];
    
    if (xvar[ns]==1){pf<-pf+c(0,0,beta1*h1var[ns],0)}
    if(xvar[ns]==2){pf<-pf+c(0,0,0,beta1*h1var[ns])} 
    
    js<-samp(pf)
    
    n_index[s,trin]<-ns;
    j_index[s,trin]<-js;
    i_index[s,trin]<-is;
    
    xvar[ns]<-js;
    
    R<-R-Rnvar[ns];
   Rnvar[ns]<-di[xvar[ns]]+(xvar[ns]==1 | xvar[ns]==2)*h1var[ns]*beta1;
    R<-R+Rnvar[ns];
    
    
    #update_____________________________________
    
   
        
   
    #-----------------------------------------    
if ((is==1 || is==2)&& (js==3 || js==4) && I1[ns]!=0){#s1 to I1
      
        for(ind in I1[ns] : I2[ns]){
          
          
          h1var[Nei1[1,ind]]<-h1var[Nei1[1,ind]]+1;
          
          
          if(xvar[Nei1[1,ind]]==1 || xvar[Nei1[1,ind]]==2){Rnvar[Nei1[1,ind]]<-di[xvar[Nei1[1,ind]]]+beta1*h1var[Nei1[1,ind]];
          R<-R+beta1;}
          
          
          
        }
      
      
   
  } 
      
    #-----------------------------------------    
 if ((is==3 || is==4) && (js==1 || js==2) && I1[ns]!=0){#I1 to s1
      
      
        for(c in I1[ns]:I2[ns]){
          
          
          h1var[Nei1[1,c]]<-h1var[Nei1[1,c]]-1
          
          
          if(xvar[Nei1[1,c]]==1 | xvar[Nei1[1,c]]==2){Rnvar[Nei1[1,c]]<-di[xvar[Nei1[1,c]]]+beta1*h1var[Nei1[1,c]];
          R<-R-beta1;
          }
          
          
       }
      

  
      
    } 
    
    
    
    
    
    
    
    Tf<-Tf+tper;

    
    
    
     };
  lasteventnumbervec[trin]<-s;
  Tfvec[trin]<-Tf;
  if((Tf<Runtime)&& (s==maxNumevent)){print("increase maxNumevent inorder to match Runtime");};
  lst=list(ts,n_index,i_index,j_index,Tfvec,lasteventnumbervec);








############ 





    
tsa<-lst[[1]];
n_indexa<-lst[[2]];
i_indexa<-lst[[3]];
j_indexa<-lst[[4]];
Tf<-lst[[5]];
lasteventnumbera<-lst[[6]];
    
    
 
    
    
    #calculating population of each comportmant through time
    runforplot=1;
    lst2<-Post_Population(x0,M,N,tsa[,runforplot],i_indexa[,runforplot]
                          ,j_indexa[,runforplot],lasteventnumbera[runforplot]);
    T<-lst2[[1]];
    StateCount<-lst2[[2]];
    
    I_1<-StateCount[3,];
    I_2<-StateCount[4,];
    
plot(T,I_1+I_2)



    

