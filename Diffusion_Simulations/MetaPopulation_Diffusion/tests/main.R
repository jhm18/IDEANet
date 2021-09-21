source("NeighborhoodDataWD.R");
source("Net_Import.R");
source("Post_Population.R");

source("GEMF_SIM_several_switching_exact.R");

source("Para_switching_exact.R");



#nets--------------------------------------------------------------

N=500;
Net1=Net_Import("a1wdedges.txt",N);

Net2=Net_Import("a2wdedges.txt",N);


v1=Net1[[2]]; v2=Net1[[3]];
IN1=rbind(v1,v2); Nei1=Net1[[1]][[1]];

v1=Net2[[2]]; v2=Net2[[3]];
IN2=rbind(v1,v2); Nei2a=Net2[[1]][[1]];

Nei2b=Nei2a*0;
for (n in 1:N){
  if(IN2[1,n]!=0){
    for(ind in IN2[1,n]:IN2[2,n]){
      neigbor<-Nei2a[1,ind];
      
      neisofneighbor<-Nei2a[1,IN2[1,neigbor]:IN2[2,neigbor]];
      corind=(which(neisofneighbor == n)[[1]])-1+IN2[1,neigbor];
      Nei2b[1,ind]<-corind;
    }
    
  }
  
}




Nei2=rbind(Nei2a,Nei2b)
rm(v1,v2,Net1,Net2,Nei2a,Nei2b)

#initial conditions__________________________________________________________________________
M=4;

x0=matrix(as.integer(4),1,N); 

#simulation terminator 
maxNumevent=1000000;Runtime=150;

p0=0.5 
gam=0.1 
bet=0.4 
      
      
      
      Para=Para_switching_exact(delta=1,gamma_1=gam,gamma_2=gam,beta1=bet,beta2=bet,p0);
      
      
     lst<- GEMF_SIM_several_switching_exact(Para,IN1,IN2,Nei1,Nei2,x0,maxNumevent,Runtime,N,numrun=2);  
      
      
      




ts<-lst[[1]];
n_index<-lst[[2]];
i_index<-lst[[3]];
j_index<-lst[[4]];
Tf<-lst[[5]];
lasteventnumber<-lst[[6]];


#calculating population of each comportmant through time
runforplot=1;
lst2<-Post_Population(x0,M,N,ts[,runforplot],i_index[,runforplot]
                      ,j_index[,runforplot],lasteventnumber[runforplot]);
T<-lst2[[1]];
StateCount<-lst2[[2]];

I_1<-StateCount[3,];
I_2<-StateCount[4,];
plot(T,I_1+I_2);



