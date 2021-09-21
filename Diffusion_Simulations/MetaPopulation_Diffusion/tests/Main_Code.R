source("NeighborhoodDataWD.R");
source("Net_Import.R");
source("Post_Population.R");
source("Post_NodesState.R");
source("GEMF_SIM_several_switching_exact.R");
source("NetGen_Geo.R");
source("NetGen_ER.R");
source("Para_switching_exact.R");
source("Initial_Condition_Generator.R")

require(lubridate)
require(stringdist)
library(readxl)

#nets--------------------------------------------------------------

#N=500;
#Net1=Net_Import("a1wdedges.txt",N);

#Net2=Net_Import("a2wdedges.txt",N);

Data <- read_excel("Uganda Ebola Network.xlsx", sheet=1)


Population<- ceiling(Data[5]/1000);

N<- 11065


seed_number<-1     # How many initial infected
seed_location<-10    # Kasese district
outbreak_type<- 'single'

inducer<-4



x0 <- Initial_Condition_Generator (N, 'single', seed_location, seed_number, inducer)

Net1=Net_Import("Permanent_Layer.txt",N);

Net2=Net_Import("Temporal_Layer",N);


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
      #corind=(which(neisofneighbor == n)[[1]])-1+IN2[1,neigbor];
      corind=IN2[1,neigbor];
      Nei2b[1,ind]<-corind;
      
    }
    
  }
  
}




Nei2=rbind(Nei2a,Nei2b)
rm(v1,v2,Net1,Net2,Nei2a,Nei2b)

#initial conditions__________________________________________________________________________

M=5;


#simulation terminator 
maxNumevent=1000000;Runtime=150;


# p0=0.7;
# gam=0.1;
# bet=1.5;
numrun=10;

#Para=Para_switching_exact(delta=1,gamma_1=gam,gamma_2=gam,beta1=bet,beta2=bet,p0);

#lst<- GEMF_SIM_several_switching_exact(Para,IN1,IN2,Nei1,Nei2,x0,maxNumevent,Runtime,N,numrun); 

for(p0 in c( 0.7)){
  
  for(gam in c(0.5)){
    
    for(bet in c(0.2, 0.5, 1.7, 2.5)){
      
      
      
      Para=Para_switching_exact(delta=0.05,gamma_1=gam,gamma_2=gam,beta1=bet,beta2=bet,p0);
      
      
      #lst<- GEMF_SIM_several_switching_exact(Para,IN1,IN2,Nei1,Nei2,x0,maxNumevent,Runtime,N,numrun=2);  
      
      
      for(palrun in 1:5){
        library(parallel);
        numcor=detectCores(); 
        cl<-makeCluster(numcor)
        result<-clusterCall(cl, GEMF_SIM_several_switching_exact, Para,IN1,IN2,Nei1,Nei2,x0,maxNumevent,Runtime,N,numrun) 
        stopCluster(cl);
        
        filename=paste("C:/Users/mahbubriad/Desktop/Temporal Network paper/Temporal Network result/test1",
                       "/p0",toString(p0),"gamma",toString(gam),"bet",toString(bet),
                       "palrun",toString(palrun),".RData",sep="")
        
        
        #filename=paste("C:/Users/mahbubriad/Desktop/Temporal Network paper/Temporal Network result/test",palrun,".RData",sep="")
        
        save(result,file=filename)
        rm(result);
      };
    };
  };
};



####################################################################################


###########################################################

##             Time Regularization   for a simgle parameter, run the previous part first and then this part for time regularization  ####################



timstp=0.5;

Tr<-seq(0,Runtime,timstp);


#x0=matrix(as.integer(4),1,N);

#M=5;

EI=vector("integer",length(Tr));
ER=vector("integer",length(Tr));
count=matrix(0,  N );
Timed=matrix(0,  N );
numtry=0;

for(palrun in 1:5){
  
  filename=paste("C:/Users/mahbubriad/Desktop/Temporal Network paper/Temporal Network result/test1",
                 "/p0",toString(p0),"gamma",toString(gam),"bet",toString(bet),
                 "palrun",toString(palrun),".RData",sep="")
  load(filename)
  
  for(corr in 1:4){
    lst=result[[corr]];
    numrun=length(lst[[5]]);
    for(runn in 1:numrun){
      ts<-lst[[1]][,runn];
      n_index<-lst[[2]][,runn];
      i_index<-lst[[3]][,runn];
      j_index<-lst[[4]][,runn];
      Tf<-lst[[5]][runn];
      lasteventnumber<-lst[[6]][runn]; 
      
      
      lst3<-Post_NodesState(x0,M,N,ts,n_index,j_index,lasteventnumber,timstp,Runtime)
      nodstt=lst3[[2]];
      
      for (i in 1:N){
        if(length(which(nodstt[,i]>=3))){
          count[i]=count[i]+1;
          Timed[i]=Timed[i]+ which(nodstt[,i]>=3)[1];
        }
      }
      
      Inf1=rowSums(nodstt==4);
      Inf2=rowSums(nodstt==3);
      Rec=rowSums(nodstt==5);
      
      Infected=Inf1+Inf2;
      
      EI=EI+Infected;
      ER=ER+Rec
      numtry=numtry+1;
    };  
  };
  
};

EI=EI/numtry;
ER=ER/numtry;

Prob<- count/numtry;     # Probability of nodes being infected
Nodes<-1:N;
Speed<-(Timed)/numtry;   # Average speed of the virus spread to specific node


par(mfrow=c(4, 1))
plot(Tr,ER); 
plot(Tr,EI);
plot(Nodes, Prob);

plot(Nodes, Speed)

filename2=paste("C:/Users/mahbubriad/Desktop/Temporal Network paper/Temporal Network result/test1",
                "/p0",toString(p0),"gamma",toString(gam),"bet",toString(bet),
                "finalresult",".RData",sep="")

dumb=list(Tr,EI, ER, Prob, Speed);
save(dumb,file=filename2)

# };


