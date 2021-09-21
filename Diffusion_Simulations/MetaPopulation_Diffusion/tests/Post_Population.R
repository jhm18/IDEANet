Post_Population<-function(x0,M,N,ts,i_index,j_index,lasteventnumber){

   StateCount=matrix(0,M,lasteventnumber+1);
 
 for(i in 1:N){StateCount[x0[i],1]<-StateCount[x0[i],1]+1};
 
 for(e in 1:lasteventnumber){
     StateCount[,e+1]<-StateCount[,e];
     StateCount[i_index[e],e+1]<-StateCount[i_index[e],e+1]-1;                              
     StateCount[j_index[e],e+1]<-StateCount[j_index[e],e+1]+1;
     };
T<-c(0,cumsum(ts[1:lasteventnumber]));
list(T,StateCount);
}