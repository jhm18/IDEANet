undirNet_to_edglis<-function(Net){
  Neighv=Net[[1]][[1]][1,];
  en=length(Neighv)/2;
  I1=Net[[2]][1,];
  I2=Net[[3]][1,];
  edl=matrix(0,en,2)
  l=0;
  for(n in 1:length(I1)){
    for (i in I1[n]:I2[n]){
      if(Neighv[i]>n){
        l=l+1;
        edl[l,]=c(n,Neighv[i]);
      
      }
      
    }
  }
  
  edl;
}