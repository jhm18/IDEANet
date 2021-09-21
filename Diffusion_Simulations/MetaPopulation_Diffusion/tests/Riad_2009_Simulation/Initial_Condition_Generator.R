Initial_Condition_Generator<-function(N, method, seed_location, seed_number, inducer){
  
  x0=matrix(as.integer(1),1,N); 
  
  
  if (amatch(method, 'single')){
    
    if (seed_location==1)
      
    {samp<- sample(seq(0, cumsum(Population)[seed_location,]), seed_number);
    x0[samp]<-inducer;
    }
    
    
    if (seed_location != 1)
      
    {samp<- sample(seq(cumsum(Population)[seed_location-1,], cumsum(Population)[seed_location,]), seed_number);
    x0[samp]<-inducer;
    }
    
  }
  

  
  
  return(x0);
  
}


