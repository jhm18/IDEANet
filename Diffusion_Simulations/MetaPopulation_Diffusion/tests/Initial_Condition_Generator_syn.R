Initial_Condition_Generator_syn<-function(N, method, seed_number, inducer){
  
  require(stringdist)
  
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


