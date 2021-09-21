samp<-function(p){






a=c(0,cumsum(p[1:length(p)-1]))/sum(p);
b=cumsum(p)/sum(p);
toss=runif(1);
k=which(a<toss & b>=toss);
k}