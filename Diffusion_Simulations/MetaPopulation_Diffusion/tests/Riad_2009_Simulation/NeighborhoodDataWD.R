NeighborhoodDataWD<-function(N,L1,L2,L3){
 pr=sort.list(L1);
 NeighVec=L2[pr];
 NeighWei=L3[pr];
 Nodes=L1[pr];
l=length(L1);
Neigh<-matrix(0,2,l);
Neigh[1,]<-NeighVec;
Neigh[2,]<-NeighWei;
d=vector("numeric",N);
I1=vector("numeric",N);
I2=vector("numeric",N);
i=1;
while(i<l){
 node=Nodes[i];
 I1[node]=i;
     while (Nodes[i+1]==Nodes[i]){
       d[node]=d[node]+1;
       i=i+1;
       if(i==l){break;};
       };
 i=i+1; 
 };
if (i==l){node=Nodes[i];I1[node]=i;d[node]=0;}
I2=I1+d;
d=I2-I1+(I1!=0);

list(Neigh,I1,I2,d);
}