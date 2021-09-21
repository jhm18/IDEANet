

/* Meta-population diffusion model macro.
	This macro produces a cell-based/macro-level diffusion simulation for use on agregate 
	poopulations.  The idea is that given some probability of contact between units & known
	infection process we can simulate the growth & contagion across settings.

	Input is a node list that contains information on the population characteristics of each
     place and a weighted edgelist, where the edgeweights represent the probability that a 
     person from place i would interact with a person from place j.

	Revision of the safegraph_... files.

	Author: Moody
	Date: 7.12.2021
	*/
libname in1 'Z:\workspace\Meta-Population Simulations\Data & Scripts';
options compress=yes;

/* load the macro file */
%include "Z:\workspace\Meta-Population Simulations\Data & Scripts\SAS_Scripts\metadiff2.mac";


/*get a node-level dataset with just info needed for transmission, assume we'll use
  attributes and such to set these values in real sim */
/* read in teh raw datafiles */

DATA rawnodes;
 keep nodeid totalpop avgdeg t pct_black;
 set in1.durham_geo_base;
  avgdeg=6;  /* assume number of close contacts per person in block, constant for testing */
RUN;
 
data rawedges;
 set in1.Durham_geo_edgeset;
  rename imlvalue=totprob_a; /* not sure why...just thought "imlvalue" was silly */
run;


/* now pull an instance to simulate on.  Here the key issue is picking a seed set.  
     In this case I hardcode it, will show how to randomly select in seperate code */

  data nodes_i;
   set rawnodes;
   if nodeid in(1,10,15,20,50) then do; /* picking seeds */
    num_i=1; /* just test seeds - will likely want to randomize this. */
   end;
   else do;
    num_i=0;
   end;
   num_r=0; 
   num_s=totalpop-num_i;
   t=0;
  run;

%MetaDiff(nodedata=nodes_i,
				nodeid=nodeid,
			   	nodepop=totalpop,
			   	n_i=num_i,
				n_s=num_s,
				n_r=num_r,
			    degvar=avgdeg,
				days=50,
			    inftime=3,
				beta=0.1,
			    edges=work.rawedges,
				nodeid1=nodeid1, 
				nodeid2=nodeid2,
				tranprob=totprob_a);

/* history of the total epidemic */
proc sgplot data=_trackt;
 series x=t y=num_i;
 series x=t y=num_r;
 series x=t y=num_s /y2axis;
run;

/* history of new transmissions by ID */
proc sgplot data=_infstack;
 where nodeid in(29, 13, 35, 132); /* just plot a handful */
 series x=t y=totnumtrans /group=nodeid;
run;


