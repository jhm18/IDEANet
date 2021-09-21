#Script 1: Direct Implementation of meta_Diff2.sas
#Jonathan H. Morgan
#3 August 2021

# Design Steps
# 1) Direct Implementation of meta_Diff2.sas
# 2) Incorporate Gillepsie Algorythm Scheduler
# 3) Implement Dynamic Mulit-Layer SIS Build
# 4) Implement Geospatial Visualizations (Templates Built in Earlier Versions)
# 5) Implement SIAS Model


# Clear Out Console Script
  cat("\014")

# Setting Work Directory
  setwd("Z:/workspace/Meta-Population Simulations/Data & Scripts/R_Scripts")
  getwd()

# Options
  options(stringsAsFactors = FALSE)
  options(mc.cores = parallel::detectCores())
  
################
#   PACKAGES   #
################
  library("cowplot")  #Used to transform Base R Graphics in ggplot2 objects 
  library('gmt')      #Used to plot Geo-spatial Data Using Generic Mapping Tools (Stable Open-Source Tool Libraries)
  
#################
#   FUNCTIONS   #
#################

######################################
#   IMPORTING & PLOTTING TEST DATA   #
######################################
  
  # Importing Raw Edges & Node_i for the Purposes of Comparison
    rawedges <- readr::read_csv("Z:/workspace/Meta-Population Simulations/Data & Scripts/rawedges.csv")
    rawnodes <- readr::read_csv("Z:/workspace/Meta-Population Simulations/Data & Scripts/rawnodes.csv")
    nodes_i <- readr::read_csv("Z:/workspace/Meta-Population Simulations/Data & Scripts/nodes_i.csv")
    
  # Looking at edges file, we see all pair combinations (i.e. all potentialities)
    
#############################
#   SPECIFYING SIMULATION   #
#############################
    
    
#############################
#   DIFFUSION EXPERIMENTS   #
#############################
    
  
# Single SQL run of the transmission prob bit.   *
# proc sql;
#   create table _newtrans as select
#   a.&nodeid1 as &nodeid, a.&nodeid2 as infsource, &t as t,
#   (1-((c.&n_i/c.&nodepop)*&beta))**(b.&degvar*a.&tranprob) as pNotTrans_ij, b.&n_s as &n_s
#   from &edges as a, &nodedata as b, &nodedata as c
#   where a.&tranprob>0 & b.&n_i>0 & c.&n_s>0 & a.&nodeid1=b.&nodeid & a.&nodeid2=c.&nodeid
#   order by (a.&nodeid1);
# quit;
    
# Parameters
  beta <- 0.1
  seedsize <- 4
  numinf <- 1
  nodedata <- nodes_i
  nodeid <- c('nodeid')
  n_i <- c('num_i')
  edges <- rawedges
  inftime <- 3
  days <- 50
  
  # Note: Walk through to standardize names.

  # Looping Over Time Periods
    for(i in 1:days) {
      # Iteration i
        t <- i
      
      # Randomly Identifying 1 Compartment with 4 Infections
        if (i == 1) {
          # Iteration Values
            trackt <- as.data.frame(cbind(seedsize, 0, (sum(nodedata$totalpop) - seedsize), 0))
            colnames(trackt) <- c('num_i', 'num_r', 'num_s', 't')
        
          # Infection History
            variables <- c('t', nodeid, n_i)
            infstack <- nodedata[variables]
            infstack <- infstack[(infstack[n_i] > 0), ]
            colnames(infstack)[[3]] <- 'totnumtrans'
        } else {
          # Using trackt & infection history from last iteration
            trackt <- trackt
            infstack <- infstack
        }
    
      # Constructing Transmission Table (Assuming Symmetric Ties)
        nodes_i <- edges[(edges[[1]] %in% infstack$nodeid), ]
        colnames(nodes_i)[[1]] <- c('nodeid')
        nodes_i <- dplyr::left_join(nodes_i, nodedata, by=c('nodeid'))
        
        nodes_j <- edges[(edges[[2]] %in% infstack$nodeid), ]
        colnames(nodes_j)[[1]] <- c('nodeid')
        nodes_j <- dplyr::left_join(nodes_j, nodedata, by=c('nodeid'))
        
      # Stacking & Eliminating Duplicates
        newtrans <- rbind(nodes_i, nodes_j)
        newtrans <- newtrans[!duplicated(newtrans[c(1,2)]),]
        colnames(newtrans)[[1]] <- c('nodeid1')
        rm(nodes_i, nodes_j)
        
      # Filtering Rows
        newtrans <- newtrans[(newtrans$totprob_a > 0 & newtrans$num_s > 0 & newtrans$num_i > 0), ]
        
      # Calculating Inverse Transmission Probability
        newtrans$pNotTrans_ij = (1 - (newtrans$num_i/newtrans$totalpop * beta))**(newtrans$avgdeg*newtrans$totprob_a)
        newtrans <- newtrans[c('nodeid1', 'nodeid2', 't', 'pNotTrans_ij', 'num_s')]
        newtrans$t <- t
        
      # Generating New Infections by Taking the Cumulative Joint Product
        sus_set <- sort(unique(newtrans$nodeid1))
        newinf <- vector('list', length(sus_set))
        for (j in seq_along(sus_set)){
          j_values <- newtrans[(newtrans$nodeid1 %in% sus_set[[j]]), ]
          notprod <- prod(j_values$pNotTrans_ij)
          newinf[[j]] <- cbind(notprod, j_values[nrow(j_values), ])
          newinf[[j]]$totnumtrans <- newinf[[j]]$num_s*(1-notprod)
          rm(j_values, notprod)
        }
        newinf <- as.data.frame(do.call("rbind", newinf))
        rm(sus_set)
        
      # Appending t and t+1 infections
        colnames(newinf)[[2]] <- c('nodeid')
        infstack <- rbind(infstack[,c('t', 'nodeid', 'totnumtrans')], newinf[,c('t', 'nodeid', 'totnumtrans')])
        
      # Control Loop
        if (t <= inftime) {
            # Nobody at Risk of Recovery
              nodedata <- nodedata[,c('nodeid', 'num_i', 'num_r', 'num_s', 'totalpop', 'avgdeg', 't')]
            
            # Replacing NA values with 0 Where Necessary
              newinf[c('totnumtrans')][is.na(newinf[c("totnumtrans")])] <- 0
            
            # Updating with New Values
              node_set <- newinf$nodeid
              for (j in seq_along(node_set)){
                nodedata[(nodedata$nodeid == node_set[[j]]),c('nodeid','num_s', 't')] <- newinf[(newinf$nodeid == node_set[[j]]), c('nodeid', 'num_s', 't')]
              }
            
            # Adjusting Number Infected & Number Susceptible
              for (j in seq_along(node_set)) {
                j_values <- nodedata[(nodedata$nodeid == node_set[[j]]),c('num_i', 'num_s')]
                totnumtrans_j <- newinf[(newinf$nodeid == node_set[[j]]), c('totnumtrans')]
                j_values$num_i <- j_values$num_i + totnumtrans_j
                j_values$num_s <- j_values$num_s - totnumtrans_j
                nodedata[(nodedata$nodeid == node_set[[j]]),c('num_i', 'num_s')] <- j_values
                rm(j_values, totnumtrans_j)
              }
          }else{
            # Some at Risk
              newrec <- infstack[(infstack$t == t - (inftime + 1)), ]
              colnames(newrec)[[3]] <- c('recnum')
            
            # Sub-Setting Data            
              nodedata <- nodedata[,c('nodeid', 'num_i', 'num_r', 'num_s', 'totalpop', 'avgdeg', 't')]
            
            # Replacing NA values with 0 Where Necessary
              newinf[c('totnumtrans')][is.na(newinf[c("totnumtrans")])] <- 0
              newrec[c('recnum')][is.na(newrec[c("recnum")])] <- 0
            
            # Updating with New Values
              node_set <- newinf$nodeid
              for (j in seq_along(node_set)){
                nodedata[(nodedata$nodeid == node_set[[j]]),c('nodeid','num_s', 't')] <- newinf[(newinf$nodeid == node_set[[j]]), c('nodeid', 'num_s', 't')]
              }
            
            # Adjusting Number Infected & Number Susceptible
              for (j in seq_along(node_set)) {
                j_values <- nodedata[(nodedata$nodeid == node_set[[j]]),c('num_i', 'num_s', 'num_r')]
                totnumtrans_j <- newinf[(newinf$nodeid == node_set[[j]]), c('totnumtrans')]
                recnum_j <- newrec[(newrec$nodeid == node_set[[j]]), c('recnum')][[1]]
              
                j_values$num_i <- j_values$num_i + totnumtrans_j - recnum_j
                j_values$num_s <- j_values$num_s - totnumtrans_j
                j_values$num_r <- j_values$num_r + recnum_j
                nodedata[(nodedata$nodeid == node_set[[j]]),c('num_i', 'num_s', 'num_r')] <- j_values
                rm(j_values, totnumtrans_j, recnum_j)
              }
          }
        rm(node_set)
        
      # Calculating Totals for this Iteration
        tinf <- c(sum(nodedata$num_i), sum(nodedata$num_r), sum(nodedata$num_s), t)
      
      # Appending to trackt
        trackt <- rbind(trackt, tinf)
        
      # Conditional Ending the Loop if Number of the Infections or Susceptible Agents are 0
        if(round(tinf[[1]], 1) == 0 | round(tinf[[3]], 1) == 0){
          t <- days
        }
    }
        
  # Adding Constants to Log
    trackt$beta <- beta
    trackt$inftime <- inftime
  
  # Plotting Diffusion Curves
    diff_plot <- function() {
      # Plotting Base Plot
        layout(rbind(1,2), heights=c(9.5,0.5))
        par(mar = c(4.5, 5.5, 2, 5.5),  family='HersheySerif')
        plot(0, type='n', xlab=' ', ylab=' ', xlim=c(min(trackt$t), max(trackt$t)), 
             ylim=c(min(trackt$num_s), max(trackt$num_s)), family='HersheySerif', 
             las=1, main=' ', bty='n', axes=FALSE)
        grid(lwd = 2)
        
      # Adding Traces
        lines(x=trackt$t, y=trackt$num_s, col="blue")
        
      # Adding Number Susceptible
        mtext("Number Susceptible",side=2,line=4.5, family='HersheySerif') 
        axis(2, ylim=c(min(trackt$num_s), max(trackt$num_s)),las=1, family='HersheySerif')
      
      # Layering Recovered
        par(new=TRUE)
        
      # Add Second Plot Dimensions
        plot(0, type='n', xlab=' ', ylab=' ', xlim=c(min(trackt$t), max(trackt$t)), 
             ylim=c(min(trackt$num_r), max(trackt$num_r)),cex.axis=1.3, family='HersheySerif', 
             las=1, main=' ', bty='n', axes=FALSE)
        
      # Adding Right Y Axis
        mtext("Number Recovered",side=4,line=4, family='HersheySerif') 
        axis(4, ylim=c(min(trackt$num_r), max(trackt$num_r)),las=1, family='HersheySerif')
        
      # Adding Recovered Trace
        lines(x=trackt$t, y=trackt$num_r, col='brown')
        lines(x=trackt$t, y=trackt$num_i, col='forestgreen')
        
      # Adding x-axis 
        mtext("Time",side=1, line=2.5, family='HersheySerif') 
        axis(1, ylim=c(min(trackt$t), max(trackt$t)),las=1, family='HersheySerif')
        
      # Adding Title
        title('Meta-Populations Model: Diffusion Curves', family='serif', cex.main=1.5)
        
      # Adding Legend
        par(mar=c(0, 0, 0, 0), xpd=NA)
        plot(0, type= 'n', bty='n', xlab='', ylab='', axes=FALSE)
        text(x=0.75, y = 0, "Susceptible")
        text(x=0.95, y=0, "Recovered")
        text(x=1.20, y=0, "Infected")
        segments(0.81, 0, 0.86, 0, col=c('blue'))
        segments(1, 0, 1.05, 0, col=c('brown'))
        segments(1.25, 0, 1.30, 0, col=c('forestgreen'))
    }      
    
    g <- cowplot::as_grob(diff_plot)
    p_1 <- cowplot::ggdraw(g)
    
    p_1
    
        
        
        
        
        
      
        
        
        
        
        
        
        lines(x=trackt$t, y=trackt$num_r, col="brown")
    
        
    
    
    
  
  
  
