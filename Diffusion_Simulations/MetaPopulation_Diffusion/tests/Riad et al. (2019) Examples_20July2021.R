# Riad et al. (2019) Model Examples
# Jonathan H. Morgan based on Examples 1 and 2 in: https://www.nature.com/articles/s41598-019-52501-1#data-availability
# 20 July 2021

# Clear Out Console Script
  cat("\014")

# Setting Work Directory
  setwd("~/Dropbox/My Mac (Jonathan’s MacBook Pro)/Desktop/DNAC/IDEANet/Data_Scripts/Riad et al, 2019 Ebola Transmission Model")
  getwd()

# Options
  options(stringsAsFactors = FALSE)
  numcor <- 4

################
#   PACKAGES   #
################
  library(parallel)
  
#################
#   FUNCTIONS   #
#################
  `%notin%` <- Negate(`%in%`)
  
# SIS Epidemic Over 1 Layer Example
  source("NeighborhoodDataWD.R")
  source("Net_Import.R")
  source("Para_SIS.R")
  source("Post_Population.R")
  source("GEMF_SIM.R")
  
# Susceptible-Alert-Infected-Susceptible Over 2 Layer Example
  source("NeighborhoodDataWD.R")
  source("NetCmbn.R")
  source("Para_SAIS_2layer.R")
  source("GEMF_SIM_Prob.R")
  source("NetGen_Geo.R")
  source("NetGen_ER.R")

##################################
#   CONSTRUCTING DATA ELEMENTS   #
##################################
  
  # Nodeslist
    nodes <- seq(1, 379, 1)
  
  # Creating All Possible Pairs
    nodes_i <- nodes
    nodes_j <- nodes
    pairs <- expand.grid(nodes_i, nodes_j)
    pairs <- pairs[(pairs$Var1 != pairs$Var2), ]
    
  # Eliminating Edges to Make into a Plausible Network
    index <- seq(1, nrow(pairs), 380)
    pairs <- cbind(seq(1, nrow(pairs), 1), pairs)
    names(pairs) <- c('Obs_ID', 'i', 'j')
    base <- pairs[pairs$Obs_ID %in% index, ]
    length(sort(unique(c(base$i, base$j))))
    
  # Adding a Few Additional Edges for Realism
    chains <- sort(sample(base$j, 100))
    extra <- pairs[pairs$Obs_ID %notin% base$Obs_ID, ]
    candidates <- extra[extra$i %in% chains, ]
    c_index <- sample(candidates$Obs_ID, 200)
    candidates <- candidates[(candidates$Obs_ID %in% c_index), ]
    candidates <- candidates[order(candidates$i, candidates$j), ]
    edges <- rbind(base, candidates)

  # Adding Weights & Tidying
    edges <- edges[c(2,3)]
    edges <- edges[order(edges$i, edges$j), ]
    edges$weight <- 1
    
  # Visualizing to Confirm Test Network
    g <- network::network.initialize(length(nodes), directed =TRUE)
    
  # Adding Edges to the Network
    el <- edges[c(1,2)]
    el[] <- lapply(el, as.character)
    g <- network::add.edges(g, el[[1]],el[[2]])
    
  # Adding Weights to the Network
    network::set.edge.value(g,"weight", edges$weight)
    
  # Visualizing to Confirm
    plot(g)
    
  # Saving to Run in Simulation
    readr::write_delim(edges,file="edgewd.txt", delim = " " )
    rm(base, candidates, el, extra, g, pairs, c_index, chains, index, nodes_i, nodes_j)
    
#########################################
#   SIS Epidemic over 1-Layer Network   #
#########################################
    
# Parameters
  File <- "edgewd.txt"
  N <- 379 
  Net <- Net_Import(File, N)  # Reading-In Network to Generate Simulation Network Object
  Para=Para_SIS(1,0.2)        # Setting-Up an SIS Epidemic Model
  
  x0 <- matrix(2,1,N)         # Generating an initial condition where all the nodes are initially infected
  i1 <- sort(sample(seq_len(length(x0)), 100, replace = FALSE))
  x0_non <- rep(1, N)
  x0[i1] <- x0_non[i1]
  
  maxNumevent=35000           # Values that Specifies when the Simulation Terminates
  Runtime=30
  
# Now the input Arguments Are Defined and We Can Run GMF_SIM
  lst <- GEMF_SIM(Para, Net, x0, maxNumevent, Runtime, N)
  
# Gathering Elements for Analysis
  ts <- lst[[1]]
  n_index <- lst[[2]]
  i_index <- lst[[3]]
  j_index <- lst[[4]]
  Tf <- lst[[5]]
  lasteventnumber <- lst[[6]]
  
# Using the Post Population Function So We can Find the Population of the Compartments Over Time
  M <- Para[[1]]
  lst2 <- Post_Population(x0, M, N, ts, i_index, j_index, lasteventnumber)
  
# Plotting Elements
  T <- lst2[[1]];
  StateCount <- lst2[[2]]
  infectedpopulation <- StateCount[2,]
  susceptiblepopulation <- StateCount[1,]
  
  x_axis <- pretty(T)
  y_axis <- pretty(c(infectedpopulation, susceptiblepopulation))
  
# Plotting
  sis_plot <- function() {
    # Plotting Base Plot
      layout(rbind(1,2), heights=c(9.5,0.5))
      par(mar = c(4,5.5,2,2),  family='HersheySerif')
      plot(0, type='n', xlab=' ', ylab=' ', xlim=c(min(x_axis), max(x_axis)), ylim=c(min(y_axis), max(y_axis)),cex.axis=1.3, family='HersheySerif', las=1, main=' ', bty='n')
      grid(lwd = 2)
      
    # Adding Traces
      lines(x=T, y=susceptiblepopulation, col="blue")
      lines(x=T, y=infectedpopulation, col="brown")
      
    # Adding Labels
      mtext(side = 1, text = "T", col = "black", line = 3, cex = 1.3, family='HersheySerif')
      mtext(side = 2, text = "Population", col = "black", line = 4, cex = 1.3, family='HersheySerif')
      
    # Adding Title
      title('SIS Epidemic Over 1-Layer Network', family='serif', cex.main=1.5)
      
    # Adding Legend
      par(mar=c(0, 0, 0, 0), xpd=NA)
      plot(0, type= 'n', bty='n', xlab='', ylab='', axes=FALSE)
      text(x=0.75, y = 0, "Susceptible")
      text(x=0.95, y=0, "Infected")
      segments(0.8, 0, 0.85, 0, col=c('blue'))
      segments(1, 0, 1.05, 0, col=c('brown'))
  }  
  
  g <- cowplot::as_grob(sis_plot)
  p_1 <- cowplot::ggdraw(g)
  
  p_1
  
################################################
#   SAIS Epidemic Model Over 2-layer Network   #
################################################
  
# Prepare the Input Arguments for the GEMF_SIM_Prob Function
  
# Generating a 2-layer Random Network with 300 nodes
  N=300 
  Net1 <- NetGen_Geo(N,0.1)
  Net2 <- NetGen_ER(N,0.03)
  NetSet <- list(Net1,Net2)
  Net <- NetCmbn(NetSet, N)
  
# Setting-Up an SAIS 2-Layer Epidemic Model
  Para <- Para_SAIS_2layer(1,0.2,0.1,0,0.1); 
  M <- Para[[1]]
  
# Generating an Initial Condition where for each Node the Probabilities of Being
# Susceptible, Infected, or Alert Initially Are 0.25, 0.5, and 0.25
  P0 <- matrix(0,M,N)
  P0[1,] <- 0.25
  P0[2,] <- 0.5
  P0[3,] <- 0.25
  
# Values that Specifies When the Simulaiton Terminates
  maxNumevent <- 100000
  Runtime <- 5
  
# Each Core of the Computer will make 25 realization of the epidemic process
  numrun <- 25;  
  comp <- c(1,2,3)
  
# We can run GEMF SIM Prob function in parallel on several cores
  cl <- makeCluster(numcor)
  
# The cluster function makes the list ”result” that has the output of GEMF_SIM_Prob, from different cores, as its elements
  timstp <- 0.1
  result <- clusterCall(cl, GEMF_SIM_Prob, Para, Net, X0=NA, maxNumevent, Runtime, N, numrun, timstp,
                        comp, drawfromprobdis=TRUE, P0)
  stopCluster(cl)
  
# We can incorporate the simulations from different cores
  Tp <- result[[1]][[1]]
  compcu <- result[[1]][[2]]
  for (j in 2:numcor){compcu <- compcu + result[[j]][[2]]}
  
# Total Number of the RealizationS for the Epidemic Process
  s <- numcor*numrun
  comppr <- compcu/s
  dim(comppr)
  susceptible <- colSums(comppr[1, , ])
  infected=colSums(comppr[2, , ])
  alert=colSums(comppr[3, , ])
  
# Plotting Elements
  x_axis <- pretty(Tp)
  y_axis <- pretty(c(susceptible, infected, alert))
  
  sais_plot <- function() {
    # Plotting Base Plot
      layout(rbind(1,2), heights=c(9.5,0.5))
      par(mar = c(4,5.5,2,2),  family='HersheySerif')
      plot(0, type='n', xlab=' ', ylab=' ', xlim=c(min(x_axis), max(x_axis)), ylim=c(min(y_axis), max(y_axis)),cex.axis=1.3, 
           family='HersheySerif', las=1, main=' ', bty='n')
      grid(lwd = 2)
      
    # Plotting Traces
      lines(x=Tp, y=susceptible, col="blue")
      lines(x=Tp, y=infected, col="brown")
      lines(x=Tp, y=alert, col="forestgreen")
      
    # Adding Labels 
      mtext(side = 1, text = "T", col = "black", line = 3, cex = 1.5, family='HersheySerif')
      mtext(side = 2, text = "Average Population", col = "black", line = 4, cex = 1.5, family='HersheySerif')
      
    # Adding Title
      title('SAIS Epidemic Model Over 2-Layer Network', family='serif', cex.main=1.5)
      
    # Adding Legends
      par(mar=c(0, 0, 0, 0), xpd=NA)
      plot(0, type= 'n', bty='n', xlab='', ylab='', axes=FALSE)
      text(x=0.8, y = 0, "Susceptible")
      text(x=1, y=0, "Infected")
      text(x=1.2, y=0, "Alert")
      segments(0.85, 0, 0.9, 0, col=c('blue'))
      segments(1.05, 0, 1.1, 0, col=c('brown'))
      segments(1.25, 0, 1.3, 0, col=c('forestgreen'))
  }
  
  g <- cowplot::as_grob(sais_plot)
  p_2 <- cowplot::ggdraw(g)
  
  p_2
  
############
#   TEST   #
############

# Examining Elements
  net_index <- Net[[1]][[1]]