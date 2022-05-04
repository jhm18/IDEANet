# Try this instead of X11, better viz without need for installing anything
# quartz(title, width, height, pointsize, family, antialias, type,
# file = NULL, bg, canvas, dpi)

# Better of splitting node-level summary into two separate plots

# Options
  options(stringsAsFactors = FALSE)
  options(scipen = 999)

# Creating Utilities Environment
  IDEANet_Utilities = new.env()

  IDEANet_Utilities$netwrite <- function(data_type = c('edgelist'), adjacency_matrix=FALSE, 
                                         adjacency_list=FALSE, nodelist=FALSE, i_elements=FALSE, 
                                         j_elements=FALSE, weights=FALSE, type=FALSE,
                                         package='igraph', missing_code=99999, 
                                         weight_type='frequency', directed=FALSE, 
                                         net_name='network') {
  
  # Installing Necessary Packages 
    list.of.packages <- c('dplyr', 'igraph', 'network', 'ggplot2', 'cowplot', 'moments')
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)
    rm(list.of.packages, new.packages)
  
  ###########################
  #    C L O S E N E S S    #
  ###########################
  
  # One issue with closeness:
  # We need indegree and outdegree closeness. At the moment we have it based on total degree.
  
  # Create an alternate closeness function
    closeness_igraph <- function(g){ 
      geo <- 1/igraph::distances(g, mode='out')
      diag(geo) <- 0 # Define self-ties as 0
      apply(geo, 1, sum) # Return sum(1/geodist) for each vertex
    }
  
  #############################
  #    B E T W E E N E S S    #
  #############################
  
  betweenness <- function(g, weights){
    # Binarizing Network if Weights Used
      if (as.logical(weights) != FALSE){
        # Binarizing Graph: Getting Nodes
          nodes <- as.data.frame(1:length(igraph::V(g)))
          colnames(nodes) <- c('id')
        
        # Binarizing Graph: Getting edges
          edges <- as.data.frame(igraph::as_edgelist(g, names=FALSE))
        
        # Creating Binarized Graph
          b_g <- igraph::graph_from_data_frame(d = edges[c(1,2)], directed = as.logical(directed), vertices = nodes$id) 
        
        # Calculating Betweenness on Binarized Graph
          b_betweenness <- igraph::betweenness(b_g, directed=as.logical(directed))
      }else{
        g <- g
      }
      
    # Calculating Betweenness 
      if(as.logical(weights) == FALSE){
        betweenness <- igraph::betweenness(g, directed=as.logical(directed), weights=NULL)
      }else{
        betweenness <- igraph::betweenness(g, directed=as.logical(directed), 
                                           weights = igraph::get.edge.attribute(g, "weight"))
      }
    
    # Creating Output Matrix
      if(as.logical(weights) != FALSE){
        betweenness <- as.data.frame(cbind(b_betweenness, betweenness))
        colnames(betweenness) <- c('binarized_betweenness', 'betweenness')
      }else{
        betweenness <- betweenness
      }
     
    # Returning Betweenness Scores
      return(betweenness)
  }
  
  #################################
  #    R E A C H A B I L I T Y    #
  #################################
  
  reachable_igraph <- function(g){
    # Isolating the node's ego-network, the number of reachable nodes, and calculating 
    # the proportion of the total
      proportion_reachable <- vector('numeric', nrow(nodes))
      if(directed == TRUE){
        for(i in seq_along(proportion_reachable)){
          # Isolating connected vertices
            ego_net <- igraph::subcomponent(g, v=igraph::V(g)[[i]], mode = c("out"))
        
          # Eliminating self-loops
            ego_net <- ego_net[ego_net != igraph::V(g)[[i]]]
        
          # Calculating the proportion reachable
            proportion_reachable[[i]] <- length(ego_net)/nrow(nodes)
            rm(ego_net)
        }
      }else{
        for(i in seq_along(proportion_reachable)){
          # Isolating connected vertices
            ego_net <- igraph::subcomponent(g, v=igraph::V(g)[[i]], mode = c("all"))
        
          # Eliminating self-loops
            ego_net <- ego_net[ego_net != igraph::V(g)[[i]]]
        
          # Calculating the proportion reachable
            proportion_reachable[[i]] <- length(ego_net)/nrow(nodes)
          rm(ego_net)
        }
      }
    
    # Writing to global environment
      assign(x = 'reachability', value = proportion_reachable,.GlobalEnv)  
  }
  
  #####################################################
  #    L A R G E S T   W E A K   C O M P O N E N T    #
  #####################################################
  
  largest_weak_component_igraph <- function(g){
    # Isolating the graph's components
      components <- igraph::clusters(g, mode="weak")
      biggest_cluster_id <- which.max(components$csize)
    
    # Extracting the ids of the largest component
      largest_component_ids <- igraph::V(g)[components$membership == biggest_cluster_id]
    
    # Extracting Subgraph
      largest_component <- igraph::induced_subgraph(g, largest_component_ids)
    
    # Assigning the ID list and Subgraph to the Global Environment
      assign(x = 'largest_component_ids', value = largest_component_ids,.GlobalEnv) 
      assign(x = 'largest_component', value = largest_component,.GlobalEnv) 
  }
  
  ###############################################
  #    L A R G E S T   B I C O M P O N E N T    #
  ###############################################
  
  largest_bicomponent_igraph <- function(g) {
    # Extracting bi-components
      bi_components <- igraph::biconnected_components(g)
      bi_component_list <- as.list(bi_components$components)
      bi_lengths <- unlist(lapply(bi_component_list, function(x) length(x)))
      bi_lengths <- cbind(as.data.frame(seq(1, length(bi_lengths), 1)), bi_lengths)
      colnames(bi_lengths) <- c('list_id', 'length')
      largest_id <- bi_lengths[(bi_lengths$length == max(bi_lengths$length)), 1]
      largest_bicomp_ids <- sort(bi_component_list[[largest_id]])
      rm(bi_components, bi_component_list, bi_lengths, largest_id)
    
    # Extracting Subgraph
      largest_bi_component <- igraph::induced_subgraph(g, largest_bicomp_ids)
    
    # Assigning the ID list and Subgraph to the Global Environment
      assign(x = 'largest_bicomponent_ids', value = largest_bicomp_ids,.GlobalEnv) 
      assign(x = 'largest_bi_component', value = largest_bi_component,.GlobalEnv) 
  }
  
  ###########################################
  #    T R A N S I T I V I T Y   R A T E    #
  ###########################################
  
  trans_rate_igraph <- function(g) {
    # Isolating One-Step Paths
      one_step_paths <- vector('list', nrow(nodes))
      names(one_step_paths) <- nodes$id
      for(i in seq_along(one_step_paths)){
        if(length(names(igraph::V(g))) == length(igraph::V(g)) ){
          one_step_paths[[i]] <- as.integer(names(igraph::neighborhood(g, order=1, mindist = 1, igraph::V(g)[[i]], mode='all')[[1]]))
        }else{
          one_step_paths[[i]] <- as.integer(igraph::neighborhood(g, order=1, mindist = 1, igraph::V(g)[[i]], mode='all')[[1]])
        }
      }
      
    # Dropping one-step paths of length 0
      one_step_paths <- one_step_paths[(lengths(one_step_paths) > 0)]
    
    # Isolating Two-Step Paths
      two_step_paths <- vector('list', length(one_step_paths))
      names(two_step_paths) <- names(one_step_paths)
      for(i in seq_along(two_step_paths)){
        paths <- vector('list', length(one_step_paths[[i]]))
        # If a named nodelist else an unnamed list
          if(length(names(igraph::V(g))) == length(igraph::V(g))){
            if(names(igraph::V(g))[[1]] == "0"){
              for(j in seq_along(paths)) {
                paths[[j]] <- as.integer(names(igraph::neighbors(g, (one_step_paths[[i]][[j]] + 1), mode=c('total'))))
              }
            }else{
              for(j in seq_along(paths)) {
                paths[[j]] <- as.integer(names(igraph::neighbors(g, (one_step_paths[[i]][[j]]), mode=c('total'))))
              }
            }
            }else{
              for(j in seq_along(paths)) {
                paths[[j]] <- as.integer(igraph::neighbors(g, (one_step_paths[[i]][[j]]), mode=c('total')))
              }
            }
        
        # Handling when there are and are not two-step paths
          if(length(paths) > 0){
            two_step_paths[[i]] <- sort(unique(unlist(paths)))
          }else{
            two_step_paths[[i]] <- NULL
          }
          rm(paths)
      }
    
    # Eliminating Any NULL Cases (Instances where there is a one-step path but not a two-step)
      two_step_paths[sapply(two_step_paths, is.null)] <- NULL
      
    # Calculating the Proportion of Two-Step Path that Are Also One-Step Paths
      proportion_two_step <- vector('numeric', length(two_step_paths))
      for(i in seq_along(proportion_two_step)) {
        # Identifying node of interest
          vertex_id <- names(two_step_paths)[[i]]
          one_step_path <- one_step_paths[names(one_step_paths) == vertex_id][[1]]
          two_step_path <- two_step_paths[names(two_step_paths) == vertex_id][[1]]
        
        # Identifying Nodes that Occur in Both Two and One-Step Paths
          shared_paths <- sort(intersect(one_step_path, two_step_path))
      
        # Identifying the proportion of nodes that occur on both paths to the number of one-step paths
          proportion_two_step[[i]] <- length(shared_paths)/length(two_step_path)
          rm(vertex_id, one_step_path, two_step_path)
      }
      
    # Creating index for the final calculation to re-incorporate the effect of isolates
      index <- as.data.frame(nodes$id)
      colnames(index) <- c('id')
      results <- cbind(names(two_step_paths), as.data.frame(proportion_two_step))
      colnames(results)[[1]] <- c('id')
      results$id <- as.numeric(results$id)
      index <- dplyr::left_join(index, results, by=c('id'))
      index$proportion_two_step[is.na(index$proportion_two_step)] <- 0
      rm(results)
    
    # Transitivity Rate
      transitivity_rate <- sum(index$proportion_two_step)/length(index$proportion_two_step)
    
    # Assigning transitivity_rate to the global environment
      assign(x = 'transitivity_rate', value = transitivity_rate,.GlobalEnv) 
      rm(one_step_paths, two_step_paths, proportion_two_step)
  } 
  
  #######################################
  #    W E I G H T E D   D E G R E E    #
  #######################################
  
  total_weighted_degree <- function(nodes, edges){
    # Isolating node_ids
      node_ids <- sort(unique(nodes$id))
      node_weights <- vector('numeric', length(node_ids))
    
    # Isolating node acting as ego and as an alter
      for(i in seq_along(node_weights)){
        ego <- edges[(edges[,3] == node_ids[[i]]), ]
        alter <- edges[(edges[,5] == node_ids[[i]]), ]
        node_edges <- rbind(ego, alter)
        node_weights[[i]] <- sum(node_edges[,6])
        rm(ego, alter, node_edges)
      }
    
    # Return node_weights
      return(node_weights)
  }
  
  #############################
  #    C O N S T R A I N T    #
  #############################
  
  constraint.orig <- function(g) {
    # Sub-setting Adjacency Matrix
      idx <- sna::degree(g, diag=FALSE, gmode=gmode, cmode='freeman', ignore.eval=TRUE) != 0
      A <- network::as.matrix.network.adjacency(g)
      A <- A[idx, idx]
      n <- sum(idx)
    
    # Calculating constraint meatures
      one <- c(rep(1,n))
      CZ <- A + t(A)
      cs <- CZ %*% one                      # degree of vertices
      ics <- 1/cs
      CS <- ics %*% t(one)                  # 1/degree of vertices
      P <- CZ * CS                          # intermediate result: proportionate tie strengths
      PSQ <- P%*%P                          # sum paths of length two
      P.bi <- as.numeric(P>0)               # exclude paths to non-contacts (& reflexive):
      PC <- (P + (PSQ*P.bi))^2              # dyadic constraint
      ci <- PC %*% one                      # overall constraint
      dim(ci) <- NULL
    
    # Assigning scores to node ids
      ci2 <- nodes$id
      ci2[idx] <- ci
      ci2[!idx] <- NaN
    
    # Assigning final scores to global environment
      assign(x = 'constraint_score', value = ci2,.GlobalEnv)  
  }
  
  #################################################
  #    D E G R E E   A S S O R T A T I V I T Y    #
  #################################################
  
  assortativity_degree <- function(edges, g) {
    # Extracting the graph's edgelist
      edges <- as.data.frame(edges)
    
    # Calculating the total degree for each node
      node_degree <- sna::degree(g, gmode=gmode, cmode='freeman', ignore.eval=TRUE)
      node_degree <- as.data.frame(cbind(seq(1, length(node_degree), 1), node_degree))
    
    # Joining i & j ids
      colnames(node_degree)[[1]] <- colnames(edges)[[3]]
      edges <- dplyr::left_join(edges, node_degree, by=colnames(edges)[[3]])
      colnames(edges)[[7]] <- c('i_degree')
    
      colnames(node_degree)[[1]] <- colnames(edges)[[5]]
      edges <- dplyr::left_join(edges, node_degree, by=colnames(edges)[[5]])
      colnames(edges)[[8]] <- c('j_degree')
      rm(node_degree)
    
    # Calculating the Pearson Correlation of i and j degree variables
      degree_assortatvity <- stats::cor(edges$i_degree, edges$j_degree, method='pearson')
    
    # Assigning correlation value to the global environment
      assign(x = 'degree_assortatvity', value = degree_assortatvity,.GlobalEnv)
  }
  
  #########################################
  #    A V E R A G E   G E O D E S I C    #
  #########################################
  
  average_geodesic <- function(g) {
    # Generating the number and lengths of all geodesics between all nodes
      gd <- sna::geodist(g, count.paths = FALSE)
    
    # Extracting the distances
      geodesics <- gd$gdist
      geodesics <- geodesics[(lower.tri(geodesics))]
    
    # Replacing infinite values with 0 for the purposes of calculating the average
      geodesics <- geodesics[!is.infinite(geodesics)]
    
    # Calculating the average shortest path length
      average_path_length <- mean(geodesics)
    
    # Assgining to the global environment       
      assign(x = 'average_path_length', value = average_path_length,.GlobalEnv) 
      rm(gd, geodesics)
  }
  
  #############################################################
  #    M U L T I P L E X   E D G E   C O R R E L A T I O N    #
  #############################################################
  
  multiplex_edge_corr_igraph <- function(edgelist, directed) {
    if('type' %in% colnames(edgelist)){
      # Creating edgelist to manipulate internally
        edges <- as.data.frame(edgelist[,])
      
      # Moving back to One-Index for Comparison Purposes
        edges[,3] <- edges[,3] + 1
        edges[,5] <- edges[,5] + 1
      
      # Recovering original weight for the purposes of comparison
        if(weight_type == 'frequency') {
          edges[,6] <- as.numeric(1/edges[,6])
        }else{
          edges[,6] <- edges[,6]
        }
      
      # Generating Correlations Either as Directed or Undirected
        if(as.logical(directed) == TRUE) {
          # Generating Sub-Networks Based on Type
            types <- sort(unique(type))
            subnets <- vector('list', length(types))
            names(subnets) <- types
            for(i in seq_along(types)){
              subnets[[i]] <- as.data.frame(edges[(type == types[[i]]), ])
              subnets[[i]] <- subnets[[i]][,c('i_id', 'j_id', 'type', 'weight')]
              colnames(subnets[[i]])[[3]] <- names(subnets)[[i]]
              colnames(subnets[[i]])[[4]] <- paste0(colnames(subnets[[i]])[[3]],'_',colnames(subnets[[i]])[[4]])
            }
        
          # Creating a Wide Data-Set to Generate Correlations
            ties <- unique(as.data.frame(edges[ ,c("i_id", "j_id")]))
            for(i in seq_along(types)){
              ties <- dplyr::left_join(ties, subnets[[i]], by=c('i_id', 'j_id'))
              ties[is.na(ties)] <- 0
            }
        
          # Calculating the Correlation for Unique Combination of Types 
            pairs <- t(combn(paste0(types,'_','weight'), 2))
            for(i in nrow(pairs)) {
              column_set <- pairs[i,]
              tie_set <- ties[,column_set]
              multiplex_edge_correlation <- paste0('Edge Correlation for ', paste(column_set, collapse= ' and '), ': ', round(stats::cor(tie_set)[1,2], digits=2))
              rm(column_set, tie_set)
            }
            rm(pairs, types, subnets, ties)
        }else{
          # Creating a separate edgelist (Symmetric Edges) to Perform Operations
            s_edges <- edges[,c('i_id', 'j_id', 'type', 'weight')]
        
          # Eliminating Duplicate Pairs
            s_edges <- s_edges[!duplicated(t(apply(s_edges[,c(1:2)], 1, sort))),]
        
          # Creating Edge Groups & Glossary
            edges_1 <- cbind(s_edges[,c(1,2)], seq(1, dim(s_edges)[[1]], 1))
            colnames(edges_1)[[3]] <- c('edge_group')
        
            edges_2 <- cbind(s_edges[,c(2,1)], seq(1, dim(s_edges)[[1]], 1))
            colnames(edges_2) <- c('i_id','j_id','edge_group')
        
            edges_glossary <- rbind(edges_1, edges_2)
            edges_glossary <- edges_glossary[order(edges_glossary$edge_group), ]
            rm(edges_1, edges_2, s_edges)
        
          # Joining edge_groups to edges
            if('Obs_ID' %in% colnames(edgelist)){
              edges <- edges
            }else{
              edges <- cbind(seq(1, dim(edges)[[1]], 1), edges)
              names(edges)[[1]] <- c('Obs_ID')
            }
            edges <- dplyr::left_join(as.data.frame(edges), edges_glossary, by=c('i_id', 'j_id'))
        
          # Eliminating Duplicates Caused by Self-Loops
            edges <- edges[!(duplicated(edges$Obs_ID)), ]
            rm(edges_glossary)
        
          # Collapsing Ties and Summing Weights
            edge_groups <- unique(edges$edge_group)
            ties <- vector('list', length(edge_groups))
            names(ties) <- edge_groups
            for(i in seq_along(edge_groups)) {
              e_group <- edges[(edges$edge_group == edge_groups[[i]]), ]
              row.names(e_group) <- seq(1, nrow(e_group), 1)
              e_types <- unique(e_group$type)
              ties[[i]] <- as.data.frame(e_group$type)
              ties[[i]]$weight <- sum(e_group$weight)
              ties[[i]]$i_id <- e_group[1,3]
              ties[[i]]$j_id <- e_group[1,5]
              colnames(ties[[i]])[[1]] <- c('type')
              ties[[i]] <- ties[[i]][,c(3,4,1,2)]
              rm(e_group, e_types)
            }
        
            ties <- do.call("rbind", ties)
        
          # Generating Sub-Networks Based on Type
            types <- sort(unique(type))
            subnets <- vector('list', length(types))
            names(subnets) <- types
            for(i in seq_along(types)){
              subnets[[i]] <- ties[(ties$type == types[[i]]), ]
              colnames(subnets[[i]])[[3]] <- names(subnets)[[i]]
              colnames(subnets[[i]])[[4]] <- paste0(colnames(subnets[[i]])[[3]],'_',colnames(subnets[[i]])[[4]])
            }
        
          # Creating a Wide Data-Set to Generate Correlations
            ties <- unique(ties[ ,c("i_id", "j_id")])
            for(i in seq_along(types)){
              ties <- dplyr::left_join(ties, subnets[[i]], by=c('i_id', 'j_id'))
              ties[is.na(ties)] <- 0
            }
        
          # Calculating the Correlation for Unique Combination of Types 
            pairs <- t(combn(paste0(types,'_','weight'), 2))
            for(i in nrow(pairs)) {
              column_set <- pairs[i,]
              tie_set <- ties[,column_set]
              multiplex_edge_correlation <- paste0('Edge Correlation for ', paste(column_set, collapse= ' and '), ': ', round(stats::cor(tie_set)[1,2], digits=2))
              rm(column_set, tie_set)
            }
            rm(pairs, types, subnets, ties)
        }
      
      # Assigning final scores to global environment
        assign(x = 'multiplex_edge_correlation', value = multiplex_edge_correlation,.GlobalEnv)  
    }else{
      edgelist <- edgelist[,]
      multiplex_edge_correlation <- 'Simplex Network'
    }
  }
  
  #########################
  #    B O N A C I C H    #
  #########################
  
  # Custom Bonacich Function
    bonacich <- function(matrix, bpct = .75) {
      # Calculate eigenvalues
        evs <- eigen(matrix)
    
      # The ones we want are in the first column
      # Something's weird here -- it's combining the two columns that SAS outputs into a single value
      # This value is a complex sum. For now, it seems like coercing the complex sum using `as.numeric`
      # gets us the values we want.
    
      # Get maximum eigenvalue
        maxev <- max(as.numeric(evs$value))
    
      # Get values that go into computation
        b <- bpct*(1/maxev) # Diameter of power weight
        n <- nrow(matrix) # Size
        i <- diag(n) # Identity matrix
        w <- matrix(rep(1, n), ncol = 1) # Column of 1s
    
      # For some reason, the output of the `solve` function is the transpose 
      # of the `INV` function in SAS. I'm going to transpose the output of `solve` here
      # so that results are consistent with what we get in SAS, but this is something
      # we need to confirm and possibly be wary of.
    
      # Key equation, this is the centrality score
        C <- (t(solve(i-b*matrix)))%*%t(matrix)%*%w 
      
      # This is Bonacich normalizing value alpha
        A <- sqrt(n/(t(C) %*% C))
      
      # This is the power centrality score
        cent <- c(A) * c(C)
    
      # This is the centralization score
        NBCNT <- sum(max(C) - C)/((n-1)*(n-2))
      
      # Vector of centralization scores
        NBCD <- c(matrix(NBCNT, ncol = n))
    
      # Collect power centrality scores and centralization scores into single dataframe
        bonacich_output <- data.frame(bonacich = cent, bon_centralization = NBCD)
    
      return(bonacich_output)
    }
  
  # Bonacich igraph
    bonacich_igraph <- function(g, directed) {
      # Store complete nodelist for merging later on
        nodelist <- data.frame(id = igraph::V(g)$name)
    
      # Detect if isolates are present in the network
        if (0 %in% igraph::degree(g, mode = "all")) {
          # If isolates are present, indicate that isolates are going to be removed
            message("(Bonacich power centrality) Isolates detected in network. Isolates will be removed from network when calculating power centrality measure, and will be assigned NA values in final output.")
      
          # Remove isolates
            g <- igraph::delete.vertices(g, v = igraph::degree(g, mode = "all", loops = F) == 0)
        }
    
      # Convert igraph object into adjacency matrix
        bon_adjmat <- as.matrix(igraph::get.adjacency(g, type = "both", names = TRUE))
    
      # We need to ensure that the adjacency matrix used in calculating eigenvectors/values
      # is not singular. The following checks for this. If the adjacency matrix is found
      # to be singular, network will be treated as undirected when calculating EVs
        if (directed == TRUE) {
          # Get generalized inverse of matrix
            inv_adj <- MASS::ginv(bon_adjmat)
      
            singular_check <- singular_check <- round(sum(diag(inv_adj %*% bon_adjmat)))
      
          if (singular_check < nrow(nodelist)) {
            message("(Bonacich power centrality) Adjacency matrix for network is singular. Network will be treated as undirected in order to calculate measures\n")
            directed <- FALSE
            g <- igraph::as.undirected(g)
          }
        }
  
      # When we have a directed network, there are three power centrality scores we can get:
      # An "indegree" one based on the original adjacency matrix, and "outdegree" one based on the
      # transpose of the original adjacency matrix, and a third one based on a symmetrized version
      # of the adjacency matrix. The following conditional flow generates all three measures
      # if netwrite is working with a directed network:
        if (directed == T) {
          # Create symmetrized (undirected) version of network
            undir_net <- igraph::as.undirected(g)
            
          # Convert undirected network into adjacency matrix
            bon_sym_mat <- as.matrix(igraph::get.adjacency(undir_net, type = "both", names = TRUE))
      
          # We now have everyting we need to get the three versions of Bonacich power centrality
          # First let's get the indegree version
            bonacich_in <- bonacich(matrix = bon_adjmat)
            
          # Update column names
            colnames(bonacich_in) <- paste(colnames(bonacich_in), "_in", sep = "")
      
          # Next we get the outdegree version (note that we're transposing `bon_adjmat` in the `matrix` argument)
            bonacich_out <- bonacich(matrix = t(bon_adjmat))
            
          # Update column names
            colnames(bonacich_out) <- paste(colnames(bonacich_out), "_out", sep = "")
      
          # Finally, we get the undirected version from the symmetrized adjacency matrix
            bonacich_sym <- bonacich(matrix = bon_sym_mat)
            
          # Update column names
            colnames(bonacich_sym) <- paste(colnames(bonacich_sym), "_sym", sep = "")
            
          # Combine into single data frame
            bon_scores <- cbind(bonacich_in, bonacich_out, bonacich_sym)
      }else{
        # If the network is undirected, proceed to get Bonacich power centrality scores on
        # just the base adjacency matrix
          bon_scores <- bonacich(bon_adjmat)
      }
    
      # Add ID variable for merging back into nodelist
        bon_scores$id <- igraph::V(g)$name
        
      # Merge scores back into nodelist
        nodelist <- dplyr::left_join(nodelist, bon_scores, by = "id")
        
      # Remove ID variable
        nodelist$id <- NULL
    
      return(nodelist)
    }
  
  #####################################################
  #    E I G E N V E C T O R   C E N T R A L I T Y    #
  #####################################################
  
  # Custom Function for Calculating Eigenvector Centrality
    eigen_custom <- function(matrix) {
      # To replicate output from SAS, we first need to transpose the adjacency matrix
        eigen_calc <- eigen(t(matrix))
    
      # Eigenvector centrality is the eigenvector associated with the largest
      # eigenvalue. I think by convention this is always the first one, but going to find
      # it explicitly to be sure:
    
      # Remove complex sums from eigenvalues
        evs <- as.numeric(eigen_calc$values) 
      
      # Indicate which is the maximum eigenvalue
        max_eval <- which(evs == max(evs)) 
    
      # If multiple eigenvectors have the maximum eigenvalue,
      # take only the first one
        if (length(max_eval) > 1) {
          max_eval <- max_eval[1]
        }
    
      # Now get that column from the eigenvector matrix;
      # these are the eigenvector centrality scores.
      # You sometimes get negative values, but that's
      # not an error. Just take the absoluate value of
      # negative values
        eigencent <- abs(as.numeric(eigen_calc$vectors[,max_eval]))
    
      # Prepare data frame for output
        eigen_df <- data.frame(eigen_centrality = eigencent)
        return(eigen_df)
    }
  
  # IGRAPH
    eigen_igraph <- function(g, directed){
      # Store complete nodelist for merging later on
        nodelist <- data.frame(id = igraph::V(g)$name)
    
      # Detect if isolates are present in the network
        if (0 %in% igraph::degree(g, mode = "all")) {
          # If isolates are present, indicate that isolates are going to be removed
            message("(Eigenvector centrality) Isolates detected in network. Isolates will be removed from network when calculating eigenvector centrality measure, and will be assigned NA values in final output.\n")
      
          # Remove isolates
            g <- igraph::delete.vertices(g, v = igraph::degree(g, mode = "all", loops = F) == 0)
        }
    
      # Assign component membership as a vertex attribute
        igraph::V(g)$component <- igraph::components(g)$membership
    
      # Calculate Eigenvector centrality for each component
      # Unique component ids
        unique_components <- unique(igraph::V(g)$component)
    
      # We need to ensure that the adjacency matrix used in calculating eigenvectors/values
      # is not singular. The following checks for this. If the adjacency matrix is found
      # to be singular, network will be treated as undirected when calculating EVs
        if (directed == TRUE) {
          # Get adjacency matrix
            check_adj <- as.matrix(igraph::as_adjacency_matrix(g, type = "both"))
      
          # Get generalized inverse of matrix
            inv_adj <- MASS::ginv(check_adj)
            singular_check <- singular_check <- round(sum(diag(inv_adj %*% check_adj)))
            
            if (singular_check < nrow(nodelist)) {
              message("(Eigenvector centrality) Adjacency matrix for network is singular. Network will be treated as undirected in order to calculate measures\n")
              directed <- FALSE
              g <- igraph::as.undirected(g)
            }
        }
        
      # Detect if multiple components exist in the network
        if (length(unique_components) > 1) {
          # Outputting message to the user
            message("(Eigenvector centrality) Network consists of 2+ unconnected components. Eigenvector centrality scores will be calculated for nodes based on their position within their respective components.\n")
      
          # Initialize data frame for storing eigen centrality measures
            eigen_scores <- data.frame()
      
          # If the network is a directed network
            if (directed == TRUE) {
              # For each component...
                for (i in 1:length(unique_components)) {
                  # Make subgraph of component
                    subgraph <- igraph::delete.vertices(g, v = igraph::V(g)$component != i)
                  
                  # Convert subgraph of component into an adjacency matrix
                    sub_adj <- as.matrix(igraph::as_adjacency_matrix(subgraph, type = "both"))
                  
                  # Make transpose of subgraph adjmat
                    sub_adj_t <- t(sub_adj)
          
                  # Make undirected version of component
                    subgraph_undir <- igraph::as.undirected(subgraph)
                  
                  # Convert into adjacency matrix
                    undir_mat <- as.matrix(igraph::as_adjacency_matrix(subgraph_undir, type = "both"))
          
                  # Get eigenvector centrality measures for indegree
                    eigen_in <- eigen_custom(sub_adj)
                  
                  # Update names to indicate indegree
                    colnames(eigen_in) <- paste(colnames(eigen_in), "_in", sep = "")
          
                  # Get eigenvector centrality measures for outdegree
                    eigen_out <- eigen_custom(sub_adj_t)
                  
                  # Update names to outdicate outdegree
                    colnames(eigen_out) <- paste(colnames(eigen_out), "_out", sep = "")
                  
                  # On symmetric matrix
                    eigen_sym <- eigen_custom(undir_mat)
                    colnames(eigen_sym) <- paste(colnames(eigen_sym), "_sym", sep = "")
          
                  # Combine into single dataframe
                    subgraph_scores <- cbind(eigen_in, eigen_out, eigen_sym)
                  
                  # Add component indicator
                    subgraph_scores$component <- i
                  
                  # Add ID variable
                    subgraph_scores$id <- igraph::V(subgraph)$name
          
                  # Bind to `eigen_scores` data frame
                    eigen_scores <- rbind(eigen_scores, subgraph_scores)
                }
            
            } else {
              # For each component...
                for (i in 1:length(unique_components)) {
                  # Make subgraph of component
                    subgraph <- igraph::delete.vertices(g, v = igraph::V(g)$component != i)
              
                  # Convert subgraph of component into an adjacency matrix
                    sub_adj <- as.matrix(igraph::as_adjacency_matrix(g, type = "both"))
          
                  # Get eigenvector centrality measures
                    subgraph_scores <- eigen_custom(sub_adj)
              
                  # Add component indicator
                    subgraph_scores$component <- i
              
                  # Add ID variable
                    subgraph_scores$id <- igraph::V(subgraph)$name
          
                  # Bind to `eigen_scores` data frame
                    eigen_scores <- rbind(eigen_scores, subgraph_scores)
                }
            }
        }else{
          # If the network is a directed network
            if (directed == TRUE) {
              # Convert graph to adjacency matrix
                eigen_adj <- as.matrix(igraph::as_adjacency_matrix(g, type = "both"))
                  
              # Make transpose of subgraph adjmat
                eigen_adj_t <- t(eigen_adj)
        
              # Make undirected version of component
                eigen_undir <- igraph::as.undirected(g)
                  
              # Convert into adjacency matrix
                undir_mat <- as.matrix(igraph::as_adjacency_matrix(g, type = "both"))
        
              # Get eigenvector centrality measures for indegree
                eigen_in <- eigen_custom(eigen_adj)
                  
              # Update names to indicate indegree
                colnames(eigen_in) <- paste(colnames(eigen_in), "_in", sep = "")
                  
              # Get eigenvector centrality measures for outdegree
                eigen_out <- eigen_custom(eigen_adj_t)
                  
              # Update names to outdicate outdegree
                colnames(eigen_out) <- paste(colnames(eigen_out), "_out", sep = "")
        
              # On symmetric matrix
                eigen_sym <- eigen_custom(undir_mat)
                colnames(eigen_sym) <- paste(colnames(eigen_sym), "_sym", sep = "")
        
              # Combine into single dataframe
                eigen_scores <- cbind(eigen_in, eigen_out, eigen_sym)
                  
              # Add ID variable
                eigen_scores$id <- igraph::V(g)$name
            } else {
              # If the network is undirected...
              # Convert graph to adjacency matrix
                eigen_adj <- as.matrix(igraph::as_adjacency_matrix(g, type = "both"))
                  
              # Get Eigenvector centrality measures
                eigen_scores <- eigen_custom(eigen_adj)
                  
              # Add ID variable
                eigen_scores$id <- igraph::V(g)$name 
            }
        }
        
      # Merge `eigen_scores` back into nodelist
        nodelist <- dplyr::left_join(nodelist, eigen_scores, by = "id")
          
      # Remove `id` variable
        nodelist$id <- NULL
    
      # Return Result
        return(nodelist)
    }
  
  ###############################################
  #    N O D E - L E V E L   M E A S U R E S    #
  ###############################################
  
  node_level_igraph <- function(nodes, g, directed) {
    total_degree <- igraph::degree(g, mode='all', loops=FALSE)
    weighted_degree <- igraph::strength(g, mode='all', loops=FALSE)
    in_degree <- igraph::degree(g, mode='in', loops=FALSE)
    out_degree <- igraph::degree(g, mode='out', loops=FALSE)
    closeness <- closeness_igraph(g)
    betweenness_scores <- betweenness(g, weights)
    # bonpow <- igraph::bonpow(g, loops=FALSE, exponent = 0.75)
    bonpow <- bonacich_igraph(g, directed=as.logical(directed))
    #eigen_cen <- igraph::eigen_centrality(g, directed=as.logical(directed), scale=FALSE)$vector
    eigen_cen <- eigen_igraph(g, directed = as.logical(directed))
    constraint <- igraph::constraint(g)
    reachability <- reachable_igraph(g)
    
    nodes <- as.data.frame(cbind(nodes, total_degree, weighted_degree, in_degree, out_degree, 
                                 closeness, betweenness_scores, bonpow, eigen_cen, constraint, reachability))
    
    return(nodes)
  }
  
  #########################################
  #    C O N D I T I O N A L   F L O W    #
  #########################################
  
  # Setting Data Type: Adjacency Matrix, Adjacency List, or Edgelist
  
  # ADJACENCY MATRIX
    if(data_type == 'adjacency_matrix'){
      # Checking for ID Column
        if (dim(adjacency_matrix)[[1]] != dim(adjacency_matrix)[[2]]){
          adjacency_matrix <- adjacency_matrix[,c(2:ncol(adjacency_matrix))]
        }else{
          adjacency_matrix <- adjacency_matrix[,]
        }    
    
        if(as.logical(directed) == TRUE){
          # Generating directed graph
            g <- igraph::graph_from_adjacency_matrix(adjacency_matrix, mode=c('directed'), diag = TRUE)
      
          # Create Nodes file with Node-level measures
            edges <- as.data.frame(igraph::as_edgelist(g, names=FALSE))
          
          # nodes <- as.data.frame(sort(unique(c(edges$V1, edges$V2))))
          # colnames(nodes) <- c('id')
          ### The above two lines were dropping isolates from the nodelist.
          ### I think this is a better alternative
            nodes <- as.data.frame(1:length(igraph::V(g)))
            colnames(nodes) <- "id"
            nodes$id <- nodes$id - 1
            nodes$label <- igraph::V(g)$name
          
          # To keep things consistent across code, we're going to reassign node names
          # to a `label` vertex attribute in `igraph` and replace the `name` attribute
          # with numeric IDs
            igraph::V(g)$name <- nodes$id
            igraph::V(g)$label <- nodes$label
          
          # Create alternate closeness function
          # Reachability function (eliminating loops)
      
          # Adding Node-level measures
            nodes <- node_level_igraph(nodes = nodes, g = g, directed = directed)
      
          # Extracting largest weakly-connected component
          # Extracting largest bicomponent
          # Calculating proportion of two-step paths that are also one-step paths (trans_rate)
      
          # Calculating system-level measures
            largest_weak_component_igraph(g)
            largest_bicomponent_igraph(g)
            degree_assortatvity <- igraph::assortativity.degree(g, directed=as.logical(directed))
            reciprocity_rate <- igraph::reciprocity(g, ignore.loops = TRUE, mode='ratio')
            trans_rate_igraph(g)
            global_clustering_coefficient <- igraph::transitivity(g, type='global')
            average_path_length <- igraph::average.path.length(g, directed=as.logical(directed))
      
          # UNDIRECTED
        } else {
          # Generating undirected graph
            g <- igraph::graph_from_adjacency_matrix(adjacency_matrix, mode=c('undirected'), diag = FALSE)
      
          # Create Nodes file with Node-level measures
            edges <- as.data.frame(igraph::as_edgelist(g, names=FALSE))
          
          # nodes <- as.data.frame(sort(unique(c(edges$V1, edges$V2))))
          # colnames(nodes) <- c('id')
          ### The above two lines were dropping isolates from the nodelist.
          ### I think this is a better alternative
            nodes <- as.data.frame(1:length(igraph::V(g)))
            colnames(nodes) <- "id"
            nodes$id <- nodes$id - 1
            nodes$label <- igraph::V(g)$name
          
          # To keep things consistent across code, we're going to reassign node names
          # to a `label` vertex attribute in `igraph` and replace the `name` attribute
          # with numeric IDs
            igraph::V(g)$name <- nodes$id
            igraph::V(g)$label <- nodes$label
      
          # Create an alternate closeness function
          # Reachablility function (Eliminate Loops, reaching yourself isn't that useful)
      
          # Adding Node-Level Measures
            nodes <- node_level_igraph(nodes = nodes, g = g, directed = directed)
      
          # Extracting the largest weakly connected component
          # Extracting the largest bi-component
          # Calculating the Proportion of Two-Step Path that Are Also One-Step Paths
      
          # Calculating System-Level Measures
            largest_weak_component_igraph(g)
            largest_bicomponent_igraph(g)
            degree_assortativity <- igraph::assortativity.degree(g, directed=as.logical(directed))
            reciprocity_rate <- igraph::reciprocity(g, ignore.loops = TRUE, mode='ratio')
            trans_rate_igraph(g)
            global_clustering_coefficient <- igraph::transitivity(g, type='global')
            average_path_length <- igraph::average.path.length(g, directed=as.logical(directed))
        } 
      
  # ADJACENCY LIST
    }else if (data_type == 'adjacency_list') {
      # Is the adjacency list a list
        if (class(adjacency_list) == 'list') {
          g <- igraph::graph_from_adj_list(adjacency_list, mode="out")
        } else {
          # IF NOT, Converting to a list
            adj_list <- vector('list', dim(adjacency_list)[[1]])
            names(adj_list) <- as.character(adjacency_list[,1])
            for(i in seq_along(adj_list)){
              adj_row <- unique(as.integer(strsplit(adjacency_list[i,2], ' ')[[1]]))
              adj_list[[i]] <- vector('list', length(adj_row))
              for(j in seq_along(adj_row)) {
                adj_list[[i]][[j]] <- adj_row[[j]]
              }
              rm(adj_row)
            }
      
          # Generating network from adjacency list
            g <- igraph::graph_from_adj_list(adj_list, mode="out")
        }
    
      # Copying igraph object
        g <- g
    
      # Creating Nodes File with Node-Level Measures
        edges <- as.data.frame(igraph::as_edgelist(g, names=FALSE))
        nodes <- as.data.frame(sort(unique(c(edges$V1, edges$V2))))
        colnames(nodes) <- c('id')
        nodes$id <- nodes$id - 1
    
      # Create an alternate closeness function
      # Reachablility function (Eliminate Loops, reaching yourself isn't that useful)
      # Adding Node-Level Measures
        nodes <- node_level_igraph(nodes = nodes, g = g, directed = directed)
        
      # Extracting the largest weakly connected component
      # Extracting the largest bi-component
      # Calculating the Proportion of Two-Step Path that Are Also One-Step Paths
      # Calculating System-Level Measures
        largest_weak_component_igraph(g)
        largest_bicomponent_igraph(g)
        degree_assortatvity <- igraph::assortativity.degree(g, directed=as.logical(directed))
        reciprocity_rate <- igraph::reciprocity(g, ignore.loops = TRUE, mode='ratio')
        trans_rate_igraph(g)
        global_clustering_coefficient <- igraph::transitivity(g, type='global')
        average_path_length <- igraph::average.path.length(g, directed=as.logical(directed))
   
  # EDGELIST   
    } else {
      # Creating Canonical Node and Edgelists
        if(weights[[1]]==FALSE){
          edgelist <-as.matrix(cbind(i_elements, j_elements))
          edgelist <-cbind(edgelist, rep(1,nrow(edgelist)))
          colnames(edgelist)[[3]] <- c('weight')
        }else{
          edgelist <-as.matrix(cbind(i_elements, j_elements, weights))
          colnames(edgelist)[[3]] <- c('weight')
        }
    
      # Checking for Edge Type
        if(type[[1]] == FALSE){
          edgelist <- edgelist
        }else if (length(type) == length(i_elements)){
          edgelist <- cbind(edgelist, type)
        }else {
          writeLines("The type indicator variable is not the same length as the network's edgelist.\nTo calculate the network's multilevel edge correlation, please supply a vector of the same length.")
        }
    
        edgelist <- edgelist[!(rowSums(is.na(edgelist))), ]
        edgelist <- edgelist[edgelist[,1] != missing_code & edgelist[,2] != missing_code, ] 
        edgelist <- cbind(seq(1,nrow(edgelist), 1), edgelist)
        colnames(edgelist)[[1]] <- c('Obs_ID')
    
      # Adding Nodes
        if(nodelist[[1]] == FALSE) {
          nodes <- as.data.frame(sort(unique(c(edgelist[,2], edgelist[,3]))))
          nodes <- cbind(seq(1,nrow(nodes),1), nodes)
          colnames(nodes) <- c('id', 'label')
      
          senders <- as.data.frame(edgelist[,c(1:2)])
          colnames(senders)[[2]] <- c('label')
          senders <- dplyr::left_join(senders, nodes, by='label')
          colnames(senders)[c(2,3)] <- c('i_elements', 'i_id')
      
          if(type[[1]] == FALSE){
            targets <- as.data.frame(edgelist[,c(1,3,4)])
          }else {
            targets <- as.data.frame(edgelist[,c(1,3,4, 5)])
          }
          
          colnames(targets)[[2]] <- c('label')
          targets <- dplyr::left_join(targets, nodes, by='label')
          if(type[[1]] == FALSE){
            colnames(targets)[c(2,4)] <- c('j_elements', 'j_id')
            targets <- targets[c(1,2,4,3)]
          } else {
            colnames(targets)[c(2,5)] <- c('j_elements', 'j_id')
            targets <- targets[c(1,2,5,3,4)]
          }
      
          edgelist <- dplyr::left_join(senders, targets, by='Obs_ID')
          edgelist <- edgelist[order(edgelist$i_id, edgelist$j_id), ]
          edgelist <- as.matrix(edgelist)
          rm(senders, targets)
      }else {
        nodes <- nodelist
      
        nodes <- cbind(as.data.frame(seq(1, length(nodes), 1)), nodes)
        colnames(nodes) <- c('id', 'label')
      
        senders <- as.data.frame(edgelist[,c(1:2)])
        colnames(senders)[[2]] <- c('label')
        senders <- dplyr::left_join(senders, nodes, by='label')
        colnames(senders)[c(2,3)] <- c('i_elements', 'i_id')
      
        if(type[[1]] == FALSE){
          targets <- as.data.frame(edgelist[,c(1,3,4)])
        }else{
          targets <- as.data.frame(edgelist[,c(1,3,4, 5)])
        }
        
        colnames(targets)[[2]] <- c('label')
        targets <- dplyr::left_join(targets, nodes, by='label')
        if(type[[1]] == FALSE){
          colnames(targets)[c(2,4)] <- c('j_elements', 'j_id')
          targets <- targets[c(1,2,4,3)]
        }else{
          colnames(targets)[c(2,5)] <- c('j_elements', 'j_id')
          targets <- targets[c(1,2,5,3,4)]
        }
      
        edgelist <- dplyr::left_join(senders, targets, by='Obs_ID')
        edgelist <- edgelist[order(edgelist$i_id, edgelist$j_id), ]
        edgelist <- as.matrix(edgelist)
        rm(senders, targets)
    }
    
    
    # Convert edgelist to data frame class and make `Obs_ID`, `i_id`, `j_id`, and `weight` numeric
      edgelist <- as.data.frame(edgelist)
      edgelist$Obs_ID <- as.numeric(edgelist$Obs_ID)
      edgelist$i_id <- as.numeric(edgelist$i_id)
      edgelist$j_id <- as.numeric(edgelist$j_id)
      edgelist$weight <- as.numeric(edgelist$weight)
    
    # Create graph objects
    
    # Make Zero-Indexed
      nodes$id <- nodes$id - 1
      edgelist[,3] <- edgelist[,3] - 1
      edgelist[,5] <- edgelist[,5] - 1
    
    # Make Weights Reflect Frequency Rather than Distance
      if(weight_type == 'frequency') {
        edgelist[,6] <- as.numeric(1/edgelist[,6])
      }else{
        edgelist[,6] <- edgelist[,6]
      }
    
    # Creating igraph object
      colnames(nodes)[[2]] <- c('attr')
      g <- igraph::graph_from_data_frame(d = edgelist[,c(3,5)], directed = as.logical(directed), vertices = nodes) 
    
    # Adding edge weights
    igraph::edge.attributes(g)$weight <- edgelist[,6]
    
    # Create an alternate closeness function
    # Reachablility function (Eliminate Loops, reaching yourself isn't that useful)
    # Adding Node-Level Measures
    nodes <- node_level_igraph(nodes = nodes, g = g, directed = directed)
    
    # Extracting the largest weakly connected component
    # Extracting the largest bi-component
    # Calculating the Proportion of Two-Step Path that Are Also One-Step Paths
    # Calculating Multiplex Edge Correlation
    # Calculating System-Level Measures
    largest_weak_component_igraph(g)
    largest_bicomponent_igraph(g)
    degree_assortatvity <- igraph::assortativity.degree(g, directed=as.logical(directed))
    reciprocity_rate <- igraph::reciprocity(g, ignore.loops = TRUE, mode='ratio')
    trans_rate_igraph(g)
    global_clustering_coefficient <- igraph::transitivity(g, type='global')
    average_path_length <- igraph::average.path.length(g, directed=as.logical(directed))
    
    multiplex_edge_corr_igraph(edgelist = edgelist, directed = as.logical(directed))
    
    # Outputting Network Objects
    assign(x = 'edgelist', value = edgelist,.GlobalEnv)  
    assign(x = 'nodelist', value = nodes,.GlobalEnv)  
    assign(x = net_name, value = g,.GlobalEnv)
  } # End edgelist condition
  
  ###########################################
  #    G E N E R A T I N G   R E P O R T    #
  ###########################################  
    
  # System-Level Data Object
    num_clusters <- igraph::clusters(g, mode="weak")[[3]]
    proportion_largest <- max(igraph::clusters(g, mode="weak")[[2]])/nrow(nodes)
  
  # Creating system-level data object
    multiplex_edge_correlation <- ifelse(type==FALSE, 'Singleplex Network', multiplex_edge_correlation)
    multiplex_edge_correlation <- multiplex_edge_correlation[[1]]
  
    measure_labels <- c('Number of Components', 'Proportion in the Largest Component',
                        'Degree Assortativity', 'Reciprocity Rate', 'Transitivity Rate', 
                        'Global Clustering Coefficient', 'Average Geodesic',
                        'Multi-Level Edge Correlation')
    
    measure_descriptions <- c('The number of weak components in the graph', 
                              'The proportion of nodes in the largest weak component of the graph',
                              'Edgewise correlation of degree', 'The proportion of directed ties that are reciprocated',
                              'The proportion of two-step paths that are also one-step paths',
                              'The proportion of closed triangles to all triangles', 'The average shortest path length',
                              'Multiplex networks edgwise correlation of relations')
    
    measures <- c(num_clusters, proportion_largest, degree_assortatvity, reciprocity_rate,
                  transitivity_rate, global_clustering_coefficient, average_path_length,
                  multiplex_edge_correlation)
    
    system_level_measures <- cbind(as.data.frame(measure_labels), measure_descriptions, measures)
  
  # Removing node-level and system-level data objects for clarity
    rm(measure_labels, measure_descriptions, num_clusters, proportion_largest, degree_assortatvity,
       reciprocity_rate, global_clustering_coefficient, average_path_length,
       multiplex_edge_correlation, measures)
  
    rm(transitivity_rate, reachability, envir = .GlobalEnv)
  
  # System & Node-Level Visualizations
    #quartz(width = 10.6806, height = 7.30556)
    system_plot <- function() {
      # Creating Layout
        viz_matrix <- matrix(c(10,10,10,10,10,10,10,10,10,
                               2,2,2,3,3,3,0,0,0,
                               1,1,1,1,1,1,4,4,4,
                               1,1,1,1,1,1,0,0,0,
                               1,1,1,1,1,1,5,5,5,
                               1,1,1,1,1,1,6,6,6,
                               1,1,1,1,1,1,0,0,0,
                               1,1,1,1,1,1,9,9,9,
                               7,7,7,8,8,8,0,0,0), 
                            ncol  = 9, byrow = TRUE)
        layout(viz_matrix)
    
      # Defining degree distribution coordinates
        y_axis <- density(nodes$total_degree)$y
        x_axis <- density(nodes$total_degree)$x
        coordinates <- cbind(as.data.frame(x_axis), y_axis)
        coordinates <- coordinates[(coordinates$x_axis >= 0), ]
        x_axis <- pretty(coordinates$x_axis)
        y_axis <- pretty(coordinates$y_axis)
        x_spacer <- x_axis[c(length(x_axis))] - x_axis[c(length(x_axis)-1)]
        x_spacer <- x_spacer*0.5
        y_spacer <- y_axis[c(length(y_axis))] - y_axis[c(length(y_axis)-1)]
        y_spacer <- y_spacer*0.5
    
      # Defining Base Degree Plot
        par(mar = c(5,6,2,2),  family='HersheySerif')
        plot(0, type='n', xlab=' ', ylab=' ', xlim=c(min(x_axis), max(x_axis)), 
             ylim=c(min(y_axis), max(y_axis)), cex.axis=1.3, family='HersheySerif', 
             las=1, main=' ', bty='n')
        grid(lwd = 2)
    
      # Adding Margin Text
        mtext(side = 1, text = 'Total Degree', col = "black", line = 3, cex = 1.5, family='HersheySerif')
        mtext(side = 2, text = 'Density', col = "black", line = 4.5, cex = 1.5, family='HersheySerif')
    
      # Plotting Degree
        lines(coordinates$x_axis, coordinates$y_axis, col='brown', lwd=1.5)
    
      # Adding Skew and Kurtosis
        skewness <- moments::skewness(nodes$total_degree)
        kurtosis <- moments::kurtosis(nodes$total_degree)
        text(x = (max(x_axis)-x_spacer), y = (max(y_axis)-y_spacer), paste('Skewness',round(skewness, digits=2)), cex=1.3)
        text(x = (max(x_axis)-x_spacer), y = (max(y_axis)-(y_spacer*2)), paste('Kurtosis',round(kurtosis, digits=2)), cex=1.3)
    
      # Adding Title
        title(c("Total Degree Distribution"), family='serif', cex.main=2)
    
      # Populating Subplots
        for(i in seq_along(system_level_measures$measure_labels)) {
          plot_measure <- system_level_measures[i,3]
      
          plot_measure <- ifelse(i < 8, as.numeric(plot_measure), plot_measure)
          plot_measure <- ifelse(i < 8, round(plot_measure, digits=2), plot_measure)
          plot_measure <- ifelse(i == 8, trimws(gsub('Edge', '', plot_measure)), plot_measure)
      
          par(mar=c(0,0,0,0), family='serif')
          plot(0, type='n', xlab=' ', ylab=' ', xlim=c(1,10), 
               ylim=c(1,10), axes=FALSE, main='', bty='n')
      
          text(x=5, y=9, system_level_measures[i,1], family='serif', font=2, cex=1.3)
          text(x=5, y=6.5, plot_measure, family='serif', cex=1.5)
          rm(plot_measure)
        }
    
      # Adding Plot Title
        par(mar=c(0,0,0,0), family='serif')
        plot(0, type='n', xlab=' ', ylab=' ', xlim=c(1,10), 
             ylim=c(1,10), axes=FALSE, main='', bty='n')
        text(x=5.5, y=5, 'System-Level Measures', family='serif', font=2, cex=3)
    } 
  
    g <- cowplot::as_grob(system_plot)
    p_1 <- cowplot::ggdraw(g)
  
    p_1 
  
    node_measures_plot <- function() {
      # Specifying nicer labels
        if(directed == TRUE){
          plot_labels <- c('Weighted Degree', 'In-Degree', 'Out-Degree', 'Closeness', 
                           'Betweenness', 
                           'Bonacich Power Centrality (In)', 'Bonacich Power Centrality (Out)', 'Bonacich Power Centrality (Sym)',
                           'Eigenvector Centrality (In)', 'Eigenvector Centrality (Out)', 'Eigenvector Centrality (Sym)', 
                           'Constraint', 'Reachability')
        }else{
          plot_labels <- c('Weighted Degree', 'Closeness', 'Betweenness', 'Bonacich Power Centrality',
                           'Eigenvector Centrality', 'Constraint', 'Reachability')
        }
    
      # Isolating the measure being visualized based on whether it's directed or not
        if(directed == TRUE){
          plot_measures <- c("weighted_degree", "in_degree", "out_degree",     
                             "closeness", "betweenness", 
                             "bonacich_in", "bonacich_out", "bonacich_sym",
                             "eigen_centrality_in", "eigen_centrality_out", "eigen_centrality_sym",
                             "constraint", "reachability")
        }else{
          plot_measures <- c("weighted_degree", "closeness", "betweenness", 
                             "bonpow", "eigen_cen", "constraint", "reachability")
        }
    
      # Defining the layout used 
        if(directed == TRUE){
          viz_matrix <- matrix(c(14,14,14,14,14,14,14,14,14,
                                 1,1,1,2,2,2,3,3,3,
                                 1,1,1,2,2,2,3,3,3,
                                 1,1,1,2,2,2,3,3,3,
                                 4,4,4,5,5,5,0,0,0,
                                 4,4,4,5,5,5,0,0,0,
                                 4,4,4,5,5,5,0,0,0,
                                 6,6,6,7,7,7,8,8,8,
                                 6,6,6,7,7,7,8,8,8,
                                 6,6,6,7,7,7,8,8,8,
                                 9,9,9,10,10,10,11,11,11,
                                 9,9,9,10,10,10,11,11,11,
                                 9,9,9,10,10,10,11,11,11,
                                 12,12,12,13,13,13,0,0,0,
                                 12,12,12,13,13,13,0,0,0,
                                 12,12,12,13,13,13,0,0,0,
                                 0,0,0,0,0,0,0,0,0), 
                            ncol  = 9, byrow = TRUE)
        layout(viz_matrix)
      }else{
        viz_matrix <- matrix(c(8,8,8,8,8,8,8,8,8,
                             1,1,1,2,2,2,3,3,3,
                             1,1,1,2,2,2,3,3,3,
                             1,1,1,2,2,2,3,3,3,
                             4,4,4,5,5,5,6,6,6,
                             4,4,4,5,5,5,6,6,6,
                             4,4,4,5,5,5,6,6,6,
                             7,7,7,0,0,0,0,0,0,
                             7,7,7,0,0,0,0,0,0,
                             7,7,7,0,0,0,0,0,0), 
                           ncol  = 9, byrow = TRUE)
          layout(viz_matrix)
      }
    
      # Generating Subplot
        for(i in seq_along(plot_measures)){
          # Eliminating NA Values
            plot_measure <- nodes[,plot_measures[[i]]]
            plot_measure <- plot_measure[!is.na(plot_measure)]
      
          # Defining degree distribution coordinates
            y_axis <- density(plot_measure)$y
            x_axis <- density(plot_measure)$x
            coordinates <- cbind(as.data.frame(x_axis), y_axis)
            coordinates <- coordinates[(coordinates$x_axis >= min(plot_measure)), ]
            x_axis <- pretty(coordinates$x_axis)
            y_axis <- pretty(coordinates$y_axis)
    
          # Defining Base Degree Plot
            if(directed == TRUE) {
              par(mar = c(2.6,7,2.2,2), family='HersheySerif')
            } else {
              par(mar = c(5,6,2,2), family='HersheySerif')
            }
        
            plot(0, type='n', xlab=' ', ylab=' ', xlim=c(min(x_axis), max(x_axis)), 
                 ylim=c(min(y_axis), max(y_axis)), cex.axis=1.3, family='HersheySerif', 
                 las=1, main=' ', bty='n')
            grid(lwd = 2)
      
          # Adding Margin Text
            mtext(side = 1, text = plot_labels[[i]], col = "black", line = 3, cex = 1.3, family='HersheySerif')
      
          # Plotting Degree
            lines(coordinates$x_axis, coordinates$y_axis, col='brown', lwd=1.5)
        } 
    
      # Adding Title
        par(mar=c(0,0,0,0), family='serif')
        plot(0, type='n', xlab=' ', ylab=' ', xlim=c(1,10), 
             ylim=c(1,10), axes=FALSE, main='', bty='n')
        text(x=5.5, y=5, 'Node-Level Measures', family='serif', font=2, cex=3)
    }
  
    g <- cowplot::as_grob(node_measures_plot)
    p_2 <- cowplot::ggdraw(g)
  
    p_2
  
  # Assigning Report Elements to the Global Environment
    assign(x = 'system_measure_plot', value = p_1,.GlobalEnv)  
    assign(x = 'node_measure_plot', value = p_2,.GlobalEnv)
    assign(x = 'system_level_measures', value = system_level_measures, .GlobalEnv)
  
  }

###################################
#    T E S T   E D G E L I S T    #
###################################

test_edgelist <- data.frame(ego = c("Tom", "Tom", "Tom", "Jim", "Gabe", "Gabe", "Gabe", "Jon", "Jon", "Dan"),
                            alter = c("Jim", "Dan", "Gabe", "Jon", "Tom", "Dan", "Jon", "Tom", "Jim", "Gabe"))

data_type = "edgelist"
adjacency_matrix <- F
adjacency_list <- F
nodelist <- FALSE
i_elements =  test_edgelist$ego
j_elements = test_edgelist$alter
#weights = FALSE
weights <- c(0.50, 0.75, 1.00, 0.25, 0.55, 0.45, 0.35, 0.99, 0.19, 0.05)
directed =  TRUE
missing_code = 99999
type = FALSE
net_name = 'net_1'
package = 'network'
weight_type='distance'

IDEANet_Utilities$netwrite(data_type = "edgelist",
                           adjacency_matrix <- F,
                           adjacency_list <- F,
                           nodelist <- FALSE,
                           i_elements =  test_edgelist$ego,
                           j_elements = test_edgelist$alter,
                           weights = weights,
                           directed =  TRUE,
                           missing_code = 99999,
                           type = FALSE,
                           net_name = 'net_1',
                          # package = 'network',
                           weight_type='frequency')

# NetExplorer Test Data (Network with Isolates)
  test_edgelist <- read_csv("Dropbox/My Mac (Jonathan’s MacBook Pro)/Desktop/DNAC/IDEANet/Data_Scripts/NetExploration_App/test_edgelist.csv")
  test_nodedata <- read_csv("Dropbox/My Mac (Jonathan’s MacBook Pro)/Desktop/DNAC/IDEANet/Data_Scripts/NetExploration_App/test_nodedata.csv")
  
  IDEANet_Utilities$netwrite(data_type = "edgelist",
                             adjacency_matrix <- F,
                             adjacency_list <- F,
                             nodelist <- test_nodedata$Node,
                             i_elements =  test_edgelist$In,
                             j_elements = test_edgelist$Out,
                             weights = test_edgelist$Weight,
                             directed =  TRUE,
                             missing_code = 99999,
                             type = FALSE,
                             net_name = 'net_1',
                             # package = 'network',
                             weight_type='frequency')
  
# Quickly Checking the Network
  weights <- igraph::get.edge.attribute(net_1, "weight", igraph::E(net_1))
  labels <- igraph::get.vertex.attribute(net_1, "name", igraph::E(net_1))
  par(mar=c(0,0,0,0))
  plot(net_1)
  
###################################################
#    T E S T   A D J A C E N C Y   M A T R I X    #
###################################################

adjvec=c(0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1,
         0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1,
         0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1,
         0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1,
         0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0)
test_adjmat <- matrix(adjvec, nrow = 10, ncol = 10)
rownames(test_adjmat) <- c("Tom", "Jim", "Dana", "Gabe", "Jon", "Peter", "Ethan", "Dan", "Scott", "Liann")
colnames(test_adjmat) <- c("Tom", "Jim", "Dana", "Gabe", "Jon", "Peter", "Ethan", "Dan", "Scott", "Liann")
rm(adjvec)

test_adjmat[,1] <- 0
test_adjmat[1,] <- 0

# # net_igraph <- igraph::graph_from_adjacency_matrix(test_adjmat, mode = "directed")
# # net_network <- network::network(test_adjmat)
# 

IDEANet_Utilities$netwrite(data_type = "adjacency_matrix",
                           adjacency_matrix <- test_adjmat,
                           adjacency_list <- F,
                           nodelist <- FALSE,
                           i_elements =  NULL,
                           j_elements = NULL,
                           weights = FALSE,
                           directed =  TRUE,
                           missing_code = 99999,
                           type = FALSE,
                           net_name = 'test_net',
                           package = 'network',
                           weight_type='frequency')



###############
