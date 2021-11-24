# IDEANet's netwrite and netread Functions User Script
# Jonathan H. Morgan
# 23 November 2021

# Clear Out Console Script
  cat("\014")
  
# Options
  options(stringsAsFactors = FALSE)
  
# Sourcing IDEANet Functions
  source("/Users/jonathan.h.morgan/Desktop/DNAC/IDEANet/Data_Scripts/R_Scripts/IDEANet_Network Base Functions.R")

###############################
#   EXAMPLE CASE: AHS_WPVAR   #
###############################

# Setting Data Directory & Importing Example Data
  setwd('/Users/jonathan.h.morgan/Desktop/DNAC/IDEANet/Data_Scripts')
  import_data <- function(file_csv) {
    # Installing Necessary Packages 
      list.of.packages <- c('readr')
      new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
      if(length(new.packages)) install.packages(new.packages)
      rm(list.of.packages, new.packages)
  
    # Reading-In CSV with readr
      file_csv <- paste0(file_csv, '.csv')
      base_data <- readr::read_csv(file_csv, col_names = TRUE)
  
    # Creating Stacked Edgelist
      communities <- vector('list', length(unique(base_data$commcnt)))
      community_ids <- sort(unique(base_data$commcnt))
      names(communities) <- community_ids
      
    # Creating community-level list elements
      for (i in seq_along(communities)) {
        community <- base_data[base_data$commcnt == community_ids[[i]], ]
        egos <- sort(unique(community$ego_nid))
    
        alters <- vector('list', length(egos))
        names(alters) <- egos
    
        for (j in seq_along(egos)) {
          ego <- community[community$ego_nid == egos[[j]], ]
          m_friends <- as.data.frame(as.integer(ego[,c(4:8)]))
          m_friends <- cbind('Male', m_friends)
          colnames(m_friends) <- c('gender','alter_id')
      
          f_friends <- as.data.frame(as.integer(ego[,c(14:18)]))
          f_friends <- cbind('Female', f_friends)
          colnames(f_friends) <- c('gender','alter_id')
      
          alter_net <- rbind(m_friends, f_friends)
          alter_net <- cbind(egos[[j]], alter_net)
          alter_net <- cbind(community_ids[[i]], alter_net)
          colnames(alter_net)[c(1:2)] <- c('commcnt', 'ego_nid')
      
          alters[[j]] <- alter_net
      
          rm(alter_net, ego, m_friends, f_friends)
      }
    
        communities[[i]] <- do.call("rbind", alters)
    
        rm(community, egos, alters)
      }
      
    # Assigning communities edglist to the Global Environment
      assign(x = 'communities_edgelist', value = communities,.GlobalEnv)  
  
      return(base_data)
  }

  base_data <- import_data('ahs_wpvar')

######################################
#   EXAMPLE APPLICATIONS OF NETWRITE # 
######################################

# Evoking a Custom Plot Window
  x11(width=10.6806, height=7.30556)

# EDGELIST EXAMPLE: MULTIPLE COMMUNITIES
  
# Arguments for edge_examples()
# 1)  net_package: Specifies whether iGraph or network graph objects should be used.
# 2)  directed_graph: Specifies whether the graph is directed or undirected.
# 3)  community_count: Specifies how many communities to example, with a max of 85 communities
  edge_examples <- function(net_package='igraph', directed_graph=TRUE, community_count=3) {
    # Creating Output Lists: Communities 1-3
      if(directed_graph==TRUE){
        directed_example_networks <- vector('list', community_count)
        names(directed_example_networks) <- paste0('net_',seq(1, community_count, 1))
      }else{
        undirected_example_networks <- vector('list', community_count)
        names(undirected_example_networks) <- paste0('net_',seq(1, community_count, 1))
      }
    
      node_lists <- vector('list', community_count)
      names(node_lists) <- paste0('net_',seq(1, community_count, 1))
    
      edge_lists <- vector('list', community_count)
      names(edge_lists) <- paste0('net_',seq(1, community_count, 1))
    
      largest_components <- vector('list', community_count)
      names(largest_components) <- paste0('net_',seq(1, community_count, 1))
    
      largest_comp_ids <- vector('list', community_count)
      names(largest_comp_ids) <- paste0('net_',seq(1, community_count, 1))
    
      largest_bicomponents <- vector('list', community_count)
      names(largest_bicomponents) <- paste0('net_',seq(1, community_count, 1))
    
      largest_bicomp_ids <- vector('list', community_count)
      names(largest_bicomp_ids) <- paste0('net_',seq(1, community_count, 1))
    
      system_measures <- vector('list', community_count)
      names(system_measures) <- paste0('net_',seq(1, community_count, 1))
    
      system_plots <- vector('list', community_count)
      names(system_plots) <- paste0('net_',seq(1, community_count, 1))
    
      node_level_plots <- vector('list', community_count)
      names(node_level_plots) <- paste0('net_',seq(1, community_count, 1))
    
    # Generating objects and measures 
      for(i in seq_along(node_lists)) {
        # Evoking netwrite to generate network objects and metrics
          IDEANet_Utilities$netwrite(data_type = c('edgelist'), adjacency_matrix=FALSE, 
                                     adjacency_list=FALSE,nodelist=FALSE, 
                                     i_elements=communities_edgelist[[i]]$ego_nid, 
                                     j_elements=communities_edgelist[[i]]$alter_id,
                                     weights=FALSE, type=FALSE, package=net_package, 
                                     missing_code=99999, weight_type='frequency', 
                                     directed=directed_graph, net_name='net')
    
        # Populating Output Objects
          if(directed_graph == TRUE){
            directed_example_networks[[i]] <- net
          }else{
            undirected_example_networks[[i]] <- net
          }
          
          node_lists[[i]] <- nodelist
          edge_lists[[i]] <- edgelist
          largest_components[[i]] <- largest_component
          largest_bicomponents[[i]] <- largest_bi_component
      
          if(net_package == 'igraph'){
            largest_comp_ids[[i]] <- cbind(as.data.frame(names(largest_component_ids)), 
                                           as.integer(largest_component_ids))
            colnames(largest_comp_ids[[i]]) <- c('igraph_label', 'igraph_id')
          
            largest_bicomp_ids[[i]] <- cbind(as.data.frame(names(largest_bicomponent_ids)), 
                                             as.integer(largest_bicomponent_ids))
            colnames(largest_bicomp_ids[[i]]) <- c('igraph_label', 'igraph_id')
          }else{
            largest_comp_ids[[i]] <- largest_component_ids
            largest_bicomp_ids[[i]] <- largest_bicomponent_ids
          }
      
          system_measures[[i]] <- system_level_measures
          system_plots[[i]] <- system_measure_plot
          node_level_plots[[i]] <- node_measure_plot
    
          rm(net, nodelist, edgelist, largest_component, largest_component_ids, 
             largest_bi_component, largest_bicomponent_ids, 
             system_level_measures, system_measure_plot, node_measure_plot, envir=.GlobalEnv)
      }
      
    # Assigning Objects
      if(directed_graph==TRUE){
        assign(x = 'directed_example_networks', value = directed_example_networks,.GlobalEnv)  
      }else{
        assign(x = 'undirected_example_networks', value = undirected_example_networks,.GlobalEnv) 
      }
        
      assign(x = 'node_lists', value = node_lists,.GlobalEnv)
      assign(x = 'edge_lists', value = edge_lists,.GlobalEnv)
      assign(x = 'largest_components', value = largest_components,.GlobalEnv)
      assign(x = 'largest_comp_ids', value = largest_comp_ids,.GlobalEnv)
      assign(x = 'largest_bicomponents', value = largest_bicomponents,.GlobalEnv)
      assign(x = 'largest_bicomp_ids', value = largest_bicomp_ids,.GlobalEnv)
      assign(x = 'system_measures', value = system_measures,.GlobalEnv)
      assign(x = 'system_plots', value = system_plots,.GlobalEnv)
      assign(x = 'node_level_plots', value = node_level_plots,.GlobalEnv)
  }

  edge_examples(net_package='igraph', directed_graph=FALSE, community_count = 3)
  
# Edgelist Example: Plotting Directed Networks
  par(mar=c(0,0,0,0))
  plot(directed_example_networks$net_1)
  plot(directed_example_networks$net_2)
  plot(directed_example_networks$net_3)
  
# Edgelist Example: Plotting Undirected Networks
  par(mar=c(0,0,0,0))
  plot(undirected_example_networks$net_1)
  plot(undirected_example_networks$net_2)
  plot(undirected_example_networks$net_3)
  
# Edgelist Example: Plotting Largest Weak Components
  par(mar=c(0,0,0,0))
  plot(largest_components$net_1)
  plot(largest_components$net_2)
  plot(largest_components$net_3)
  
# Edgelist Example: Plotting Largest Bi-Components
  par(mar=c(0,0,0,0))
  plot(largest_bicomponents$net_1)
  plot(largest_bicomponents$net_2)
  plot(largest_bicomponents$net_3)
  
# Edgelist Example: Examining System-Level Metrics
  system_plots$net_1
  system_plots$net_2
  system_plots$net_3
  
# Edgelist Example: Examining Node-Level Metrics
  node_level_plots$net_1
  node_level_plots$net_2
  node_level_plots$net_3

# MULTIPLEX NETWORK USING AN EDGELIST

# Creating an Example Tie Types, Say Friendship and Advice
  types <- c(1, 2)
  tie_type <- sample(types, dim(communities_edgelist[[1]])[[1]], replace=TRUE)
  rm(types)
  
# Also generating some weights to indicate the strength of the tie
  t_weights <- c(1, 2, 3, 4, 5, 6, 7)
  tie_weight <- sample(t_weights, dim(communities_edgelist[[1]])[[1]], replace=TRUE)
  rm(t_weights)
  
# Looking at Multiplex Edge Correlation of the Types for Network 1
  IDEANet_Utilities$netwrite(data_type = c('edgelist'), adjacency_matrix=FALSE, 
                             adjacency_list=FALSE, nodelist=FALSE, 
                             i_elements=communities_edgelist[[1]]$ego_nid, 
                             j_elements=communities_edgelist[[1]]$alter_id, 
                             weights=tie_weight, type=tie_type, package='igraph', 
                             missing_code=99999, weight_type='frequency', 
                             directed='TRUE', net_name='net_1')

  # Plotting Graph
    edge_colors <- as.character(ifelse(edgelist[,7] == 1, 'black', 'red'))
    igraph::E(net_1)$color <- edge_colors
    igraph::E(net_1)$width <- as.numeric(edgelist[,6])
    
    par(mar=c(0,0,0,0))
    plot(net_1)
    
  # Examining System-Level Measures
    system_measure_plot

  # Examining Node-Level Measures
    node_measure_plot

  rm(largest_component, largest_component_ids, largest_bi_component, largest_bicomponent_ids,
     system_level_measures, system_measure_plot, node_measure_plot)
  
# Supplying a nodelist as well as weights and tie types
  nodes <- nodelist$attr
  
  IDEANet_Utilities$netwrite(data_type = c('edgelist'), adjacency_matrix=FALSE, 
                             adjacency_list=FALSE, nodelist=nodes, 
                             i_elements=communities_edgelist[[1]]$ego_nid, 
                             j_elements=communities_edgelist[[1]]$alter_id, 
                             weights=tie_weight, type=tie_type, package='network', 
                             missing_code=99999, weight_type='frequency', 
                             directed='TRUE', net_name='net_1')

  # Plotting Graph
    par(mar=c(0,0,0,0))
    plot(net_1, edge.col=edge_colors, edge.lwd=(as.numeric(edgelist[,6])*0.10))
    
  # Examining System-Level Measures
    system_measure_plot
  
  # Examining Node-Level Measures
    node_measure_plot
  
  rm(largest_component, largest_component_ids, largest_bi_component, largest_bicomponent_ids,
     system_level_measures, system_measure_plot, node_measure_plot)
  
# GENERATING NETWORKS WITH AN ADJACENCY MATRIX
  
# Creating an Adjacency Matrix from net_1
  edges <- as.data.frame(as.matrix(net_1, matrix.type="edgelist"))
  edges <- cbind(seq(1, nrow(edges), 1), edges)
  colnames(edges) <- c('Obs_ID', 'i_id', 'j_id')

  adj_mat <- matrix(nrow=nrow(nodelist), ncol=nrow(nodelist))
  for(i in seq_along(edges$Obs_ID)){
    edge <- edges[i, ]
    adj_mat[edge$j_id, edge$i_id] <- 1
    rm(edge)
  }
  adj_mat[is.na(adj_mat)] <- 0
  colnames(adj_mat) <- seq(1, nrow(nodelist), 1)
  rm(edges)

  IDEANet_Utilities$netwrite(data_type = c('adjacency_matrix'), adjacency_matrix=adj_mat, 
                             adjacency_list=FALSE, nodelist=FALSE, 
                             i_elements=FALSE, j_elements=FALSE, 
                             weights=FALSE, type=FALSE,
                             package='igraph', missing_code=99999, 
                             weight_type='frequency', 
                             directed='TRUE', net_name='net_2')

  # Plotting Graph
    par(mar=c(0,0,0,0))
    plot(net_2)
    
  # Examining System-Level Measures
    system_measure_plot

  # Examining Node-Level Measures
    node_measure_plot

  rm(largest_component, largest_component_ids, largest_bi_component, largest_bicomponent_ids,
      system_level_measures, system_measure_plot, node_measure_plot)

# GENERATING NETWORKS WITH AN ADJACENCY LIST  
  
# Adjacency List Example
  adjacency_list <- igraph::as_adj_list(net_2, mode='out')

  adj_list <- vector('list', length(adjacency_list))
  for (i in seq_along(adjacency_list)) {
    adj_row <- unique(as.integer(adjacency_list[[i]]))
    adj_row <- paste(adj_row, collapse = ' ')
    adj_list[[i]] <- adj_row
    rm(adj_row)
  }

  adj_list <- as.data.frame(do.call("rbind", adj_list))
  adj_list <- cbind(names(adjacency_list), adj_list)
  colnames(adj_list) <- c('ego', 'alters')

  IDEANet_Utilities$netwrite(data_type = c('adjacency_list'), adjacency_matrix=FALSE, 
                             adjacency_list=adj_list, nodelist=FALSE, 
                             i_elements=FALSE, j_elements=FALSE, 
                             weights=FALSE, type=FALSE, package='igraph', 
                             missing_code=99999, weight_type='frequency', 
                             directed='TRUE', net_name='net_3')

  # Plotting Graph
    par(mar=c(0,0,0,0))
    plot(net_3)
    
  # Examining System-Level Measures
    system_measure_plot

  #  Examining Node-Level Measures
     node_measure_plot

  rm(largest_component, largest_component_ids, largest_bi_component, largest_bicomponent_ids,
     system_level_measures, system_measure_plot, node_measure_plot)
  
# GENERATING A MULTIPLEX NETWORK FROM ADJACENCY MATRICES

  
#####################################
#   EXAMPLE APPLICATIONS OF NETREAD # 
#####################################

# network
  IDEANet_Utilities$netread(package='network', network_object=net_1)
  
# iGraph
  IDEANet_Utilities$netread(package='igraph', network_object=net_2)
  