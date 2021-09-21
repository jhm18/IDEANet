#Transforming Data into Pajek Dataset
#Jonathan H. Morgan
#15 January 2020

#Clear Out Console Script
cat("\014")


#Setting Work Directory
setwd("/Users/jonathan.h.morgan/Desktop")
getwd()

#Options
options(stringsAsFactors = FALSE)

################
#   PACKAGES   #
################

library(network)  

#################
#   FUNCTIONS   #
#################

source("/Applications/Pajek64/Pajek R Tools/R_Pajek Functions_15July2019.R")

######################
#   IMPORTING DATA   #
######################

load('primate_parasite_network3.rdata')

##############################
#   PULLING-OUT ATTRIBUTES   #
##############################

# Looking at Data
  summary(net)

# Getting Edgelist
  edges <- as.data.frame(network::as.edgelist(net))
  
# Getting Edge Values: I see none in the object (assuming 1)
  edges$value <- rep(1, nrow(edges))
  
# Getting Nodes
  nodes <- (unique(append(edges$V1, edges$V2)))
  
# Getting Vertext attributes
  v_attributes <- c('Parasite.Type', 'taxgroup', 'Threat', 'other.mammals', 'vertex.names')
  attribute_list <- vector('list', length(v_attributes))
  names(attribute_list) <- v_attributes
  for (i in seq_along(v_attributes)) {
    attribute_list[[i]] <- get.vertex.attribute(net, v_attributes[[i]])
  }
  
# Making Nodelist with Attributes
  nodelist <- as.data.frame(cbind(nodes, attribute_list[[5]], attribute_list[[1]], attribute_list[[2]], attribute_list[[3]], attribute_list[[4]]))
  colnames(nodelist) <- c('id', 'label', 'parasite.type', 'taxgroup', 'threat', 'other.mammals' )
  nodelist$id <- as.integer(nodelist$id)
  
#########################################
#   TRANSFORM ATTRIBUTES INTO INTEGERS  #
#########################################
  
parasites <- unique(nodelist$parasite.type)
parasites <- as.data.frame(cbind(parasites,seq(1,length(parasites) ,1)))
colnames(parasites) <- c('parasite.type', 'parasite_id')
parasites$parasite_id <- as.numeric(parasites$parasite_id)

taxgroup <- unique(nodelist$taxgroup)
taxgroup  <- as.data.frame(cbind(taxgroup,seq(1,length(taxgroup) ,1)))
colnames(taxgroup) <- c('taxgroup', 'taxgroup_id')
taxgroup$taxgroup_id <- as.numeric(taxgroup$taxgroup_id)

other.mammals <- unique(nodelist$other.mammals)
other.mammals <- as.data.frame(cbind(other.mammals,seq(1,length(other.mammals) ,1)))
colnames(other.mammals) <- c('other.mammals', 'other.mammals_id')
other.mammals$other.mammals_id <- as.numeric(other.mammals$other.mammals_id)

##############################
#   JOINING NUMERIC VALUES   #
##############################

nodelist <- dplyr::left_join(nodelist, parasites, by=c('parasite.type'))
nodelist <- dplyr::left_join(nodelist, taxgroup, by=c('taxgroup'))
nodelist <- dplyr::left_join(nodelist, other.mammals, by=c('other.mammals'))

##############################
#  OUPUTTING PAJEK OBJECTS   #
##############################

# .net file
  write_net('edges', nodelist$label, 
            ' ', ' ', ' ', 
           'blue' , 'black',
           edges$V1, edges$V2, edges$value, ' ', 
           'primate_net_15Jan2020', TRUE)


# Partition Files (.clu)
  write_clu(nodelist$parasite_id, 'parasite')
  write_clu(nodelist$taxgroup_id, 'taxgroup')
  write_clu(nodelist$other.mammals_id, 'other.mammals')
  
# Outputting CSV of node file for reference
  readr::write_csv(nodelist, file.path(getwd(), "primate_nodelist_15Jan2020.csv"))
  
  
  




