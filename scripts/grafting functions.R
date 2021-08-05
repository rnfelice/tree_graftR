#grafting functions

paste_tips<-function(tree, tipname, node_number, tip_age){
  target_node_date <- as.numeric(suppressMessages(dateNodes(tree)[node_number]))
  if(!any(which(tree$edge[,2]==node_number))){
    tree$root.edge<-3
    required_edge_length<-c((target_node_date-tip_age) + .5, (target_node_date-tip_age) + 1.0, (target_node_date-tip_age) + 2.0)
    steps <- c(.5, 1.0, 2.0)
    newtrees <- lapply(1:length(steps), function(k) phytools::bind.tip(tree = tree, tip.label = tipname, edge.length = required_edge_length[k],where = node_number,position = steps[k]))
    class(newtrees)<-"multiPhylo"
     } else{

  if (target_node_date>tip_age){
    edgelength<-tree$edge.length[which(tree$edge[,2]==node_number)]
    steps<-seq((edgelength/4), (edgelength-(edgelength/4)), length.out = 3)
    required_edge_length <- target_node_date+steps-tip_age
    tree$root.edge<-10
    #i think i can comment out this line
    #newtrees <- lapply(1:length(steps), function(k) phytools::bind.tip(tree = tree, tip.label = tipname, edge.length = required_edge_length[k],where = node_number,position = steps[k]))
    #class(newtrees)<-"multiPhylo"

  }
  else {
    parent_of_target <- tree$edge[which(tree$edge[,2]==node_number),1]
    parentnode_date <- as.numeric(suppressMessages(dateNodes(tree)[parent_of_target]))
    offset <- tip_age - target_node_date
    edgelength <- parentnode_date-tip_age
    required_edge_length<-seq((edgelength/4), (edgelength-(edgelength/4)), length.out = 3)
    steps <- required_edge_length+offset
  }

  newtrees <- lapply(1:length(steps), function(k) phytools::bind.tip(tree = tree, tip.label = tipname, edge.length = required_edge_length[k],where = node_number,position = steps[k]))
  class(newtrees)<-"multiPhylo"
  }
  return(newtrees)
}



tree_grafter <- function(tree_w_fossils, target_tree, species_data,max_trees=100,verbose=TRUE){
  composite_tree<-target_tree
  time_travelers<-c()
  tree_to_annimate<-composite_tree
  nodedates <- suppressMessages(dateNodes(tree_w_fossils))

  nodetable <- tibble(Taxon = tree_w_fossils$tip.label , Tip_Date = round(nodedates, 4)[1:length(tree_w_fossils$tip.label)])

  #good_molec_tree<-keep.tip(tree_badnames,tree_badnames$tip.label[which(tree_badnames$tip.label %in% nodetable$Taxon)])

  full_species_data <- left_join(species_data, nodetable, by = c("Tip_Label" = "Taxon" ))

  #inputs:
  #old
  fossil_species_data<- full_species_data %>%
    filter(., !Tip_Label %in% target_tree$tip.label) %>%
    select(., c(Tip_Label, Tip_Date)) %>%
    arrange(., Tip_Date)

  phylodists <- as.matrix(distTips(tree_w_fossils))

  writeLines("Abandon all hope, ye who enter here.")
  for (i in 1:nrow(fossil_species_data)){


    if( class(composite_tree) == "phylo"){
      one_topology <- TRUE
      composite_tree <- c(composite_tree,composite_tree)
    } else {
      one_topology <- FALSE
    }

    nodenum_a<-which(tree_w_fossils$tip.label == fossil_species_data$Tip_Label[i])
    tip_name_to_paste <- fossil_species_data$Tip_Label[i]
    #get sister clade to fossil
    tips_nearby<-"none"
    sister_clade_to_fossil<-ips::sister(tree_w_fossils,tip_name_to_paste,label = FALSE)
    if(any(composite_tree[[1]]$tip.label %in% tree_w_fossils$tip.label[sister_clade_to_fossil[which(sister_clade_to_fossil<=Ntip(tree_w_fossils))]])){
      tips_nearby <- composite_tree[[1]]$tip.label[which(composite_tree[[1]]$tip.label %in% tree_w_fossils$tip.label[sister_clade_to_fossil[which(sister_clade_to_fossil<=Ntip(tree_w_fossils))]])]
    }
    while(tips_nearby[[1]]=="none"){
      nodenum_b <- getMRCA(tree_w_fossils, c(nodenum_a, sister_clade_to_fossil))
      sister_clade_to_fossil<-ips::sister(tree_w_fossils,nodenum_b,label = FALSE)
      if(any(composite_tree[[1]]$tip.label %in% tree_w_fossils$tip.label[sister_clade_to_fossil[which(sister_clade_to_fossil<=Ntip(tree_w_fossils))]])){
        tips_nearby <- composite_tree[[1]]$tip.label[which(composite_tree[[1]]$tip.label %in% tree_w_fossils$tip.label[sister_clade_to_fossil[which(sister_clade_to_fossil<=Ntip(tree_w_fossils))]])]
      }
    }

    anc_node_in_fossil_tree <- tree_w_fossils$edge[which(tree_w_fossils$edge[,2] == nodenum_a),1]

    fossil_divergence<- nodedates[tree_w_fossils$edge[which(tree_w_fossils$edge[,2] == anc_node_in_fossil_tree),1]]

    #are any of these in the tree already?
    if(length(tips_nearby)==1){
      target_node <-  which(composite_tree[[1]]$tip.label==tips_nearby)
    } else {
      target_node <- getMRCA(composite_tree[[1]], tips_nearby)
    }


    #####THIS IS WHERE WE GET A LIST OF NODE DATES!
    #date the current tree:
    current_nodedates <- suppressMessages(lapply(composite_tree, dateNodes))


    #get the age of the target node on the current tree:
    target_node_ages <- sapply(current_nodedates, "[[", target_node)

    anc_to_target <- composite_tree[[1]]$edge[which(composite_tree[[1]]$edge[,2]==target_node),1]

    ###GET ANCESTOR NODE AGES
    if(length(anc_to_target) == 0){
     ancestor_node_ages <- target_node_ages + 10
    } else {
    ancestor_node_ages<- sapply(current_nodedates, "[[", anc_to_target)
    }
    ##check if target tip age is older than its ancestors node
   # if (min(ancestor_node_ages) < fossil_species_data$Tip_Date[i]){
   #   time_travelers <- c(time_travelers, fossil_species_data$Tip_Label[i])
   # }

    #check which topologies would produce invalid trees with this paste
    #if the tip date is older than the ancestral node you will get a negative branch so we will skip it

   impossible_topoology_numbers<-which((ancestor_node_ages-fossil_species_data$Tip_Date[i])<=0)

   #do i need a singleton checker here?

   if (length(impossible_topoology_numbers)>0){
         for (jj in 1:length(impossible_topoology_numbers)){
           num<-impossible_topoology_numbers[jj]
           composite_tree[[num]]$edge.length[which(composite_tree[[num]]$edge[,1]==anc_to_target)] <- composite_tree[[num]]$edge.length[which(composite_tree[[num]]$edge[,1]==anc_to_target)] + (fossil_divergence-ancestor_node_ages[num])
           composite_tree[[num]]$edge.length[which(composite_tree[[num]]$edge[,2]==anc_to_target)] <- composite_tree[[num]]$edge.length[which(composite_tree[[num]]$edge[,2]==anc_to_target)] - (fossil_divergence-ancestor_node_ages[num])
           #time_travelers <- c(time_travelers, fossil_species_data$Tip_Label[i])
         }
       }
#  if (length(impossible_topoology_numbers)==length(composite_tree)){
#    for (jj in 1:length(composite_tree)){
#      composite_tree[[jj]]$edge.length[which(composite_tree[[jj]]$edge[,1]==anc_to_target)] <- composite_tree[[jj]]$edge.length[which(composite_tree[[jj]]$edge[,1]==anc_to_target)] + (fossil_divergence-ancestor_node_ages[jj])
#      composite_tree[[jj]]$edge.length[which(composite_tree[[jj]]$edge[,2]==anc_to_target)] <- composite_tree[[jj]]$edge.length[which(composite_tree[[jj]]$edge[,2]==anc_to_target)] - (fossil_divergence-ancestor_node_ages[jj])
#      #time_travelers <- c(time_travelers, fossil_species_data$Tip_Label[i])
#    }
#  }
     impossible_topoology_numbers<-which((ancestor_node_ages-fossil_species_data$Tip_Date[i])<=0)

   if(length(impossible_topoology_numbers) > 0 & length(impossible_topoology_numbers) < length(composite_tree)){
     composite_tree<-composite_tree[-impossible_topoology_numbers]
   }
#
   if(length(composite_tree)>max_trees){
     composite_tree<-sample(composite_tree, size = max_trees)
   }

   for (nn in 1:length(composite_tree)){
    if(any(composite_tree[[nn]]$edge.length<0)){
     current_nodedates <- suppressMessages(dateNodes(composite_tree[[nn]]))
     datediff<-tibble(anc_node = composite_tree[[nn]]$edge[,1],
                      desc_node =  composite_tree[[nn]]$edge[,2],
                      age_anc_node = as.integer(current_nodedates[composite_tree[[nn]]$edge[,1]]),
                      age_desc_node = as.integer(current_nodedates[composite_tree[[nn]]$edge[,2]])) %>%
       mutate(., age_offset = age_anc_node - age_desc_node)
     problem_nodes<-datediff$anc_node[which(datediff$age_offset<0)]
     problem_node_offsets<-datediff$age_offset[which(datediff$age_offset<0)]

     for (ii in 1:length(problem_nodes)){

       composite_tree[[nn]]$edge.length[which(composite_tree[[nn]]$edge[,1]==problem_nodes[ii])] <- composite_tree[[nn]]$edge.length[which(composite_tree[[nn]]$edge[,1]==problem_nodes[ii])] + (1-problem_node_offsets[ii])
       composite_tree[[nn]]$edge.length[which(composite_tree[[nn]]$edge[,2]==problem_nodes[ii])] <- composite_tree[[nn]]$edge.length[which(composite_tree[[nn]]$edge[,2]==problem_nodes[ii])] - (1-problem_node_offsets[ii])
       #time_travelers <- c(time_travelers, fossil_species_data$Tip_Label[i])

     }
   }
   }
    ####MAYBE CAN SPEED THIS UP WITH PARALLEL FOREACH OR A MCLAPPLY?
    if(one_topology){
      temp.trees<-paste_tips(tree=composite_tree[[1]],tipname =  tip_name_to_paste , node_number = target_node,tip_age = as.numeric(fossil_species_data$Tip_Date[i]))
      partially.finished.trees <- temp.trees
    } else {
      for(j in 1:length(composite_tree)){
        target_node_date <- suppressMessages(dateNodes(composite_tree[[j]])[target_node])
       #removed
        #if (target_node_date>fossil_species_data$Tip_Date[i]){
          temp.trees<-paste_tips(tree=composite_tree[[j]],tipname =  fossil_species_data$Tip_Label[i], node_number = target_node,tip_age = as.numeric(fossil_species_data$Tip_Date[i]))
          if (j==1){
            partially.finished.trees<-temp.trees
          } else {
            partially.finished.trees<-c(partially.finished.trees,temp.trees)
          }
        #removed with if above }
      }
    }
    composite_tree<-partially.finished.trees
    #tree_to_animate<-c(tree_to_annimate,composite_tree[[1]])
   if(verbose){
     writeLines(paste0(fossil_species_data$Tip_Label[i], " added, now ", Ntip(composite_tree[[1]]), " tips"))
   }
    }
  return(composite_tree)
}
