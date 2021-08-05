library(paleotree)
library(geiger)
library(tidyverse)
library(treeio)
library(ggtree)
library(adephylo)
library(abind)
library(phangorn)
library(ips)

#load functions:

source("./Grafting_trees/grafting functions.R")

#load trees
#molecular tree from Upham et al 2019
good_molec_tree<-read.nexus("./Raw_Data/trees/Upham_tree_FBD_extantsub.nex")
#tree that contains all the fossils in the dataset
mammaltree1 <- read.nexus("./Raw_Data/trees/Bayesresults2b.con.tre")
#remove "Balaena_mysticetus"  becuase its not in the dataset
mammaltree1 <- drop.tip(mammaltree1, "Balaena_mysticetus")
my_species_data <- read_csv("./Raw_Data/full_species_data.csv")

#posterior distribution of trees from Upham et al downloaded from treelife.org
treefolder <- "/Users/rnf/Downloads/Completed_5911sp_topoCons_NDexp/"
treefiles <- dir(treefolder)
mytrees <- lapply(1:10000, function(x) read.tree(paste0(treefolder,treefiles[x])))
mytrees2 <- lapply(1:length(mytrees), function(x) keep.tip(mytrees[[x]], good_molec_tree$tip.label))

distances <- lapply(1:length(mytrees2), function(x) RF.dist(good_molec_tree, mytrees2[[x]]))

mytrees3<-mytrees2[which(distances==0)]
class(mytrees3) <- "multiPhylo"

dates<-lapply(1:length(mytrees3), function(x) max(dateNodes(mytrees3[[x]])))
dates2<-tibble(Root_Age=unlist(dates))
ggplot(data=dates2, aes(dates2$Root_Age)) +
  geom_histogram(breaks = seq(70,115,by=5))

mammaltree2 <- drop.tip(mammaltree1, c("Palaeoprionodon_lamandini"))
correct_node<-getMRCA(mammaltree2, c("Acinonyx_jubatus" , "Nandinia_binotata"))
correct_node
mammaltree2<-paste_tips(mammaltree2, "Palaeoprionodon_lamandini", correct_node, 31.7)

trees_per_bin<-100
upham_70_to_75<-which(dates2$Root_Age>70 & dates2$Root_Age<75)
upham_75_to_80<-which(dates2$Root_Age>75 & dates2$Root_Age<80)
upham_80_to_85<-which(dates2$Root_Age>80 & dates2$Root_Age<85) %>% sample(., size=trees_per_bin)
upham_85_to_90<-which(dates2$Root_Age>85 & dates2$Root_Age<90) %>% sample(., size=trees_per_bin)
upham_90_to_95<-which(dates2$Root_Age>90 & dates2$Root_Age<95) %>% sample(., size=trees_per_bin)
upham_95_to_100<-which(dates2$Root_Age>95 & dates2$Root_Age<100) %>% sample(., size=trees_per_bin)

grafted_trees_70_to_75<-lapply(1:length(upham_70_to_75), function(x) tree_grafter(tree_w_fossils=mammaltree2,     target_tree=mytrees3[[upham_70_to_75[x]]], species_data=my_species_data, max_trees=trees_per_bin))
grafted_trees_75_to_80<-lapply(1:length(upham_75_to_80), function(x) tree_grafter(tree_w_fossils=mammaltree2,     target_tree=mytrees3[[upham_75_to_80[x]]], species_data=my_species_data, max_trees=trees_per_bin))
grafted_trees_80_to_85<-lapply(1:length(upham_80_to_85), function(x) tree_grafter(tree_w_fossils=mammaltree2,     target_tree=mytrees3[[upham_80_to_85[x]]], species_data=my_species_data, max_trees=trees_per_bin))
grafted_trees_85_to_90<-lapply(1:length(upham_85_to_90), function(x) tree_grafter(tree_w_fossils=mammaltree2,     target_tree=mytrees3[[upham_85_to_90[x]]], species_data=my_species_data, max_trees=trees_per_bin))
grafted_trees_90_to_95<-lapply(1:length(upham_90_to_95), function(x) tree_grafter(tree_w_fossils=mammaltree2,     target_tree=mytrees3[[upham_90_to_95[x]]], species_data=my_species_data, max_trees=trees_per_bin))
grafted_trees_95_to_100<-lapply(1:length(upham_95_to_100), function(x) tree_grafter(tree_w_fossils=mammaltree2,   target_tree=mytrees3[[upham_95_to_100[x]]], species_data=my_species_data, max_trees=trees_per_bin))
