# Handling and Analysing Phylogenetic Tree's in R

#SET UP

library(ape)

#MAKING PRETEND IDEALISED TREES

star_tree <- stree(n=10, type = "star")
#Creates a tree with 10 sp/groups in the type star, see notes of for the other types
#of tree structure and what they mean. 

plot(star_tree)
#tree where ALL the tips come from one common ancestor, all speciatied at the same
#time ??
#R uses plot to work out the correct specific function to use to plot it, e.g. in this
#scenario it finds the plot.phylofunction from APE. 

my_tip_labels <- paste("Species", 1:8)
#each , creates a new word, if want them separated by a different character use 
#sep = "character"

right_imbalanced_tree <- stree(n=8, type="right", tip.label= my_tip_labels)
#tip.label allows you to assign the names of the tips
plot(right_imbalanced_tree)

balanced_tree <- stree(n=8, type="balanced", tip.label= my_tip_labels)
plot(balanced_tree)

#HOW R ENCODES INFO

balanced_tree #gives the summary of properties R assigned this object e.g. 
# a phylo tree with 8 tips and 7 internal nodes, tip labels = ..., rooted with no branch lengths

str(balanced_tree)
#shows what the object is actually made up of in terms of R + coding shit

balanced_tree$edge #specifically printing a summary of the topology (that is what edge looks at)
#in this case edges = the branches, numbers are those assigned to the tip and 
#internal node of each branch. 
#Tips are numbered 1-8, the internal nodes (each branching bit) are 9-15.
#first column gives ancestor and second gives descendant

plot(balanced_tree) 
nodelabels() #adds labels to the internal nodes, but wouldnt usually do this.

#Pulling entries out of a matrix
#uses the format [row,column]

i <- 1
balanced_tree$edge[i,]
#extracts branch i ( in this case 1) from the table, has the descendant 10 and 
#the ancestor 9, is the first branch down on the right of the tree.
balanced_tree$edge[i,1]
#This tells us specifically the ancestor of branch i (1)
balanced_tree$edge[i,2]
#This tells us specifically the descendant of branch i (1)

balanced_tree$edge[which(balanced_tree$edge[,2]==i),1]
#find the branch that leads to node/tip i
#gives you the ancestor ,1 in the column where the descendant is i

balanced_tree$edge[which(balanced_tree$edge[,2]==which(balanced_tree$tip.label=="Species 5")),1]
#gives you the ancestor of a specifically named node

#ROTATING INTERNAL NODES

#display order of the two branches does not convey info when topology and
#branch lengths remain the same so can rearranged stuff to make patterns more obvious.
# Rotate the internal nodes using rotate

for(each.node in c(9:15)) {
  balanced_tree <- rotate(balanced_tree, node=each.node)
}

plot(balanced_tree)
nodelabels()
#reorganises so they are in species numerical order

#RANDOM TREE'S

random_coalescent <- rcoal(n=8, tip.label=my_tip_labels)
plot(random_coalescent)
#creates a different random tree each time, is an ultrametric tree though so all
#sp have been evolving for the same amount of time. 

random_coalescent$edge
random_coalescent$tip.label
#because it is a random tree need to also work out which row is which species tip.
str(random_coalescent)

#PREPARING AND MANIPULATING TREES

random_tree <- rtree(n=8, tip.label = my_tip_labels)
plot(random_tree)
#rtree makes a random NON+ULTRAMETRIC tree

is.ultrametric(random_tree)
is.ultrametric(random_coalescent)
#checks if a tree is ultrametric

ultrametric_tree <- chronos(random_tree)
#chronos is going to make an non-ultra tree into an ultra tree
plot(ultrametric_tree)

#also ways to include fossil or other calibrations and do more sophisticated transforming
#but this is a simple way to re-estimate branch lengths.

#Can extract branch lengths from ultrametric tree's using the APE function branching.times(tree)
branching.times(ultrametric_tree)

#Athough ideal tree is fully bifurcating in reality often we are not able to resolve
#to this level and end up with polytomies. 
#
#A hard polytomy is a specified node and will 3 or more brances, while a soft polytomy
#may appear topologically bifurcating but some some branch lengths are so clsoe to zero
#it looks like a polytomy when plotted. Many functions assume a fully resolved tree.

is.binary(star_tree) #Used to check if a tree is fully dichotomous.

#We can resolve trees with multichomy at the root, like the start tree.
new_star_tree <- multi2di(star_tree)
plot(new_star_tree)

#COSMETICS FOR PLOTTING TREES

par(mfrow=c(1,2))
ladderized_tree <- ladderize(random_coalescent, right = TRUE) #ladderize reorganises the internal
#tree structure to get a ladderise affect when plotted. Right = specifices whether the
#smallest clade is on the left or right.
ladderized_tree2 <- ladderize(random_coalescent, right = FALSE)
plot(ladderized_tree)
plot(ladderized_tree2)
axisPhylo() #adds a scaled axis to the phylogenetic trees if they are ultrametric

#if non-ultrametric can add scale bar using add.scale.bar()

#plot.phylo can also be used to plot other types of tree's like unrooted, radial and cladograms.
#look into help(plot.phylo) for more info on how to customise how the trees look.

#READING IN REAL TREES AND MATCHING DATA

#Standard phylo software will save trees in the Newick format, uses () to denote topology
#of the tree, and branch lengths are specified using :length after taxon name for terminal
#branches. 
#Can read trees like this directly into R

tr <- read.tree(text = "(((A:1.5,B:1.5):1,C:2.5):2,D:4.5);")
plot(tr) #start with the innermost/most recently diverged taxon, e.g. A is most recent 
# in this tree. 

#Can usually read in trees from files though. 
fusarium_tree <- read.tree(file = "Ffuji.tre")
plot(fusarium_tree)
axisPhylo()

#Some tree files might have multiple trees in them, use [[number]], to access certain 
#trees on the list. Some trees might also be saved as nexus files so you read
# in using read.nexus. 

#Matching tree with data

#Sometimes you will want to map phylogenetic data onto a tree. 

fusarium_data <- read.csv("Ffuji.data.csv", header = TRUE)
head(fusarium_data)
#The names are written in a different way to the tree, we can fix this though. 

new_names <- gsub("_", "", fusarium_data$Taxon) #This is telling R to substitute _ for 
#a blank space, this lets us get rid of the underscores in the names in Fusarium_data. 

new_names [1:5] <- substring(new_names[1:5],1,4) #Make the first five names only 4 letters long
#and input them into the new names variable.

fusarium_data$Taxon <- new_names

#Also have the problem that the names are not in the same order as on the tree. 
#This doesnt matter for all function, but better to be safe than sorry to stop stuff
#being read in in the wrong order. 

pmatch(fusarium_tree$tip.label, fusarium_data$Taxon)
#Finds matches between characters in one vector and another.
#Output tell us what position each element of the first vector occupies in the second, 
#use this to to sort the table into the same order as the tree. 

sorted_fusarium_data <- fusarium_data[pmatch(fusarium_tree$tip.label, fusarium_data$Taxon),]

#ANALYSING THE PROPERTIES OF TREES

#Total branch length

sum(fusarium_tree$edge.length) #Gives you the total branch length
#Total branch length is 76.40981

#Can even remove species and recalculate

fusarium_tree_trimmed <- drop.tip(fusarium_tree, c("Fver", "Fpro"))
plot(fusarium_tree_trimmed)
axisPhylo()

sum(fusarium_tree_trimmed$edge.length)
# New total branch length is 56.81702.

#Can also remove tips randomly
fusarium_tree_random <- drop.tip(fusarium_tree, c(fusarium_tree$tip.label[sample (1:11, 3)]))
plot(fusarium_tree_random)

#TO estimate average loss would make a for loop.

#To calculate the amount of unique evo history of a given tip is the length of the terminal 
#branch leading to it.

terminal_branch <- which(fusarium_tree$edge [,2]== i)
fusarium_tree$edge.length[terminal_branch]
#8.028423

#How imbalanced is the tree ?

#Imbalance can suggest certain sub-clades are more diverse, and also that there might
#be progressive branching. This is when the next branching even is most likely to 
#be on the descendant of the most recent branch rather than any other branch in tree.

#balance function counts the number of descendants of each of the 2 clades coming
#off the original split.

balance_table <- balance(random_coalescent)

#now we want to find the average balance across the tree.

#use a for loop to sort each row so the most divers clade is in column 1

for(i in (1:nrow(balance_table))) {
  if (balance_table [i,1] < balance_table [i,2]) {
    balance_table[i,1:2] <- balance_table[i,2:1]}
}

