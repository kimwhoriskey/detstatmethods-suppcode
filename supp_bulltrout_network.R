# Supporting information for Whoriskey et al. Statistical methods for detection data.

# Analyzing the bull trout dataset with network analysis. 
# Author: Kim Whoriskey

######## 
# packages
require(igraph)
require(vegan)


########
# load the detection data
dets <- read.csv("Jan_2011_BT.csv")
# get the unique fish ids
ids <- unique(dets$FISHID)


########
# load the receiver metadata and get the receiver locations
metadata <- read.csv("Kinbasket_Receiver_groups.csv")
# take only the locations of receivers where we have detections
receiver.locations <- data.frame(table(dets$Receiver))
names(receiver.locations) <- c("receiver", "freq")
posx <- numeric()
posy <- numeric()
for(i in 1:dim(receiver.locations)[1]){
  posx[i] <- metadata[metadata$receiver_sn == receiver.locations$receiver[i],]$x_plane
  posy[i] <- metadata[metadata$receiver_sn == receiver.locations$receiver[i],]$y_plane
}
receiver.locations <- data.frame(receiver.locations, posx, posy)
l <- cbind(receiver.locations$posx, receiver.locations$posy)

########
# function to get the network movements between the receivers
make_network_data <- function(dat){
  
  ids <- unique(dat$FISHID)
  
  from <- numeric() #the receiver the fish moved FROM
  to <- numeric() #the receiver the fish moved TO
  fish <- numeric()
  
  for(i in 1:length(ids)){
    subs <- dat[dat$FISHID==unique(dat$FISHID)[i],]
    # here rle calculates the runs of detections from one fish at a receiver
    # separate into from and to directions
    from <- append(from, rle(subs$Receiver)$values[-length(rle(subs$Receiver)$values)])
    to <- append(to, rle(subs$Receiver)$values[-1])
    # the fish id 
    fish <- append(fish, rep(unique(subs$FISHID), length(rle(subs$Receiver)$values)-1))
  }
  
  # movements between receivers
  individual.moves <- data.frame(from, to, fish) 
  
  # gives the pairwise counts of all fish moving from one receiver to another, in matrix form
  moves.matrix <- table(individual.moves[,1:2])
  
  # summarizes the pairwise counts
  moves <- data.frame(moves.matrix) 
  moves <- moves[moves$Freq!=0,]
  
  return(list(moves=moves, moves.matrix=moves.matrix, individual.moves=individual.moves))
}


# get the network data for the full dataset
full.network <- make_network_data(dat=dets)
full.graph <- graph_from_data_frame(full.network$moves, directed=TRUE,
                               vertices=receiver.locations)

# plot the full network, where the node size is proportional to node degree, 
# and the thickness of the edges is proportional to the number of moves between receivers.
plot(full.graph, edge.arrow.size=.05, arrow.mode=1, edge.width=full.network$moves$Freq,
     vertex.size=degree(full.graph, mode="all"), vertex.label=NA, vertex.color="lightseagreen",
     layout=l)

# get the network data for females and males
males.dets <- dets[dets$Sex=="m",]
males.network <- make_network_data(dat=males.dets)
males.graph <- graph_from_data_frame(males.network$moves, directed=TRUE,
                                 vertices=receiver.locations$receiver)

females.dets <- dets[dets$Sex=="f",]
females.network <- make_network_data(dat=females.dets)
females.graph <- graph_from_data_frame(females.network$moves, directed=TRUE,
                                 vertices=receiver.locations$receiver)

# plot the two networks side by side, again where the node size is proportional to node degree,
# and the thickness of the edges is proportional to the number of moves between receivers. 
par(mfrow=c(1,2))
plot(males.graph, edge.arrow.size=.2, arrow.mode=2, edge.width=males.network$moves$Freq,
     vertex.size=degree(males.graph, mode="all"), vertex.label=NA, vertex.color="lightseagreen",
     layout=l)
plot(females.graph, edge.arrow.size=.2, arrow.mode=2, edge.width=females.network$moves$Freq,
     vertex.size=degree(females.graph, mode="all"), vertex.label=NA, vertex.color="tomato2",
     layout=l)



########
# use a Mantel test to see if there is a significant correlation between male and female 
# movement/use in the reservoir. 

# have to create two adjacency matrices 
# the adjacency matrix should be symmetric (dissimilarity matrix)
# so where we had directed movements before, now we have undirected
m = males.network$individual.moves
f = females.network$individual.moves

# all receivers that had detections on them
all.receivers = unique(c(f$from, f$to, m$from, m$to))

# function to get the matrices
get_adj_matx = function(dat, all.receivers){
  adj.matx = as.data.frame(matrix(0, nrow=length(all.receivers), ncol=length(all.receivers), dimnames=list(all.receivers, all.receivers)))
  for(j in 1:length(all.receivers)){
    for(k in j:length(all.receivers)){
      sub1 = dat[(dat$from==all.receivers[j]) & (dat$to==all.receivers[k]),]
      sub2 = dat[(dat$from==all.receivers[k]) & (dat$to==all.receivers[j]),]
      adj.matx[as.character(all.receivers[j]),as.character(all.receivers[k])] = dim(sub1)[1] + dim(sub2)[1]
      adj.matx[as.character(all.receivers[k]),as.character(all.receivers[j])] = dim(sub1)[1] + dim(sub2)[1]
      # index based on receiver name to be sure it's the right spot
    }
  }
  return(adj.matx)
}

# now get the matrices
m.adj.matx = get_adj_matx(dat=m, all.receivers)
f.adj.matx = get_adj_matx(dat=f, all.receivers)

# now perform the test using package vegan, with both pearson's and spearman's correlation. 
mantel(xdis=m.adj.matx, ydis=f.adj.matx, method="pearson", permutations=999) #pearson's correlation
mantel(xdis=m.adj.matx, ydis=f.adj.matx, method="spearman", permutations=999) #spearman's correlation












