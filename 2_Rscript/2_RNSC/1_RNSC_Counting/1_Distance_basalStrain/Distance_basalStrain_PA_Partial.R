############
# PACKAGES #
############
#install.packages("devtools")
#devtools::install_github("thibautjombart/OutbreakTools")

#require(OutbreakTools)
library(treeio)
library(ape)
library(maps)
library(phytools)

##########
# SCRIPT #
##########
#####################
# 1. fetch the data #
#####################
baseDir <- "H3N2/1_Data/2_RNSC/2_Treelog/1_Trees_basalStrain"
#baseDir <- "~/Documents/vaccineLibrary/"
setwd(baseDir)
segments <- c("PA")
segment.index <- 1
currentSegment <- segments[segment.index]
# A. the trees.   ATTENTION: this assumes that the burnin phase trees are not present any more! 
#                 The naming should end at '.noBI.trees'. 
# A test with 50 HA trees shows that it takes ~14 seconds to read this file:
# > system.time(read.beast(file = treeFiles[currentSegment.treeFile.index]))
# user  system elapsed 
# 14.187   0.054  14.258 
# hence, it is easiest to just read all trees into memory and subset afterwards if required.
treeFiles <- list.files(".", pattern="\\_noBI.trees$")
currentSegment.treeFile.index <- grep(pattern = paste0("_", currentSegment, "_"), x = treeFiles)
trees <- read.mcmctree(file = treeFiles[currentSegment.treeFile.index])
nrLookedAtTrees <- length(trees)

# collect tip info in list of data frames 
vaccinePeriods <- 1:2
##########
# for PA #
##########
inFileNames <- c("9_PA_HK15_195seq.fasta","10_PA_Kan17_18seq.fas")
vaccineStrains <- c("EPI686117_A_Hong_Kong_15611_2015",
                    "MG974447_A_Kansas_14_2017")
vaccineStrains <- sapply(X = vaccineStrains, function(x) paste(strsplit(x = x, split = "_")[[1]][1:9], collapse = "_"))

####################################
# 2. loop over the vaccine periods #
####################################
for (current.vaccinePeriod in 1:2){
  
  vaccineStrain <- vaccineStrains[current.vaccinePeriod]
  #the tips to be contrasted:
  
  setwd("/1_Data/2_RNSC/3_DistanceLog/1_Before2020_601/PA")
  fasta = scan(inFileNames[current.vaccinePeriod], what="", sep="\n", quiet=TRUE)
  seq_IDs = fasta[which(grepl(">",fasta))]
  seq_IDs <- gsub('>', '', seq_IDs)
  #avoid contrasting the vaccine strain against itself:
  tipNameVector <- seq_IDs[which(seq_IDs != vaccineStrain)]
  # remove accession number from tip labels 
  tipNameVector <- sapply(X = tipNameVector, function(x) paste(strsplit(x = x, split = "_")[[1]][1:9], collapse = "_"))
  setwd(baseDir)
  
  #############################
  # 3. prepare saving results #
  #############################
  loggedVariables <- c("tree.nr", "N", "S", "vaccineStrain", "compareStrain", "tipDate")
  distancePerSiteDF <- as.data.frame(matrix(NA,  ncol = length(loggedVariables), #also record vaccine period 
                                            nrow = nrLookedAtTrees*length(tipNameVector) ))
  names(distancePerSiteDF) <- loggedVariables
  treeNrColumnIndex <- length(distancePerSiteDF[1,]) - 5
  nonSynonColumnIndex <- length(distancePerSiteDF[1,]) - 4
  synonColumnIndex <- length(distancePerSiteDF[1,]) - 3
  vaccineStrainColumnIndex <- length(distancePerSiteDF[1,]) - 2
  compareStrainColumnIndex <- length(distancePerSiteDF[1,]) - 1
  tipDateColumnIndex <- length(distancePerSiteDF[1,])
  
  focus <- paste0(currentSegment, "_G", current.vaccinePeriod, "_distances.txt")
  outFile2 <- paste("../",focus, sep="")
  write.table(x = as.matrix(t(loggedVariables)), file = outFile2, append = FALSE, sep = "\t", 
              row.names = FALSE, col.names=FALSE)
  
  ############################
  # 4. the actual looking up #
  ############################
  
  for (t in 1:length(trees)){
    current.tree <- trees[[t]]
    rowIndexPerSite <- (t - 1) * length(tipNameVector)
    # remove accession number from tip labels 
    current.tree@phylo$tip.label <- sapply(X = current.tree@phylo$tip.label, function(x) unname(paste(strsplit(x = x, split = "_")[[1]][1:9], collapse = "_")))
    
    for (i in 1:length(tipNameVector)){ 
      currentTip <- as.character(tipNameVector[i])
      #print(currentTip)
      rowIndexPerSite <- rowIndexPerSite + 1
      #tipDate is always last element after '_'
      tipdate <- strsplit(currentTip, "_")[[1]][length(strsplit(currentTip, "_")[[1]])]
      if (tipdate == "NA"){
        tipdate <- strsplit(currentTip, "_")[[1]][length(strsplit(currentTip, "_")[[1]])-1]
      }
      mrca <- phytools::findMRCA(tree = current.tree@phylo, tips = c(currentTip, vaccineStrain), type = "node")
      
      #follow path from mrca to tip and mrca to vaccineStrain in edge matrix
      ### path from tip to mrca ###
      current.parent.node <- which(current.tree@phylo$tip.label == currentTip)
      nodePath.mrca.to.tip <- vector()
      mrcaEncountered <- FALSE
      while (mrcaEncountered == FALSE){
        current.child.node <- current.parent.node
        nodePath.mrca.to.tip <- append(x = nodePath.mrca.to.tip, values = current.child.node, after = length(nodePath.mrca.to.tip))
        current.parent.node <- current.tree@phylo$edge[which(current.tree@phylo$edge[,2] == current.child.node), 1]
        if ( current.parent.node == mrca ){
          mrcaEncountered <- TRUE
        }  
      } #END: while
      ### path from tip to vaccineStrain ###
      mrcaEncountered <- FALSE
      current.parent.node <- which(current.tree@phylo$tip.label == vaccineStrain)
      nodePath.mrca.to.vaccinestrain <- vector()
      while (mrcaEncountered == FALSE){
        current.child.node <- current.parent.node
        nodePath.mrca.to.vaccinestrain <- append(x = nodePath.mrca.to.vaccinestrain, values = current.child.node, after = length(nodePath.mrca.to.vaccinestrain)) 
        current.parent.node <- current.tree@phylo$edge[which(current.tree@phylo$edge[,2] == current.child.node), 1]
        if ( current.parent.node == mrca ){
          mrcaEncountered <- TRUE
        }  
      } #END: while
      
      ####################
      # fetch trait info #
      ####################
      ### non synon distance ###
      mrcaToTip <- 0
      for (j in 1:length(nodePath.mrca.to.tip)){
        current.child.node <- nodePath.mrca.to.tip[j] 
        traitIndex <- which(current.tree@data$node == current.child.node) 
        mrcaToTip <- mrcaToTip + unname(current.tree@data$N[traitIndex])
      } #END: j
      mrcaToVaccine <- 0
      for (j in 1:length(nodePath.mrca.to.vaccinestrain)){
        current.child.node <- nodePath.mrca.to.vaccinestrain[j] 
        traitIndex <- which(current.tree@data$node == current.child.node) 
        mrcaToVaccine <- mrcaToVaccine + unname(current.tree@data$N[traitIndex])
      } #END: j 
      total.non.synon.dist <- mrcaToTip + mrcaToVaccine   
      ### synon distance ###
      mrcaToTip <- 0
      for (j in 1:length(nodePath.mrca.to.tip)){
        current.child.node <- nodePath.mrca.to.tip[j] 
        traitIndex <- which(current.tree@data$node == current.child.node) 
        mrcaToTip <- mrcaToTip + unname(current.tree@data$S[traitIndex])
      } #END: j
      mrcaToVaccine <- 0
      for (j in 1:length(nodePath.mrca.to.vaccinestrain)){
        current.child.node <- nodePath.mrca.to.vaccinestrain[j] 
        traitIndex <- which(current.tree@data$node == current.child.node) 
        mrcaToVaccine <- mrcaToVaccine + unname(current.tree@data$S[traitIndex])
      } #END: j 
      total.synon.dist <- mrcaToTip + mrcaToVaccine   
      
      #log all info:
      distancePerSiteDF[rowIndexPerSite,treeNrColumnIndex] <- t
      distancePerSiteDF[rowIndexPerSite,nonSynonColumnIndex]  <- total.non.synon.dist
      distancePerSiteDF[rowIndexPerSite,synonColumnIndex]  <- total.synon.dist
      distancePerSiteDF[rowIndexPerSite,vaccineStrainColumnIndex] <- vaccineStrain
      distancePerSiteDF[rowIndexPerSite,compareStrainColumnIndex] <- currentTip
      distancePerSiteDF[rowIndexPerSite,tipDateColumnIndex] <- tipdate
    } #end i - for each tip 
  } # end loop t over trees
  write.table(x = distancePerSiteDF[1:rowIndexPerSite,], file = outFile2, 
              append = TRUE, sep = "\t", row.names = FALSE, col.names=FALSE)
}



