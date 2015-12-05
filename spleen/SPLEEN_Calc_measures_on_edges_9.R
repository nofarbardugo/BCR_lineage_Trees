# set directory : change path depending on work environment
source.dir <-'/home/bardugn1/spleen' # pep4

setwd(source.dir)

getEdgeName <-function (father,son)
  # this function get a father and son nodes and return its edges' name
  # Args: 
  #   father- a father node in the tree
  #   son - the son node of the father in the tree
  # Returns: sutable edge name 
{

  fatherIsotype <- isotypesNames.df[which(isotypesNames.df[,"head"] == father),"IsotypeName"]
  sonIsotype <- isotypesNames.df[which(isotypesNames.df[,"head"] == son),"IsotypeName"]
  
  combine <- paste0(fatherIsotype,sonIsotype)
  return(switch(combine,"igMminusigMminus" = "minus-minus","igMplusigMminus" = "plus-minus", "igMplusigMplus" = "plus-plus"))
                
}

GoOverTree<-function(node)
  # this function goes over the tree and save for each node type how much directes kids it has
  # Args: 
  #   node- an edge of the tree   
  # Returns: 
{
  
  # if condition-   what we work on[the condition, if the condition is true give me
  #                                                                  the Relevant information]
  # get all the suns
  suns <- edges.df [edges.df [,"from"]==node,"to"]
  
  # save suns num 
  sunsNum <- length(suns)
  
  # if its not a leaf
  if(sunsNum > 0)
  {
    list  <- allEdgeDetails
    
    p<-1
    
    while(p <=sunsNum)
    {
      GoOverTree(suns[p])
      
      # get edges typ 
      type <-getEdgeName(node,suns[p])
      
      if(suns[p] !="VJ_1")
      {
        syn<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms"]
        nonS<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms"]
        
        frSyn<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms_FR"]
        frNonS<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_FR"]
        
        cdrSyn<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms_CDR"]
        cdrNonS<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_CDR"]

        list[[type]] <- rbind(allEdgeDetails[[type]],c(syn,nonS,frSyn,frNonS,cdrSyn, cdrNonS)) 
       
      }
      p <- p+1

    }
    
    # save "list" in the outer "allEdgeDetails" data frame
    assign("allEdgeDetails",list, envir = .GlobalEnv)
  }
  
}  


############## main ###########################

EDGE_NUM <- 3
# path for trees folder
tree.dir <- paste0(source.dir,"/spleenHumanTree/")

# get list of fold
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)

# FOR ALL PATIANTs TOGTHER ------ data fream that save details about each tree edge mean Synonymus mutataion - 
fd_ALL.TotalEdgeSynonymus<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c(),"place" = c())

# FOR ALL PATIANTs TOGTHER ------data fream that save details about each tree edge mean NON Synonymus mutataion - 
fd_ALL.TotalEdgeNonSynonymus<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c(),"place" = c())

# FOR ALL PATIANTs TOGTHER ------ data fream that save details about each tree edge mean Synonymus CDR mutataion - 
fd_ALL.TotalEdgeSynonymus_CDR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c(),"place" = c())

# FOR ALL PATIANTs TOGTHER ------ data fream that save details about each tree edge mean Synonymus FR mutataion - 
fd_ALL.TotalEdgeSynonymus_FR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c(),"place" = c())

# FOR ALL PATIANTs TOGTHER ------data fream that save details about each tree edge mean NON Synonymus CDR mutataion - 
fd_ALL.TotalEdgeNonSynonymus_CDR<-data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c(),"place" = c())

# FOR ALL PATIANTs TOGTHER ------data fream that save details about each tree edge mean NON Synonymus FR mutataion - 
fd_ALL.TotalEdgeNonSynonymus_FR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c(),"place" = c())

# FOR ALL PATIANTs TOGTHER ------data fream that save details about each tree edge mean CDR mutataion - 
fd_ALL.TotalEdgeCDR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c(),"place" = c())

# FOR ALL PATIANTs TOGTHER ------data fream that save details about each tree edge mean FR mutataion - 
fd_ALL.TotalEdgeFR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c(),"place" = c())


# FOR ALL PATIANTs TOGTHER ------ data fream that save details about each tree edge mean depth mutataion -(s+ns)
fd_ALL.TotalEdgeDepth<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c(),"place" = c())

# FOR ALL PATIANTs TOGTHER ------data fream that save details about each tree edge mean ratio mutataion - (ns/(s+ns))
fd_ALL.TotalEdgeRatio<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c(),"place" = c())

# FOR ALL PATIANTs TOGTHER ------data fream that save details about each tree edge mean ratio between CDR to FR - (CDR/(FR+CDR))
fd_ALL.TotalEdgeRatio_CDRvsFR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c(),"place" = c())


# FOR ALL PATIANTs TOGTHER ------data fream that save details about each tree edge mean ratio between CDR (ns/(s+ns))
fd_ALL.TotalEdgeRatio_CDR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c(),"place" = c())

# FOR ALL PATIANTs TOGTHER ------data fream that save details about each tree edge mean ratio between FR (ns/(s+ns))
fd_ALL.TotalEdgeRatio_FR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c(),"place" = c())


NumberOfTreesCount <- 0

# going over each patient
for(k in 2:length(patient.dirs))
{
  # decleare her to get detailes for each patient, can move out in order to get details for all patient togther
  
  # data fream that save details about each tree edge mean Synonymus mutataion - 
  fd.TotalEdgeSynonymus<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c())
  
  # data fream that save details about each tree edge mean NON Synonymus mutataion - 
  fd.TotalEdgeNonSynonymus<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c())
  
  # data fream that save details about each tree edge mean Synonymus CDR mutataion - 
  fd.TotalEdgeSynonymus_CDR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c())
  
  # data fream that save details about each tree edge mean Synonymus FR mutataion - 
  fd.TotalEdgeSynonymus_FR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c())
  
  # data fream that save details about each tree edge mean NON Synonymus CDR mutataion - 
  fd.TotalEdgeNonSynonymus_CDR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c())
  
  # data fream that save details about each tree edge mean NON Synonymus FR mutataion - 
  fd.TotalEdgeNonSynonymus_FR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c())
  
  # data fream that save details about each tree edge mean CDR mutataion - 
  fd.TotalEdgeCDR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c())
  
  # data fream that save details about each tree edge mean FR mutataion - 
  fd.TotalEdgeFR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c())
  
  # data fream that save details about each tree edge mean depth mutataion -(s+ns)
  fd.TotalEdgeDepth<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c())
  
  # data fream that save details about each tree edge mean ratio mutataion - (ns/(s+ns))
  fd.TotalEdgeRatio<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c())
  
  # data fream that save details about each tree edge mean ratio between CDR to FR - (CDR/(FR+CDR))
  fd.TotalEdgeRatio_CDRvsFR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c())
  
  
  # data fream that save details about each tree edge mean ratio between CDR (ns/(s+ns))
  fd.TotalEdgeRatio_CDR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c())
  
  # data fream that save details about each tree edge mean ratio between FR (ns/(s+ns))
  fd.TotalEdgeRatio_FR<- data.frame("minus-minus"=c() ,"plus-minus" = c(),"plus-plus" =c())
  
  
  
  # get list of folders
  dirs <- list.files(paste0(tree.dir,patient.dirs[k]),full.names = F)
  
  # going over each tree folder`
  for(i in 1:length(dirs))
  {
    dirPath <- paste0(tree.dir,patient.dirs[k],'/',dirs[i],'/',dirs[i])
    filePath <-paste0(dirPath,"_edgesCut.tab")
    
    # case file exist (if there is not enough or too much sequences, wile be only fasta file)
    if(file.exists(filePath))
    {
      
      # read the edges table into data frame  
      edges.df <- tryCatch(read.table(file = filePath ,header = T,stringsAsFactors = F,sep ="\t"), error=function(e) NULL)
      
      if(!is.null(edges.df))
      {
        # get details about current tree
        totalLeavesType.fd <-tryCatch(read.table(file = paste0(dirPath,"_totalLeavesType.tab")
                                                 ,header = T,stringsAsFactors = F,sep ="\t"),
                                      error=function(e) NULL)
        
        # get table with all node's isotypes names
        isotypesNames.df <- read.table(file= paste0(dirPath,"_isotypeNames.tab")
                                       ,header = T,stringsAsFactors = F,sep ="\t")
        
        # get a root for each clone
        roots <- setdiff(edges.df[,"from"],edges.df[,"to"])
        
        # goes over each root and send the root to "GoOverTree" function 
        # sapply(roots,GoOverTree)
        j<-1
        
        while(j<=length(roots))
        { 
          
          # count modes only on trees that have at last 10 leave  
          if(!(is.null(totalLeavesType.fd))  && (totalLeavesType.fd[j,"sum_Leaves"]>10) )
          {
              
            # for count how many treee ther is
            NumberOfTreesCount <- NumberOfTreesCount +1
            print("in")
            # create datafream for each edge that will save: Synonymus, and nonSynonymus
            # save all data freams into a list for save all togther
            allEdgeDetails <- list()
        
            allEdgeDetails[["minus-minus"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),
                                                          cdrSyn = c(),cdrNonSyn = c(),frSyn = c(),frNonSyn = c())
            allEdgeDetails[["plus-minus"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),
                                                         cdrSyn = c(),cdrNonSyn = c(),frSyn = c(),frNonSyn = c())
            allEdgeDetails[["plus-plus"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),
                                                        cdrSyn = c(),cdrNonSyn = c(),frSyn = c(),frNonSyn = c())
            
            GoOverTree(roots[j]) 
            
            AvgSYN.vector <- vector()
            AvgNON.vector <- vector()
            AvgDepth.vector <- vector()
            AvgRatio.vector <- vector()
            
            AvgSYN_CDR.vector <- vector()
            AvgSYN_FR.vector <- vector()
            AvgNON_CDR.vector <- vector()
            AvgNON_FR.vector <- vector()
            AvgCDR.vector <- vector()
            AvgFR.vector <- vector()
            AvgRatioNON_CDR.vector <- vector() # save the ratio between CDR to FR mutations
            AvgRatio_CDR.vector <- vector() # save the ratio between ns/s only on CDR mutations
            AvgRatio_FR.vector <- vector() # save the ratio between ns/s only on FR mutations
            
            # calc avg of each dataframe and add to the big dataFrame
            for(p in 1:EDGE_NUM)
            {
              if(length(allEdgeDetails[[p]]) > 0)
              {
                #syn/nonsyn
                AvgSYN.vector[p] <-  mean((allEdgeDetails[[p]])[,1]) # avg of the synonymus in this edge in the tree
                AvgNON.vector[p] <-  mean((allEdgeDetails[[p]])[,2]) # avg of the nonsynonymus in this edge in the tree
                #cdr/fr
                AvgSYN_CDR.vector[p] <-  mean((allEdgeDetails[[p]])[,3]) # avg of the synonymus in this edge in the tree
                AvgNON_CDR.vector[p] <-  mean((allEdgeDetails[[p]])[,4]) # avg of the nonsynonymus in this edge in the tree
                AvgSYN_FR.vector[p] <-  mean((allEdgeDetails[[p]])[,5]) # avg of the synonymus in this edge in the tree
                AvgNON_FR.vector[p] <-  mean((allEdgeDetails[[p]])[,6]) # avg of the nonsynonymus in this edge in the tree
                
                CDR <-  (allEdgeDetails[[p]])[,3]+ (allEdgeDetails[[p]])[,4] # avg of the CDR in this edge in the tree
                AvgCDR.vector[p] <- mean(CDR)
                
                FR <-  (allEdgeDetails[[p]])[,5]+ (allEdgeDetails[[p]])[,6] # avg of the FR in this edge in the tree
                AvgFR.vector[p] <- mean(FR)
                
                Depth <- (allEdgeDetails[[p]])[,1]+ (allEdgeDetails[[p]])[,2] # avg of the depth in this edge in the tree
                AvgDepth.vector[p] <-  mean(Depth) # avg of Depth

                # ratio 
                div <- c()
                divPlace <- c()
                divCDR <- c()
                divFR <- c()
                
                for(d in 1: length(Depth))
                {
                  if(Depth[d]!=0)
                  {
                    div[d] <- (allEdgeDetails[[p]])[d,2]/Depth[d]
                    divPlace[d] <- CDR[d]/Depth[d]
                    
                  }
                  else
                  {
                    div[d] <- 0
                    divPlace[d] <- 0
                  }
                  
                  if(CDR[d]!=0)
                  {
                    divCDR[d] <- (allEdgeDetails[[p]])[d,4]/CDR[d]
                  }
                  else
                  {
                    divCDR[d] <- 0
                  }
                  
                  if(FR[d]!=0)
                  {
                    divFR[d] <- (allEdgeDetails[[p]])[d,6]/FR[d]
                  }
                  else
                  {
                    divFR[d] <- 0
                  }

                }
                  
                AvgRatio.vector[p] <-  mean(div) # avg of Depth
                AvgRatioNON_CDR.vector[p] <- mean(divPlace)
                AvgRatio_CDR.vector[p] <- mean(divCDR)
                AvgRatio_FR.vector[p] <- mean(divFR)
              }
              else
              {
                AvgSYN.vector[p]=NA
                AvgNON.vector[p]= NA
                AvgDepth.vector[p]=NA
                AvgRatio.vector[p]= NA
                AvgSYN_CDR.vector[p] = NA
                AvgSYN_FR.vector[p] = NA
                AvgNON_CDR.vector[p] = NA 
                AvgNON_FR.vector[p] = NA
                AvgCDR.vector[p] = NA 
                AvgFR.vector[p] = NA
                AvgRatioNON_CDR.vector[p] = NA
                AvgRatio_CDR.vector[p] = NA
                AvgRatio_FR.vector[p] = NA
              }
           
            }
            
            # add new line to the big data freams
            fd.TotalEdgeSynonymus <- rbind(fd.TotalEdgeSynonymus,AvgSYN.vector)
            fd.TotalEdgeNonSynonymus <- rbind(fd.TotalEdgeNonSynonymus,AvgNON.vector)
            fd.TotalEdgeDepth <- rbind(fd.TotalEdgeDepth,  AvgDepth.vector)
            fd.TotalEdgeRatio <- rbind(fd.TotalEdgeRatio,AvgRatio.vector)

            fd.TotalEdgeSynonymus_CDR<- rbind(fd.TotalEdgeSynonymus_CDR,AvgSYN_CDR.vector)
            fd.TotalEdgeSynonymus_FR<-  rbind( fd.TotalEdgeSynonymus_FR,AvgSYN_FR.vector)
            fd.TotalEdgeNonSynonymus_CDR<- rbind(fd.TotalEdgeNonSynonymus_CDR,AvgNON_CDR.vector)
            fd.TotalEdgeNonSynonymus_FR<-  rbind(fd.TotalEdgeNonSynonymus_FR,AvgNON_FR.vector)
            fd.TotalEdgeCDR<-  rbind(fd.TotalEdgeCDR,AvgCDR.vector)
            fd.TotalEdgeFR<-  rbind(fd.TotalEdgeFR,AvgFR.vector)
            fd.TotalEdgeRatio_CDRvsFR<-  rbind(fd.TotalEdgeRatio_CDRvsFR,AvgRatioNON_CDR.vector)
            fd.TotalEdgeRatio_CDR <- rbind(fd.TotalEdgeRatio_CDR,AvgRatio_CDR.vector)
            fd.TotalEdgeRatio_FR <- rbind(fd.TotalEdgeRatio_FR,AvgRatio_FR.vector)

          }
          j<-j + 1
        }
      }
    }
  }


  fd.TotalEdgeSynonymus$place <- patient.dirs[k]
  fd.TotalEdgeNonSynonymus$place <- patient.dirs[k]
  fd.TotalEdgeDepth$place <- patient.dirs[k]
  fd.TotalEdgeRatio$place <- patient.dirs[k]
  fd.TotalEdgeSynonymus_CDR$place <- patient.dirs[k]
  fd.TotalEdgeSynonymus_FR$place <- patient.dirs[k]
  fd.TotalEdgeNonSynonymus_CDR$place <- patient.dirs[k]
  fd.TotalEdgeNonSynonymus_FR$place <- patient.dirs[k]
  fd.TotalEdgeCDR$place <- patient.dirs[k]
  fd.TotalEdgeFR$place <- patient.dirs[k]
  fd.TotalEdgeRatio_CDRvsFR$place <- patient.dirs[k]
  fd.TotalEdgeRatio_CDR$place <-  patient.dirs[k]
  fd.TotalEdgeRatio_FR$place <-  patient.dirs[k]

  names(fd.TotalEdgeSynonymus) <-c("minus-minus" ,"plus-minus" ,"plus-plus","place")
  names(fd.TotalEdgeNonSynonymus) <- c("minus-minus" ,"plus-minus" ,"plus-plus","place")
  names(fd.TotalEdgeDepth) <- c("minus-minus" ,"plus-minus" ,"plus-plus","place")
  names(fd.TotalEdgeRatio) <- c("minus-minus" ,"plus-minus" ,"plus-plus","place")
  names(fd.TotalEdgeSynonymus_CDR)<- c("minus-minus" ,"plus-minus" ,"plus-plus","place")
  names(fd.TotalEdgeSynonymus_FR)<- c("minus-minus" ,"plus-minus" ,"plus-plus","place")
  names(fd.TotalEdgeNonSynonymus_CDR)<-c("minus-minus" ,"plus-minus" ,"plus-plus","place")
  names(fd.TotalEdgeNonSynonymus_FR)<-  c("minus-minus" ,"plus-minus" ,"plus-plus","place")
  names(fd.TotalEdgeCDR) <-c("minus-minus" ,"plus-minus" ,"plus-plus","place")
  names(fd.TotalEdgeFR) <-c("minus-minus" ,"plus-minus" ,"plus-plus","place")
  names(fd.TotalEdgeRatio_CDRvsFR)<-c("minus-minus" ,"plus-minus" ,"plus-plus","place")
  names(fd.TotalEdgeRatio_CDR) <- c("minus-minus" ,"plus-minus" ,"plus-plus","place")
  names(fd.TotalEdgeRatio_FR) <- c("minus-minus" ,"plus-minus" ,"plus-plus","place")
  
  ### add the current patiant dataframe into the data frame that saves all patiants togther
  fd_ALL.TotalEdgeSynonymus <- rbind(fd_ALL.TotalEdgeSynonymus,fd.TotalEdgeSynonymus)
  fd_ALL.TotalEdgeNonSynonymus <- rbind(fd_ALL.TotalEdgeNonSynonymus,fd.TotalEdgeNonSynonymus)
  fd_ALL.TotalEdgeDepth <- rbind(fd_ALL.TotalEdgeDepth,fd.TotalEdgeDepth)
  fd_ALL.TotalEdgeRatio <- rbind(fd_ALL.TotalEdgeRatio,fd.TotalEdgeRatio)
  
  fd_ALL.TotalEdgeSynonymus_CDR<- rbind(fd_ALL.TotalEdgeSynonymus_CDR,fd.TotalEdgeSynonymus_CDR)
  fd_ALL.TotalEdgeSynonymus_FR<-  rbind(fd_ALL.TotalEdgeSynonymus_FR,fd.TotalEdgeSynonymus_FR)
  fd_ALL.TotalEdgeNonSynonymus_CDR<- rbind(fd_ALL.TotalEdgeNonSynonymus_CDR,fd.TotalEdgeNonSynonymus_CDR)
  fd_ALL.TotalEdgeNonSynonymus_FR<-  rbind(fd_ALL.TotalEdgeNonSynonymus_FR,fd.TotalEdgeNonSynonymus_FR)
  fd_ALL.TotalEdgeCDR<-  rbind(fd_ALL.TotalEdgeCDR,fd.TotalEdgeCDR)
  fd_ALL.TotalEdgeFR<-  rbind(fd_ALL.TotalEdgeFR,fd.TotalEdgeFR)
  fd_ALL.TotalEdgeRatio_CDRvsFR<-  rbind(fd_ALL.TotalEdgeRatio_CDRvsFR,fd.TotalEdgeRatio_CDRvsFR)
  fd_ALL.TotalEdgeRatio_CDR <- rbind(fd_ALL.TotalEdgeRatio_CDR,fd.TotalEdgeRatio_CDR)
  fd_ALL.TotalEdgeRatio_FR <- rbind(fd_ALL.TotalEdgeRatio_FR,fd.TotalEdgeRatio_FR)
}  

write.csv(fd_ALL.TotalEdgeSynonymus,file = "SPLEEN_Edges_SYN.csv")
write.csv(fd_ALL.TotalEdgeNonSynonymus,file = "SPLEEN_Edges_NonSYN.csv")
write.csv(fd_ALL.TotalEdgeDepth,file = "SPLEEN_Edges_Depth.csv")
write.csv(fd_ALL.TotalEdgeRatio,file = "SPLEEN_Edges_Ratio.csv")
write.csv(fd_ALL.TotalEdgeSynonymus_CDR,file = "SPLEEN_Edges_SYN_CDR.csv")
write.csv(fd_ALL.TotalEdgeSynonymus_FR,file = "SPLEEN_Edges_SYN_FR.csv")
write.csv(fd_ALL.TotalEdgeNonSynonymus_CDR,file = "SPLEEN_Edges_NonSYN_CDR.csv")
write.csv(fd_ALL.TotalEdgeNonSynonymus_FR,file = "SPLEEN_Edges_NonSYN_FR.csv")
write.csv(fd_ALL.TotalEdgeCDR,file = "SPLEEN_Edges_CDR.csv")
write.csv(fd_ALL.TotalEdgeFR,file = "SPLEEN_Edges_FR.csv")
write.csv(fd_ALL.TotalEdgeRatio_CDRvsFR,file = "SPLEEN_Edges_cdrVSfr_Ratio.csv")
write.csv(fd_ALL.TotalEdgeRatio_CDR,file = "SPLEEN_Edges_Ratio_CDR.csv")
write.csv(fd_ALL.TotalEdgeRatio_FR,file = "SPLEEN_Edges_Ratio_FR.csv")

