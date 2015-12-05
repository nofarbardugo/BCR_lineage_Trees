# set directory : change path depending on work environment
#source.dir <- 'G:/project/files/pfizer' # window
source.dir <-'/home/bardugn1/pfizer' # pep4
#source.dir<- '/media/nofar/F6F9-0D76/project/files/pfizer' #ubuntuo
#source.dir <-'C:/Users/Nofar/Desktop/trees' 
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
  return(switch(combine,
   "IGHGIGHG"= "G-G" , "IGHAIGHA"= "A-A","IGHEIGHE" ="E-E","IGHMIGHM"="M-M",
   "naive_IGHDnaive_IGHD" = "ND-ND", "naive_IGHMnaive_IGHM" ="NM-NM",
   "naive_IGHDIGHM" = "ND-M","naive_IGHDIGHG" = "ND-G","naive_IGHDIGHA" = "ND-A","naive_IGHDIGHE" = "ND-E",
   "naive_IGHMIGHM" = "NM-M","naive_IGHMIGHG" = "NM-G", "naive_IGHMIGHA" = "NM-A","naive_IGHMIGHE" = "NM-E",
   "IGHMIGHG" = "M-G","IGHMIGHA" = "M-A", "IGHMIGHE" = "M-E",
   "IGHGIGHA" ="G-A","IGHGIGHE" ="G-E"))
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
      # get edges type
      type <-getEdgeName(node,suns[p])
      
      if(suns[p] !="VJ_1")
      {
        syn<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms"]
        nonS<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms"]
        list[[type]] <- rbind(allEdgeDetails[[type]],c(syn,nonS)) 
       
      }
      p <- p+1

    }
    
    # save "list" in the outer "allEdgeDetails" data frame
    assign("allEdgeDetails",list, envir = .GlobalEnv)
  }
  
}  


############## main ###########################

EDGE_NUM <- 19
# path for trees folder
tree.dir <- paste0(source.dir,"/pfizerHumanTree/")

# get list of fold
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)

# FOR ALL PATIANTs TOGTHER ------ data fream that save details about each tree edge mean Synonymus mutataion - 
fd_ALL.TotalEdgeSynonymus<- data.frame("G-G"=c() ,"A-A" = c(),"E-E" =c(),"M-M" = c(),"ND-ND" = c(), "NM-NM" = c(),
                                   "ND-M"=c() ,"ND-G" = c(),"ND-A" =c(),"ND-E" = c(),
                                   "NM-M"=c() ,"NM-G" = c(),"NM-A" =c(),"NM-E" = c(),
                                   "M-G" = c(),"M-A" = c(), "M-E" = c(),
                                   "G-A"=c() ,"G-E" = c())

# FOR ALL PATIANTs TOGTHER ------data fream that save details about each tree edge mean NON Synonymus mutataion - 
fd_ALL.TotalEdgeNonSynonymus<- data.frame("G-G"=c() ,"A-A" = c(),"E-E" =c(),"M-M" = c(),"ND-ND" = c(), "NM-NM" = c(),
                                      "ND-M"=c() ,"ND-G" = c(),"ND-A" =c(),"ND-E" = c(),
                                      "NM-M"=c() ,"NM-G" = c(),"NM-A" =c(),"NM-E" = c(),
                                      "M-G" = c(),"M-A" = c(), "M-E" = c(),
                                      "G-A"=c() ,"G-E" = c())

# FOR ALL PATIANTs TOGTHER ------ data fream that save details about each tree edge mean depth mutataion -(s+ns)
fd_ALL.TotalEdgeDepth<- data.frame("G-G"=c() ,"A-A" = c(),"E-E" =c(),"M-M" = c(),"ND-ND" = c(), "NM-NM" = c(),
                                       "ND-M"=c() ,"ND-G" = c(),"ND-A" =c(),"ND-E" = c(),
                                       "NM-M"=c() ,"NM-G" = c(),"NM-A" =c(),"NM-E" = c(),
                                       "M-G" = c(),"M-A" = c(), "M-E" = c(),
                                       "G-A"=c() ,"G-E" = c())

# FOR ALL PATIANTs TOGTHER ------data fream that save details about each tree edge mean ratio mutataion - (ns/(s+ns))
fd_ALL.TotalEdgeRatio<- data.frame("G-G"=c() ,"A-A" = c(),"E-E" =c(),"M-M" = c(),"ND-ND" = c(), "NM-NM" = c(),
                                          "ND-M"=c() ,"ND-G" = c(),"ND-A" =c(),"ND-E" = c(),
                                          "NM-M"=c() ,"NM-G" = c(),"NM-A" =c(),"NM-E" = c(),
                                          "M-G" = c(),"M-A" = c(), "M-E" = c(),
                                          "G-A"=c() ,"G-E" = c())

NumberOfTreesCount <- 0

# going over each patient
for(k in 1:length(patient.dirs))
{
  # decleare her to get detailes for each patient, can move out in order to get details for all patient togther
  
  # data fream that save details about each tree edge mean Synonymus mutataion - 
  fd.TotalEdgeSynonymus<- data.frame("G-G"=c() ,"A-A" = c(),"E-E" =c(),"M-M" = c(),"ND-ND" = c(), "NM-NM" = c(),
                                     "ND-M"=c() ,"ND-G" = c(),"ND-A" =c(),"ND-E" = c(),
                                     "NM-M"=c() ,"NM-G" = c(),"NM-A" =c(),"NM-E" = c(),
                                     "M-G" = c(),"M-A" = c(), "M-E" = c(),
                                     "G-A"=c() ,"G-E" = c())
  
  # data fream that save details about each tree edge mean NON Synonymus mutataion - 
  fd.TotalEdgeNonSynonymus<- data.frame("G-G"=c() ,"A-A" = c(),"E-E" =c(),"M-M" = c(),"ND-ND" = c(), "NM-NM" = c(),
                                        "ND-M"=c() ,"ND-G" = c(),"ND-A" =c(),"ND-E" = c(),
                                        "NM-M"=c() ,"NM-G" = c(),"NM-A" =c(),"NM-E" = c(),
                                        "M-G" = c(),"M-A" = c(), "M-E" = c(),
                                        "G-A"=c() ,"G-E" = c())
  
  # data fream that save details about each tree edge mean depth mutataion -(s+ns)
  fd.TotalEdgeDepth<- data.frame("G-G"=c() ,"A-A" = c(),"E-E" =c(),"M-M" = c(),"ND-ND" = c(), "NM-NM" = c(),
                                 "ND-M"=c() ,"ND-G" = c(),"ND-A" =c(),"ND-E" = c(),
                                 "NM-M"=c() ,"NM-G" = c(),"NM-A" =c(),"NM-E" = c(),
                                 "M-G" = c(),"M-A" = c(), "M-E" = c(),
                                 "G-A"=c() ,"G-E" = c())
  
  # data fream that save details about each tree edge mean ratio mutataion - (ns/(s+ns))
  fd.TotalEdgeRatio<- data.frame("G-G"=c() ,"A-A" = c(),"E-E" =c(),"M-M" = c(),"ND-ND" = c(), "NM-NM" = c(),
                                 "ND-M"=c() ,"ND-G" = c(),"ND-A" =c(),"ND-E" = c(),
                                 "NM-M"=c() ,"NM-G" = c(),"NM-A" =c(),"NM-E" = c(),
                                 "M-G" = c(),"M-A" = c(), "M-E" = c(),
                                 "G-A"=c() ,"G-E" = c())
  
  
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
          
          # count modes only on trees that have at last 10 leave and 2 types of isotipes 
        #  if(!(is.null(totalLeavesType.fd)) && (totalLeavesType.fd[j,"type_num"]>1)
        #     && (totalLeavesType.fd[j,"sum_Leaves"]>10) )
        #  {
             
            # count modes only on tree plot only trees that have at last 10 leave and only 1 type of isotipes 
            if(!(is.null(totalLeavesType.fd)) && (totalLeavesType.fd[j,"type_num"]==1)
               && (totalLeavesType.fd[j,"sum_Leaves"]>10) )
            {  
              
            # for count how many treee ther is
            NumberOfTreesCount <- NumberOfTreesCount +1
            print("in")
            # create datafream for each edge that will save: Synonymus, and nonSynonymus
            # save all data freams into a list for save all togther
            allEdgeDetails <- list()
            allEdgeDetails[["G-G"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["A-A"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["E-E"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["M-M"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["ND-ND"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["NM-NM"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["ND-M"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["ND-G"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["ND-A"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["ND-E"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["NM-M"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["NM-G"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["NM-A"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["NM-E"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["M-G"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["M-A"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["M-E"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["G-A"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            allEdgeDetails[["G-E"]] <- data.frame(Synonymus =c(),nonSynonymus =c())
            
            GoOverTree(roots[j]) 
            
            AvgSYN.vector <- vector()
            AvgNON.vector <- vector()
            AvgDepth.vector <- vector()
            AvgRatio.vector <- vector()
            
            # calc avg of each dataframe and add to the big dataFrame
            for(p in 1:EDGE_NUM)
            {
              if(length(allEdgeDetails[[p]]) > 0)
              {
                AvgSYN.vector[p] <-  mean((allEdgeDetails[[p]])[,1]) # avg of the synonymus in this edge in the tree
                AvgNON.vector[p] <-  mean((allEdgeDetails[[p]])[,2]) # avg of the nonsynonymus in this edge in the tree
                Depth <- (allEdgeDetails[[p]])[,1]+ (allEdgeDetails[[p]])[,2] # avg of the depth in this edge in the tree
                AvgDepth.vector[p] <-  mean(Depth) # avg of Ratio

                div <- c()
                for(d in 1: length(Depth))
                {
                  if(Depth[d]!=0)
                  {
                    div[d] <- (allEdgeDetails[[p]])[d,2]/Depth[d]
                  }
                  else
                  {
                    div[d] <- 0
                  }
                }
                  
                AvgRatio.vector[p] <-  mean(div) # avg of Depth
              }
              else
              {
                AvgSYN.vector[p]=NA
                AvgNON.vector[p]= NA
                AvgDepth.vector[p]=NA
                AvgRatio.vector[p]= NA
              }
           
            }
            
            # add new line to the big data freams
            fd.TotalEdgeSynonymus <- rbind(fd.TotalEdgeSynonymus,AvgSYN.vector)
            fd.TotalEdgeNonSynonymus <- rbind(fd.TotalEdgeNonSynonymus,AvgNON.vector)
            fd.TotalEdgeDepth <- rbind(fd.TotalEdgeDepth,  AvgDepth.vector)
            fd.TotalEdgeRatio <- rbind(fd.TotalEdgeRatio,AvgRatio.vector)
          }
          j<-j + 1
        }
      }
    }
  }

  names(fd.TotalEdgeSynonymus) <- c("G-G" ,"A-A" ,"E-E","M-M" ,"ND-ND", "NM-NM",
                                    "ND-M","ND-G","ND-A","ND-E",
                                    "NM-M","NM-G","NM-A","NM-E",
                                    "M-G","M-A", "M-E",
                                    "G-A" ,"G-E")
  
  names(fd.TotalEdgeNonSynonymus) <- c("G-G" ,"A-A" ,"E-E","M-M" ,"ND-ND", "NM-NM",
                                       "ND-M","ND-G","ND-A","ND-E",
                                       "NM-M","NM-G","NM-A","NM-E",
                                       "M-G","M-A", "M-E",
                                       "G-A" ,"G-E")
  
  names(fd.TotalEdgeDepth) <- c("G-G" ,"A-A" ,"E-E","M-M" ,"ND-ND", "NM-NM",
                                    "ND-M","ND-G","ND-A","ND-E",
                                    "NM-M","NM-G","NM-A","NM-E",
                                    "M-G","M-A", "M-E",
                                    "G-A" ,"G-E")
  
  names(fd.TotalEdgeRatio) <- c("G-G" ,"A-A" ,"E-E","M-M" ,"ND-ND", "NM-NM",
                                       "ND-M","ND-G","ND-A","ND-E",
                                       "NM-M","NM-G","NM-A","NM-E",
                                       "M-G","M-A", "M-E",
                                       "G-A" ,"G-E")
  
  ### add the current patiant dataframe into the data frame that saves all patiants togther
  fd_ALL.TotalEdgeSynonymus <- rbind(fd_ALL.TotalEdgeSynonymus,fd.TotalEdgeSynonymus)
  fd_ALL.TotalEdgeNonSynonymus <- rbind(fd_ALL.TotalEdgeNonSynonymus,fd.TotalEdgeNonSynonymus)
  fd_ALL.TotalEdgeDepth <- rbind(fd_ALL.TotalEdgeDepth,fd.TotalEdgeDepth)
  fd_ALL.TotalEdgeRatio <- rbind(fd_ALL.TotalEdgeRatio,fd.TotalEdgeRatio)

}  

#write into file

write.csv(fd_ALL.TotalEdgeSynonymus,file = "PFIZER_Edges_SYN.csv")
write.csv(fd_ALL.TotalEdgeNonSynonymus,file = "PFIZER_Edges_NonSYN.csv")
write.csv(fd_ALL.TotalEdgeDepth,file = "PFIZER_Edges_Depth.csv")
write.csv(fd_ALL.TotalEdgeRatio,file = "PFIZER_Edges_Ratio.csv")\

