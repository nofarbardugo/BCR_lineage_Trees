# set directory : change path depending on work environment

source.dir <-'/home/bardugn1/pfizer' # pep4
setwd(source.dir)


getIsotypeName <-function (node)
  # this function get a node and return its isotype name
  # Args: 
  #   node- an node of the tree   
  # Returns: isotype nane of the node
{
  return (isotypesNames.df[which(isotypesNames.df[,"head"] == node),"IsotypeName"])
}



GoOverTree<-function(node,mutationNum)
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
  
  # will savw the data dor the father antil end to update all children
  df.tmp <-data.frame(IsotypeName = c(0), Son =c(0),theSame =c(0),depth = c(0))
  
  
  # if its not a leaf
  if(sunsNum > 0)
  {
    p<-1  
    while(p <=sunsNum)
    {
      
      # get isotype name 
      type <-getIsotypeName(node)
      
      if(suns[p] !="VJ_1")
      {
        # get mutation number antil child level
        depth<- edges.df[(edges.df [,"to"]==suns[p]),"distance"] + mutationNum
        
        GoOverTree(suns[p],depth) # calc child 
        
        # check if the child type is like the father
        sontype <- getIsotypeName(suns[p])
        
        # the same - add one to same count
        if(sontype==type)
        {
          df.tmp[p,] <- c(type,sontype,TRUE,mutationNum)
          
        }
        # not the same - add one to the diffrent count
        else
        {
          df.tmp[p,] <- c(type,sontype,FALSE,mutationNum)
        }
        
      }
      
      p <- p+1    
    }  #if(sunsNum > 0)
    
    data  <- df.allIsotypeDetails
    data <- rbind(data,df.tmp) 
    
    # save "list" in the outer "allEdgeDetails" data frame
    assign("df.allIsotypeDetails",data, envir = .GlobalEnv)
  }
}  



############## main ###########################

# path for trees folder
tree.dir <- paste0(source.dir,"/pfizerHumanTreeNew/")

# get list of folders
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)

# create datafream that will save for each node 
df.allIsotypeDetails <- data.frame(IsotypeName = c(), Son =c(),theSame =c(),depth = c())

# going over each patient
for(k in 1:length(patient.dirs))
{
  # get list of folders
  dirs <- list.files(paste0(tree.dir,patient.dirs[k]),full.names = F)
  
  # going over each tree folder
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
        j<-1
        
        while(j<=length(roots))
        { 
          # count modes only on trees that have at last 10 leave and 2 types of isotipes 
          if(!(is.null(totalLeavesType.fd)) && (totalLeavesType.fd[j,"type_num"]>1)
             && (totalLeavesType.fd[j,"sum_Leaves"]>10) )
          {
            
            # count modes only on tree plot only trees that have at last 10 leave and only 1 type of isotipes 
            #   if(!(is.null(totalLeavesType.fd)) && (totalLeavesType.fd[j,"type_num"]==1)
            #       && (totalLeavesType.fd[j,"sum_Leaves"]>10) )
            #    {  
            
            print ("in")
            GoOverTree(roots[j],0) 
            
          }
          j<-j + 1
        }
      }
    }
  }
  
}

# write to file 


names(df.allIsotypeDetails)<- c ("Isotype Name","Son type","Is son like father", "Depth")

write.csv(df.allIsotypeDetails,file = "PFIZER_Node_son_divide_per_isotype.csv")

# if the script run on NOCDR3 mut mode
