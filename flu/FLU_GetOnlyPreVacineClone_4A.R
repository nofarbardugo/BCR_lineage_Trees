
TREE.THRESHOLD <- 4

# set directory :change path depending on work environment
source.dir <- '/home/bardugn1/flu2' #pep 4

setwd(source.dir)

tree.dir <- paste0(source.dir,"/fluHumanTree15oct2015/")



GoOverTree<-function(father,child)
  # this function goes over the tree and delete singltons
  # Args:
  #   father- the father name of the child node
  #   child- an edge of the tree
  #   
{
  
  # if condition-   what we work on[the condition, if the condition is true give me
  #                                                                  the Relevant information]
  # get all the suns of the child 
  suns <- df.edgesCut [df.edgesCut [,"from"]==child,"to"]
  
  # if we reached to an edge with suns
  if(length(suns)>0)
  {
    
    # get to singlton - delete  
    if(length(suns)==1)
    {  
      
      df <-df.edgesCut 

      # if it's a root - only delete row, child become the root like above
      # else: father of the child become father of the sun 
      if(!(is.na(father)))
      {
        
        # delete child to sun from data
        df <-df.edgesCut[!(df.edgesCut[,"from"]==child & df.edgesCut[,"to"]==suns[1]),]
        
        df[which(df[,"to"]==child),"to"] <- suns[1]  
        
        child <- father
      }
      
      # save changes in the original data cut
      assign("df.edgesCut",df, envir = .GlobalEnv)
      
      # (a for loop) apply the "GoOverTree" function on all the suns
      # beacuse we delete child , the father of the sun is father 
      # if the child was a root, the sun is the new root and the father is NA anyway
      sapply(suns,GoOverTree,father = father)
      
    }
    else
    {
      # (a for loop) apply the "GoOverTree" function on all the suns
      sapply(suns,GoOverTree,father = child)
    }
    
  }
}


######################main#####################
# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)
patient.dirs <- patient.dirs[-grep('IB',patient.dirs)]

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
      df.edges <- read.table(file = filePath ,header = T,stringsAsFactors = F,sep ="\t")
      
      if (nrow(df.edges)!=0)
      {
        # get Isotypes names/time point
        isotypesNames.df <- read.table(file= paste0(dirPath,"_isotypeNames.tab")
                                       ,header = T,stringsAsFactors = F,sep ="\t")
        
        # get all leave that thier time point is pre vacine
        isotypesNames <- c(isotypesNames.df[grep('RL013', isotypesNames.df$TimePoint),"head"], 
                           isotypesNames.df[grep('RL014', isotypesNames.df$TimePoint),"head"] ,
                           isotypesNames.df[grep('RL015', isotypesNames.df$TimePoint),"head"] )
        
        #save only the edges that theier leaves are in pre vacine time
        df.edgesCut <- df.edges[!is.na(as.numeric(df.edges[,"to"])) | df.edges[,"to"]%in% isotypesNames,]
        
        # delete internal node that became leaves as a result of cut:
        repeat
        {
          # get all leaves (nodes that found only in the "to" column)
          leaves <- setdiff(df.edgesCut[,"to"],df.edgesCut[,"from"])
          
          # and save only leaves that needs to be delete
          #(Distinction: internal node has a numeric name) 
          leaves <- leaves[!is.na(as.numeric(leaves))]
          
          # if internal leaves not found, exit loop
          if(length(leaves)==0)
          {
            break
          }
          
          # save only edge to real leaves  
          v <-vector()
          for(j in 1:length(df.edgesCut[,"to"]))
          {
            v[j]<-!(df.edgesCut[j,"to"] %in% leaves)
          }
          df.edgesCut<- df.edgesCut[v,]
          
        }
        
        # delete singlton places
        # get a root for each clone
        roots <- setdiff(df.edgesCut[,"from"],df.edgesCut[,"to"])
        
        # goes over each root and send the root to "GoOverTree" function 
        sapply(roots,GoOverTree,father = NA)
        
        # save the new table into a file - if file exit - dystroyed old file
        write.table(df.edgesCut, file = paste0(dirPath, "_edgesCut_preVacine.tab"), quote=F, sep='\t', col.names=T, row.names=F, append = F)
        print(paste0(dirPath, "_edgesCut_preVacine.tab"))
      }
      
    }
    
  }
  
  
}

print('finish')