 
#const
#WITH_NAIVE<-7

# set directory : change path depending on work environment
#source.dir <- 'G:/project/files/pfizer' #window
#source.dir <- '/u/peptibase2/bardugn1' # peptibase
#source.dir <-'C:/Users/Nofar/Desktop/trees' #window
source.dir <- "/home/bardugn1/spleen" #pep 4
#source.dir <- "/home/nofar/Desktop/spleen"
setwd(source.dir)



ChangeIsotypeToNumber <- function(isotypeCharacter)
{
  # the function get Isotype and retern is numeric representation
  #
  # Args:
  #   isotypeCharacter: character  
  # Returns:
  #   reutrn the isotypes numeric representation  
  return (switch(isotypeCharacter,"igMplus"=2, "igMminus"=1, NA))
}

ChangeNumberToIsotype <- function(isotypeNumber)
  # the function get Isotype number and retern is charactric representation
  #
  # Args:
  #   isotypeNumber: number
  # Returns:
  #   reutrn the isotypes string  
{
  return (switch(as.character(isotypeNumber), "2"="igMplus", "1"="igMminus"))
}
  
IsotypeRules <- function(isotypeVector)
  # the function get a set of isotypes and return thier father
  #
  # Args:
  #   isotypeVector: vector with isotypes (suns of the node)
  # Returns:
  #   reutrn the isotypes father 
{
  # represent the isotype with numeric values
  numberIsotype.vector <- sapply(isotypeVector,ChangeIsotypeToNumber)
  
  # get the maximum number (the highest isotype in the hierarchy)
  maxNumber <- max(numberIsotype.vector,na.rm=TRUE)
  
  # return the father isotype
  return (ChangeNumberToIsotype(maxNumber))
}

GoOverTree<-function(node)
  # this function goes over the tree and add the isotype name to the names file
  # Args: 
  #   node- an edge of the tree   
  #   leavesTypeData - a data frame that save for isotype how much leaves it has
  # Returns:the node isotype name
{

  
  # if condition-   what we work on[the condition, if the condition is true give me
  #                                                                  the Relevant information]
  # get all the suns 
  suns <- df.edges[df.edges[,"from"]==node,"to"]
  
  # if we reached to the leafs
  if(length(suns)==0)
  {
    # gets the full name of the leaf
    fullName <- df.names[which(df.names[,"head"] == node),"head2"]
    
    # search for igMplus 
    pluseNumber <- gregexpr(pattern ='IgMplus',fullName)
    
    #if it not pluse is negative igM 
    
    # if the full name of the leaf contains the word "naive" in it we will includ it in the
    # "IsotypeName" column, otherwise we will write there only the i isotype name
    
    if(pluseNumber[[1]][1]==-1)
    {
      isotypeName <- "igMminus"
    
    }
    else
    {
      isotypeName <- "igMplus"      
    }
    
    # create local data frame
    df <-df.names
    
    df[which(df[,"head"] == node),"IsotypeName"] <- isotypeName
    
    # save isotype in the outer name data frame
    assign("df.names",df, envir = .GlobalEnv)
    
    if(!is.na(isotypeName))
    {
      # add isotypeName of the leave to "fd.leavesType" count
      # create local data frame
      df <- fd.leavesType
      colToadd <- ChangeIsotypeToNumber(isotypeName)+1
      df[1,colToadd] <- df[1,colToadd]+1
      
      # save "df_TypeData" in the outer "fd.leavesType" data frame
      assign("fd.leavesType",df, envir = .GlobalEnv)
    }
      
    return (isotypeName)
  }
  # if we reached to an edge with suns
  else
  {
    # (a for loop) apply the "GoOverTree" function on all the suns to get their isotype
    # (before we will get the fathers isotype) 
    sunsIsotype.vector <- sapply(suns,GoOverTree)
 
     print(sunsIsotype.vector)
     isotypeNode <- IsotypeRules(sunsIsotype.vector)
    
    # create local data frame
    df <-df.names
    df[which(df[,"head"] == node),"IsotypeName"] <- isotypeNode
    
    # save isotype in the outer name data frame
    assign("df.names",df, envir = .GlobalEnv)

    return (isotypeNode)
  }
}


########### MAIN #########################
tree.dir <- paste0(source.dir,"/spleenHumanTree/")

# get list of fold
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)

# going over each patient
for(k in 1:length(patient.dirs))
{
  # get list of folders
  dirs <- list.files(paste0(tree.dir,patient.dirs[k]),full.names = F)
  
  # data fream that save details about each tree leaves - just for get details
  fd.TotalLeaves<- data.frame(treeNUM=c() ,igmminus =c(),igMplus = c(),sum_Leaves =c(),type_num = c())
  
    #going over each tree folder
  for(i in 1:length(dirs))
  {
    dirPath <- paste0(tree.dir,patient.dirs[k],'/',dirs[i],'/',dirs[i])
    filePath <-paste0(dirPath,"_edgesCut.tab")
    
    # case file exist (if there is not enough or too much sequences, wile be only fasta file)
    if(file.exists(filePath))
    {
  
      # get all the edges of the current tree
      df.edges <- tryCatch(read.table(file = filePath ,header = T,stringsAsFactors = F,sep ="\t"), error=function(e) NULL)
              
      if(!is.null(df.edges))
      {
        # get the names of the edges of the tree and add a column do the relavent isotype name
        df.names <- read.table(file= paste0(dirPath,"_names.tab"),header = T,stringsAsFactors = F,sep ="\t")
        
        # add column that will save the isotype type for each node
        df.names$IsotypeName <- '-'
        
        # get a root for each clone
        roots <- setdiff(df.edges[,"from"],df.edges[,"to"])
        
        #create dataFrame that will save for each clone how much leave from each type we have
        fd.TotalLeavesType<- data.frame(treeNUM=c() ,igMminus=c(),igMplus= c(),sum_Leaves =c(),type_num = c())
        p<-1
        while(p <=length(roots))
        {
          
          fd.leavesType <- data.frame(treeNUM=c(paste(dirs[i],roots[p],sep = '_')) 
                                      ,igMminus=c(0),igMplus = c(0)
                                      ,sum_Leaves =c(0),type_num = c(0))
          
          
          GoOverTree(roots[p])
          
          # calc total leaves num 
          fd.leavesType[1,"sum_Leaves"] <- fd.leavesType[1,"igMplus"]+fd.leavesType[1,"igMminus"] 
          
          # get type isotype num 
          num <- 0
          for (t in 2:3)
          {
            if (fd.leavesType[1,t] >0)
            {
              num <- num + 1
            }
          }
          
          fd.leavesType[1,"type_num"] <-num
          
          # add fd.leavesType to total
          fd.TotalLeavesType <- rbind(fd.TotalLeavesType, fd.leavesType)
          
          p<-p+1
        }
        
        
        #fd.leavesType = data.frame(treeNUM,naive_IGHM, naive_IGHD, IGHM, IGHG, IGHE, IGHA)
        # goes over each root and send the root to "GoOverTree" function 
        #sapply(roots,GoOverTree,leavesTypeData =fd.leavesType[,])
        
        # save the TotalLeavesType table into a file
        write.table(fd.TotalLeavesType, file = paste0(dirPath,"_totalLeavesType.tab"),
                    quote=F, sep='\t', col.names=T, row.names=T,append = F)
       
        # save the new isotypeNames table into a file
        write.table(df.names, file = paste0(dirPath,"_isotypeNames.tab"),
                    quote=F, sep='\t', col.names=T, row.names=F,append = F)
        
        # to remove after run - just for get details
        # add fd.leavesType to total
        fd.TotalLeaves <- rbind(fd.TotalLeaves,fd.TotalLeavesType)
      }
    }
  }
  
  # save all tree with more than 1 isotype
  fd.TotalOver2 <- fd.TotalLeaves[which(fd.TotalLeaves[,"type_num"]>1),]
  write.csv(fd.TotalOver2, file = paste0(source.dir,"/TotalOver1_",dirs[k],".csv"))
  fd.TotalOver3 <- fd.TotalLeaves[which(fd.TotalLeaves[,"type_num"]>2),]
  fd.TotalOver4 <- fd.TotalLeaves[which(fd.TotalLeaves[,"type_num"]>3),]
}

