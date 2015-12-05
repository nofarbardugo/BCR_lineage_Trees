
#const
WITH_NAIVE<-5

# set directory : change path depending on work environment
#source.dir <- 'G:/project/files/pfizer'
#source.dir <- '/u/peptibase2/bardugn1'
#source.dir <-'C:/Users/Nofar/Desktop/trees'
source.dir <- "/home/bardugn1/pfizer"
setwd(source.dir)



ChangeIsotypeToNumber <- function(isotypeCharacter)
{
  # the function get Isotype and retern is numeric representation
  #
  # Args:
  #   isotypeCharacter: character  
  # Returns:
  #   reutrn the isotypes numeric representation  
  return (switch(isotypeCharacter, "naive_IGHM_B"=12,"naive_IGHM_A"=11,"naive_IGHD_B"=10, "naive_IGHD_A"=9, "IGHM_B"=8, 
                                    "IGHM_A"=7,"IGHG_B"=6,"IGHG_A"=5,  "IGHE_B"=4,"IGHE_A"=3 ,"IGHA_B"=2 , "IGHA_A"=1, NA))
}

GoOverTree<-function(node)
  # this function goes over the tree and count isotype leaves in the tree
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
  if(length(suns)==0  )
  {
    if(node != "VJ_1")
    {
      
      # gets the full name of the leaf
      fullName <- df.names[which(df.names[,"head"] == node),"head2"]
      
      # the number of "_" in the full name of the leaf
      makafNumber <- gregexpr(pattern ='_',fullName)
      
      # if the full name of the leaf contains the word "naive" in it we will includ it in the
      # "IsotypeName" column, otherwise we will write there only the i isotype name
      
      if(length(makafNumber[[1]])==WITH_NAIVE)
      {
        sample <- substr(fullName,makafNumber[[1]][4] -1,makafNumber[[1]][4]-1)
      }
      else
      {
        sample <- substr(fullName,makafNumber[[1]][3] -1,makafNumber[[1]][3]-1)
      }
      
      isotypeName <-  paste(df.names[which(df.names[,"head"] == node),"IsotypeName"],sample,sep = "_")
      
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
      
      
    }
  }
  # if we reached to an edge with suns
  else
  {
    # (a for loop) apply the "GoOverTree" function on all the suns to get their isotypes
    sapply(suns,GoOverTree)
  }
}


########### MAIN #########################
tree.dir <- paste0(source.dir,"/pfizerHumanTree/")

# get list of fold
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)

# data fream that save details about each tree leaves - just for get details
fd.TotalLeavesForALL<- data.frame(treeNUM=c() ,IGHA_A = c(),IGHA_B = c(),IGHE_A =c(),IGHE_B = c(),
                                      IGHG_A = c(),IGHG_G = c() ,IGHM_A = c(),IGHM_B = c(),
                                      naive_IGHD_A = c(),naive_IGHD_B = c(),
                                      naive_IGHM_A = c(),naive_IGHM_B = c(),sum_Leaves =c(),type_num = c())
# going over each patient
for(k in 1:length(patient.dirs))
{
  # get list of folders
  dirs <- list.files(paste0(tree.dir,patient.dirs[k]),full.names = F)
  
  # data fream that save details about each tree leaves - just for get details
  fd.TotalLeavesForPatient<- data.frame(treeNUM=c() ,IGHA_A = c(),IGHA_B = c(),IGHE_A =c(),IGHE_B = c(),
                              IGHG_A = c(),IGHG_B = c() ,IGHM_A = c(),IGHM_B = c(),
                              naive_IGHD_A = c(),naive_IGHD_B = c(),
                              naive_IGHM_A = c(),naive_IGHM_B = c(),sum_Leaves =c(),type_num = c())
  
  
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
        df.names <- read.table(file= paste0(dirPath,"_isotypeNames.tab"),header = T,stringsAsFactors = F,sep ="\t")
        
        # get a root for each clone
        roots <- setdiff(df.edges[,"from"],df.edges[,"to"])
        
        #create dataFrame that will save for each clone how much leave from each type we have
        fd.TotalLeavesTypeForClone<- data.frame(treeNUM=c() ,IGHA_A = c(),IGHA_B = c(),IGHE_A =c(),IGHE_B = c(),
                                                IGHG_A = c(),IGHG_B = c() ,IGHM_A = c(),IGHM_B = c(),
                                                naive_IGHD_A = c(),naive_IGHD_B = c(),
                                                naive_IGHM_A = c(),naive_IGHM_B = c(),sum_Leaves =c(),type_num = c())
        p<-1
        while(p <=length(roots))
        {
          
          fd.leavesType <- data.frame(treeNUM=c(paste(dirs[i],roots[p],sep = '_')) ,
                                      IGHA_A = c(0),IGHA_B = c(0),IGHE_A =c(0),IGHE_B = c(0),
                                      IGHG_A = c(0),IGHG_B = c(0) ,IGHM_A = c(0),IGHM_B = c(0),
                                      naive_IGHD_A = c(0),naive_IGHD_B = c(0),
                                      naive_IGHM_A = c(0),naive_IGHM_B = c(0)
                                      ,sum_Leaves =c(0),type_num = c(0))
          
          GoOverTree(roots[p])
          
          # calc total leaves num 
          fd.leavesType[1,"sum_Leaves"] <- fd.leavesType[1,"IGHA_A"] + fd.leavesType[1,"IGHA_B"] +
                                           fd.leavesType[1,"IGHE_A"] + fd.leavesType[1,"IGHE_B"] +
                                           fd.leavesType[1,"IGHG_A"] + fd.leavesType[1,"IGHG_B"] +
                                           fd.leavesType[1,"IGHM_A"] + fd.leavesType[1,"IGHM_A"] +
                                           fd.leavesType[1,"naive_IGHD_A"] +
                                           fd.leavesType[1,"naive_IGHD_B"] +
                                           fd.leavesType[1,"naive_IGHM_A"] +
                                           fd.leavesType[1,"naive_IGHM_B"]
          
          # get type isotype num 
          num <- 0
          for (t in 2:13)
          {
            if (fd.leavesType[1,t] >0)
            {
              num <- num + 1
            }
          }
          
          fd.leavesType[1,"type_num"] <-num
          
          # add fd.leavesType to total
          fd.TotalLeavesTypeForClone <- rbind( fd.TotalLeavesTypeForClone, fd.leavesType)
          
          p<-p+1
        }
        
        
        #fd.leavesType = data.frame(treeNUM,naive_IGHM, naive_IGHD, IGHM, IGHG, IGHE, IGHA)
        # goes over each root and send the root to "GoOverTree" function 
        #sapply(roots,GoOverTree,leavesTypeData =fd.leavesType[,])
        
        # save the TotalLeavesType table into a file
        write.table(fd.TotalLeavesTypeForClone, file = paste0(dirPath,"_totalLeavesType_A_B.csv"),
                    quote=F, sep='\t', col.names=T, row.names=T,append = F)
      
        # add fd.leavesType to total
        fd.TotalLeavesForPatient <- rbind(fd.TotalLeavesForPatient, fd.TotalLeavesTypeForClone)
      }
    }
  }
  
  # save all petient togther
  write.table(fd.TotalLeavesForPatient, file = paste0(source.dir,"/",patient.dirs[k],"-ALL_totalLeavesType_A_B.csv"),
              quote=F, sep='\t', col.names=T, row.names=T,append = F)
  
  fd.TotalLeavesForALL
  
  # add fd.leavesType for petient  to total All
  fd.TotalLeavesForALL <- rbind(fd.TotalLeavesForALL,fd.TotalLeavesForPatient)
  
  # save all tree with more than 1 isotype
  #fd.TotalOver2 <- fd.TotalLeaves[which(fd.TotalLeaves[,"type_num"]>1),]
  # write.csv(fd.TotalOver2, file = paste0(source.dir,"/TotalOver1_",dirs[k],".csv"))
  #fd.TotalOver3 <- fd.TotalLeaves[which(fd.TotalLeaves[,"type_num"]>2),]
  #fd.TotalOver4 <- fd.TotalLeaves[which(fd.TotalLeaves[,"type_num"]>3),]
}

# save all petient togther
write.table(fd.TotalLeavesForALL, file = paste0(source.dir,"/ALL_totalLeavesType_A_B.csv"),
            quote=F, sep='\t', col.names=T, row.names=T,append = F)

