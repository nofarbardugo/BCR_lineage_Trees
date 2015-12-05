# set directory : change path depending on work environment
#source.dir <- 'G:/project/files/pfizer'
#source.dir<- '/media/nofar/F6F9-0D76/project/files/pfizer'
#source.dir <-'C:/Users/Nofar/Desktop/trees'
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
  
  # if we reached to the leafs - return 1 to indicate one child
  if(sunsNum==0)
  {
    return (1);
  }
  else
  {
    # save all childern num (include childern of children..)
    childNum <- 0  
    p<-1
    
    while(p <=sunsNum)
    {
      childNum <- childNum + GoOverTree(suns[p])
      p <- p+1
    }
  
    firstVSall <- sunsNum/childNum 
    # save how mach childs the cur nodes have
    type <-getIsotypeName(node)
    list  <- allIsotypeDetails
    list[[type]] <- rbind(allIsotypeDetails[[type]],c(sunsNum,childNum,firstVSall)) 
    
    # save "df_TypeData" in the outer "fd.leavesType" data frame
    assign("allIsotypeDetails",list, envir = .GlobalEnv)
    
    
    return (childNum)
  }
 
}  


############## main ###########################

# path for trees folder
tree.dir <- paste0(source.dir,"/pfizerHumanTree/")

# get list of fold
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)

# data fream that save details about each tree leaves - just for get details
fd.TotalLeaves<- data.frame(treeNUM=c() ,IGHA = c(),IGHE =c(),IGHG = c()
                            ,IGHM = c(), naive_IGHD = c(),naive_IGHM = c(),sum_Leaves =c(),type_num = c())

# create datafream for each isotype that will save: number of suns (first childs only), all child, suns/children
# save all data freams into a list for save all togther
TOTAL.allIsotypeDetails <- list()
TOTAL.allIsotypeDetails[["IGHA"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c())
TOTAL.allIsotypeDetails[["IGHE"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c())
TOTAL.allIsotypeDetails[["IGHG"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c())
TOTAL.allIsotypeDetails[["IGHM"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c())
TOTAL.allIsotypeDetails[["naive_IGHD"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c())
TOTAL.allIsotypeDetails[["naive_IGHM"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c())

# going over each patient
for(k in 1:length(patient.dirs))
{
  # get list of folders
  dirs <- list.files(paste0(tree.dir,patient.dirs[k]),full.names = F)
  
  # create datafream for each isotype that will save: number of suns (first childs only), all child, suns/children
  # save all data freams into a list for save all togther
  allIsotypeDetails <- list()
  allIsotypeDetails[["IGHA"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c())
  allIsotypeDetails[["IGHE"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c())
  allIsotypeDetails[["IGHG"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c())
  allIsotypeDetails[["IGHM"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c())
  allIsotypeDetails[["naive_IGHD"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c())
  allIsotypeDetails[["naive_IGHM"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c())
  
  
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
        
        # get tabale with all node's isotypes names
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
       #   if(!(is.null(totalLeavesType.fd)) && (totalLeavesType.fd[j,"type_num"]>1)
        #     && (totalLeavesType.fd[j,"sum_Leaves"]>10) )
        #  {
            
          # count modes only on tree plot only trees that have at last 10 leave and only 1 type of isotipes 
          if(!(is.null(totalLeavesType.fd)) && (totalLeavesType.fd[j,"type_num"]==1)
             && (totalLeavesType.fd[j,"sum_Leaves"]>10) )
          {     
            
            print ("in")
            GoOverTree(roots[j]) 
            
          }
          j<-j + 1
        }
      }
    }
  }
   
  
  names(allIsotypeDetails[["IGHA"]]) =names(TOTAL.allIsotypeDetails[["IGHA"]])
  names(allIsotypeDetails[["IGHE"]]) =names(TOTAL.allIsotypeDetails[["IGHE"]])
  names(allIsotypeDetails[["IGHG"]]) =names(TOTAL.allIsotypeDetails[["IGHG"]])
  names(allIsotypeDetails[["IGHM"]]) =names(TOTAL.allIsotypeDetails[["IGHM"]])
  names(allIsotypeDetails[["naive_IGHM"]]) =names(TOTAL.allIsotypeDetails[["naive_IGHM"]])
  names(allIsotypeDetails[["naive_IGHD"]]) =names(TOTAL.allIsotypeDetails[["naive_IGHD"]])
  ####   add the current patiant list into the list that saves all patiants togther
  TOTAL.allIsotypeDetails[["IGHA"]] <- rbind(TOTAL.allIsotypeDetails[["IGHA"]],allIsotypeDetails[["IGHA"]])
  TOTAL.allIsotypeDetails[["IGHE"]] <- rbind(TOTAL.allIsotypeDetails[["IGHE"]],allIsotypeDetails[["IGHE"]])
  TOTAL.allIsotypeDetails[["IGHG"]] <- rbind(TOTAL.allIsotypeDetails[["IGHG"]],allIsotypeDetails[["IGHG"]])
  TOTAL.allIsotypeDetails[["IGHM"]] <- rbind(TOTAL.allIsotypeDetails[["IGHM"]],allIsotypeDetails[["IGHM"]])
  TOTAL.allIsotypeDetails[["naive_IGHD"]] <- rbind(TOTAL.allIsotypeDetails[["naive_IGHD"]],allIsotypeDetails[["naive_IGHD"]])
  TOTAL.allIsotypeDetails[["naive_IGHM"]] <- rbind(TOTAL.allIsotypeDetails[["naive_IGHM"]],allIsotypeDetails[["naive_IGHM"]])
  
} 


# write to file 

names(TOTAL.allIsotypeDetails[["IGHA"]])<- c ("first_generation_child_num","All_generation_child_num","firstVSallChilden")
names(TOTAL.allIsotypeDetails[["IGHG"]])<- c ("first_generation_child_num","All_generation_child_num","firstVSallChilden")
names(TOTAL.allIsotypeDetails[["IGHM"]])<- c ("first_generation_child_num","All_generation_child_num","firstVSallChilden")
names(TOTAL.allIsotypeDetails[["naive_IGHM"]])<- c ("first_generation_child_num","All_generation_child_num","firstVSallChilden")

write.csv(TOTAL.allIsotypeDetails[["IGHA"]],file = "PFIZER_node_IGHA.csv")
write.csv(TOTAL.allIsotypeDetails[["IGHG"]],file = "PFIZER_node_IGHG.csv")
write.csv(TOTAL.allIsotypeDetails[["IGHM"]],file = "PFIZER_node_IGHM.csv")
write.csv(TOTAL.allIsotypeDetails[["naive_IGHM"]],file = "PFIZER_node_naive_IGHM.csv")


