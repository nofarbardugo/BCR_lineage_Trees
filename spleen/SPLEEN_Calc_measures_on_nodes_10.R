# set directory : change path depending on work environment
source.dir <-'/home/bardugn1/spleen' # pep4
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
tree.dir <- paste0(source.dir,"/spleenHumanTree/")

# get list of fold
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)

# data fream that save details about each tree leaves - just for get details
#fd.TotalLeaves<- data.frame(treeNUM=c() igMminus = c(),igMplus = c(),sum_Leaves =c(),type_num = c())

# create datafream for each isotype that will save: number of suns (first childs only), all child, suns/children
# save all data freams into a list for save all togther
TOTAL.allIsotypeDetails <- list()
TOTAL.allIsotypeDetails[["igMminus"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c(),place = c())
TOTAL.allIsotypeDetails[["igMplus"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c(),place = c())

# going over each patient
for(k in 1:length(patient.dirs))
{
  # get list of folders
  dirs <- list.files(paste0(tree.dir,patient.dirs[k]),full.names = F)
  
  # create datafream for each isotype that will save: number of suns (first childs only), all child, suns/children
  # save all data freams into a list for save all togther
  allIsotypeDetails <- list()
  allIsotypeDetails[["igMminus"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c())
  allIsotypeDetails[["igMplus"]] <- data.frame(First_generation_child_num =c(),All_generation_child_num =c(),firstVSallChilden = c())
  
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
          # count modes only on trees that have at last 10 leave 
          if(!(is.null(totalLeavesType.fd)) && (totalLeavesType.fd[j,"sum_Leaves"]>10))
          {      
            
            print ("in")
            GoOverTree(roots[j]) 
            
          }
          j<-j + 1
        }
      }
    }
  }
  allIsotypeDetails[["igMminus"]] <- cbind(allIsotypeDetails[["igMminus"]], place = patient.dirs[k])
  allIsotypeDetails[["igMplus"]] <-cbind(allIsotypeDetails[["igMplus"]], place = patient.dirs[k])
  
  names(allIsotypeDetails[["igMminus"]]) =names(TOTAL.allIsotypeDetails[["igMminus"]])
  names(allIsotypeDetails[["igMplus"]]) =names(TOTAL.allIsotypeDetails[["igMplus"]])
  
  
  ####   add the current patiant list into the list that saves all patiants togther
  TOTAL.allIsotypeDetails[["igMminus"]] <- rbind(TOTAL.allIsotypeDetails[["igMminus"]],allIsotypeDetails[["igMminus"]])
  TOTAL.allIsotypeDetails[["igMplus"]] <- rbind(TOTAL.allIsotypeDetails[["igMplus"]],allIsotypeDetails[["igMplus"]])
  
} 

# write to file 
names(TOTAL.allIsotypeDetails[["igMminus"]])<- c("first_generation_child_num","All_generation_child_num","firstVSallChilden", "place")
names(TOTAL.allIsotypeDetails[["igMplus"]])<- c("first_generation_child_num","All_generation_child_num","firstVSallChilden","place")

write.csv(TOTAL.allIsotypeDetails[["igMminus"]],file = "SPLEEN_Nodes_igMminus.csv")
write.csv(TOTAL.allIsotypeDetails[["igMplus"]],file = "SPLEEN_Nodes_igMplus.csv")

