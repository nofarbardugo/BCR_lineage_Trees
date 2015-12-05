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

isIsotypeLeaf <-function (node)
  # this function get a node and return if its a leaf
  # Args: 
  #   node- an node of the tree   
  # Returns: True id its a leaf
{
  s <- isotypesNames.df[which(isotypesNames.df[,"head"] == node),"head"]
  return (startsWith(s, 'L', trim=TRUE)==1)
}

GoOverTree<-function(sameAsFather,node)
  # this function goes over the tree and save for each node type how much directes kids it has
  # Args: 
  #   node- an edge of the tree   
  # Returns: 
{
  
  # get node isotype
  type <-getIsotypeName(node)
  
  # if condition-   what we work on[the condition, if the condition is true give me
  #                                                                  the Relevant information]
  # get all the suns
  suns <- edges.df [edges.df [,"from"]==node,"to"]
  sunsNum <-length(suns)
  
  countLeafChildFromTheSameIsotype <- 0 
  countBrunchFromTheSameIsotype <- 0 

  syn <- 0
  nonSyn <- 0
  syn_cdr <- 0
  syn_fr <-  0
  nonSyn_cdr <- 0
  nonSyn_fr <- 0 
  
  # calc mesurmeant for the suns
  p <-1
  while(p <=sunsNum)
  {
    if(suns[p]!= "VJ_1")
    {
      
      # if this is the same isotype type:
      if(type == getIsotypeName(suns[p])) 
      {
        # add 1 to the edge count
        countBrunchFromTheSameIsotype <- countBrunchFromTheSameIsotype + 1
        
        # check if its a leaf -  # add 1 to the leaf count
        if(isIsotypeLeaf(suns[p])==T)
        {
          countLeafChildFromTheSameIsotype <- countLeafChildFromTheSameIsotype + 1
        }
        # not a leaf - calc sun node details and it to the count
        else
        {
          
          sonDetails.vector <- GoOverTree(T,suns[p])
          
          # add the details of the sun:
          countLeafChildFromTheSameIsotype <- countLeafChildFromTheSameIsotype + sonDetails.vector[1] 
          countBrunchFromTheSameIsotype <- countBrunchFromTheSameIsotype + sonDetails.vector[2] 
          syn <- syn + sonDetails.vector[3] 
          nonSyn <- nonSyn +  sonDetails.vector[4] 
          syn_cdr <- syn_cdr +  sonDetails.vector[5] 
          syn_fr <- syn_fr +  sonDetails.vector[6] 
          nonSyn_cdr <- nonSyn_cdr +  sonDetails.vector[7] 
          nonSyn_fr <- nonSyn_fr +  sonDetails.vector[8] 
          
        }
        
        # get mutations information from cure node to each of its child 
        syn <- syn +  edges.df [(edges.df [,"to"]==suns[p]),"synonyms"]
        nonSyn <- nonSyn +  edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms"]
        syn_cdr <- syn_cdr +  edges.df [(edges.df [,"to"]==suns[p]),"synonyms_CDR"]
        syn_fr <- syn_fr +  edges.df [(edges.df [,"to"]==suns[p]),"synonyms_FR"]
        nonSyn_cdr <- nonSyn_cdr +edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_CDR"]
        nonSyn_fr <- nonSyn_fr + edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_FR"]
        
      }
      # its not the same type and if not a leaf - not enter to the node calc (has it one calc)
      # if its a leaf ignore it
      else if(isIsotypeLeaf(suns[p])==F)
      {
        GoOverTree(F,suns[p])     
      }
    }
   
    p <- p+1 
    
  }# end:while(p <=sunsNum)
  
  # - add a new line to "df.allIsotypeDetails" datafream (to internal or not)
  data <- df.allIsotypeDetails
  depth <- syn+nonSyn
  ratio <-0
  if(depth >0)
  {
    ratio <-nonSyn/depth
  }
  
  d <- data.frame(type,sameAsFather,countLeafChildFromTheSameIsotype,countBrunchFromTheSameIsotype,syn,nonSyn,depth,ratio,syn_cdr,syn_fr,nonSyn_cdr,nonSyn_fr, stringsAsFactors=F)
  data <- rbind(data,d) 
  
  # save "data" in the outer "df.allIsotypeDetails" data frame
  assign("df.allIsotypeDetails",data, envir = .GlobalEnv)
  
  # if this node is the same as its father - send it all data to its father
  if(sameAsFather==T)
  {
    return(c(countLeafChildFromTheSameIsotype,countBrunchFromTheSameIsotype,syn,nonSyn,syn_cdr,syn_fr,nonSyn_cdr,nonSyn_fr))
  }

} 


############## main ###########################

require(gdata)
# path for trees folder
tree.dir <- paste0(source.dir,"/spleenHumanTree/")

# get list of folders
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)

# create datafream that will save for each isotype :children from with the same isotype number, Synonymus mutataion,
#                                                   NON Synonymus mutataion,Depth = s+ ns, Ratio = ns/(ns + s)   
#                                                   synonyms_CDR mutation, synonyms_FR mutation,
#                                                   nonSynonyms_CDR mutations and nonSynonyms_FR mutations
df.allIsotypeDetails <- data.frame(IsotypeName = c(),IsInternal = c(),ChildrenNum = c(),edgesNum = c(), Synonymus =c(),nonSynonymus =c(),
                                   depth = c(),Ratio = c(), synonyms_CDR= c(),synonyms_FR = c(),
                                   nonSynonyms_CDR = c(),nonSynonyms_FR = c())

# going over each patient
for(k in 2:2)#length(patient.dirs))
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
          # count modes only on trees that have at last 10 leave
         if(!(is.null(totalLeavesType.fd)) && (totalLeavesType.fd[j,"sum_Leaves"]>10) )
          {
          
          # count modes only on tree plot only trees that have at last 10 leave and only 1 type of isotipes 
       #   if(!(is.null(totalLeavesType.fd)) && (totalLeavesType.fd[j,"type_num"]==1)
      #       && (totalLeavesType.fd[j,"sum_Leaves"]>10) )
      #    {  
            
            print ("in")
            GoOverTree(F,roots[j]) 
            
          }
          j<-j + 1
        }
      }
    }
  }
  
}

# write to file 

names(df.allIsotypeDetails)<- c ("Isotype Name","Is internal ","Children Number","Number of edges from node to leaf",
                                 "Synonymus","nonSynonymus","depth","Ratio", "synonyms_CDR","synonyms_FR",
                                  "nonSynonyms_CDR","nonSynonyms_FR")

write.csv(df.allIsotypeDetails,file = "spleen_SPL_mutation_num_per_isotype_with_certain_children_numberWithInternal.csv")

# if the script run on NOCDR3 mut mode
