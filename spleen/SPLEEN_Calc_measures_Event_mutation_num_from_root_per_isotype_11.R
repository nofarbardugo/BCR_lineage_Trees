# set directory : change path depending on work environment
#source.dir <- 'G:/project/files/pfizer'
#source.dir<- '/media/nofar/F6F9-0D76/project/files/pfizer'
#source.dir <-'C:/Users/Nofar/Desktop/trees'
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

GoOverTree<-function(node,syn, nonSyn,syn_cdr,syn_fr,nonSyn_cdr,nonSyn_fr)
  # this function goes over the tree and save for each node type how much directes kids it has
  # Args: 
  #   node- an edge of the tree   
  # Returns: 
{
  
  # save mutations for this isotype
  type <-getIsotypeName(node)
  list  <- allIsotypeDetails
  depth <- syn+nonSyn
  ratio <-0
  if(depth >0)
  {
    ratio <-nonSyn/depth
  }

  list[[type]] <- rbind(allIsotypeDetails[[type]],c(syn,nonSyn,depth,ratio,syn_cdr,syn_fr,nonSyn_cdr,nonSyn_fr)) 
  
  # save "df_TypeData" in the outer "fd.leavesType" data frame
  assign("allIsotypeDetails",list, envir = .GlobalEnv)
  
  # if condition-   what we work on[the condition, if the condition is true give me
  #                                                                  the Relevant information]
  # get all the suns
  suns <- edges.df [edges.df [,"from"]==node,"to"]
  
  sunsNum <-length(suns)
  p <-1
  # if this is not a leaf - calc the children mut num
  while(p <=sunsNum)
  {
     # get syn and non syn from cure node to each of its child 
     directSyn<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms"]
     directNonS<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms"]
    
     directSyn_CDR<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms_CDR"]
     directNonS_CDR<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_CDR"]
     directSyn_FR<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms_FR"]
     directNonS_FR<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_FR"]
      
    GoOverTree(suns[p],syn + directSyn,nonSyn + directNonS,syn_cdr +directSyn_CDR ,
               syn_fr + directSyn_FR,nonSyn_cdr+ directNonS_CDR, nonSyn_fr + directNonS_FR)
    p <- p+1 
   }
}  


############## main ###########################

# path for trees folder
tree.dir <- paste0(source.dir,"/spleenHumanTree/")

# get list of folders
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)

# data fream that save details about each tree leaves - just for get details
#fd.TotalLeaves<- data.frame(treeNUM=c() igMminus = c(),igMplus = c(),sum_Leaves =c(),type_num = c())

# create datafream for each isotype that will save: Synonymus mutataion,NON Synonymus mutataion,
#                                                   Depth = s+ ns, Ratio = ns/(ns + s)   
# save all data freams into a list for save all togther
TOTAL.allIsotypeDetails <- list()
TOTAL.allIsotypeDetails[["igMminus"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                                synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                                nonSynonyms_FR = c(),place = c())
TOTAL.allIsotypeDetails[["igMplus"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                                synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(), 
                                                nonSynonyms_FR = c(),place = c())


# going over each patient
for(k in 1:length(patient.dirs))
{

  # get list of folders
  dirs <- list.files(paste0(tree.dir,patient.dirs[k]),full.names = F)
  
  # create datafream for each isotype that will save: Synonymus mutataion,NON Synonymus mutataion,
  #                                                   Depth = s+ ns, Ratio = ns/(ns + s)   
  # save all data freams into a list for save all togther
  allIsotypeDetails <- list()
  allIsotypeDetails[["igMminus"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                            synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                            nonSynonyms_FR = c())
  
  allIsotypeDetails[["igMplus"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                            synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                            nonSynonyms_FR = c())
  
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
         if(!(is.null(totalLeavesType.fd)) && (totalLeavesType.fd[j,"sum_Leaves"]>10))
         {        
            print ("in")
            GoOverTree(roots[j],0,0,0,0,0,0) 
            
          }
          j<-j + 1
        }
      }
    }
  }
  
  # add "pbl" or "spl" place into the list
  allIsotypeDetails[["igMminus"]] <- cbind(allIsotypeDetails[["igMminus"]], place = patient.dirs[k])
  allIsotypeDetails[["igMplus"]] <- cbind(allIsotypeDetails[["igMplus"]], place = patient.dirs[k])
  
  names(allIsotypeDetails[["igMminus"]]) =names(TOTAL.allIsotypeDetails[["igMminus"]])
  names(allIsotypeDetails[["igMplus"]]) =names(TOTAL.allIsotypeDetails[["igMplus"]])
 
  ####   add the current patiant list into the list that saves all patiants togther
  TOTAL.allIsotypeDetails[["igMminus"]] <- rbind(TOTAL.allIsotypeDetails[["igMminus"]],allIsotypeDetails[["igMminus"]])
  TOTAL.allIsotypeDetails[["igMplus"]] <- rbind(TOTAL.allIsotypeDetails[["igMplus"]],allIsotypeDetails[["igMplus"]])

}

# write to file 

names(TOTAL.allIsotypeDetails[["igMminus"]])<- c ("Synonymus","nonSynonymus","depth","Ratio",
                                                  "synonyms_CDR","synonyms_FR" ,"nonSynonyms_CDR",
                                                  "nonSynonyms_FR","place")

names(TOTAL.allIsotypeDetails[["igMplus"]])<- c ("Synonymus","nonSynonymus","depth","Ratio",
                                                 "synonyms_CDR","synonyms_FR" ,"nonSynonyms_CDR",
                                                 "nonSynonyms_FR","place")

write.csv(TOTAL.allIsotypeDetails[["igMminus"]],file = "SPLEEN_event_igMminus.csv")
write.csv(TOTAL.allIsotypeDetails[["igMplus"]],file = "SPLEEN_event_igMplus.csv")
