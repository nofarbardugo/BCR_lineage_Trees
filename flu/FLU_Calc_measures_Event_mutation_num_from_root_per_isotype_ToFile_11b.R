# set directory : change path depending on work environment

source.dir <-'/home/bardugn1/flu2' # pep4
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
    
    directSyn_FR<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms_FR"]
    directNonS_FR<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_FR"]
    
    # for check of mutation cdr with mutation in cdr 3:
      directSyn_CDR<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms_CDR"]
      directNonS_CDR<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_CDR"]
    
    # for check of mutation cdr without mutation in cdr 3
    #directSyn_CDR3<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms_CDR3"]
    #directNonS_CDR3<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_CDR3"]
    #directSyn<- directSyn - directSyn_CDR3
    #directNonS<-directNonS - directNonS_CDR3
    
   # directSyn_CDR<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms_CDR_NO3"]
  #  directNonS_CDR<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_CDR_NO3"]
    
    
    GoOverTree(suns[p],syn + directSyn,nonSyn + directNonS,syn_cdr +directSyn_CDR ,
               syn_fr + directSyn_FR,nonSyn_cdr+ directNonS_CDR, nonSyn_fr + directNonS_FR)
    p <- p+1 
  }
} 



############## main ###########################
#### for each  isotype get a file with all its node details
###############################################

# path for trees folder
tree.dir <- paste0(source.dir,"/fluHumanTree15oct2015/")

# get list of folders
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)
patient.dirs <- patient.dirs[-grep('IB',patient.dirs)]
# data fream that save details about each tree leaves - just for get details
fd.TotalLeaves<- data.frame(treeNUM=c() ,IGHA = c(),IGHE =c(),"IGHG-1" = c(),"IGHG-2" = c()
                            ,IGHM = c(), IGHD = c(),naive_IGHM = c(),sum_Leaves =c(),type_num = c())

# create datafream for each isotype that will save: Synonymus mutataion,NON Synonymus mutataion,
#                                                   Depth = s+ ns, Ratio = ns/(ns + s)   
# save all data freams into a list for save all togther
TOTAL.allIsotypeDetails <- list()
TOTAL.allIsotypeDetails[["IGHA"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                                synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                                nonSynonyms_FR = c())
TOTAL.allIsotypeDetails[["IGHE"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                                synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                                nonSynonyms_FR = c())
TOTAL.allIsotypeDetails[["IGHG-1"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                                  synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                                  nonSynonyms_FR = c())
TOTAL.allIsotypeDetails[["IGHG-2"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                                  synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                                  nonSynonyms_FR = c())
TOTAL.allIsotypeDetails[["IGHM"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                                synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                                nonSynonyms_FR = c())
TOTAL.allIsotypeDetails[["IGHD"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                                      synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                                      nonSynonyms_FR = c())
TOTAL.allIsotypeDetails[["naive_IGHM"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                                      synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                                      nonSynonyms_FR = c())


# going over each patient
for(k in 1:length(patient.dirs))
{
  # get list of folders
  dirs <- list.files(paste0(tree.dir,patient.dirs[k]),full.names = F)
  
  # create datafream for each isotype that will save: Synonymus mutataion,NON Synonymus mutataion,
  #                                                   Depth = s+ ns, Ratio = ns/(ns + s)   
  # save all data freams into a list for save all togther
  allIsotypeDetails <- list()
  allIsotypeDetails[["IGHA"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                            synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                            nonSynonyms_FR = c())
  allIsotypeDetails[["IGHE"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                            synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                            nonSynonyms_FR = c())
  allIsotypeDetails[["IGHG-1"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                              synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                              nonSynonyms_FR = c())
  allIsotypeDetails[["IGHG-2"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                              synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                              nonSynonyms_FR = c())
  allIsotypeDetails[["IGHM"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                            synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                            nonSynonyms_FR = c())
  allIsotypeDetails[["IGHD"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
                                                  synonyms_CDR= c(),synonyms_FR = c(),nonSynonyms_CDR = c(),
                                                  nonSynonyms_FR = c())
  allIsotypeDetails[["naive_IGHM"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),
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
          # count modes only on trees that have at last 10 leave and 2 types of isotipes 
          if(!(is.null(totalLeavesType.fd)) && (totalLeavesType.fd[j,"type_num"]>1)
             && (totalLeavesType.fd[j,"sum_Leaves"]>10) )
          {
            print ("in")
            GoOverTree(roots[j],0,0,0,0,0,0) 
            
          }
          j<-j + 1
        }
      }
    }
  }
  
  
  names(allIsotypeDetails[["IGHA"]]) =names(TOTAL.allIsotypeDetails[["IGHA"]])
  names(allIsotypeDetails[["IGHE"]]) =names(TOTAL.allIsotypeDetails[["IGHE"]])
  names(allIsotypeDetails[["IGHG-1"]]) =names(TOTAL.allIsotypeDetails[["IGHG-1"]])
  names(allIsotypeDetails[["IGHG-2"]]) =names(TOTAL.allIsotypeDetails[["IGHG-2"]])
  names(allIsotypeDetails[["IGHM"]]) =names(TOTAL.allIsotypeDetails[["IGHM"]])
 # names(allIsotypeDetails[["naive_IGHM"]]) =names(TOTAL.allIsotypeDetails[["naive_IGHM"]])
  names(allIsotypeDetails[["IGHD"]]) =names(TOTAL.allIsotypeDetails[["IGHD"]])
  ####   add the current patiant list into the list that saves all patiants togther
  TOTAL.allIsotypeDetails[["IGHA"]] <- rbind(TOTAL.allIsotypeDetails[["IGHA"]],allIsotypeDetails[["IGHA"]])
  TOTAL.allIsotypeDetails[["IGHE"]] <- rbind(TOTAL.allIsotypeDetails[["IGHE"]],allIsotypeDetails[["IGHE"]])
  TOTAL.allIsotypeDetails[["IGHG-1"]] <- rbind(TOTAL.allIsotypeDetails[["IGHG-1"]],allIsotypeDetails[["IGHG-1"]])
  TOTAL.allIsotypeDetails[["IGHG-2"]] <- rbind(TOTAL.allIsotypeDetails[["IGHG-2"]],allIsotypeDetails[["IGHG-2"]])
  TOTAL.allIsotypeDetails[["IGHM"]] <- rbind(TOTAL.allIsotypeDetails[["IGHM"]],allIsotypeDetails[["IGHM"]])
  TOTAL.allIsotypeDetails[["IGHD"]] <- rbind(TOTAL.allIsotypeDetails[["IGHD"]],allIsotypeDetails[["IGHD"]])
  #TOTAL.allIsotypeDetails[["naive_IGHM"]] <- rbind(TOTAL.allIsotypeDetails[["naive_IGHM"]],allIsotypeDetails[["naive_IGHM"]])
  
}


# write to file 

names(TOTAL.allIsotypeDetails[["IGHA"]])<- c ("Synonymus","nonSynonymus","depth","Ratio",
                                              "synonyms_CDR","synonyms_FR" ,"nonSynonyms_CDR",
                                              "nonSynonyms_FR")

names(TOTAL.allIsotypeDetails[["IGHE"]])<- c ("Synonymus","nonSynonymus","depth","Ratio",
                                              "synonyms_CDR","synonyms_FR" ,"nonSynonyms_CDR",
                                              "nonSynonyms_FR")

names(TOTAL.allIsotypeDetails[["IGHG-1"]])<- c ("Synonymus","nonSynonymus","depth","Ratio",
                                              "synonyms_CDR","synonyms_FR" ,"nonSynonyms_CDR",
                                              "nonSynonyms_FR")

names(TOTAL.allIsotypeDetails[["IGHG-2"]])<- c ("Synonymus","nonSynonymus","depth","Ratio",
                                                    "synonyms_CDR","synonyms_FR" ,"nonSynonyms_CDR",
                                                    "nonSynonyms_FR")

names(TOTAL.allIsotypeDetails[["IGHM"]])<- c ("Synonymus","nonSynonymus","depth","Ratio",
                                              "synonyms_CDR","synonyms_FR" ,"nonSynonyms_CDR",
                                              "nonSynonyms_FR")

names(TOTAL.allIsotypeDetails[["IGHD"]])<- c ("Synonymus","nonSynonymus","depth","Ratio",
                                              "synonyms_CDR","synonyms_FR" ,"nonSynonyms_CDR",
                                              "nonSynonyms_FR")

#names(TOTAL.allIsotypeDetails[["naive_IGHM"]])<- c ("Synonymus","nonSynonymus","depth","Ratio",
 #                                                   "synonyms_CDR","synonyms_FR" ,"nonSynonyms_CDR",
#                                                    "nonSynonyms_FR")


#write.csv(TOTAL.allIsotypeDetails[["IGHA"]],file = "FLU_event_all_IGHA.csv")
#write.csv(TOTAL.allIsotypeDetails[["IGHE"]],file = "FLU_event_all_IGHE.csv")
#write.csv(TOTAL.allIsotypeDetails[["IGHG-1"]],file = "FLU_event_all_IGHG1.csv")
#write.csv(TOTAL.allIsotypeDetails[["IGHG-2"]],file = "FLU_event_all_IGHG2.csv")
#write.csv(TOTAL.allIsotypeDetails[["IGHM"]],file = "FLU_event_all_IGHM.csv")
#write.csv(TOTAL.allIsotypeDetails[["IGHD"]],file = "FLU_event_all_naive_IGHD.csv")

write.csv(TOTAL.allIsotypeDetails[["IGHA"]],file = "FLU_event_IGHA.csv")
write.csv(TOTAL.allIsotypeDetails[["IGHE"]],file = "FLU_event_IGHE.csv")
write.csv(TOTAL.allIsotypeDetails[["IGHG-1"]],file = "FLU_event_IGHG1.csv")
write.csv(TOTAL.allIsotypeDetails[["IGHG-2"]],file = "FLU_event_IGHG2.csv")
write.csv(TOTAL.allIsotypeDetails[["IGHM"]],file = "FLU_event_IGHM.csv")
write.csv(TOTAL.allIsotypeDetails[["IGHD"]],file = "FLU_event_naive_IGHD.csv")
write.csv(TOTAL.allIsotypeDetails[["naive_IGHM"]],file = "FLU_event_naive_IGHM.csv")



