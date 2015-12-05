# set directory : change path depending on work environment

source.dir <-'/home/bardugn1/pfizer' # pep4
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

    p<-1
    
    while(p <=sunsNum)
    {
      
      # get edges typ 
      type <-getEdgeName(node,suns[p])
      
      if(suns[p] !="VJ_1")
      {
        syn<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms"]
        nonS<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms"]
        
        frSyn<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms_FR"]
        frNonS<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_FR"]
        
        cdrSyn<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms_CDR"]
        cdrNonS<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_CDR"]
        
        depth <- syn+nonS
        ratio <-0
        if(depth >0)
        {
          ratio <-nonS/depth
        }
        
        data  <- df.allIsotypeDetails
        d <- data.frame(type,syn,nonS,depth,ratio,cdrSyn,frSyn,cdrNonS,frNonS, stringsAsFactors=F)
        data <- rbind(data,d) 
        
        # save "list" in the outer "allEdgeDetails" data frame
        assign("df.allIsotypeDetails",data, envir = .GlobalEnv)
        
        GoOverTree(suns[p])
        
      }
      
      p <- p+1
      
    }
    

  }
  
}  


############## main ###########################

# path for trees folder
tree.dir <- paste0(source.dir,"/pfizerHumanTree15oct2015/")

# get list of folders
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)

# create datafream that will save for each edge :   Synonymus mutataion,
#                                                   NON Synonymus mutataion,Depth = s+ ns, Ratio = ns/(ns + s)   
#                                                   synonyms_CDR mutation, synonyms_FR mutation,
#                                                   nonSynonyms_CDR mutations and nonSynonyms_FR mutations
df.allIsotypeDetails <- data.frame(EdgeName = c(), Synonymus =c(),nonSynonymus =c(),
                                   depth = c(),Ratio = c(), synonyms_CDR= c(),synonyms_FR = c(),
                                   nonSynonyms_CDR = c(),nonSynonyms_FR = c())



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
            GoOverTree(roots[j]) 
            
          }
          j<-j + 1
        }
      }
    }
  }
  
}

# write to file 

names(df.allIsotypeDetails)<- c ("Edge Name",
                                 "Synonymus","nonSynonymus","depth","Ratio", "synonyms_CDR","synonyms_FR",
                                  "nonSynonyms_CDR","nonSynonyms_FR")

write.csv(df.allIsotypeDetails,file = "PFIZER_Edges_mutation_num_per_edge3.csv")

# if the script run on NOCDR3 mut mode
