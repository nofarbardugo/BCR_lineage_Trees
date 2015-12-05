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

getTimePoint <-function (node)
  # this function get a node and return its timePoint
  # Args: 
  #   node- an node of the tree   
  # Returns: time point of the node
{
  x <-isotypesNames.df[which(isotypesNames.df[,"head"] == node),"TimePoint"]
  
  # if there is no time point
  if(x == '-')
  {
    return (0)
  }
  
  # get only the number
  return(as.numeric(substr(x,start = 4, stop = 5)))
}

GoOverTree<-function(node,syn, nonSyn)
  # this function goes over the tree and save for each node syn non syn from root
  # Args: 
  #   node- an edge of the tree   
  # Returns: 
{
  
  # if condition-   what we work on[the condition, if the condition is true give me
  #                                                                  the Relevant information]
  # get all the suns
  suns <- edges.df [edges.df [,"from"]==node,"to"]
  
  sunsNum <-length(suns)
  
  # if its a leaf
  if(sunsNum==0)
  {
    # save mutations for this isotype
    type <-getIsotypeName(node)
    time <- getTimePoint(node)
    list  <- allIsotypeDetails
    depth <- syn+nonSyn
    ratio <-0
    if(depth >0)
    {
      ratio <-nonSyn/depth
    }
    
    list[[type]] <- rbind(allIsotypeDetails[[type]],c(syn,nonSyn,depth,ratio,time)) 
    
    # save "df_TypeData" in the outer "fd.leavesType" data frame
    assign("allIsotypeDetails",list, envir = .GlobalEnv)
    
  }
 
  p <-1
  # if this is not a leaf - calc the children mut num
  while(p <=sunsNum)
  {
     # get syn and non syn from cure node to each of its child 
     directSyn<- edges.df [(edges.df [,"to"]==suns[p]),"synonyms"]
     directNonS<- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms"]
      
    GoOverTree(suns[p],syn + directSyn,nonSyn + directNonS)
    p <- p+1
   }
}  


############## main ###########################

# path for trees folder
tree.dir <- paste0(source.dir,"/pfizerHumanTree/")

# get list of folders
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)

# data fream that save details about each tree leaves - just for get details
fd.TotalLeaves<- data.frame(treeNUM=c() ,IGHA = c(),IGHE =c(),"IGHG-1" = c(),"IGHG-2" = c()
                            ,IGHM = c(), IGHD = c(),naive_IGHM = c(),sum_Leaves =c(),type_num = c())

# create datafream for each isotype that will save: Synonymus mutataion,NON Synonymus mutataion,
#                                                   Depth = s+ ns, Ratio = ns/(ns + s)   
# save all data freams into a list for save all togther
TOTAL.allIsotypeDetails <- list()
TOTAL.allIsotypeDetails[["IGHA"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),timePoint = c(),stringsAsFactors=FALSE )
TOTAL.allIsotypeDetails[["IGHE"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),timePoint = c(),stringsAsFactors=FALSE )
TOTAL.allIsotypeDetails[["IGHG-1"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),timePoint = c(),stringsAsFactors=FALSE )
TOTAL.allIsotypeDetails[["IGHG-2"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),timePoint = c(),stringsAsFactors=FALSE )
TOTAL.allIsotypeDetails[["IGHM"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),timePoint = c(),stringsAsFactors=FALSE )
TOTAL.allIsotypeDetails[["IGHD"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),timePoint = c(),stringsAsFactors=FALSE )
TOTAL.allIsotypeDetails[["naive_IGHM"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),timePoint = c(),stringsAsFactors=FALSE )


# going over each patient
for(k in 1:length(patient.dirs))
{
  # get list of folders
  dirs <- list.files(paste0(tree.dir,patient.dirs[k]),full.names = F)
  
  # create datafream for each isotype that will save: Synonymus mutataion,NON Synonymus mutataion,
  #                                                   Depth = s+ ns, Ratio = ns/(ns + s)   
  # save all data freams into a list for save all togther
  allIsotypeDetails <- list()
  allIsotypeDetails[["IGHA"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),timePoint = c(),stringsAsFactors=FALSE )
  allIsotypeDetails[["IGHE"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),timePoint = c(),stringsAsFactors=FALSE )
  allIsotypeDetails[["IGHG-1"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),timePoint = c(),stringsAsFactors=FALSE )
  allIsotypeDetails[["IGHG-2"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),timePoint = c(),stringsAsFactors=FALSE )
  allIsotypeDetails[["IGHM"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),timePoint = c(),stringsAsFactors=FALSE )
  allIsotypeDetails[["IGHD"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),timePoint = c(),stringsAsFactors=FALSE )
  allIsotypeDetails[["naive_IGHM"]] <- data.frame(Synonymus =c(),nonSynonymus =c(),depth = c(),Ratio = c(),timePoint = c(),stringsAsFactors=FALSE )
  
  
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
            GoOverTree(roots[j],0,0) 
            
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
  names(allIsotypeDetails[["naive_IGHM"]]) =names(TOTAL.allIsotypeDetails[["naive_IGHM"]])
  names(allIsotypeDetails[["IGHD"]]) =names(TOTAL.allIsotypeDetails[["IGHD"]])
  ####   add the current patiant list into the list that saves all patiants togther
  TOTAL.allIsotypeDetails[["IGHA"]] <- rbind(TOTAL.allIsotypeDetails[["IGHA"]],allIsotypeDetails[["IGHA"]])
  TOTAL.allIsotypeDetails[["IGHE"]] <- rbind(TOTAL.allIsotypeDetails[["IGHE"]],allIsotypeDetails[["IGHE"]])
  TOTAL.allIsotypeDetails[["IGHG-1"]] <- rbind(TOTAL.allIsotypeDetails[["IGHG-1"]],allIsotypeDetails[["IGHG-1"]])
  TOTAL.allIsotypeDetails[["IGHG-2"]] <- rbind(TOTAL.allIsotypeDetails[["IGHG-2"]],allIsotypeDetails[["IGHG-2"]])
  TOTAL.allIsotypeDetails[["IGHM"]] <- rbind(TOTAL.allIsotypeDetails[["IGHM"]],allIsotypeDetails[["IGHM"]])
  TOTAL.allIsotypeDetails[["IGHD"]] <- rbind(TOTAL.allIsotypeDetails[["IGHD"]],allIsotypeDetails[["IGHD"]])
  TOTAL.allIsotypeDetails[["naive_IGHM"]] <- rbind(TOTAL.allIsotypeDetails[["naive_IGHM"]],allIsotypeDetails[["naive_IGHM"]])
  
}

jpeg(filename=  paste0('/home/bardugn1/flu2/Measurements/event-ratioVSdapthPlot-AllPatientsLEAF.jpeg'),width = 1400, height = 800)
colorsL <- vector()
namesL <- vector()
c <-1
for(m in 1:7)
{
  if (length(TOTAL.allIsotypeDetails[[m]])!= 0)
  {  
    names (TOTAL.allIsotypeDetails[[m]]) <-vectorMesNames
    maxDepth <- max(TOTAL.allIsotypeDetails[[m]][,"Depth"])
    plotMeanRatio.vector <-c()
    plotMeanDepth.vector <-c()
    plotSD <-c()
    for(j in 1:maxDepth)
    {
      plotMeanRatio.vector[j]<-mean(TOTAL.allIsotypeDetails[[m]][,"Ratio"][TOTAL.allIsotypeDetails[[m]][,"Depth"]==j] )
      plotMeanDepth.vector[j]<-j
      sqrtOFlength<- sqrt(length(TOTAL.allIsotypeDetails[[m]][,"Ratio"][TOTAL.allIsotypeDetails[[m]][,"Depth"]==j]))
      plotSD <- sd(TOTAL.allIsotypeDetails[[m]][,"Ratio"][TOTAL.allIsotypeDetails[[m]][,"Depth"]==j])/sqrtOFlength
    }
    
    if(m==1)
    { 
      plot(x = plotMeanDepth.vector,y = plotMeanRatio.vector,type = 'l',xlim = c(0,100),ylim = c(0,1),lwd = 2.5,xlab = "Depth: NS + S",ylab = "Ratio: NS/(NS + S)")
      segments(plotMeanDepth.vector, plotMeanRatio.vector-plotSD.vector,plotMeanDepth.vector, plotMeanRatio.vector+plotSD.vector,col=m)
      epsilon = 0.1
      segments(plotMeanDepth.vector-epsilon,plotMeanRatio.vector-plotSD.vector,plotMeanDepth.vector+epsilon,plotMeanRatio.vector-plotSD.vector,col=m)
      segments(plotMeanDepth.vector-epsilon,plotMeanRatio.vector+plotSD.vector,plotMeanDepth.vector+epsilon,plotMeanRatio.vector+plotSD.vector,col=m)

    }
    else
    {
      points(x = plotMeanDepth.vector,y = plotMeanRatio.vector,type = 'l',col=m,lwd = 2.5)
      segments(plotMeanDepth.vector, plotMeanRatio.vector-plotSD.vector,plotMeanDepth.vector, plotMeanRatio.vector+plotSD.vector,col=m)
      epsilon = 0.1
      segments(plotMeanDepth.vector-epsilon,plotMeanRatio.vector-plotSD.vector,plotMeanDepth.vector+epsilon,plotMeanRatio.vector-plotSD.vector,col=m)
      segments(plotMeanDepth.vector-epsilon,plotMeanRatio.vector+plotSD.vector,plotMeanDepth.vector+epsilon,plotMeanRatio.vector+plotSD.vector,col=m)
      
    } 
    colorsL[c]<-m
    namesL[c]<-names(allIsotypeDetails[m])
    c <- c+1
  }
}

legend ( "topright" ,4 , legend =namesL, lwd = 4:4,
         col = colorsL , lty =1:1 )
dev.off()