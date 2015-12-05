# set directory : change path depending on work environment
#source.dir <- 'G:/project/files/pfizer'
#source.dir<- '/media/nofar/F6F9-0D76/project/files/pfizer'
#source.dir <-'C:/Users/Nofar/Desktop/trees'
source.dir <-'/home/bardugn1/pfizer' # pep4
#source.dir <- "/home/nofar/Desktop/spleen"
setwd(source.dir)
require("ape")

#givEdgeShape <-function(type)
#{
#  return (switch(type, "same" =1,"dif"= 3))

#}

givEdgeShape <-function(type)
{
  return (switch(type, "equal" =1,"cdr"=2,"fr"= 10))
  
}

givEdgeColor <-function(type)
{
  return (switch(type, "equal" ="black","syn"= "tomato",nonSyn = "purple" )) 
}

giveColor <- function(isotypeName)
  # the function get Isotype  and retern is color
  # Args: 
  #   isotypeName- isotype to get it's color   
  # Returns: the isotype color
{
  return (switch(isotypeName, "IGHM"="blue", "IGHG"="red", "IGHD"="purple","IGHA" ="orange","naive_IGHM"="green","plum2"))
}

GoOverTree<-function(node)
  # this function goes over the tree and create a newick format
  # Args: 
  #   node- an edge of the tree   
  # Returns: the nweic format of the node and it's suns
{
  
  # if condition-   what we work on[the condition, if the condition is true give me
  #                                                                  the Relevant information]
  # get all the suns 
  suns <- edges.df [edges.df [,"from"]==node,"to"]
  
  # if we reached to the leafs
  if(length(suns)==0)
  {
    # return the  name and the distance of the leaf
    return(paste0(getIsotypeName(node),":", edges.df[which(edges.df[,"to"]==node) ,"distance"] ))
    
  }
  
  # if we reached to an edge with suns
  else
  {
    # (a for loop) apply the "GoOverTree" function on all the suns to get their newick
    # (before we will get the fathers isotype) 
    sunsNewick.vector <- sapply(suns,GoOverTree)
    
    # add an extra node (with garbidge value) for resolve the singlton problem
    if(length(suns)==1)
    {
      sunsNewick.vector <- sapply(suns,GoOverTree)
      sunsNewickString <- paste0(sunsNewick.vector[1], ",demo:0")
    }
    else
    {
      sunsNewickString <- paste0(sunsNewick.vector, collapse = ',')
    }
    
    # build the node newick format
    nodeNewickString <-paste0(getIsotypeName(node),":", edges.df[which(edges.df[,"to"]==node) ,"distance"] )
    
    # return the whole newick format : the node with all is suns
    return (paste0("(",sunsNewickString,")",nodeNewickString))
  }
}

GetEdgsDetails<-function(node)
  # this function goes over the tree and create a newick format
  # Args: 
  #   node- an edge of the tree   
  # Returns: the nweic format of the node and it's suns
{
  
  # if condition-   what we work on[the condition, if the condition is true give me
  #                                                                  the Relevant information]
  # get all the suns 
  suns <- edges.df [edges.df [,"from"]==node,"to"]
  
  # if we reached to the leafs
  if(length(suns)==0)
  {
    vec <-c (sameIsotype.vector,"same")
    
    if(node!= "VJ_1")
      # save if the is an isotype switch or not
      if(getIsotypeName(node) != getIsotypeName(edges.df[which(edges.df[,"to"]==node) ,"from"]))
      {
        vec <- c(sameIsotype.vector,"dif")
        
      }
    
    # save "vec" in the outer "sameIsotype.vector" vector
    assign("sameIsotype.vector",vec, envir = .GlobalEnv) 
    
    vec <-c(majorMut.vector,"equal")
    
    syn = edges.df[which(edges.df[,"to"]==node) ,"synonyms"]
    nonSyn = edges.df[which(edges.df[,"to"]==node) ,"nonSynonyms"]
    
    # save if the is syn or non syn
    if(syn > nonSyn)
    {
      vec <- c(majorMut.vector,"syn")
    }
    else if(syn < nonSyn)
      
    {
      vec <- c(majorMut.vector,"nonSyn")
    }
    
    # save "vec" in the outer "sameIsotype.vector" vector
    assign("majorMut.vector",vec, envir = .GlobalEnv) 
    
    # 
    
    vec <-c(placeMut.vector,"equal")
    
    cdr = edges.df[which(edges.df[,"to"]==node) ,"synonyms_CDR"] + edges.df[which(edges.df[,"to"]==node) ,"nonSynonyms_CDR"]
    fr =  edges.df[which(edges.df[,"to"]==node) ,"synonyms_FR"] + edges.df[which(edges.df[,"to"]==node) ,"nonSynonyms_FR"]
    
    # save if the mutaion take place in the cdr or in the fr
    if(cdr > fr)
    {
      vec <- c(placeMut.vector,"cdr")
    }
    else if(cdr < fr)
      
    {
      vec <- c(placeMut.vector,"fr")
    }
    
    # save "vec" in the outer "sameIsotype.vector" vector
    assign("placeMut.vector",vec, envir = .GlobalEnv) 
    
  }
  
  # if we reached to an edge with suns
  else
  {
    
    # if it as a father
    f <- edges.df[which(edges.df[,"to"]==node) ,"from"]
    
    if(length(f) != 0)
    {
      vec <-c (sameIsotype.vector,"same")
      
      # save if the is an isotype switch or not
      if(getIsotypeName(node) != getIsotypeName(f))
      {
        vec <- c(sameIsotype.vector,"dif")
      }
      
      # save "vec" in the outer "sameIsotype.vector" vector
      assign("sameIsotype.vector",vec, envir = .GlobalEnv)
    }
    
    
    syn = edges.df[which(edges.df[,"to"]==node) ,"synonyms"]
    nonSyn = edges.df[which(edges.df[,"to"]==node) ,"nonSynonyms"]
       if(length(syn) != 0)
        {
    vec <-c ( majorMut.vector,"equal")
    
    # save if the is an isotype switch or not
    if(syn > nonSyn)
    {
      vec <- c(majorMut.vector,"syn")
    }
    else if(syn < nonSyn)
      
    {
      vec <- c(majorMut.vector,"nonSyn")
    }
    
    # save "vec" in the outer "sameIsotype.vector" vector
    assign("majorMut.vector",vec, envir = .GlobalEnv) 
    
     }
    
    # 
    cdr = edges.df[which(edges.df[,"to"]==node) ,"synonyms_CDR"] + edges.df[which(edges.df[,"to"]==node) ,"nonSynonyms_CDR"]
    fr =  edges.df[which(edges.df[,"to"]==node) ,"synonyms_FR"] + edges.df[which(edges.df[,"to"]==node) ,"nonSynonyms_FR"]
      if(length(cdr) != 0)
      {
    vec <-c(placeMut.vector,"equal")
    # save if the mutaion take place in the cdr or in the fr
    if(cdr > fr)
    {
      vec <- c(placeMut.vector,"cdr")
    }
    else if(cdr < fr)
      
    {
      vec <- c(placeMut.vector,"fr")
    }
    
    # save "vec" in the outer "sameIsotype.vector" vector
    assign("placeMut.vector",vec, envir = .GlobalEnv) 
      }
    
    
    sunsNewick.vector <- sapply(suns,GetEdgsDetails)
    
  }
}



getIsotypeName <-function (node)
  # this function get a node and return its isotype name
  # Args: 
  #   node- an node of the tree   
  # Returns: isotype nane of the node
{
  return (isotypesNames.df[which(isotypesNames.df[,"head"] == node),"IsotypeName"])
}


############## main ###########################
tree.dir <- paste0(source.dir,"/pfizerHumanTree15oct2015/")

# get list of fold
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)


# going over each patient
for(k in 1:length(patient.dirs))
{
  pdf(paper = "a4",file =paste0("/home/bardugn1/pfizer/PlotNew/",patient.dirs[k],"_treeImage.pdf"),onefile = T)
  
  # get list of folders
  dirs <- list.files(paste0(tree.dir,patient.dirs[k]),full.names = F)
  
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
        
        isotypesNames.df <- read.table(file= paste0(dirPath,"_isotypeNames.tab")
                                       ,header = T,stringsAsFactors = F,sep ="\t")
        
        # get a root for each clone
        roots <- setdiff(edges.df[,"from"],edges.df[,"to"])
        
        # goes over each root and send the root to "GoOverTree" function 
        # sapply(roots,GoOverTree)
        newickString <-vector()
        j<-1
        
        while(j<=length(roots))
        {
          # create a tree plot only trees that have at last 10 leave and 2 types of isotipes 
           if(!(is.null(totalLeavesType.fd)) && (totalLeavesType.fd[j,"type_num"]>3)
                                             && (totalLeavesType.fd[j,"sum_Leaves"]>10) )
           {
          
          # create a tree plot only trees that have at last 10 leave  
        #  if(!(is.null(totalLeavesType.fd)))
            #  && (totalLeavesType.fd[j,"sum_Leaves"]>10) )
        #  {      
            
            sameIsotype.vector <- c();
            majorMut.vector <- c();
            placeMut.vector <- c();
            
            # create a phylo pormat from the newick format we build                      
            tree.phylo <- read.tree(text =  paste0(GoOverTree(roots[j]),"0.0;"))
            
            GetEdgsDetails(roots[j])
            
            # get the leaves colors
            colors.tip =  sapply(tree.phylo[["tip.label"]],giveColor)
            
            # get the internal node colors
            colors.node = sapply(tree.phylo[["node.label"]],giveColor)
            
            # get edges info 
          #  shape.edge = sapply(sameIsotype.vector,givEdgeShape)
            shape.edge = sapply(placeMut.vector,givEdgeShape)
            col.edge = sapply(majorMut.vector,givEdgeColor)
            # plot the tree
            
            #png(filename= paste0('C:/Users/Nofar/Desktop/',dirs[k],'/',dirs[i],"_",j,"_treeImage.jpeg"),width = 1400, height = 800)
     

            
       #   plot.phylo(main = paste0(dirs[i],"_",j) ,tree.phylo,tip.color = colors.tip
      #                ,edge.width = 3,use.edge.length = F,edge.color = col.edge,edge.lty = shape.edge
      #                , font = 1,label.offset = 0.055)
        
           plot.phylo(main = paste0(dirs[i],"_",j) ,tree.phylo,tip.color = colors.tip,show.node.label = F,no.margin = F,
                   ,edge.width = 3,use.edge.length = T,edge.color = col.edge,edge.lty = shape.edge
                   , font = 1,label.offset = 0.055)
            
            # nodelabels(text = tree.phylo[["node.label"]],bg = colors.node,frame = "none",col =  colors.node)
            nodelabels(bg = colors.node,col =  colors.node)
      
      
          }
          
          j<-j + 1
        }
      }
    }
  }
  
  dev.off()
}  

