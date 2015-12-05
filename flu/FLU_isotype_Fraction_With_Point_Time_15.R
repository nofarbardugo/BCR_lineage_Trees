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

getIsotypeTimePoint <-function (node)
  # this function get a node and return its Time point
  # Args: 
  #   node- an node of the tree   
  # Returns: Time Point of the node
  # Note : only leaves have TimePoint (internal nodes are not real samples)
{
  return (isotypesNames.df[which(isotypesNames.df[,"head"] == node),"TimePoint"])
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

GoOverTree<-function(node)
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
  
  # temp datafream that will sava all the new line to add 
  d <-data.frame(IsotypeName = c(),ChildName = c(),Synonymus =c(),nonSynonymus =c(),
                                   depth = c(),Ratio = c(), synonyms_CDR= c(),synonyms_FR = c(),
                                   nonSynonyms_CDR = c(),nonSynonyms_FR = c(),TimePoint = c(),stringsAsFactors=F)
  # calc mesurmeant for the suns
  p <-1
  while(p <=sunsNum)
  {
    if(suns[p]!= "VJ_1")
    {  
        # check if its a leaf -  # add new line to data fream
        if(isIsotypeLeaf(suns[p])==T)
        {
           # get mutations information from cure node to each of its child         
           syn <- edges.df [(edges.df [,"to"]==suns[p]),"synonyms"]
           nonSyn <- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms"]
           syn_cdr <- edges.df [(edges.df [,"to"]==suns[p]),"synonyms_CDR"]
           syn_fr <-  edges.df [(edges.df [,"to"]==suns[p]),"synonyms_FR"]
           nonSyn_cdr <- edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_CDR"]
           nonSyn_fr <-  edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_FR"]
	         depth <- syn+nonSyn
           ratio <-0
           if(depth >0)
           {
              ratio <-nonSyn/depth
           }
           
           tmp <-data.frame(IsotypeName = c(type),ChildName = c(getIsotypeName(suns[p])),Synonymus =c(syn),nonSynonymus =c(nonSyn),
                          depth = c(depth),Ratio = c(ratio), synonyms_CDR= c(syn_cdr),synonyms_FR = c(syn_fr),
                          nonSynonyms_CDR = c(nonSyn_cdr),nonSynonyms_FR = c(nonSyn_fr),TimePoint = c(getIsotypeTimePoint(suns[p])),
                          stringsAsFactors=F)
           
          d <- rbind(d,tmp)
        }
        # not a leaf - calc son node fraction details and add it to the father data.frame (with require adjustment)
      	# the childern of the son is also childern of the father - 
        else
        {

          sonDetails.df <- GoOverTree(suns[p])
          
          # adjustment the details of the son:
       	  sonDetails.df$IsotypeName <-type # change the son name into the father 
	        sonDetails.df$Synonymus <- sonDetails.df$Synonymus + edges.df [(edges.df [,"to"]==suns[p]),"synonyms"]
	        sonDetails.df$nonSynonymus <- sonDetails.df$nonSynonymus + edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms"]
          sonDetails.df$synonyms_CDR <- sonDetails.df$synonyms_CDR + edges.df [(edges.df [,"to"]==suns[p]),"synonyms_CDR"]
          sonDetails.df$synonyms_FR <- sonDetails.df$synonyms_FR + edges.df [(edges.df [,"to"]==suns[p]),"synonyms_FR"]
 	        sonDetails.df$nonSynonyms_CDR <- sonDetails.df$nonSynonyms_CDR + edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_CDR"]
	        sonDetails.df$nonSynonyms_FR <- sonDetails.df$nonSynonyms_FR + edges.df [(edges.df [,"to"]==suns[p]),"nonSynonyms_FR"]
	        sonDetails.df$depth <- sonDetails.df$Synonymus + sonDetails.df$nonSynonymus

	        for(mm in 1:nrow(sonDetails.df))
      	  {
		         if(sonDetails.df[mm,"depth"]!= 0)
		         {
			          sonDetails.df[mm,"Ratio"] <- sonDetails.df[mm,"nonSynonymus"]/ sonDetails.df[mm,"depth"]
		         }
	 
          }

	        # add all new line to the father lines
          d <- rbind(d,sonDetails.df)

        }        
      
    }#if(suns[p]!= "VJ_1")
   
    p <- p+1 
    
  }# end:while(p <=sunsNum)
  
  # - add new lines to "df.allIsotypeDetails" datafream 
  data <- df.allIsotypeDetails
  data <- rbind(data,d) 
  
  # save "data" in the outer "df.allIsotypeDetails" data frame
  assign("df.allIsotypeDetails",data, envir = .GlobalEnv)
  
  return(d) # return the new lines into the father

} 


############## main #####################

require(gdata)
# path for trees folder
tree.dir <- paste0(source.dir,"/fluHumanTree15oct2015/")  # change here

# get list of folders
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)
patient.dirs <- patient.dirs [-grep('IB',patient.dirs)]

# create datafream that will save for each isotype :IsotypeName,kidsName, Synonymus mutataion,
#                                                   NON Synonymus mutataion,Depth = s+ ns, Ratio = ns/(ns + s)   
#                                                   synonyms_CDR mutation, synonyms_FR mutation,
#                                                   nonSynonyms_CDR mutations and nonSynonyms_FR mutations, time point
df.allIsotypeDetails <- data.frame(IsotypeName = c(),ChildName = c(),Synonymus =c(),nonSynonymus =c(),
                                   depth = c(),Ratio = c(), synonyms_CDR= c(),synonyms_FR = c(),
                                   nonSynonyms_CDR = c(),nonSynonyms_FR = c(),TimePoint = c())

# going over each patient
for(k in 1:length(patient.dirs))
{
  # get list of folders
  dirs <- list.files(paste0(tree.dir,patient.dirs[k]),full.names = F)
  
  # going over each tree folder
  for(i in 1:length(dirs))
  {
    dirPath <- paste0(tree.dir,patient.dirs[k],'/',dirs[i],'/',dirs[i])
    filePath <-paste0(dirPath,"_edgesCut_postVacine.tab") # change here
    
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

names(df.allIsotypeDetails)<- c ("Isotype Name","Child Name ","Synonymus","nonSynonymus","depth","Ratio",
				 "synonyms_CDR","synonyms_FR","nonSynonyms_CDR","nonSynonyms_FR",
				 "Time Point")



write.csv(df.allIsotypeDetails,file = "FLU_isotype_Fraction_With_Point_Time.csv")

