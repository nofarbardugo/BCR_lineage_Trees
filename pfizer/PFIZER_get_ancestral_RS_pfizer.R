library(seqinr)

get.RF <- function(Vlen, Vgerm){
  
  # germline V sequence lengths
  germ.len <- nchar(Vgerm)
  
  # trim beginning of germline V sequences to match Vlen
  trim.germ <- lapply(Vgerm, function(x) substr(x, nchar(x)-Vlen+1, nchar(x)))
  
  # modulu of full germline V end
  rem <- germ.len%%3
  
  # modulu of trimmed germline V end (without extra NTs at the end of the full sequences)
  rem2 <- nchar(mapply(function(x,y) substr(x, 1, nchar(x)-y), x=trim.germ, y=rem))%%3
  
  # reading frame of trimmed germline V sequences
  V.RF <- rem2+1;
  
  # if there are germline sequences shorter than Vlen - set RF to be NA
  V.RF[germ.len<Vlen] <- NA
  
  return(V.RF)
}


get.regions <- function(regions, Vlen, V.RF){
  
  # for too short sequences - regions will be automatically NAs
  # calculate variable region length according to Vlen
  diff <- regions$V_END-Vlen+V.RF-1
  
  # fix region indices
  regions[,-1] <- regions[,-1]-diff
  
  # for Vlen which does not include all regions, set first indices to -1
  # set first index to 1
  for(i in 1:nrow(regions)){
    neg.ind <- which(regions[i,-1]<0)
    if(length(neg.ind)==1)
      regions[i,-1][neg.ind] <- 1
    else{
      regions[i,-1][neg.ind[1:length(neg.ind)-1]] <- -1
      regions[i,-1][neg.ind[length(neg.ind)]] <- 1
    }    
  }
  
  return(regions)
}




# Read FASTA file and conert to data frame
#
# Params:   FASTA file name
#           file directory
#
# Returns:  data frame of sequences
read.FASTA <- function(file.name, in.dir){
  
  fasta.list <- read.fasta(paste0(in.dir, file.name), seqtype = "DNA", as.string = T, forceDNAtolower = F,
                           set.attributes = F, legacy.mode = T, seqonly = F, strip.desc = F)
  fasta.df <- data.frame(head = names(fasta.list), stringsAsFactors=F)
  fasta.df$seq <- unlist(fasta.list)
  return(fasta.df)
}

Vlen <- 130

w.dir = "/home/bardugn1/" #pep 4
setwd(w.dir)

# load germline
germline.path <- '/home/bardugn1/'
Vgerm <- read.FASTA(file.name = "/Vgermline_human_VH.fasta", in.dir =w.dir)
Jgerm <- read.FASTA(file.name = "/Jgermline_human_VH.fasta", in.dir =w.dir)
# load region file
regions <- read.csv(file = "CDR_FR_regions_IgBLAST_human_VH.csv",header= T, stringsAsFactors= F)

#Vgerm <- read.fasta("Vgermline_human_VH.fasta", 
#                    seqtype = "DNA", as.string = T, forceDNAtolower = F, set.attributes = F, legacy.mode = T, seqonly = F, strip.desc = F) 

df.dir <- 'pfizer/Formatted/'
tree.dir <- 'pfizer/pfizerHumanTree15oct2015/'
out.dir <- 'pfizer/mut/'

dem <-c()
index <-1
res <- data.frame(sample = c(),isotype = c(),tree = c(), rootNum = c(),distance = c(),
		              synonyms =c(),nonSynonyms =c(),
                  synonyms_CDR= c(),synonyms_FR = c(),
                  nonSynonyms_CDR = c(),nonSynonyms_FR = c(),
                  synonyms_CDR1 = c(),synonyms_CDR2 = c(), synonyms_CDR3 = c(),
           	      synonyms_FR1 = c(), synonyms_FR2 = c(), synonyms_FR3 = c(),
           	      nonSynonyms_CDR1 = c(), nonSynonyms_CDR2 = c(), nonSynonyms_CDR3 = c(),
           	      nonSynonyms_FR1 = c(), nonSynonyms_FR2 = c(), nonSynonyms_FR3 = c(),
  	   	          synonyms_CDR_NO3 = c(), nonSynonyms_CDR_NO3 = c(), stringsAsFactors=F
                  )

numberOfcloneConut <-0
numberNotIF <-0
numberJ_NotMod3 <-0
numberDistMinus <- 0
numberAll <-0
  
# get list of samples - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)
patient.dirs <- patient.dirs[-grep('136',patient.dirs)]
patient.dirs <- patient.dirs[-grep('273',patient.dirs)]
# going over each patient
for(k in 1:length(patient.dirs))
{
  # get all tree folders in sample 'k'
  dirs <- list.files(paste0(tree.dir,patient.dirs[k]),full.names = F)

  # going over each tree folder
  for(i in 1:length(dirs))
  {
    dirPath <- paste0(tree.dir,patient.dirs[k],'/',dirs[i],'/',dirs[i])
    filePath <-paste0(dirPath,"_edgesCut.tab")
   
    # case file exist (if there is not enough or too much sequences, wile be only fasta file)
    if(file.exists(filePath))
    {
      
      # get all the edges of the current tree
      trim.edge.df <- tryCatch(read.table(file = filePath ,header = T,stringsAsFactors = F,sep ="\t"), error=function(e) NULL)
      
      if(!is.null(trim.edge.df))
      {
              
        # get details about current tree
        totalLeavesType.fd <-tryCatch(read.table(file = paste0(dirPath,"_totalLeavesType.tab")
                                                 ,header = T,stringsAsFactors = F,sep ="\t"),
                                      error=function(e) NULL)
        
        nameSeq.df <- read.table(file= paste0(dirPath,"_isotypeNames.tab"),header = T,stringsAsFactors = F,sep ="\t")
        fasta.df <- read.FASTA(file.name = ".fasta", in.dir =dirPath)
        
        # get a root for each clone
        sub.roots <- setdiff(trim.edge.df[,"from"],trim.edge.df[,"to"])
        
        # case there is no roots
        if (length(sub.roots)==0)
        {
          next
        }
        
      	# check if clone is in-frame
        J <- as.numeric(substr(dirs[i],start = 5,stop = 7))
        if(J<0)
        {
          dem[index] <- dirs[i]

          numberOfcloneConut <-0
          numberNotIF <-numberNotIF +1

          next
        }
      	Jseq.len <- nchar(Jgerm[J,"seq"])-1
      	seq.len <- nchar(fasta.df[4,"seq"])-20 + Jseq.len
	

	      if(seq.len%%3 !=0)
	      {
	        numberJ_NotMod3 <-numberJ_NotMod3 +1

    		#  next
    	  }
       
         # check that the distanse between V to J is higher then -30
        if (as.numeric(substr(dirs[i],start = 9,stop = 11)) < -25)
        {
          numberDistMinus <- numberDistMinus +1
          next
        }
        
       
      	# get germline sequence
      	V <- as.numeric(substr(dirs[i],start = 1,stop = 3))
      	Vend <- as.numeric(substr(dirs[i],start = nchar(dirs[i])-2,stop = nchar(dirs[i])))
      	Vgene <- Vgerm[V,"head"]
      	germ.seq  <- Vgerm[V,"seq"]
      	germ.len <- nchar(germ.seq)
      	germ.seq <- substr(germ.seq, germ.len-Vend+1, germ.len-10) # minus 10 becuase we dont want the cdr3  
      	germ.char <- unlist(strsplit(germ.seq,''))
      
     	 # split tree into clones - get all nodes in current clone
    	 # for each clone root - compute R/S ratio in CDR and FR between root and germline V gene
      	for(j in 1:length(sub.roots))
      	{
      	 
      	  numberAll <- numberAll+1
      	  # count modes only on trees that have at last 10 leave and 2 types of isotipes 
      	  if(!(is.null(totalLeavesType.fd)) && (totalLeavesType.fd[j,"type_num"]>=1)
           && (totalLeavesType.fd[j,"sum_Leaves"]>1) )
      	  {
      	    numberOfcloneConut <-numberOfcloneConut +1
      	    root <- sub.roots[j]
      	    
      	    # get root sequence
      	    seq <- fasta.df[which(fasta.df$head%in%root), "seq"]
      	    #seq <- substr(seq, 1+Vlen-Vend, Vlen-10)
      	    seq <- substr(seq, 1, Vend-10)
      	    seq.char <- unlist(strsplit(seq,''))
      	    
      	    # get mutation posiiton
      	    
      	    
      	    V.RF <- get.RF(Vend, Vgerm)
      	    region <- get.regions(regions, Vend, V.RF)
      	    seq.region <- region[V,]
      	    
      	    germ.len <- length(germ.char)
      	    seq.len <- length(seq.char)
      	    
      	    if(germ.len != seq.len)
      	    {
      	      print("error - germ len soppuse to be equal to seq len")
      	      print(dirs[i])
      	      stop()
      	    }
      	    
      	    distance <- length(which(germ.char!=seq.char))
      	    synonyms <- 0
      	    nonSynonyms <- 0
      	    synonyms_CDR <- 0
      	    nonSynonyms_CDR <- 0
      	    synonyms_FR <- 0
      	    nonSynonyms_FR <- 0
      	    synonyms_CDR1 <- 0
      	    synonyms_CDR2 <- 0
      	    synonyms_CDR3 <- 0
      	    synonyms_FR1 <- 0
      	    synonyms_FR2 <- 0
      	    synonyms_FR3 <- 0
      	    
      	    nonSynonyms_CDR1 <- 0
      	    nonSynonyms_CDR2 <- 0
      	    nonSynonyms_CDR3 <- 0
      	    nonSynonyms_FR1 <- 0
      	    nonSynonyms_FR2 <- 0
      	    nonSynonyms_FR3 <- 0
      	    
      	    synonyms_CDR_NO3 <- 0
      	    nonSynonyms_CDR_NO3 <- 0
      	    
      	    from<-1
      	    to <- 3
      	    
      	    while(from <=germ.len)
      	    {
      	      # get num of mutations in a threesome
      	      
      	      mutNum <-length(which(seq.char[from:to] != germ.char[from:to]))
      	      if(mutNum>0)
      	      {
      	        to.amino <- seqinr::translate(seq.char[from:to], numcode = 1, NAstring = "X", ambiguous = FALSE)
      	        
      	        from.amino <-seqinr::translate(germ.char[from:to], numcode = 1, NAstring = "X", ambiguous = FALSE)
      	        
      	        # if the amino acid identical : add "mutNum" to synonyms
      	        if(to.amino==from.amino)
      	        {
      	          synonyms <- mutNum + synonyms
      	          
      	          # found mutaion place
      	          if(from<region[1,"CDR1"])
      	          {
      	            synonyms_FR1 <- mutNum + synonyms_FR1
      	            synonyms_FR <- mutNum + synonyms_FR
      	          }
      	          else if (from<region[1,"FR2"]) 
      	          {
      	            synonyms_CDR1 <- mutNum + synonyms_CDR1
      	            synonyms_CDR <- mutNum + synonyms_CDR
      	            synonyms_CDR_NO3 <- mutNum + synonyms_CDR_NO3
      	          }
      	          else if (from<region[1,"CDR2"]) 
      	          {
      	            synonyms_FR2 <- mutNum + synonyms_FR2
      	            synonyms_FR <- mutNum + synonyms_FR
      	          }
      	          else if (from<region[1,"FR3"]) 
      	          {
      	            synonyms_CDR2 <- mutNum + synonyms_CDR2
      	            synonyms_CDR <- mutNum + synonyms_CDR
      	            synonyms_CDR_NO3 <- mutNum + synonyms_CDR_NO3
      	          }
      	          else if (from<region[1,"CDR3"]) 
      	          {
      	            synonyms_FR3 <- mutNum + synonyms_FR3
      	            synonyms_FR <- mutNum + synonyms_FR
      	          }
      	          else
      	          {
      	            synonyms_CDR3 <- mutNum + synonyms_CDR3
      	            synonyms_CDR <- mutNum + synonyms_CDR
      	          }
      	          
      	        }#if(to.amino==from.amino)
      	        # add "mutNum" to non synonyms
      	        else
      	        {
      	          nonSynonyms <- mutNum + nonSynonyms
      	          
      	          # found mutaion place
      	          if(from<region[1,"CDR1"])
      	          {
      	            nonSynonyms_FR1 <- mutNum + nonSynonyms_FR1
      	            nonSynonyms_FR <- mutNum + nonSynonyms_FR
      	          }
      	          else if (from<region[1,"FR2"]) 
      	          {
      	            nonSynonyms_CDR1<- mutNum + nonSynonyms_CDR1
      	            nonSynonyms_CDR<- mutNum + nonSynonyms_CDR
      	            nonSynonyms_CDR_NO3<- mutNum + nonSynonyms_CDR_NO3
      	          }
      	          else if (from<region[1,"CDR2"]) 
      	          {
      	            nonSynonyms_FR2 <- mutNum + nonSynonyms_FR2
      	            nonSynonyms_FR <- mutNum + nonSynonyms_FR
      	          }
      	          else if (from<region[1,"FR3"]) 
      	          {
      	            nonSynonyms_CDR2 <- mutNum + nonSynonyms_CDR2
      	            nonSynonyms_CDR <- mutNum + nonSynonyms_CDR
      	            nonSynonyms_CDR_NO3 <- mutNum + nonSynonyms_CDR_NO3
      	          }
      	          else if (from<region[1,"CDR3"]) 
      	          {
      	            nonSynonyms_FR3 <- mutNum + nonSynonyms_FR3
      	            nonSynonyms_FR <- mutNum + nonSynonyms_FR
      	          }
      	          else
      	          {
      	            nonSynonyms_CDR3 <- mutNum + nonSynonyms_CDR3
      	            nonSynonyms_CDR <- mutNum + nonSynonyms_CDR
      	          }  
      	        }	
      	        
      	      }#if(mutNum>0)
      	      
      	      # promote to next amino acid
      	      from <- from+3
      	      to <- to+3
      	      
      	    }#while(from <=germ.len)
      	    
      	    # enter new record into res table
      	    
      	    
      	    newRow <- data.frame(sample = c(patient.dirs[k]),isotype =nameSeq.df[nameSeq.df$head==root,"IsotypeName"],tree = c(dirs[i]), rootNum = c(sub.roots[j]),distance = c(distance),
      	                         synonyms =c(synonyms),nonSynonyms =c(nonSynonyms),
      	                         synonyms_CDR= c(synonyms_CDR),synonyms_FR = c(synonyms_FR),
      	                         nonSynonyms_CDR = c(nonSynonyms_CDR),nonSynonyms_FR = c(nonSynonyms_FR),
      	                         synonyms_CDR1 = c(synonyms_CDR1),synonyms_CDR2 = c(synonyms_CDR2), synonyms_CDR3 = c(synonyms_CDR3),
      	                         synonyms_FR1 = c(synonyms_FR1), synonyms_FR2 = c(synonyms_FR2), synonyms_FR3 = c(synonyms_FR3),
      	                         nonSynonyms_CDR1 = c(nonSynonyms_CDR1), nonSynonyms_CDR2 = c(nonSynonyms_CDR2), nonSynonyms_CDR3 = c(nonSynonyms_CDR3),
      	                         nonSynonyms_FR1 = c(nonSynonyms_FR1), nonSynonyms_FR2 = c(nonSynonyms_FR2), nonSynonyms_FR3 = c(nonSynonyms_FR3),
      	                         synonyms_CDR_NO3 = c(synonyms_CDR_NO3), nonSynonyms_CDR_NO3 = c(nonSynonyms_CDR_NO3), stringsAsFactors=F
      	    )
            
      	    res <- rbind(res,newRow)
      	    
      	    
      	  }

      	}#for(j in 1:length(sub.roots))
       
      	print(patient.dirs[k])
      }#if(!is.null(trim.edge.df))
      	  
    }#if(file.exists(filePath))
  }#for(i in 1:length(dirs))
}#for(k in 1:length(patient.dirs))

print(numberOfcloneConut) 
print(numberNotIF)
print(numberJ_NotMod3)
print(numberDistMinus)
print(numberAll)
 
#save table
#write.csv(res, paste0(out.dir, 'PFIZER_Combine_isotype_RS.csv'), row.names=F)  
#write.csv(res, paste0(out.dir, 'PFIZER_all_above_10_RS.csv'), row.names=F)   
write.csv(res, paste0(out.dir, 'PFIZER_RS.csv'), row.names=F)  
print("finish")
