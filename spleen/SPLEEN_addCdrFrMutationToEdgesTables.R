
require(seqinr)
require(Biostrings)
source.dir <- "/home/bardugn1/"
#source.dir <- "/home/nofar/Desktop/spleen/"
setwd(source.dir)



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

# Change sequence names
#
# Params:   data frame of sequences
#
# Returns:  data frame of sequences with extra column of new headers
change.names <- function(fasta.df, outgroup){
  
  # Format - L*ID*_*CP ('Leave', ID number, copy number)
  fasta.df$head2 <- ""
  for(i in 1:nrow(fasta.df)){
    if(fasta.df[i,'head']==paste0(outgroup,'_1'))
      fasta.df[i,'head2'] <- fasta.df[i,'head']
    else{
      tmp <- unlist(strsplit(fasta.df[i,'head'],'_'))
      fasta.df[i,'head2'] <- paste0('L', i, '_', tmp[length(tmp)])
    }
  }
  
  return(fasta.df)
}

# Fix ambiguous nucleotides in internal sequences, according to IUPAC code
#
# Params:   nameSeq.df - data frame with headers and sequences, order by distance from root
#           edge.df - data frame with columns - parent(from), child(to), distance(weight)
#           outgroup name
#
# Returns:  - nameSeq.df - data frame with fixed sequences
fix.internal <- function(nameSeq.df, edge.df, outgroup){
  
  nucleotides <- c('A', 'C','G','T','-')
  
  # from leaves to root (except outgroup)
  # fix only internal nodes!
  intern.seq <- rev(setdiff(2:nrow(nameSeq.df),grep('^L', nameSeq.df[,'head'], perl=T, fixed=F)))
  for(i in intern.seq){
    # get ambiguous nucleotides
    curr.seq <- nameSeq.df[i,'seq']
    amb.pos <- setdiff(1:nchar(curr.seq),unlist(lapply(nucleotides, function(x) {gregexpr(pattern =x,curr.seq)})))
    if(length(amb.pos) > 0){
      # get son sequences
      sons.seq <- sapply(edge.df[which(edge.df[,'from']==nameSeq.df[i,'head']),'to'], function(x) {nameSeq.df[which(nameSeq.df[,'head']==x),'seq']})
      # do not fix outgroup sequence
      outgroup.ind <- which(names(sons.seq)==paste0(outgroup, '_1'))
      if(length(outgroup.ind)!=0) { sons.seq <- sons.seq[-outgroup.ind] }
      org.len <- nchar(curr.seq) ###
      for(j in amb.pos){
        amb.nuc <- substr(curr.seq,j,j)
        sons.nuc <- substr(sons.seq,j,j)
        # if deletion
        if(amb.nuc == 'O')
          curr.seq <- paste0(substr(curr.seq,1, j-1), '-', substr(curr.seq,j+1, nchar(curr.seq)))
        else{
          # check if at least one of the sons has one match for ambiguous letter
          if(is.element(amb.nuc, c('?', 'X')))
            matches <- which(is.element(sons.nuc, nucleotides))
          else
            matches <- which(is.element(sons.nuc, toupper(amb(amb.nuc, forceToLower = TRUE))))
          if(length(matches)==1) # if only one son had a match
            curr.seq <- paste0(substr(curr.seq,1, j-1), sons.nuc[matches], substr(curr.seq,j+1, nchar(curr.seq)))
          else if(length(matches)==0){ # choose randomly from son nucleotides
            rand <- sample(1:length(sons.nuc), 1, replace=T)
            curr.seq <- paste0(substr(curr.seq,1, j-1), sons.nuc[rand], substr(curr.seq,j+1, nchar(curr.seq)))
          }else{ # choose randomly from son nucleotides that matched
            rand <- sample(1:length(matches), 1, replace=T)
            curr.seq <- paste0(substr(curr.seq,1, j-1), sons.nuc[matches[rand]], substr(curr.seq,j+1, nchar(curr.seq)))
          }
        }
        if(org.len!=nchar(curr.seq)){
          print('HEREEEEEEEEEE')
        }
      }
      nameSeq.df[i,'seq'] <- curr.seq 
    }
  }
  return(nameSeq.df)
}

# Order sequence from root to leaves
#
# Params:   root sequence name
#           nameSeq.df - data frame with headers and sequences
#           edge.df - data frame with columns - parent(from), child(to), distance(weight)
#           outgroup name
#
# Returns:  nameSeq.df - data frame with headers and sequences, ordered by distance from root
order.nameSeq <- function(root, nameSeq.df, edge.df, outgroup){
  
  n.seq <- nrow(nameSeq.df)
  nameSeq.df2 <- nameSeq.df[which(nameSeq.df[,'head']==paste0(outgroup, '_1')),] # outgroup
  nameSeq.df2 <- rbind(nameSeq.df2, nameSeq.df[which(nameSeq.df[,'head']==root),]) # root
  for(i in 2:n.seq){ #1 - outgroup, #2 - root
    sons <- edge.df[which(edge.df[,'from']==nameSeq.df2[i,'head']),2]
    if(length(sons)>0){
      outgroup.ind <- which(sons==paste0(outgroup, '_1'))
      if(length(outgroup.ind)!=0) { sons <- sons[-outgroup.ind] }
      for(j in 1:length(sons))
        nameSeq.df2 <- rbind(nameSeq.df2, nameSeq.df[which(nameSeq.df[,'head']==sons[j]),])
    }
  }
  row.names(nameSeq.df2)<-NULL
  return(nameSeq.df2)
}

# Parse dnapars output tree file, get internal sequences and tree structure
#
# Params:   input file in phylip format and directory
#           outgroup name
#
# Returns:  - list of:
#           nameSeq.df - data frame with headers and sequences, ordered by distance from root
#           edge.df - data frame with columns - parent(from), child(to), distance(weight)
parse.dnapars <- function(file.name, in.dir, outgroup){
  
  # read output tree file
  out.tree <- scan(paste0(in.dir, file.name, '_out.txt'), what='character',sep='\n',
                   blank.lines.skip=F, strip.white=F)
  # check if tree was build
  if (any(grepl('-1 trees in all found', out.tree))) { return(NULL) }
  
  # get internal sequences
  seq.start <- min(grep('From\\s+To\\s+Any Steps\\?\\s+State at upper node', out.tree, perl=T, fixed=F))
  seq.empty <- grep('^\\s*$', out.tree[seq.start:length(out.tree)], perl=T, fixed=F)
  seq.len <- seq.empty[min(which(seq.empty[-1] == (seq.empty[-length(seq.empty)] + 1)))]
  seq.block <- paste(out.tree[(seq.start + 2):(seq.start + seq.len - 2)], collapse='\n')
  seq.df <- read.table(textConnection(seq.block), as.is=T, fill=T, blank.lines.skip=T,  row.names = NULL, header=F, stringsAsFactors=F)
  
  # fix root lines and remove empty rows
  fix.row <- which(seq.df[,3]!="yes" & seq.df[,3]!="no" & seq.df[,3]!="maybe")
  seq.df[fix.row, ] <- cbind(seq.df[fix.row, 1], seq.df[fix.row, 2],'no', seq.df[fix.row, 3:6], stringsAsFactors=F)
  
  # save full sequences as a data frame
  names <- unique(seq.df[, 2])
  seq <- sapply(names, function(x) { paste(t(as.matrix(seq.df[seq.df[, 2] == x, -c(1:3)])), collapse='') })
  nameSeq.df <- data.frame(head=names, seq=seq, stringsAsFactors=F, row.names=NULL)
  
  # get tree structure
  edge.start <- min(grep('between\\s+and\\s+length', out.tree, perl=T, fixed=F))
  edge.len <- min(grep('^\\s*$', out.tree[edge.start:length(out.tree)], perl=T, fixed=F))
  edge.block <- paste(out.tree[(edge.start + 2):(edge.start + edge.len - 2)], collapse='\n')
  edge.df <- read.table(textConnection(edge.block), col.names=c('from', 'to', 'weight'), as.is=T, stringsAsFactors=F)
  
  # order sequences by distance from root
  root <- seq.df[which(seq.df[,1]=='root')[1],2]
  nameSeq.df <- order.nameSeq(root, nameSeq.df, edge.df, outgroup)
  
  return(list(nameSeq.df, edge.df))
}


# Parse neighbor output tree file, get internal sequences and tree structure
#
# Params:   input file in phylip format and directory
#           data frame of input sequences
#           outgroup name
#,
# Returns:  - list of:
#           nameSeq.df - data frame with headers and sequences
#           edge.df - data frame with columns - parent(from), child(to), distance(weight)
parse.neighbor <- function(file.name, in.dir, fasta.df, outgroup){
  
  # read ouput tree file
  out.tree <- scan(paste0(in.dir, file.name, '_out.txt'), what='character',sep='\n',
                   blank.lines.skip=F, strip.white=F)
  
  # check if tree was build
  if (any(grepl('-1 trees in all found', out.tree))) { return(NULL) }
  
  # get tree structure
  edge.start <- min(grep('Between\\s+And\\s+Length', out.tree, perl=T, fixed=F))
  edge.len <- min(grep('^\\s*$', out.tree[edge.start:length(out.tree)], perl=T, fixed=F))
  edge.block <- paste(out.tree[(edge.start + 2):(edge.start + edge.len - 2)], collapse='\n')
  edge.df <- read.table(textConnection(edge.block), col.names=c('from', 'to', 'weight'), as.is=T)
  
  # create nameSeq.df from input sequence only (for now)
  nameSeq.df <- data.frame(head = fasta.df[,'head2'], seq = fasta.df[,'seq'], stringsAsFactors=F) 
  
  return(list(nameSeq.df, edge.df))
}

# Traverse tree in postorder (step 2 in Fitch's algorithm)
#
# Params:  father and its 2 son names
#          nameSeq.list - list of sequences with headers
#          edge.df - data frame with columns - parent(from), child(to), distance(weight)
#
# Returns: modified nameSeq.list - list of sequences with headers
traverse.up <- function(father, sons, edge.df, nameSeq.list){
  
  sons.seq <- sapply(sons,function(x) NULL)
  
  # check if each son sequence a leaf or is already reconstructed
  for(i in sons){
    if(!is.element(i, names(nameSeq.list)))
      nameSeq.list <-  traverse.up(i, edge.df[edge.df[,'from']==i,'to'], edge.df, nameSeq.list)
    sons.seq[i] <- nameSeq.list[i]
  }  
  seq.len <- length(sons.seq[[1]])
  
  father.seq <- character(seq.len)
  for(i in 1:seq.len){ # for each position in sequence
    curr.nuc <- sapply(sons.seq, function(x)  x[[i]])
    nuc <- Reduce(intersect, strsplit(curr.nuc,""))
    if(length(nuc)>0) # save intersection
      father.seq[i] <- paste(nuc, collapse='')
    else # save union
      father.seq[i] <- paste( Reduce(union, strsplit(curr.nuc,'')),collapse="")
  }
  nameSeq.list[[father]]<-father.seq # save modified sequence
  
  return(nameSeq.list)
}


# compute edge lengths
#
# Params:  nameSeq.df - data frame with headers and sequences
#          edge.df - data frame with columns - parent(from), child(to), edge weight (weight) 
#
# Returns: edge.df - data frame with columns - parent(from), child(to), edge weight (weight) , edge length(distance)
compute.edge <- function(nameSeq.df, edge.df,region){
  
  n.edge <- nrow(edge.df)
  
  # add new columns
  edge.df$distance <- rep(0, n.edge)
  edge.df$synonyms <-rep(0, n.edge)
  edge.df$nonSynonyms <-rep(0, n.edge)
  edge.df$synonymsCDR <-rep(0, n.edge)
  edge.df$nonSynonymsCDR <-rep(0, n.edge)
  edge.df$synonymsFR <-rep(0, n.edge)
  edge.df$nonSynonymsFR <-rep(0, n.edge)
  
  
  for(i in 1:n.edge)
  {
    from.seq <- unlist(strsplit(nameSeq.df[nameSeq.df[,'head']==as.character(edge.df[i,'from']),'seq'],''))
    to.seq <- unlist(strsplit(nameSeq.df[nameSeq.df[,'head']==as.character(edge.df[i,'to']),'seq'],''))
    edge.df[i,'distance'] <- length(which(from.seq!=to.seq))
    
    # find Synonyms and non Synonyms
    len.to <- length(to.seq)
    len.from <- length(from.seq)
    
    
    if((len.to + 4 < len.from)||(len.to  > len.from + 4))
    {
      print("len.to")
      print(len.to)      
      print("to.seq")
      print (to.seq)
      print("len.from")
      print(len.from)
      print("from.seq")
      print (from.seq)
      stop()
    }
    
    len <- min(len.to,len.from)
    j<-1
    k <- 3
    
    while(j <=len)
    {
      # get num of mutations in a threesome
      mutNum <-length(which(to.seq[j:k] != from.seq[j:k]))
      
      to.amino <- seqinr::translate(to.seq[j:k], numcode = 1, NAstring = "X", ambiguous = FALSE)
      
      from.amino <-seqinr::translate(from.seq[j:k], numcode = 1, NAstring = "X", ambiguous = FALSE)
      
      # if the amino acid identical : add "mutNum" to synonyms
      if(to.amino==from.amino)
      {
        edge.df[i,'synonyms']<- mutNum + edge.df[i,'synonyms']
        
      }
      # add "mutNum" to non synonyms
      else
      {
        edge.df[i,'nonSynonyms']<- mutNum + edge.df[i,'nonSynonyms']
      }
      
      # promote to next amino acid
      j <- j+3
      k <- k+3
      
    }
  }
  
  return(edge.df)
}




########### MAIN #########################

#  lengh for V and J (change according to dataset)
Vlen <- 291

MIN.SEQ <- 3
MID.SEQ <- 100
MAX.SEQ <- 7000

# get data frame of all 'v' options
Vgerm <- read.FASTA(file.name = "/Vgermline_human_VH.fasta",in.dir = source.dir)

# read the regions tables
fileAddress <- paste0(source.dir,"CDR_FR_regions_IgBLAST_human_VH.csv")
regions <- read.csv(file = fileAddress,header= T, stringsAsFactors= F)

V.RF <- get.RF(Vlen, Vgerm)
regions <- get.regions(regions, Vlen, V.RF)


tree.dir <- paste0(source.dir,"spleen/spleenHumanTree/")
fasta.dir <- paste0(source.dir,"spleen/spleenHumanFasta/")

# get list of fold
dirs = list.files(tree.dir,full.names = F)

# get list of folders - each folder represent patient
patient.dirs <- list.files(tree.dir,full.names = F)

# going over each patient
for(k in 1:length(patient.dirs))
{
  k<-1
  # get list of folders
  dirs <- list.files(paste0(tree.dir,patient.dirs[k]),full.names = F)
  
  #going over each tree folder
  for(i in 1:length(dirs))
  {
    
    dirPath <- paste0(tree.dir,patient.dirs[k],'/',dirs[i],'/')
    filePath <-paste0(dirPath,dirs[i],"_edges.tab")
    FASTA.file <- paste0(dirs[i],".fasta")
    
    # case file exist (if there is not enough or too much sequences, wile be only fasta file)
    if(file.exists(filePath))
    {

      # read each V-J-distance FASTA file
      fasta.df <- read.FASTA(FASTA.file, dirPath)
      file.name <-dirs[i]
      file.dir <- paste0(tree.dir,patient.dirs[k],"/",dirs[i],"/")
      n.seq <- nrow(fasta.df)
      
      if(n.seq >= MIN.SEQ & n.seq < MID.SEQ ){
       
        # parse dnapars output 
        tmp <- parse.dnapars(file.name , dirPath, "VJ")
        if(is.null(tmp)) { return(F) } # if tree building failed
        nameSeq.df <- tmp[[1]]
        edges.df <- tmp[[2]]
        nameSeq.df <- fix.internal(nameSeq.df,edges.df, "VJ")
        
        # if number of sequence is between MID.SEQ and MAX.SEQ - build tree with Neigbor joining
      }
   else if( n.seq >= MID.SEQ & n.seq < MAX.SEQ )
        {
          
         # change input sequence name to shorter ones
         fasta.df <- change.names(fasta.df, "VJ")
        # parse neighbor output
        tmp <- parse.neighbor( file.name , dirPath, fasta.df,"VJ")
        if(is.null(tmp)) { return(F) } # if tree building failed
        nameSeq.df <- tmp[[1]]
        edges.df <- tmp[[2]]
        nameSeq.df <- get.internal(nameSeq.df, edges.df, "VJ")
      }
      
      # add the CDR/FR mutation
      edge.df <- compute.edge(nameSeq.df,edges.df,regions)
      
      # save the new table into a file - if file exit - dystroyed old file
      write.table(edge.df, file = paste0(dirPath, "_edgesCut.tab"), quote=F, sep='\t', col.names=T, row.names=F, append = F)
      print(paste0(dirPath, "_edgesCut.tab"))
      
    }
  }
  
}

