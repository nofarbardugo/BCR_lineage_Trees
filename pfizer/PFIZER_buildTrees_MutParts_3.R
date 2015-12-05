
MIN.SEQ <- 3
MID.SEQ <- 100
MAX.SEQ <- 7000


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

# Create fasta file from clone data.frame
#
# Params:  clone.df = data.frame of clone with [taxa, seq] columns
#          out.file = the file name to write to
#
# Returns: NULL
write.FASTA <- function(file.name, out.dir, nameSeq.df) {
  sequences <- nameSeq.df$seq
  names(sequences) <- nameSeq.df$head
  writeXStringSet(DNAStringSet(sequences), file=paste0(out.dir, file.name, '.fasta'), width=1000)
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

# Create phylip output file
#
# Params:   data frame of sequences
#           output file name and directory
#
# Returns:  - 
fasta2phylip <- function(fasta.df, out.file, out.dir){
  phy.df <- rbind(data.frame(head=sprintf('%-9s', nrow(fasta.df)), seq=nchar(fasta.df$seq[1]), stringsAsFactors=F), 
                  data.frame(head=sprintf('%-9s', fasta.df$head2), seq=fasta.df$seq, stringsAsFactors=F))
  write.table(phy.df, file=paste0(out.dir, out.file, '.phy'), 
              quote=F, sep=' ', col.names=F, row.names=F)    
}

# Run dnapars (maximum parsimony)
#
# Params:   input file in phylip format and directory
#           phylip program path
#           outgroup index in file
#
# Returns:  - 
run.dnapars <- function( file.name, in.dir, phylip.path, outgroup.ind){
  
  curr.dir <- getwd()
  setwd(phylip.path) 
  
  # options from dnapars program
  # index of outgroup - last in data frame
  pars.options <- c(paste0(file.name, '.phy'), 'O', outgroup.ind, 'V', '1', '5', '.', '2', 'Y')
  
  # run dnapars
  system2('./dnapars', input=pars.options, stdout=NULL)
  
  # move .phy file and output tree files
  file.rename(from=paste0(phylip.path, 'outfile'), to=paste0(in.dir, file.name, '_out.txt'))
  file.rename(from=paste0(phylip.path, 'outtree'), to=paste0(in.dir, file.name, '_tree.txt'))
  
  setwd(curr.dir)
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

# Compute distance matrix (Kimura model) for neighbor joining using dnadist (Phylip package)
#
# Params:  input file name
#          phylip directory
#
# Returns: -
run.dnadist <- function(file.name,  phylip.path){
  
  # options from dnadist program
  dnadist.options <- c(paste0(file.name, '.phy'), 'D', '2', 'Y')
  
  # run dnadist
  system2('./dnadist', input=dnadist.options,, stdout=NULL)
  
  # move .phy file and output tree files
  file.rename(from=paste0(phylip.path, 'outfile'), to=paste0(phylip.path, file.name, '.dis'))
}


# Run neighbor (Neighbor joining)
#
# Params:  input file name and directory
#          phylip directory
#          outgroup index 
#
# Returns: -
run.neighbor <- function(file.name, in.dir, phylip.path, outgroup.ind){
  
  curr.dir <- getwd()
  setwd(phylip.path) 
  
  # run dnadist to create distance matrix
  run.dnadist(file.name,  phylip.path )
  
  # options for neighbor program
  # index of outgroup - last data frame
  neigh.options <- c(paste0(file.name, '.dis'), 'O', outgroup.ind, '2', 'Y')
  
  # run neighbor
  system2('./neighbor', input=neigh.options, stdout=NULL)
  
  # move .phy, .dis and output tree files
  file.rename(from=paste0(phylip.path, file.name, '.dis'), to=paste0(in.dir, file.name, '.dis'))
  file.rename(from=paste0(phylip.path, file.name, '.phy'), to=paste0(in.dir, file.name, '.phy'))
  file.rename(from=paste0(phylip.path, 'outfile'), to=paste0(in.dir, file.name, '_out.txt'))
  file.rename(from=paste0(phylip.path, 'outtree'), to=paste0(in.dir, file.name, '_tree.txt'))
  
  setwd(curr.dir) 
}

# Parse neighbor output tree file, get internal sequences and tree structure
#
# Params:   input file in phylip format and directory
#           data frame of input sequences
#           outgroup name
#
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

# Recursively traverse tree in preorder (step 1 in Fitch's algorithm)
#
# Params:  father and son names
#          nameSeq.list - list of sequences with headers
#          edge.df - data frame with columns - parent(from), child(to), distance(weight)
#
# Returns: modified nameSeq.list - list of sequences with headers
traverse.down <- function(father, son, edge.df, nameSeq.list, outgroup){
  
  son.seq <- nameSeq.list[[son]]
  seq.len <- length(son.seq)
  
  if(is.null(father)){ # root
    for(i in 1:seq.len){
      if(nchar(son.seq[i])>1){ # more than 1 option in intersection - choose randomly
        #rand <- floor(runif(1, min=1, max=nchar(son.seq[i])+1))
        rand <- sample(1:nchar(son.seq[i]), 1, replace=T)
        son.seq[i] <- substr(son.seq[i],rand, rand)
      }
    }
  }else{
    father.seq <- nameSeq.list[[father]]
    for(i in 1:seq.len){
      if(nchar(son.seq[i])>1){
        nuc <- Reduce(intersect, strsplit(c(father.seq[i],son.seq[i]),''))
        # if only one nucleotide in intersection - keep it
        if(length(nuc)==0){  # no intersection - choose randomly from son's options
          #rand <- floor(runif(1, min=1, max=nchar(son.seq[i])+1))
          rand <- sample(1:nchar(son.seq[i]), 1, replace=T)
          son.seq[i] <- substr(son.seq[i],rand, rand)
        }else # if intersection is not empty, but son has more than 1 option
          son.seq[i] <- nuc   
      }
    }
  }
  nameSeq.list[[son]] <- son.seq # save modified sequence
  
  # recursive call for each son's sons 
  sons <- edge.df[edge.df[,'from']==son,'to']
  if(length(sons) > 0){
    # remove outgroup from root's children
    outgroup.ind <- which(sons==paste0(outgroup, '_1'))
    if(length(outgroup.ind)!=0) { sons <- sons[-outgroup.ind] }
    nameSeq.list <- traverse.down(son, sons[1], edge.df, nameSeq.list, outgroup)
    nameSeq.list <- traverse.down(son, sons[2], edge.df, nameSeq.list, outgroup)
  }
  return(nameSeq.list)
}

# Reconstruct internal sequences with Fitch algorithm (for neighbor trees)
#
# Params:  nameSeq.df - data frame with headers and input sequences only
#          edge.df - data frame with columns - parent(from), child(to), distance(weight)
#          outgroup name
#
# Returns: nameSeq.df - data frame with headers and all sequences, ordered by distance from root
get.internal <- function(nameSeq.df, edge.df, outgroup){
  
  # find root - parent of outgroup
  root <- as.character(edge.df[edge.df[,'to']==paste0(outgroup,'_1'),'from'])
  root.sons <- edge.df[edge.df[,'from']==root,'to']
  
  # remove outgroup from root's children
  outgroup.ind <- which(root.sons==paste0(outgroup, '_1'))
  if(length(outgroup.ind)!=0) { root.sons <- root.sons[-outgroup.ind] }
  
  # convert sequences to lists 
  nameSeq.list <- sapply(nameSeq.df[,'seq'], function(x) {strsplit(x, "", fixed=FALSE)})
  names(nameSeq.list) <-  nameSeq.df[,'head']
  
  # Step 1 - preorder on tree - get intersection, otherwise - get union
  nameSeq.list <- traverse.up(root, root.sons, edge.df, nameSeq.list)
  
  # Step 2 - postorder on tree - get intersction, otherwise choose randomly
  nameSeq.list <- traverse.down(NULL, root, edge.df, nameSeq.list, outgroup)
  
  # convert sequence list to dataframe
  nameSeq.list2 <-lapply(nameSeq.list,function(x) {paste(x, collapse="")})
  nameSeq.df <-data.frame(head = names(nameSeq.list2), stringsAsFactors=F)
  nameSeq.df$seq <- unlist(nameSeq.list2)
  
  # order sequences in nameSeq.df by distance from root
  nameSeq.df <- order.nameSeq(root, nameSeq.df, edge.df, outgroup)
  return(nameSeq.df) 
}

# Match old input sequence names with new in final output sequence data frame
#
# Params:  nameSeq.df - data frame with headers and all sequences
#          fasta.df - data frame with new and old headers and input sequences only
#
# Returns: nameSeq.df - data frame with headers and sequences, with old and new names
match.names <- function(fasta.df, nameSeq.df){
  
  nameSeq.df$head2 <- rep("-", nrow(nameSeq.df))
  for(i in 1:nrow(fasta.df))
    nameSeq.df[which(fasta.df[i,'head2']==nameSeq.df[,'head']),'head2'] <- fasta.df[i,'head']
  
  return(nameSeq.df)
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
  edge.df$synonyms_CDR <-rep(0, n.edge)
  edge.df$nonSynonyms_CDR <-rep(0, n.edge)
  edge.df$synonyms_FR <-rep(0, n.edge)
  edge.df$nonSynonyms_FR <-rep(0, n.edge)
  edge.df$synonyms_CDR1 <-rep(0, n.edge)
  edge.df$synonyms_CDR2 <-rep(0, n.edge)
  edge.df$synonyms_CDR3 <-rep(0, n.edge)
  edge.df$synonyms_FR1 <-rep(0, n.edge)
  edge.df$synonyms_FR2 <-rep(0, n.edge)
  edge.df$synonyms_FR3 <-rep(0, n.edge)
  
  edge.df$nonSynonyms_CDR1 <-rep(0, n.edge)
  edge.df$nonSynonyms_CDR2 <-rep(0, n.edge)
  edge.df$nonSynonyms_CDR3 <-rep(0, n.edge)
  edge.df$nonSynonyms_FR1 <-rep(0, n.edge)
  edge.df$nonSynonyms_FR2 <-rep(0, n.edge)
  edge.df$nonSynonyms_FR3 <-rep(0, n.edge)
  
  edge.df$synonyms_CDR_NO3 <-rep(0, n.edge)
  edge.df$nonSynonyms_CDR_NO3 <-rep(0, n.edge)
  
  for(i in 1:n.edge)
  {
    #print(edge.df[i,'to'])
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
      if(mutNum>0) {
        to.amino <- seqinr::translate(to.seq[j:k], numcode = 1, NAstring = "X", ambiguous = FALSE)
        
        from.amino <-seqinr::translate(from.seq[j:k], numcode = 1, NAstring = "X", ambiguous = FALSE)
        
        # if the amino acid identical : add "mutNum" to synonyms
        if(to.amino==from.amino)
        {
          edge.df[i,'synonyms']<- mutNum + edge.df[i,'synonyms']
          
          # found mutaion place
          if(j<region[1,"CDR1"])
          {
            edge.df[i,'synonyms_FR1']<- mutNum + edge.df[i,'synonyms_FR1']
            edge.df[i,'synonyms_FR']<- mutNum + edge.df[i,'synonyms_FR']
          }
          else if (j<region[1,"FR2"]) 
          {
            edge.df[i,'synonyms_CDR1']<- mutNum + edge.df[i,'synonyms_CDR1']
            edge.df[i,'synonyms_CDR']<- mutNum + edge.df[i,'synonyms_CDR']
            edge.df[i,'synonyms_CDR_NO3']<- mutNum + edge.df[i,'synonyms_CDR_NO3']
          }
          else if (j<region[1,"CDR2"]) 
          {
            edge.df[i,'synonyms_FR2']<- mutNum + edge.df[i,'synonyms_FR2']
            edge.df[i,'synonyms_FR']<- mutNum + edge.df[i,'synonyms_FR']
          }
          else if (j<region[1,"FR3"]) 
          {
            edge.df[i,'synonyms_CDR2']<- mutNum + edge.df[i,'synonyms_CDR2']
            edge.df[i,'synonyms_CDR']<- mutNum + edge.df[i,'synonyms_CDR']
            edge.df[i,'synonyms_CDR_NO3']<- mutNum + edge.df[i,'synonyms_CDR_NO3']
          }
          else if (j<region[1,"CDR3"]) 
          {
            edge.df[i,'synonyms_FR3']<- mutNum + edge.df[i,'synonyms_FR3']
            edge.df[i,'synonyms_FR']<- mutNum + edge.df[i,'synonyms_FR']
          }
          else
          {
            edge.df[i,'synonyms_CDR3']<- mutNum + edge.df[i,'synonyms_CDR3']
            edge.df[i,'synonyms_CDR']<- mutNum + edge.df[i,'synonyms_CDR']
          }
          
        }
        # add "mutNum" to non synonyms
        else
        {
          edge.df[i,'nonSynonyms']<- mutNum + edge.df[i,'nonSynonyms']

          # found mutaion place
          if(j<region[1,"CDR1"])
          {
            edge.df[i,'nonSynonyms_FR1']<- mutNum + edge.df[i,'nonSynonyms_FR1']
            edge.df[i,'nonSynonyms_FR']<- mutNum + edge.df[i,'nonSynonyms_FR']
          }
          else if (j<region[1,"FR2"]) 
          {
            edge.df[i,'nonSynonyms_CDR1']<- mutNum + edge.df[i,'nonSynonyms_CDR1']
            edge.df[i,'nonSynonyms_CDR']<- mutNum + edge.df[i,'nonSynonyms_CDR']
            edge.df[i,'nonSynonyms_CDR_NO3']<- mutNum + edge.df[i,'nonSynonyms_CDR_NO3']
          }
          else if (j<region[1,"CDR2"]) 
          {
            edge.df[i,'nonSynonyms_FR2']<- mutNum + edge.df[i,'nonSynonyms_FR2']
            edge.df[i,'nonSynonyms_FR']<- mutNum + edge.df[i,'nonSynonyms_FR']
          }
          else if (j<region[1,"FR3"]) 
          {
            edge.df[i,'nonSynonyms_CDR2']<- mutNum + edge.df[i,'nonSynonyms_CDR2']
            edge.df[i,'nonSynonyms_CDR']<- mutNum + edge.df[i,'nonSynonyms_CDR']
            edge.df[i,'nonSynonyms_CDR_NO3']<- mutNum + edge.df[i,'nonSynonyms_CDR_NO3']
          }
          else if (j<region[1,"CDR3"]) 
          {
            edge.df[i,'nonSynonyms_FR3']<- mutNum + edge.df[i,'nonSynonyms_FR3']
            edge.df[i,'nonSynonyms_FR']<- mutNum + edge.df[i,'nonSynonyms_FR']
          }
          else
          {
            edge.df[i,'nonSynonyms_CDR3']<- mutNum + edge.df[i,'nonSynonyms_CDR3']
            edge.df[i,'nonSynonyms_CDR']<- mutNum + edge.df[i,'nonSynonyms_CDR']
          }  
        }
        
      }
      
      # promote to next amino acid
      j <- j+3
      k <- k+3
    
    }
  }
  
  return(edge.df)
}


# Main function - Build phylogenetic tree with maximum parsomony or neigbor joining 
#
# Params:  file.name - original FASTA file name
#          input directory for 
#
# Returns: -
build.trees <- function( FASTA.file, in.dir, out.dir, outgroup, phylip.path,regions){
  
  # read each V-J-distance FASTA file
  fasta.df <- read.FASTA(FASTA.file, in.dir)
  
  # skip small or big files
  n.seq <- nrow(fasta.df)
  if( n.seq > MAX.SEQ){
    # IN FUTURE - ADD OPTION FOR ALTERNATIVE PROGRAM
    print(paste0(FASTA.file, ' - Too many sequences to make tree with Phylip'))
    return(F)
  }
  
  # create output directory
  file.name <- substr(FASTA.file, 1, nchar(FASTA.file)-6)
  file.dir <- paste0(out.dir,'/' ,file.name, '/')
  dir.create(file.dir, showWarnings = FALSE)
  
  if( n.seq < MIN.SEQ) { 
    print(paste0(FASTA.file,' - Not enough sequences to make tree with Phylip'))
    # only 1 sequence - save 
    write.FASTA(file.name, file.dir, fasta.df)
    return(F) 
  }
  
  # change input sequence name to shorter ones
  fasta.df <- change.names(fasta.df, outgroup)
  
  # convert sequences to alugnment format for phylip programs
  fasta2phylip(fasta.df, file.name, phylip.path)
  
  # if number of sequence is between MIN.SEQ and MID.SEQ - build tree with Maximum Parsimony
  if( n.seq >= MIN.SEQ & n.seq < MID.SEQ ){
    print(paste0(FASTA.file, ' - maximum parsimony'))
    
    # run dnapars
    run.dnapars(file.name,in.dir = file.dir, phylip.path, n.seq)
    
    # parse dnapars output 
    tmp <- parse.dnapars(file.name, file.dir, outgroup)
    if(is.null(tmp)) { return(F) } # if tree building failed
    nameSeq.df <- tmp[[1]]
    edge.df <- tmp[[2]]
    nameSeq.df <- fix.internal(nameSeq.df,edge.df, outgroup)
  
  # if number of sequence is between MID.SEQ and MAX.SEQ - build tree with Neigbor joining
  }else if( n.seq >= MID.SEQ & n.seq < MAX.SEQ ){
    print(paste0(FASTA.file, ' - neighbor joining'))
    
    # run neighbor
    run.neighbor(file.name, file.dir, phylip.path, n.seq)
    
    # parse neighbor output
    tmp <- parse.neighbor(file.name, file.dir, fasta.df,outgroup)
    if(is.null(tmp)) { return(F) } # if tree building failed
    nameSeq.df <- tmp[[1]]
    edge.df <- tmp[[2]]
    nameSeq.df <- get.internal(nameSeq.df, edge.df, outgroup)
  }
  
  # compute edge lengths
  edge.df <- compute.edge(nameSeq.df, edge.df,regions)
  
  # retrieve old names for input sequences
  nameSeq.df <- match.names(fasta.df, nameSeq.df)
  
  # save output files
  
  # save sequence as FASTA file
  write.FASTA(file.name, file.dir, nameSeq.df)
  # save Fome/To/distances table (edges)
  write.table(edge.df, file=paste0(file.dir, file.name, '_edges.tab'), quote=F, sep='\t', col.names=T, row.names=F)
  # save old and new sequence names
  write.table(nameSeq.df[,c('head','head2')], file=paste0(file.dir, file.name, '_names.tab'), quote=F, sep='\t', col.names=T, row.names=F)
  
  return(T)
}


#############main##############################################
# RUN COMBINED TREES

# imports
require(seqinr)
require(Biostrings)

# flags
OUTGROUP <- 'VJ'
TREE.THRESHOLD <- 4

# set directory : change path depending on work environment
w.dir = "/home/bardugn1/" #pep 4
setwd(w.dir)

# directories for the "phylip" place :change path depending on work environment
PHYLIP.PATH <- 'Phylip/phylip-3.695/exe/' # for run on peptibase/ pep4
#PHYLIP.PATH <- '/phylip-3.69/exe/' # run on windows

out.dir <- "/home/bardugn1/pfizer/pfizerHumanTree15oct2015" #pep4

#  lengh for V and J
Vlen <- 130

# change path depending on work environment
#source(paste0(w.dir ,'/scripts/SPLEEN_buildTrees.R'))

# input directory of FASTA files :change path depending on work environment
in.dir <- "/home/bardugn1/pfizer/pfizerHumanFasta15oct2015/" #pep4


# get all fasta files names
VJdis.files <- list.files(path = in.dir, pattern='*.fasta', full.names = F)

phylip.path <- paste0(w.dir, PHYLIP.PATH)

print('Build phylogenetic trees')

# get data frame of all 'J' options
Jgermline <- read.FASTA(file.name = "/Jgermline_human_VH.fasta", in.dir= w.dir)

# get data frame of all 'v' options
Vgermline <- read.FASTA(file.name = "/Vgermline_human_VH.fasta", in.dir =w.dir)

# read the regions tables
fileAddress <- paste0(w.dir,"/CDR_FR_regions_IgBLAST_human_VH.csv")
regions <- read.csv(file = fileAddress,header= T, stringsAsFactors= F)

#V.RF <- get.RF(Vlen, Vgermline)
#regions <- get.regions(regions, Vlen, V.RF)

# going over all fasta files and build tree for each one
for(f in 1:length(VJdis.files))
{
  Vlen <- as.numeric(substr(x =VJdis.files[f],start = nchar(VJdis.files[f])-8,stop = nchar(VJdis.files[f])-6))
  V.RF <- get.RF(Vlen, Vgermline)
  regions <- get.regions(regions, Vlen, V.RF)
  
  # get patient number
  patientNum <- substr(VJdis.files[f],start = 13 ,stop = 15)
  tree.dir <- paste(out.dir,patientNum, sep = "/" )
  
  # create patient folder if not exist
  dir.create(tree.dir, recursive = F, showWarnings = F)
  
  V <- as.numeric(substr(x =VJdis.files[f],start = 1,stop = 3))

  # build the tree - calling function in the "buildTrees.r" script file 
  build.trees(VJdis.files[f], in.dir, tree.dir, OUTGROUP, phylip.path,regions[V,])
}