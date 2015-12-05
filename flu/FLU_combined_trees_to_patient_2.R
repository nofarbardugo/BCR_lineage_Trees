
add.outgroup <- function(dis, Vlen, Jlen, Vseq, Jseq ){
  #  this function create V-J outgroup sequence
  # arg: 
  #   dis: The distance between V and J
  #   Vlen: length for the V sequence
  #   Jlen: length for the J sequence
  #   Vseq: the V sequence 
  #   Jseq: the J sequence
  #
  # returns:
  #     V-J outgroup sequence after aditing 
  ################################################
  
  # V segment
  if(nchar(Vseq) < Vlen) {
    # if Vseq is too short - add gaps at beginning
    Vseq <- paste0(paste(rep('-',Vlen-nchar(Vseq)), collapse = ""), Vseq)
  }
  else {
   Vseq <- substr(Vseq, nchar(Vseq)-Vlen+1, nchar(Vseq))
  } 
  
  # J segment
  Jseq <- substr(Jseq, 1, Jlen)
  
  # cast dis to numeric
  dis <-as.numeric(dis)  
  
  # distance region
  if(dis < 0){  # if V-J distance is negative, removes symetrically NTs from V and J segments
    if(abs(dis) %% 2 == 0){ # distance is even
      VJseq <- paste0(substr(Vseq,1, nchar(Vseq)-ceiling(abs(dis/2))), substr(Jseq, ceiling(abs(dis/2))+1, nchar(Jseq)))
    } else { # distance is odd
      VJseq <- paste0(substr(Vseq,1, nchar(Vseq)-ceiling(abs(dis/2))), substr(Jseq, ceiling(abs(dis/2)), nchar(Jseq)))
    }
  }else
  {
    VJseq <- paste0(Vseq, paste(rep('-',dis), collapse = ""), Jseq)
  }
  
  return(VJseq)
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
PHYLIP.PATH <- '/Phylip/phylip-3.695/exe/' # for run on peptibase/ pep4


# home directory of current dataset :change path depending on work environment
out.dir <- "/home/bardugn1/flu2/fluHumanTree15oct2015" #pep4

#  length for V and J
#Vlen <- 291 # each fsta file hase its own vlen
Jlen <- 40

# change path depending on work environment
#source(paste0(w.dir ,'/scripts/FLU_buildTrees_MutParts.R'))

# input directory of FASTA files :change path depending on work environment
in.dir <- "/home/bardugn1/flu2/fluHumanFasta15oct2015/" #pep4

# get all fasta files names
VJdis.files <- list.files(path = in.dir, pattern='*.fasta', full.names = F)

phylip.path <- paste0(w.dir, PHYLIP.PATH)

# get data frame of all 'J' options
Jgermline <- read.FASTA(file.name = "/Jgermline_human_VH.fasta", in.dir= w.dir)

# get data frame of all 'v' options
Vgermline <- read.FASTA(file.name = "/Vgermline_human_VH.fasta", in.dir =w.dir)

# going over all fasta files and add VJ to fasta file # run only in first time/ 
for(j in 1:length(VJdis.files))
{
  #add VJ into the fasta file
  dis <- substr(x =VJdis.files[j],start = 9,stop = 11)
  vSeq <- Vgermline[as.numeric(substr(VJdis.files[j],start = 1,stop = 3)),"seq"]
  jSeq <- Jgermline[as.numeric(substr(VJdis.files[j],start = 5,stop = 7)),"seq"]
  Vlen <- substr(x =VJdis.files[j],start = nchar(VJdis.files[j])-8,stop = nchar(VJdis.files[j])-6)
  # get Outgroup sequence for the tree
  VJOutgroup <- add.outgroup(dis, as.double(Vlen), Jlen, vSeq, jSeq)
  
  # add into the fasta file the Outgroup 
  write.fasta(sequences = c(VJOutgroup),names = c("VJ_1"),  file.out = paste0(in.dir,VJdis.files[j]), open = "a", nbchar = 60)

}

