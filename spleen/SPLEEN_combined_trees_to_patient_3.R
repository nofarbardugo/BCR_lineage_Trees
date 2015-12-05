
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
#PHYLIP.PATH <- '/phylip-3.69/exe/' # run on windows

out.dir <- "/home/bardugn1/spleen/spleenHumanTree" #pep4

#  lengh for V and J
Vlen <- 291
Jlen <- 46

# change path depending on work environment
#source(paste0(w.dir ,'/scripts/SPLEEN_buildTrees.R'))

# input directory of FASTA files :change path depending on work environment
in.dir <- "/home/bardugn1/spleen/spleenHumanFasta/" #pep4

# get all fasta files names
VJdis.files <- list.files(path = in.dir, pattern='*.fasta', full.names = F)

phylip.path <- paste0(w.dir, PHYLIP.PATH)

print('Build phylogenetic trees')

# get data frame of all 'J' options
Jgermline <- read.FASTA(file.name = "/Jgermline_human_VH.fasta", in.dir= w.dir)

# get data frame of all 'v' options
Vgermline <- read.FASTA(file.name = "/Vgermline_human_VH.fasta", in.dir =w.dir)

# read the regions tables
#fileAddress <- paste0(source.dir,"CDR_FR_regions_IgBLAST_human_VH.csv")
#regions <- read.csv(file = fileAddress,header= T, stringsAsFactors= F)

#V.RF <- get.RF(Vlen, Vgermline)
#regions <- get.regions(regions, Vlen, V.RF)

# going over all fasta files and add VJ to fasta file # run only in first time/ can 
for(j in 1:length(VJdis.files))
{

  #add VJ into the fasta file
  dis <- substr(x =VJdis.files[j],start = 9,stop = 11)
 vSeq <- Vgermline[as.numeric(substr(VJdis.files[j],start = 1,stop = 3)),"seq"]
  jSeq <- Jgermline[as.numeric(substr(VJdis.files[j],start = 5,stop = 7)),"seq"]
  # get Outgroup sequence for the tree
  VJOutgroup <- add.outgroup(dis, Vlen, Jlen, vSeq, jSeq)
  
  # add into the fasta file the Outgroup 
 write.fasta(sequences = c(VJOutgroup),names = c("VJ_1"),  file.out = paste0(in.dir,VJdis.files[j]), open = "a", nbchar = 60)
}

# going over all fasta files and build tree for each one
#for(j in 1:length(VJdis.files))
#{
  # get patient number
#  patientNum <- substr(VJdis.files[j],start = 13 ,stop = 15)
#  tree.dir <- paste(out.dir,patientNum, sep = "/" )
  
  # create patient folder if not exist
#  dir.create(tree.dir, recursive = F, showWarnings = F)
  
  # build the tree - calling function in the "buildTrees.r" script file 
#  build.trees(VJdis.files[j], in.dir, tree.dir, OUTGROUP, phylip.path,regions)
#}

