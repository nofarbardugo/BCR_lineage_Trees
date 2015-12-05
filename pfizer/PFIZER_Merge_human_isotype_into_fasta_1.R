
# imports
require(seqinr)
require(Biostrings)

Enter_Isotype<-function(path,folder,fileAddress)
{
  #  this function reads a file and enter
  #           type of isotype into "sequence" column
  # arg: 
  #   folder name: folder to read the file from
  #   file address: name of the file
  #   path: path of the folder
  #
  # returns:
  #     data frame with "sequence", "sequence_id and "vdjdis" columns 
  ################################################
  
  # add the folder name into the path 
  fileAddress <- paste(path,fileAddress,sep="/")
  
  # read the file into data frame
  df <- read.csv(file = fileAddress ,header= T, stringsAsFactors= F)
  
  VJDis.vector <- c() 
  
  # enter the type of isotype into "SEQUENCE_ID" 
  for(i in 1:length(df[,"SEQUENCE_ID"]))
  {
    # get places of the char '_' in "SEQUENCE_ID" column
    lPlacesSequenceId <- gregexpr(pattern ='_',df[i,"SEQUENCE_ID"])
    
    # get Patient number
    lPlacesCloneId <- gregexpr(pattern ='_',df[i,"CLONE_ID"]) # get places of the char '_' in "CLONE_ID" column
    patientNumber<- substr(df[i,"CLONE_ID"],start =lPlacesCloneId[[1]][5] + 1,stop = lPlacesCloneId[[1]][6] -1)
    
    VJDis.vector[i] <- paste(substr(df[i,"SEQUENCE_ID"],start= lPlacesSequenceId[[1]][1] +1,
                                   stop=lPlacesSequenceId[[1]][2]-1),
                            substr(patientNumber,start = 1, stop = nchar(patientNumber)-1),
                            sep =  "_")
    
    
    df[i,"SEQUENCE_ID"] <- paste(substr(df[i,"SEQUENCE_ID"],start = 1, stop = lPlacesSequenceId[[1]][1] -1),
                                substr(df[i,"SEQUENCE_ID"],start =lPlacesSequenceId[[1]][1],stop = nchar(df[i,"SEQUENCE_ID"])),
                                sep =  paste0("_",folder,"_",patientNumber))
  }
  
  # return the edited data frame
  return (data.frame( df[,"SEQUENCE_ID"], VJDis.vector, df[,"SEQUENCE"], df[,"V_END"], stringsAsFactors=F))
}

########################main#############################################

# set directory - when run need to change her the direction
#setwd("G:/project/files/pfizer")

# get list of folders

path <- "/home/bardugn1/pfizer/Formatted/" # pep4
dirs <- list.files(path , full.names = F)
#dirs <- list.files( full.names = F) 

# going over each folder and get it's files
for(i in 1:length(dirs))
  {
   # get list of files 
   folderPath <- paste0(path,dirs[i]) 
   files <- list.files(folderPath, full.names = F)
   
   # going over each file  
   for(j in 1:length(files))
   {
      # get data frame of the file and add it to current dataframe
      if(i==1 & j==1)
      {
        df = Enter_Isotype(folderPath,dirs[i], files[j])
       
      }
      else
      {
        df= rbind(df,Enter_Isotype(folderPath,dirs[i], files[j]))
      }
    }
  }

# gives columns names 
names(df) <- c("SEQUENCE_ID","VJDis","SEQUENCE","V_END")

j <- 1
# create fasta files:
while(j <= length(df[,"VJDis"]))
{ 
  # give a path to the new fasta file - change according to wanted place 
 # fileName <- paste0("G:\\project\\files\\pfizer\\pfizerHumanFasta\\",df[j,"VJDis"],".fasta")
  #fileName <- paste0("/u/peptibase2/bardugn1/pfizerHumanFasta/",df[j,"VJDis"],".fasta")
  #fileName <- paste0("C:\\Users\\Nofar\\Desktop\\pfizerHumanFasta\\",df[j,"VJDis"],".fasta")
  #C:\Users\Nofar\Desktop\pfizerHumanFasta
  #fileName <- paste0(" C:\\Users\\Nofar\\Desktop\\FASTA\\",df[j,"VJDis"],".fasta")
 
  fileName <- paste0("/home/bardugn1/pfizer/pfizerHumanFasta15oct2015/",df[j,"VJDis"],"_",df[j,"V_END"],".fasta")#pep 4
  #fileName <- paste0("fasta\\",df[j,"VJDis"],".fasta")
  write.fasta(c(df[j,"SEQUENCE"]),c(df[j,"SEQUENCE_ID"]), file.out = fileName, open = "a", nbchar = 60)#, as.string = F) 
  j <- j+1
}
