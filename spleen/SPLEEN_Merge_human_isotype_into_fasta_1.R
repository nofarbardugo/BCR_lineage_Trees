
# imports
require(seqinr)
#require(Biostrings)

#Enter_Isotype<-function(path,folder,fileAddress)
Enter_Isotype<-function(path,fileAddress)
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
  
  # add the file name into the path 
  fileAddress <- paste0(path,fileAddress)

  # read the file into data frame
  df <- read.csv(file = fileAddress ,header= T, stringsAsFactors= F)
  
  VJDis.vector <- c() 
  
  # enter the type of isotype into "SEQUENCE_ID" 
  for(i in 1:length(df[,"SEQUENCE_ID"]))
  {
    # get places of the char '_' in "SEQUENCE_ID" column
    lPlacesSequenceId <- gregexpr(pattern ='_',df[i,"SEQUENCE_ID"])
    
    # get Patient number + isotype IGM plus or minus + place
    lPlacesCloneId <- gregexpr(pattern ='_',df[i,"CLONE_ID"]) # get places of the char '_' in "CLONE_ID" column
    start <- length(lPlacesCloneId[[1]]) -1
    
    sampleDetails<- substr(df[i,"CLONE_ID"],start =1,stop = lPlacesCloneId[[1]][start] -1)
  
    vdgFromClone = substr(df[i,"CLONE_ID"],start= lPlacesCloneId[[1]][start] +1,stop=lPlacesCloneId[[1]][start +1]-1)
      
    # get Patient VDJ + the patient number "000.000.000_XXX"
    VJDis.vector[i] <- paste(vdgFromClone,substr(df[i,"CLONE_ID"],start =0,stop = 3),
                          sep =  "_")
    
    
    df[i,"SEQUENCE_ID"] <- paste(substr(df[i,"SEQUENCE_ID"],start = 1, stop = lPlacesSequenceId[[1]][1] ),
                                substr(df[i,"SEQUENCE_ID"],start =lPlacesSequenceId[[1]][1] ,stop = nchar(df[i,"SEQUENCE_ID"])),
                                sep =  sampleDetails)
    
    if( VJDis.vector[i] == "S35_SPL")
    {
      #print("bkmflb")
      print(i) 
    }

  }
  
  # return the edited data frame
  return (data.frame( df[,"SEQUENCE_ID"], VJDis.vector, df[,"SEQUENCE"], stringsAsFactors=F))
}

########################main#############################################

# set directory - when run need to change her the direction
#setwd("G:/project/files/pfizer")

# get list of folders
#path <- "/media/nofar/CDD8-A371/project/files/pfizer/humenIsotype/" #ubuntu
#path <-"/home/nofar/Desktop/try/" #window
#path <-"G:/project/files/pfizer/humenIsotype/" #window
#path <-"/u/peptiqnap2/labdata/Genetics/Repertoires/Bcells/Heavy/Human/pfizer/" # peptibase
path <- "/home/bardugn1/spleen/Formatted/"  # pep 4
#dirs <- list.files(path , full.names = F)
#dirs <- list.files( full.names = F) 

# going over each folder and get it's files

print ("start")
#for(i in 1:length(dirs))
#  {
   # get list of files 
 #  folderPath <- paste0(path,dirs[i]) 
   files <- list.files(path, full.names = F)
   
   # going over each file  
   for(j in 1:length(files))
   {
      # get data frame of the file and add it to current dataframe
      if(j==1)
      {
       # df = Enter_Isotype(path,dirs[i], files[j])
        df = Enter_Isotype(path, files[j])
      }
      else
      {
       # df= rbind(df,Enter_Isotype(path,dirs[i], files[j]))
        df= rbind(df,Enter_Isotype(path, files[j]))
      }
    }
 # }

# gives columns names 
names(df) <- c("SEQUENCE_ID","VJDis","SEQUENCE")

j <- 1

# create fasta files:
while(j <= length(df[,"VJDis"]))
{ 
  # give a path to the new fasta file - change according to wanted place 

  fileName <- paste0("/home/bardugn1/spleen/spleenHumanFasta/",df[j,"VJDis"],".fasta") #pep 4
  
  print( fileName)
  write.fasta(c(df[j,"SEQUENCE"]),c(df[j,"SEQUENCE_ID"]), file.out = fileName, open = "a", nbchar = 60, as.string = F) 
  j <- j+1
}
