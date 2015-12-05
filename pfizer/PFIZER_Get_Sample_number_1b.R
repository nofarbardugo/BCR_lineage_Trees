
# imports
require(seqinr)
require(Biostrings)
require(plyr)

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
  return (data.frame( df[,"SEQUENCE_ID"], VJDis.vector, df[,"SEQUENCE"], stringsAsFactors=F))
}

########################main#############################################

# set directory - when run need to change her the direction

# get list of folders
path <- "/home/bardugn1/pfizer/Formatted/" # pep4
dirs <- list.files(path , full.names = F)
#dirs <- list.files( full.names = F) 

# going over each folder and get it's files

print ("start")
for(i in 1:length(dirs))
{
  # get list of files 
  folderPath <- paste0(path,dirs[i]) 
  files <- list.files(folderPath, full.names = F)
  
  # save data for each isotype
  df <- data.frame("SEQUENCE_ID" = c(0),"VJDis" = c(0),"SEQUENCE" = c(0))
  
  # going over each file  
  for(j in 1:length(files))
  { 
    
    if( grepl("IB", substr(files[j],start= 1, stop=2))==F)
    {
      # get data frame of the file and add it to current dataframe
      if(j==1)
      {
        df <- Enter_Isotype(folderPath,dirs[i], files[j])
        
      }
      else
      {
        df <- rbind(df,Enter_Isotype(folderPath,dirs[i], files[j]))
      }
    }
  }
  
  
  if(i ==1)
  {
    df.sampleDate <-df
    df.sampleDate$Isotype <- dirs[i]
    names(df.sampleDate) <- c("SEQUENCE_ID","VJDis","SEQUENCE","ISOTYPE")
  }
  else
  {
    tmp <-df
    tmp$Isotype <- dirs[i]
    names(tmp) <- c("SEQUENCE_ID","VJDis","SEQUENCE","ISOTYPE")
    names(df.sampleDate) <- c("SEQUENCE_ID","VJDis","SEQUENCE","ISOTYPE")
    df.sampleDate <-rbind(tmp,df.sampleDate)       
  }
  
}

x <- 5


