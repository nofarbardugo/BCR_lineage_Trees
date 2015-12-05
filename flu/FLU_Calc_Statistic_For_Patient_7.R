# set directory : change path depending on work environment
#source.dir <- 'G:/project/files/pfizer' # window
#source.dir <-'C:/Users/Nofar/Desktop/trees' #window
#source.dir<- '/media/nofar/CDD8-A371/project/files/pfizer' # ubuntu
source.dir <-'/home/bardugn1/flu2' # pep4
setwd(source.dir)

########### MAIN #########################
# change here - depend on rand or real data. 
#statistic.dir <- paste0(source.dir,"/randStatistic/")
statistic.dir <- paste0(source.dir,"/statistic/")
#path <- "randStatistic/"
path <- "statistic/"
# get list of files

files = list.files(statistic.dir,full.names = F)

ziroVec <-rep(0,length(files))

# create data fram for the histogram plot
fd.hist <- data.frame (sum_Leaves = c(),type_num = c())

# to remove after run - just for get details
fd.PatientIsotype<- data.frame(PatientNum=c(files) ,"naive_IGHD" = c(ziroVec),"IGHG.2" = c(ziroVec),"IGHG.1" = c(ziroVec),"IGHM" = c(ziroVec),"IGHE" = c(ziroVec),"IGHA" = c(ziroVec),
                               "IGHA + IGHE" = c(ziroVec),"IGHA + IGHG.1" = c(ziroVec),"IGHA + IGHM" = c(ziroVec),
                               "IGHA + N_IGHD" = c(ziroVec),"IGHA + IGHG.2" = c(ziroVec),"IGHE + IGHG.1" = c(ziroVec),
                               "IGHE + IGHM" = c(ziroVec),"IGHE + naive_IGHD" = c(ziroVec),"IGHE + IGHG.2" = c(ziroVec),
                               "IGHG.1 + IGHM" = c(ziroVec),"IGHG.1 + naive_IGHD" = c(ziroVec),"IGHG.1 + IGHG.2" = c(ziroVec),
                               "IGHM + naive_IGHD" = c(ziroVec),"IGHM + IGHG.2" = c(ziroVec),"IGHG.2 + naive_IGHD" = c(ziroVec),
                               "IGHA + IGHE + IGHG.1" = c(ziroVec),"IGHA + IGHE + IGHM" = c(ziroVec),
                               "IGHA + IGHE + naive_IGHD" = c(ziroVec),"IGHA + IGHE + IGHG.2" = c(ziroVec),  
                               "IGHA + IGHG.1 + IGHM" = c(ziroVec),"IGHA + IGHG.1 + naive_IGHD" = c(ziroVec),
                               "IGHA + IGHG.1 + IGHG.2" = c(ziroVec),"IGHA + IGHM + naive_IGHD" = c(ziroVec),
                               "IGHA + IGHG.2 + IGHM" = c(ziroVec),"IGHA + IGHG.2 + naive_IGHD" = c(ziroVec),
                               "IGHE + IGHG.1 + IGHM" = c(ziroVec), "IGHE + IGHG.1 + naive_IGHD" = c(ziroVec),
                               "IGHE + IGHG.1 + IGHG.2" = c(ziroVec), "IGHE + IGHM + naive_IGHD" = c(ziroVec),
                               "IGHE + IGHM + IGHG.2" = c(ziroVec), "IGHE + naive_IGHD + IGHG.2" = c(ziroVec),
                               "IGHG.1 + IGHM + naive_IGHD" = c(ziroVec), "IGHG.1 + IGHM + IGHG.2" = c(ziroVec),
                               "IGHG.1 + IGHG.2 + naive_IGHD " = c(ziroVec),"IGHM + naive_IGHD + IGHG.2" = c(ziroVec),
                               "IGHA + IGHE + IGHG.1 + IGHM"  = c(ziroVec),"IGHA + IGHE + IGHG.1 + naive_IGHD"  = c(ziroVec),
                               "IGHA + IGHE + IGHG.1 + IGHG.2"  = c(ziroVec),"IGHA + IGHE + IGHM + naive_IGHD"  = c(ziroVec),
                               "IGHA + IGHE + IGHM + IGHG.2"  = c(ziroVec),"IGHA + IGHE + naive_IGHD + IGHG.2"  = c(ziroVec),
                               "IGHA + IGHG.1 + naive_IGHD + IGHG.2"  = c(ziroVec),"IGHA + IGHG.1 + IGHM + naive_IGHD"  = c(ziroVec),
                               "IGHA + IGHG.1 + IGHM + IGHG.2"  = c(ziroVec),"IGHA + IGHM + naive_IGHD + IGHG.2"  = c(ziroVec),
                               "IGHE + IGHG.1 + IGHM + naive_IGHD"  = c(ziroVec),"IGHE + IGHG.1 + IGHM + IGHG.2"  = c(ziroVec),
                               "IGHE + IGHG.1 + IGHG.2 + naive_IGHD" = c(ziroVec),
                               "IGHE + IGHM + naive_IGHD + IGHG.2"  = c(ziroVec), "IGHG.1 + IGHM + naive_IGHD + IGHG.2"  = c(ziroVec),
                               "IGHA + IGHE + IGHG.1 + IGHM + naive_IGHD"  = c(ziroVec),
                               "IGHA + IGHE + IGHG.1 + IGHM + IGHG.2"  = c(ziroVec),
                               "IGHA + IGHG.1 + IGHM + naive_IGHD + IGHG.2"  = c(ziroVec),
                               "IGHA + IGHE + IGHM + naive_IGHD + IGHG.2"  = c(ziroVec),
                               "IGHA + IGHE + IGHG.1 + naive_IGHD + IGHG.2"  = c(ziroVec),
                               "IGHE + IGHG.1 + IGHM + naive_IGHD + IGHG.2"  = c(ziroVec),
                               'IGHA + IGHE + IGHG.1 + IGHM + naive_IGHD + IGHG.2'  = c(ziroVec))

# give names to data fream
names(fd.PatientIsotype) <- c("PatientNum" ,"naive_IGHD","IGHG.2","IGHG.1" ,"IGHM" ,"IGHE","IGHA",
                              "IGHA + IGHE" ,"IGHA + IGHG.1" ,"IGHA + IGHM" ,
                              "IGHA + naive_IGHD" ,"IGHA + IGHG.2","IGHE + IGHG.1",
                              "IGHE + IGHM","IGHE + naive_IGHD","IGHE + IGHG.2" ,
                              "IGHG.1 + IGHM","IGHG.1 + naive_IGHD" ,"IGHG.1 + IGHG.2",
                              "IGHM + naive_IGHD","IGHG.2 + IGHM","IGHG.2 + naive_IGHD" ,
                              "IGHA + IGHE + IGHG.1","IGHA + IGHE + IGHM",
                              "IGHA + IGHE + naive_IGHD" ,"IGHA + IGHE + IGHG.2" ,  
                              "IGHA + IGHG.1 + IGHM" ,"IGHA + IGHG.1 + naive_IGHD",
                              "IGHA + IGHG.1 + IGHG.2" ,"IGHA + IGHM + naive_IGHD",
                              "IGHA + IGHG.2 + IGHM","IGHA + IGHG.2 + naive_IGHD",
                              "IGHE + IGHG.1 + IGHM" , "IGHE + IGHG.1 + naive_IGHD",
                              "IGHE + IGHG.1 + IGHG.2", "IGHE + IGHM + naive_IGHD",
                              "IGHE + IGHM + IGHG.2", "IGHE + naive_IGHD + IGHG.2",
                              "IGHG.1 + IGHM + naive_IGHD" , "IGHG.1 + IGHG.2 + IGHM",
                              "IGHG.1 + IGHG.2 + naive_IGHD" ,"IGHG.2 + IGHM + naive_IGHD" ,
                              "IGHA + IGHE + IGHG.1 + IGHM" ,"IGHA + IGHE + IGHG.1 + naive_IGHD",
                              "IGHA + IGHE + IGHG.1 + IGHG.2","IGHA + IGHE + IGHM + naive_IGHD",
                              "IGHA + IGHE + IGHG.2 + IGHM","IGHA + IGHE + IGHG.2 + naive_IGHD " ,
                              "IGHA + IGHG.1 + IGHG.2 + naive_IGHD","IGHA + IGHG.1 + IGHM + naive_IGHD" ,
                              "IGHA + IGHG.1 + IGHG.2 + IGHM" ,"IGHA + IGHG.2 + IGHM + naive_IGHD",
                              "IGHE + IGHG.1 + IGHM + naive_IGHD" ,"IGHE + IGHG.1 + IGHG.2 + IGHM",
                              "IGHE + IGHG.1 + IGHG.2 + naive_IGHD",
                              "IGHE + IGHG.2 + IGHM + naive_IGHD", "IGHG.1 + IGHG.2 + IGHM + naive_IGHD",
                              "IGHA + IGHE + IGHG.1 + IGHM + naive_IGHD",
                              "IGHA + IGHE + IGHG.1 + IGHG.2 + IGHM",
                              "IGHA + IGHG.1 + IGHG.2 + IGHM + naive_IGHD",
                              "IGHA + IGHE + IGHG.2 + IGHM + naive_IGHD",
                              "IGHA + IGHE + IGHG.1 + IGHG.2 + naive_IGHD",
                              "IGHE + IGHG.1 + IGHG.2 + IGHM + naive_IGHD",
                              "IGHA + IGHE + IGHG.1 + IGHG.2 + IGHM + naive_IGHD")

# going over each file
for(k in 1:length(files))
{
  # reas patient tree details
  totalLeavesType.fd <-read.csv(file = paste0(path,files[k]),header = T, stringsAsFactors = F)
  
  # 
  fd.hist <- rbind (fd.hist,totalLeavesType.fd[,9:10])
  
  # going over all each tree detailes
  for(i in 1:nrow(totalLeavesType.fd))
  {
    fieldName <- ""
    for(j in 3:8)
    {
      # if there is a leave from current isotype add to fieldName
      if(totalLeavesType.fd[i,j]>0)
      {
        fieldName <- paste0(fieldName," + ",colnames(totalLeavesType.fd[j]))
      }
    }
    
    # delete remainder from tail
    fieldName <-substr(x = fieldName ,start = 4, stop = nchar(fieldName))
    
    
   # if((nchar(fieldName) != 4) && nchar(fieldName) != 10)
  #  {
      fd.PatientIsotype [k,fieldName] <- fd.PatientIsotype [k,fieldName] + 1
  #  }
    
    
  }
  
 
}
write.csv(fd.PatientIsotype, file = paste0(statistic.dir,"PatientIsotypeStatistic.csv"))

#### plot histogram#############

d <- density(fd.hist[which(fd.hist[,"type_num"]== 2),"sum_Leaves"],width =10, from =0, to =100 )
plot (d,main = "Data", type='l',col = 2,lwd = 2)
# create histogram for 2,3, 4 isotype type 
for(i in 3:4)
{
  d <- density(fd.hist[which(fd.hist[,"type_num"]== i),"sum_Leaves"],width =10, from =0, to =100 )
  points(d,type='l',col = i,lwd = 2)
  # hist(main = i, xlab = "number of leaves", ylab = "number of trees",fd.hist[which(fd.hist[,"type_num"]== i),"sum_Leaves"],breaks  = 100)
}

legend ( "topright" ,4 , legend =c("1 isotype","2 isotype","3 isotype","4 isotype"), lwd = 4:4,
         col = 1:4 , lty = 1:1 )

