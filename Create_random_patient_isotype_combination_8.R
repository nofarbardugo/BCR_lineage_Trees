require("pracma")

# set directory : change path depending on work environment
#source.dir <- 'G:/project/files/pfizer'
#source.dir <-'C:/Users/Nofar/Desktop/trees'
source.dir <-'/home/bardugn1/pfizer' # pep4
setwd(source.dir)


########### MAIN #########################
statistic.dir <- paste0(source.dir,"/statistic/")
rand.dir <- paste0(source.dir,"/randStatistic/")
path <- "statistic/"
# get list of files
files = list.files(statistic.dir,full.names = F)

# going over each file
for(k in 1:length(files))
{
  # read patient tree details
  totalLeavesType.fd <-read.csv(file = paste0(path,files[k]),header = T, stringsAsFactors = F)
  
  ziroVec <-rep(0,nrow(totalLeavesType.fd))
  # enter each tree to the tree vector "leaves sum" times:
  treeSum.vector = c();
  index <- 1;
  # going over all each tree detailes
  for(i in 1:nrow(totalLeavesType.fd))
  {
    # add tree name to list X times
    for(j in 1:totalLeavesType.fd[i,"sum_Leaves"])
    {
      treeSum.vector[index] = totalLeavesType.fd[i,"treeNUM"]
      index <- index +1
    }
  }
  
  isotypeSum.vector = c();
  index <-1
  # create isotype Vector
  for(i in 3:8)
  {
    isotypeINumber <- sum (totalLeavesType.fd[,i])
    j <-0
    while(j < isotypeINumber)
    {
      isotypeSum.vector[index] = i-1
      index <- index +1
      j <- j +1
    }
  }
  
  # use randperm funcion for mix isotypeSum veactor randomly
  RandisotypeSum.vector <- randperm(isotypeSum.vector, length(isotypeSum.vector))
  
  # create data fream with the patiant random details
  fd.TotalLeavesType<- data.frame(treeNUM=c(totalLeavesType.fd[,"treeNUM"]) 
                                  ,IGHA = c(ziroVec),IGHE =c(ziroVec),IGHG = c(ziroVec)
                                  ,IGHM = c(ziroVec), naive_IGHD = c(ziroVec),naive_IGHM = c(ziroVec)
                                  ,sum_Leaves =c(ziroVec),type_num = c(ziroVec))
  
  # file "fd.TotalLeavesType"
  for (i in 1:length(isotypeSum.vector))
  {
    fd.TotalLeavesType[which(fd.TotalLeavesType[,"treeNUM"]==treeSum.vector[i]),(colnames(fd.TotalLeavesType)[(RandisotypeSum.vector[i])])]= fd.TotalLeavesType[which(fd.TotalLeavesType[,"treeNUM"]==treeSum.vector[i]),(colnames(fd.TotalLeavesType)[(RandisotypeSum.vector[i])])] +1;
  }
  
  # calc total leaves num 
  for (i in 1:nrow(fd.TotalLeavesType))
  {
    fd.TotalLeavesType[i,"sum_Leaves"] <- fd.TotalLeavesType[i,"IGHA"]+fd.TotalLeavesType[i,"IGHE"] +
                                          fd.TotalLeavesType[i,"IGHG"] + fd.TotalLeavesType[i,"IGHM"] +
                                          fd.TotalLeavesType[i,"naive_IGHD"]+ fd.TotalLeavesType[i,"naive_IGHM"]
  
    # get type isotype num 
    num <- 0
    for (t in 2:7)
    {
      if (fd.TotalLeavesType[i,t] >0)
      {
        num <- num + 1
      }
    }
  
    fd.TotalLeavesType[i,"type_num"] <-num
  }
  
 # save the TotalLeavesType table into a file
 write.csv(fd.TotalLeavesType, file = paste0(rand.dir,files[k]))
} 