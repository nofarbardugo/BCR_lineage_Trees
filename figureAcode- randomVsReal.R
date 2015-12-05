source.dir <-'/home/nofar/Desktop/project papaers/dataForFinalGraph'
setwd(source.dir)
library(ggplot2)
library(reshape2)
library(grid)
library(scales)
library(plyr)
library(gridExtra)
library(cowplot)
library(lattice)
library(RGraphics)


# pfizer details
df.real_data<- read.csv('PatientIsotypeStatistic.csv', header = T, stringsAsFactors = F)
df.real_data$Isotype_number <- as.character(df.real_data$type_num)

figureA <-ggplot(df.real_data, aes(x=sum_Leaves,fill = Isotype_number)) + geom_histogram(aes(y=..density..),adjust = 0.1,position ="dodge") +
  scale_fill_brewer(palette="Set1")+
  xlab("Leaves number") + 
  xlim(c(0,100))+
  ylim(c(0,0.3)) +
  ylab("Density - Number of trees") +
  ggtitle("Density of number of leaves per number of isotypes") +
  theme(text = element_text(lineheight=.8, face="bold",size = 16)) + theme_gray(17) 

df.rand_data<- read.csv('Rand_PatientIsotypeStatistic.csv', header = T, stringsAsFactors = F)
df.rand_data$Isotype_number <- as.character(df.rand_data$type_num)
figureB <-ggplot(df.rand_data, aes(x=sum_Leaves,fill = Isotype_number)) + geom_histogram(aes(y=..density..),adjust = 0.1,position ="dodge") +
  scale_fill_brewer(palette="Set1")+
  xlab("Leaves number") + 
  xlim(c(0,100))+
  ylim(c(0,0.3)) +
  ylab("Density - Number of trees") +
  ggtitle("Random data - Density of number of leaves per number of isotypes") +
  theme(text = element_text(lineheight=.8, face="bold",size = 16)) + theme_gray(17) 


figure1 <- grid.arrange(figureA, figureB , ncol=1) 


#######################################################################################################3


df.flu<- read.csv('FLU_SEQ.csv', header = T, stringsAsFactors = F)
df.flu$VJDis <- substr(df.flu$VJDis,start =13 ,stop = length(df.flu$VJDis))


df.fluData <- ddply(df.flu, c("VJDis","ISOTYPE"), summarise,
                             SeqNumNotUniqe = length(SEQUENCE),
                             SeqNumUniqe = length(unique(SEQUENCE)))


df.pfizer<- read.csv('PFIZER_SEQ.csv', header = T, stringsAsFactors = F)

vectorSubject <-vector()

for(i in 1:length(df.sampleDate$VJDis))
{
  lMakaf <- gregexpr(pattern ='_', df.sampleDate[i,"SEQUENCE_ID"])
  startplace = lMakaf[[1]][length(lMakaf[[1]])-2]+1
  
  vectorSubject[i] <-substr(df.sampleDate[i,"SEQUENCE_ID"],start = startplace ,stop = startplace + 2)
}

df.sampleDate$subject <- vectorSubject

df.sampleDateToSend <- data.frame(df.sampleDate$SEQUENCE_ID,df.sampleDate$SEQUENCE,df.sampleDate$ISOTYPE,vectorSubject)

write.csv(df.sampleDateToSend,file = "PFIZER_SEQ.csv")

df.PfizerData <- ddply(df.pfizer, c("subject","ISOTYPE"), summarise,
                    SeqNumNotUniqe = length(SEQUENCE),
                    SeqNumUniqe = length(unique(SEQUENCE)))



