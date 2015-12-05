source.dir <-'/home/nofar/Desktop/LabProject/data'
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


# pfizer

# details of all trees
pfizer_all_dfIGHG <- read.csv('pfizer/PFIZER_event_ALL_IGHG.csv', header = T, stringsAsFactors = F)
pfizer_all_dfIGHG$isotype <- "IGHG"

pfizer_all_dfIGHM <- read.csv('pfizer/PFIZER_event_ALL_IGHM.csv', header = T, stringsAsFactors = F)
pfizer_all_dfIGHM$isotype <- "IGHM"

pfizer_all_dfIGHA <- read.csv('pfizer/PFIZER_event_ALL_IGHA.csv', header = T, stringsAsFactors = F)
pfizer_all_dfIGHA$isotype <- "IGHA"

pfizer_all_dfNaiveIGHM <- read.csv('pfizer/PFIZER_event_ALL_naive_IGHM.csv', header = T, stringsAsFactors = F)
pfizer_all_dfNaiveIGHM$isotype <- "Naive_IGHM"

# merge all into one dataFream
pfizer_ds_event_all <- rbind(pfizer_all_dfIGHG,pfizer_all_dfIGHM,pfizer_all_dfIGHA,pfizer_all_dfNaiveIGHM)

# get number of isotypes for each type
pfizer_isotypeSize_all <- ddply(pfizer_ds_event_all, c("isotype"), summarise,
                            repeatTimes = length(depth)
                            
)

# get number of isotype in each depth for each type
pfizer_ds_event_to_point_all <- ddply(pfizer_ds_event_all, c("isotype","depth"), summarise,
                                  repeatTimes = length(depth)                                  
)

#  enter the fracion for each row (number of isotype in a specific depth divide by total sum of the suitable isotype) 
for(i in 1: nrow(pfizer_ds_event_to_point_all))
{
  pfizer_ds_event_to_point_all[i,"totalNode"] <- pfizer_isotypeSize_all[pfizer_isotypeSize_all$isotype==pfizer_ds_event_to_point_all[i,"isotype"],"repeatTimes"]
}
pfizer_ds_event_to_point_all$fraction <- pfizer_ds_event_to_point_all$repeatTimes/pfizer_ds_event_to_point_all$totalNode
pfizer_ds_event_to_point_all$SE <- pfizer_ds_event_to_point_all$fraction *(1/sqrt(pfizer_ds_event_to_point_all$repeatTimes))


# plot a xy plot of the data
figure3a1 <-ggplot(pfizer_ds_event_to_point_all , aes(x=depth, y=fraction, colour=isotype)) +
  geom_line() + geom_point()  +
  geom_errorbar(aes(ymin=fraction-SE, ymax=fraction+SE), width=0.4) +
  xlab("Depth: S+NS") + 
  scale_y_log10(breaks = c(0.00001,0.0001,0.001,0.01,0.1)) +
  #scale_x_log10(breaks = c(0.5,1,10,20,30,40,50))+
  ylab("Density - Number of nodes") +
  ggtitle("Histogram of mutations per isotype - Pfizer data") +
  scale_colour_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","Naive_IGHM" = "green")) +
  theme_linedraw()+
  theme(panel.grid.major = element_line(colour = "blue"))
#  theme(  legend.title=element_blank(),
#          legend.position="none", # legend position
          #legend.justification=c(1,1),
#          panel.grid.minor=element_blank(), # remove grid
#          panel.grid.major=element_blank()
         
 # ) #+
  #guides(fill=guide_legend(nrow=2))

figure3a2 <-ggplot(pfizer_ds_event_to_point_all , aes(x=depth, y=fraction, colour=isotype)) +
  geom_line() + geom_point()  +
  geom_errorbar(aes(ymin=fraction-SE, ymax=fraction+SE), width=0.04) +
  xlab("Depth: S+NS") + 
  scale_y_log10(limits = c(0.001,1)) +
  scale_x_log10(breaks = seq(from = 1,to = 5,by = 1),limits = c(0.5,5))+
  ylab("Density - Number of nodes") +
  ggtitle("Histogram of mutations per isotype - Pfizer data") +
  scale_colour_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","Naive_IGHM" = "green")) +
  theme_linedraw()+
  theme(panel.grid.major = element_line(colour = "blue"))
#  theme(  legend.title=element_blank(),
#          legend.position="none", # legend position
   #       legend.position=c(1,1), # legend position
#          legend.justification=c(1,1),
#          panel.grid.minor=element_blank(), # remove grid
#          panel.grid.major=element_blank()
          
#  )# +
  #guides(fill=guide_legend(nrow=2))

#vp <- viewport(width = 0.4, height = 0.4, x = 1,
#               y = unit(20, "lines"), just = c("right","bottom"))

#figure3a <- plot(figure3a1) 
#figure3a <- figure3a + figure3a2 + viewport(width = 0.4, height = 0.4, x = 1,
#                                            y = unit(20, "lines"), just = c("right","bottom")) 



# details of only combine trees
pfizer_dfIGHG <- read.csv('pfizer/PFIZER_event_IGHG.csv', header = T, stringsAsFactors = F)
pfizer_dfIGHG$isotype <- "IGHG"

pfizer_dfIGHM <- read.csv('pfizer/PFIZER_event_IGHM.csv', header = T, stringsAsFactors = F)
pfizer_dfIGHM$isotype <- "IGHM"

pfizer_dfIGHA <- read.csv('pfizer/PFIZER_event_IGHA.csv', header = T, stringsAsFactors = F)
pfizer_dfIGHA$isotype <- "IGHA"

pfizer_dfNaiveIGHM <- read.csv('pfizer/PFIZER_event_naive_IGHM.csv', header = T, stringsAsFactors = F)
pfizer_dfNaiveIGHM$isotype <- "Naive_IGHM"

# merge all into one dataFream
pfizer_ds_event <- rbind(pfizer_dfIGHG,pfizer_dfIGHM,pfizer_dfIGHA,pfizer_dfNaiveIGHM)

# get number of isotypes for each type
pfizer_isotypeSize <- ddply(pfizer_ds_event, c("isotype"), summarise,
                            repeatTimes = length(depth)
                            
)

# get number of isotype in each depth for each type
pfizer_ds_event_to_point <- ddply(pfizer_ds_event, c("isotype","depth"), summarise,
                                  repeatTimes = length(depth)                                  
)

#  enter the fracion for each row (number of isotype in a specific depth divide by total sum of the suitable isotype) 
for(i in 1: nrow(pfizer_ds_event_to_point))
{
  pfizer_ds_event_to_point[i,"totalNode"] <- pfizer_isotypeSize[pfizer_isotypeSize$isotype==pfizer_ds_event_to_point[i,"isotype"],"repeatTimes"]
}
pfizer_ds_event_to_point$fraction <- pfizer_ds_event_to_point$repeatTimes/pfizer_ds_event_to_point$totalNode
pfizer_ds_event_to_point$SE <- pfizer_ds_event_to_point$fraction *(1/sqrt(pfizer_ds_event_to_point$repeatTimes))
# plot a xy plot of the data
figure3b1 <-ggplot(pfizer_ds_event_to_point, aes(x=depth, y=fraction, colour=isotype)) +
  geom_line() + geom_point()  +
  geom_errorbar(aes(ymin=fraction-SE, ymax=fraction+SE), width=0.4) +
  xlab("Depth: S+NS")+
  scale_y_log10(breaks = c(0.00001,0.0001,0.001,0.01,0.1))+
  #scale_x_log10(breaks = c(0.5,1,10,20,30,40,50))+
  ylab("Density - Number of nodes") +
  ggtitle("Histogram of mutations per isotype - Pfizer combine isotypes only") +
  scale_colour_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","Naive_IGHM" = "green")) +
  theme_linedraw()+
  theme(panel.grid.major = element_line(colour = "blue"))
#  theme(  legend.title=element_blank(),
#          legend.position=c(1,1), # legend position
#          legend.justification=c(1,1),
#          panel.grid.minor=element_blank(), # remove grid
#          panel.grid.major=element_blank()
          
#  ) +
#  guides(fill=guide_legend(nrow=2))

figure3b2 <-ggplot(pfizer_ds_event_to_point, aes(x=depth, y=fraction, colour=isotype)) +
  geom_line() + geom_point()  +
  geom_errorbar(aes(ymin=fraction-SE, ymax=fraction+SE), width=0.04) +
  xlab("Depth: S+NS")+
  scale_y_log10(limits = c(0.001,1)) +
  scale_x_log10(breaks = seq(from = 1,to = 5,by = 1),limits = c(0.5,5))+
  ylab("Density - Number of nodes") +
  ggtitle("Histogram of mutations per isotype - Pfizer combine isotypes only") +
  scale_colour_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","Naive_IGHM" = "green")) +
  theme_linedraw()+
  theme(panel.grid.major = element_line(colour = "blue"))
#  theme(  legend.title=element_blank(),
#          legend.position="none", # legend position
        #  legend.position=c(1,1), # legend position
#          legend.justification=c(1,1),
#          panel.grid.minor=element_blank(), # remove grid
#          panel.grid.major=element_blank()
          
#  ) #+
  #guides(fill=guide_legend(nrow=2))

#figure3b <- print(figure3b1) + print(figure3b2, vp = vp)

         
####################################################################################################
 
# flu

# details of all trees

flu_dfIGHG1_all <- read.csv('flu/FLU_event_all_IGHG1.csv', header = T, stringsAsFactors = F)
flu_dfIGHG1_all$isotype <- "IGHG"

flu_dfIGHG2_all<- read.csv('flu/FLU_event_all_IGHG2.csv', header = T, stringsAsFactors = F)
flu_dfIGHG2_all$isotype <- "IGHG"

flu_dfIGHA_all <- read.csv('flu/FLU_event_all_IGHA.csv', header = T, stringsAsFactors = F)
flu_dfIGHA_all$isotype <- "IGHA"

# too litle data - remove for now
#flu_dfIGHE_all <- read.csv('flu/FLU_event_all_IGHE.csv', header = T, stringsAsFactors = F)
#flu_dfIGHE_all$isotype <- "IGHE"

flu_dfIGHM_all <- read.csv('flu/FLU_event_all_IGHM.csv', header = T, stringsAsFactors = F)
flu_dfIGHM_all$isotype <- "IGHM"

flu_dfIGHD_all <- read.csv('flu/FLU_event_all_naive_IGHD.csv', header = T, stringsAsFactors = F)
flu_dfIGHD_all$isotype <- "IGHD"

# merge all into one dataFream
flu_ds_event_all <- rbind(flu_dfIGHG1_all,flu_dfIGHG2_all,flu_dfIGHA_all,flu_dfIGHM_all,flu_dfIGHD_all)#,flu_dfIGHE_all,)

# get number of isotypes for each type
flu_isotypeSize_all <- ddply(flu_ds_event_all, c("isotype"), summarise,
                         repeatTimes = length(depth)
                         
)

# get number of isotype in each depth for each type
flu_ds_event_to_point_all <- ddply(flu_ds_event_all, c("isotype","depth"), summarise,
                               repeatTimes = length(depth)                                  
)

#  enter the fracion for each row (number of isotype in a specific depth divide by total sum of the suitable isotype) 
for(i in 1: nrow(flu_ds_event_to_point_all))
{
  flu_ds_event_to_point_all[i,"totalNode"] <- flu_isotypeSize_all[flu_isotypeSize_all$isotype==flu_ds_event_to_point_all[i,"isotype"],"repeatTimes"]
}
flu_ds_event_to_point_all$fraction <- flu_ds_event_to_point_all$repeatTimes/flu_ds_event_to_point_all$totalNode
flu_ds_event_to_point_all$SE <- flu_ds_event_to_point_all$fraction *(1/sqrt(flu_ds_event_to_point_all$repeatTimes))
figure3c1 <-ggplot(flu_ds_event_to_point_all, aes(x=depth, y=fraction, colour=isotype)) +
  geom_line() + geom_point()  +
  geom_errorbar(aes(ymin=fraction-SE, ymax=fraction+SE), width=0.4) +
  xlab("Depth: S+NS")+
  scale_y_log10(breaks = c(0.00001,0.0001,0.001,0.01,0.1))+
 # scale_x_log10(breaks = c(0.5,1,10,20,30,40,50))+
  ylab("Density - Number of nodes") +
  ggtitle("Histogram of mutations per isotype - Flu data") +
  scale_colour_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","Naive_IGHM" = "green","IGHD" = "#FF3399")) +
  theme_linedraw()+
  theme(panel.grid.major = element_line(colour = "blue"))
#  theme(  legend.title=element_blank(),
#          legend.position="none", # legend position
    #      legend.position=c(1,1), # legend position
#          legend.justification=c(1,1),
#          panel.grid.minor=element_blank(), # remove grid
#          panel.grid.major=element_blank()
          
 # )# +
 # guides(fill=guide_legend(nrow=2))

figure3c2 <-ggplot(flu_ds_event_to_point_all, aes(x=depth, y=fraction, colour=isotype)) +
  geom_line() + geom_point()  +
  geom_errorbar(aes(ymin=fraction-SE, ymax=fraction+SE), width=0.04) +
  xlab("Depth: S+NS")+
  scale_y_log10(limits = c(0.001,1)) +
  scale_x_log10(breaks = seq(from = 1,to = 5,by = 1),limits = c(0.5,5))+
  ylab("Density - Number of nodes") +
  ggtitle("Histogram of mutations per isotype - Flu data") +
  scale_colour_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","Naive_IGHM" = "green","IGHD" = "#FF3399")) +
  theme_linedraw()+
  theme(panel.grid.major = element_line(colour = "blue"))

#  theme(  legend.title=element_blank(),
#          legend.position="none", # legend position
          #      legend.position=c(1,1), # legend position
#          legend.justification=c(1,1),
#          panel.grid.minor=element_blank(), # remove grid
#          panel.grid.major=element_blank()
          
#  )

#figure3c <- print(figure3c1) + print(figure3c2, vp = vp)

# details of only combine trees
flu_dfIGHG1 <- read.csv('flu/FLU_event_IGHG1.csv', header = T, stringsAsFactors = F)
flu_dfIGHG1$isotype <- "IGHG"

flu_dfIGHG2<- read.csv('flu/FLU_event_IGHG2.csv', header = T, stringsAsFactors = F)
flu_dfIGHG2$isotype <- "IGHG"

flu_dfIGHA <- read.csv('flu/FLU_event_IGHA.csv', header = T, stringsAsFactors = F)
flu_dfIGHA$isotype <- "IGHA"

# too litle data - remove for now
#flu_dfIGHE <- read.csv('flu/FLU_event_IGHE.csv', header = T, stringsAsFactors = F)
#flu_dfIGHE$isotype <- "IGHE"

flu_dfIGHM <- read.csv('flu/FLU_event_IGHM.csv', header = T, stringsAsFactors = F)
flu_dfIGHM$isotype <- "IGHM"

flu_dfIGHD <- read.csv('flu/FLU_event_naive_IGHD.csv', header = T, stringsAsFactors = F)
flu_dfIGHD$isotype <- "IGHD"

# merge all into one dataFream
flu_ds_event <- rbind(flu_dfIGHG1,flu_dfIGHG2,flu_dfIGHA,flu_dfIGHM,flu_dfIGHD)#,flu_dfIGHE,)

# get number of isotypes for each type
flu_isotypeSize <- ddply(flu_ds_event, c("isotype"), summarise,
                            repeatTimes = length(depth)
                           
)

# get number of isotype in each depth for each type
flu_ds_event_to_point <- ddply(flu_ds_event, c("isotype","depth"), summarise,
                                  repeatTimes = length(depth)                                  
)

#  enter the fracion for each row (number of isotype in a specific depth divide by total sum of the suitable isotype) 
for(i in 1: nrow(flu_ds_event_to_point))
{
  flu_ds_event_to_point[i,"totalNode"] <- flu_isotypeSize[flu_isotypeSize$isotype==flu_ds_event_to_point[i,"isotype"],"repeatTimes"]
}
flu_ds_event_to_point$fraction <- flu_ds_event_to_point$repeatTimes/flu_ds_event_to_point$totalNode
flu_ds_event_to_point$SE <- flu_ds_event_to_point$fraction *(1/sqrt(flu_ds_event_to_point$repeatTimes))
figure3d1 <-ggplot(flu_ds_event_to_point, aes(x=depth, y=fraction, colour=isotype)) +
  geom_line() + geom_point()  +
  geom_errorbar(aes(ymin=fraction-SE, ymax=fraction+SE), width=0.4) +
  xlab("Depth: S+NS")+
  scale_y_log10(breaks = c(0.00001,0.0001,0.001,0.01,0.1))+
  #scale_x_log10(breaks = c(0.5,1,10,20,30,40,50))+
  ylab("Density - Number of nodes") +
  ggtitle("Histogram of mutations per isotype - Flu combine isotypes only") +
  scale_colour_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","Naive_IGHM" = "green","IGHD" = "#FF3399")) +
  theme_linedraw()+
  theme(panel.grid.major = element_line(colour = "blue"))
#  theme(  legend.title=element_blank(),
#          legend.position=c(1,1), # legend position
#          legend.justification=c(1,1),
#          panel.grid.minor=element_blank(), # remove grid
#          panel.grid.major=element_blank()
          
#  ) +
#  guides(fill=guide_legend(nrow=2))

figure3d2 <-ggplot(flu_ds_event_to_point, aes(x=depth, y=fraction, colour=isotype)) +
  geom_line() + geom_point()  +
  geom_errorbar(aes(ymin=fraction-SE, ymax=fraction+SE), width=0.04) +
  xlab("Depth: S+NS")+
  scale_y_log10(limits = c(0.001,1)) +
  scale_x_log10(breaks = seq(from = 1,to = 5,by = 1),limits = c(0.5,5))+
  ylab("Density - Number of nodes") +
  ggtitle("Histogram of mutations per isotype - Flu combine isotypes only") +
  scale_colour_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","Naive_IGHM" = "green","IGHD" = "#FF3399")) +
  theme_linedraw()+
  theme(panel.grid.major = element_line(colour = "blue"))
#  theme(  legend.title=element_blank(),
#          legend.position=c(1,1), # legend position
#          legend.justification=c(1,1),
#          panel.grid.minor=element_blank(), # remove grid
#          panel.grid.major=element_blank()
          
#  ) +
#  guides(fill=guide_legend(nrow=2))

#figure3d <- print(figure3d1) + print(figure3d2, vp = vp)

figure3 <- grid.arrange(figure3a1,figure3a2,
                        figure3b1 ,figure3b2 ,
                        figure3c1,figure3c2,
                        figure3d1,figure3d2, ncol=2) 
