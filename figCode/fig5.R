
source.dir <-'/home/nofar/Desktop/project papaers/dataForFinalGraph'
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



grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

############### ~~~~~~~~~~~fig 5 ~~~~~~~~~~~~~~~~ ################################

#~~~~~ fig5A

# pfizer details

### IGHG
pfizer_dfIGHG <- read.csv('pfizer/PFIZER_event_IGHG.csv', header = T, stringsAsFactors = F)

# get the ratio beteewn each mutation type in each depth in pfizer
df.pfizer_Igg_mutation <- ddply(pfizer_dfIGHG, c("depth"), summarise,
                                mutationNumber = sum(depth),
                                meanSynFR = sum(synonyms_FR)/mutationNumber,
                                meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

# remove from her for no hols in aix x 
df.pfizer_Igg_mutation <- df.pfizer_Igg_mutation[df.pfizer_Igg_mutation$depth<35,]
# change represent of data set to fit ggplot format
df.pfizer_Igg_mutation <- subset(df.pfizer_Igg_mutation,select = -c(mutationNumber))
df.pfizer_igg_mutationRatio_depth<- melt(df.pfizer_Igg_mutation, id.vars = 'depth', variable.name='muteType')
df.pfizer_igg_mutationRatio_depth$isotype <- "IGHG"

#### IGHA
pfizer_dfIGHA <- read.csv('pfizer/PFIZER_event_IGHA.csv', header = T, stringsAsFactors = F)

# get the ratio beteewn each mutation type in each depth in pfizer
df.pfizer_Iga_mutation <- ddply(pfizer_dfIGHA, c("depth"), summarise,
                                mutationNumber = sum(depth),
                                meanSynFR = sum(synonyms_FR)/mutationNumber,
                                meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

# remove from her for no hols in aix x 
df.pfizer_Iga_mutation <- df.pfizer_Iga_mutation[df.pfizer_Iga_mutation$depth<38,]
# change represent of data set to fit ggplot format
df.pfizer_Iga_mutation <- subset(df.pfizer_Iga_mutation,select = -c(mutationNumber))
df.pfizer_iga_mutationRatio_depth<- melt(df.pfizer_Iga_mutation, id.vars = 'depth', variable.name='muteType')
df.pfizer_iga_mutationRatio_depth$isotype <- "IGHA"

#### IGHM
pfizer_dfIGHM <- read.csv('pfizer/PFIZER_event_IGHM.csv', header = T, stringsAsFactors = F)

# get the ratio beteewn each mutation type in each depth in pfizer
df.pfizer_Igm_mutation <- ddply(pfizer_dfIGHM, c("depth"), summarise,
                                mutationNumber = sum(depth),
                                meanSynFR = sum(synonyms_FR)/mutationNumber,
                                meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

# remove from her for no hols in aix x 
df.pfizer_Igm_mutation <- df.pfizer_Igm_mutation[df.pfizer_Iga_mutation$depth<51,]

# change represent of data set to fit ggplot format
df.pfizer_Igm_mutation <- subset(df.pfizer_Igm_mutation,select = -c(mutationNumber))
df.pfizer_igm_mutationRatio_depth<- melt(df.pfizer_Igm_mutation, id.vars = 'depth', variable.name='muteType')
df.pfizer_igm_mutationRatio_depth$isotype <- "IGHM"

#### naive IGHM
pfizer_dfnaiveIGM <- read.csv('pfizer/PFIZER_event_naive_IGHM.csv', header = T, stringsAsFactors = F)

# get the ratio beteewn each mutation type in each depth in pfizer
df.pfizer_naive_Igm_mutation <- ddply(pfizer_dfnaiveIGM, c("depth"), summarise,
                                      mutationNumber = sum(depth),
                                      meanSynFR = sum(synonyms_FR)/mutationNumber,
                                      meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                      meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                      meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

# remove from her for no hols in aix x 
df.pfizer_naive_Igm_mutation <- df.pfizer_naive_Igm_mutation[df.pfizer_naive_Igm_mutation$depth<51,]

# change represent of data set to fit ggplot format
df.pfizer_naive_Igm_mutation <- subset(df.pfizer_naive_Igm_mutation,select = -c(mutationNumber))
df.pfizer_naive_igm_mutationRatio_depth<- melt(df.pfizer_naive_Igm_mutation, id.vars = 'depth', variable.name='muteType')
df.pfizer_naive_igm_mutationRatio_depth$isotype <- "naive IGHM"

### combine all toogther
df.pfizer_combine_depth_mutaion <- rbind(df.pfizer_igg_mutationRatio_depth,
                                         df.pfizer_iga_mutationRatio_depth,
                                         df.pfizer_igm_mutationRatio_depth,
                                         df.pfizer_naive_igm_mutationRatio_depth)

df.pfizer_combine_depth_mutaion$type <- "Pfizer" 

# flu details

### IGHG
flu_dfIGHG <- rbind(read.csv('flu/FLU_event_IGHG1.csv', header = T, stringsAsFactors = F),
                    read.csv('flu/FLU_event_IGHG2.csv', header = T, stringsAsFactors = F))

# get the ratio beteewn each mutation type in each depth in flu
df.flu_Igg_mutation <- ddply(flu_dfIGHG, c("depth"), summarise,
                             mutationNumber = sum(depth),
                             meanSynFR = sum(synonyms_FR)/mutationNumber,
                             meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                             meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                             meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

# remove from her for no hols in aix x 
df.flu_Igg_mutation <- df.flu_Igg_mutation[df.flu_Igg_mutation$depth<41,]

# change represent of data set to fit ggplot format
df.flu_Igg_mutation <- subset(df.flu_Igg_mutation,select = -c(mutationNumber))
df.flu_igg_mutationRatio_depth<- melt(df.flu_Igg_mutation, id.vars = 'depth', variable.name='muteType')
df.flu_igg_mutationRatio_depth$isotype <- "IGHG"


#### IGHA
flu_dfIGHA <- read.csv('flu/FLU_event_IGHA.csv', header = T, stringsAsFactors = F)

# get the ratio beteewn each mutation type in each depth in flu
df.flu_Iga_mutation <- ddply(flu_dfIGHA, c("depth"), summarise,
                             mutationNumber = sum(depth),
                             meanSynFR = sum(synonyms_FR)/mutationNumber,
                             meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                             meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                             meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

# remove from her for no hols in aix x 
df.flu_Iga_mutation <- df.flu_Iga_mutation[df.flu_Iga_mutation$depth<39,]

# change represent of data set to fit ggplot format
df.flu_Iga_mutation <- subset(df.flu_Iga_mutation,select = -c(mutationNumber))
df.flu_iga_mutationRatio_depth<- melt(df.flu_Iga_mutation, id.vars = 'depth', variable.name='muteType')
df.flu_iga_mutationRatio_depth$isotype <- "IGHA"

#### IGHE
flu_dfIGHE <- read.csv('flu/FLU_event_IGHE.csv', header = T, stringsAsFactors = F)

# get the ratio beteewn each mutation type in each depth in flu
df.flu_Ige_mutation <- ddply(flu_dfIGHE, c("depth"), summarise,
                             mutationNumber = sum(depth),
                             meanSynFR = sum(synonyms_FR)/mutationNumber,
                             meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                             meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                             meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

# remove from her for no hols in aix x 
df.flu_Ige_mutation <- df.flu_Ige_mutation[df.flu_Ige_mutation$depth<9,]

# change represent of data set to fit ggplot format
df.flu_Ige_mutation <- subset(df.flu_Ige_mutation,select = -c(mutationNumber))
df.flu_ige_mutationRatio_depth<- melt(df.flu_Ige_mutation, id.vars = 'depth', variable.name='muteType')
df.flu_ige_mutationRatio_depth$isotype <- "IGHE"

#### IGHM
flu_dfIGHM <- read.csv('flu/FLU_event_IGHM.csv', header = T, stringsAsFactors = F)

# get the ratio beteewn each mutation type in each depth in flu
df.flu_Igm_mutation <- ddply(flu_dfIGHM, c("depth"), summarise,
                             mutationNumber = sum(depth),
                             meanSynFR = sum(synonyms_FR)/mutationNumber,
                             meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                             meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                             meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

df.flu_Igm_mutation <- df.flu_Igm_mutation[df.flu_Igm_mutation$depth<40,]
# change represent of data set to fit ggplot format
df.flu_Igm_mutation <- subset(df.flu_Igm_mutation,select = -c(mutationNumber))
df.flu_igm_mutationRatio_depth<- melt(df.flu_Igm_mutation, id.vars = 'depth', variable.name='muteType')
df.flu_igm_mutationRatio_depth$isotype <- "IGHM"

#### IGHD
flu_dfIGHD <- read.csv('flu/FLU_event_naive_IGHD.csv', header = T, stringsAsFactors = F)

# get the ratio beteewn each mutation type in each depth in flu
df.flu_Igd_mutation <- ddply(flu_dfIGHD, c("depth"), summarise,
                             mutationNumber = sum(depth),
                             meanSynFR = sum(synonyms_FR)/mutationNumber,
                             meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                             meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                             meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

# remove from her for no hols in aix x 
df.flu_Igd_mutation <- df.flu_Igd_mutation[df.flu_Igd_mutation$depth<35,]

# change represent of data set to fit ggplot format
df.flu_Igd_mutation <- subset(df.flu_Igd_mutation,select = -c(mutationNumber))
df.flu_igd_mutationRatio_depth<- melt(df.flu_Igd_mutation, id.vars = 'depth', variable.name='muteType')
df.flu_igd_mutationRatio_depth$isotype <- "IGHD"


### combine all toogther
df.flu_combine_depth_mutaion <- rbind(df.flu_igg_mutationRatio_depth,
                                      df.flu_iga_mutationRatio_depth,
                                  #    df.flu_ige_mutationRatio_depth,
                                      df.flu_igm_mutationRatio_depth,
                                      df.flu_igd_mutationRatio_depth)

df.flu_combine_depth_mutaion$type <- "Flu" 

# combine pfizer with flu 
df.combine_depth_mutaion <- rbind(df.pfizer_combine_depth_mutaion,df.flu_combine_depth_mutaion) 

figure5a1 <- ggplot(data=df.pfizer_combine_depth_mutaion, aes(x=depth, y=value, fill=muteType))  +
  facet_grid(isotype ~ type ,scale = "free_y") + 
  geom_bar(colour="black", stat="identity",size=.3) +   # Thinner lines
  scale_fill_brewer(palette="Set1")+
  ylab("") + # Set axis labels
  scale_y_continuous(breaks = c(0,0.5,1)) +
  xlab("") +
  scale_x_continuous(breaks = seq(from = 1,to = 20,by = 1),limits = c(0,21)) +
  
  ggtitle("Mutation distribution in each depth") +  # Set title
  background_grid(major = 'y', minor = "y")+
  scale_colour_discrete(name  ="Mutation type",
                        breaks=c("meanSynFR", "meanNonSynFR","meanSynCDR", "meanNonSynCDR"),
                        labels=c("Syn FR", "NonSyn FR","Syn CDR", "NonSyn CDR"))+
theme_linedraw()+
  theme(  panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank()
  )

figure5a2 <- ggplot(data=df.flu_combine_depth_mutaion, aes(x=depth, y=value, fill=muteType))  + 
  facet_grid(isotype ~ type ,scale = "free_y") + 
  geom_bar(colour="black", stat="identity",size=.3) +   # Thinner lines
  scale_fill_brewer(palette="Set1")+
  ylab("") + # Set axis labels
  scale_y_continuous(breaks = c(0,0.5,1)) +
  xlab("Depth") +
  scale_x_continuous(breaks = seq(from = 1,to = 20,by = 1),limits = c(0,21))+
  ggtitle("") +  # Set title
  background_grid(major = 'y', minor = "y")+
  scale_colour_discrete(name  ="Mutation type",
                        breaks=c("meanSynFR", "meanNonSynFR","meanSynCDR", "meanNonSynCDR"),
                        labels=c("Syn FR", "NonSyn FR","Syn CDR", "NonSyn CDR"))+ 
  theme_linedraw()+
  theme(  panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank()
           )

figure5a <- grid_arrange_shared_legend(figure5a1 , figure5a2)



########################################### ~~~ EDGES~~~~ #########################################

# pfizer details
df.pfizer_edges<- read.csv('pfizer/PFIZER_Edges_mutation_num_per_edge.csv', header = T, stringsAsFactors = F)

# get the ratio beteewn each mutation type in pfizer
df.pfizer_mutation_edges <- ddply(df.pfizer_edges, c("Edge.Name"), summarise,
                                  mutationNumber = sum(depth),
                                  meanSynFR = sum(synonyms_FR)/mutationNumber,
                                  meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                  meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                  meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

df.pfizer_mutation_edges <- subset(df.pfizer_mutation_edges,select = -c(mutationNumber))
df.pfizer_mutation_edges<- melt(df.pfizer_mutation_edges, id.vars = 'Edge.Name', variable.name='muteType')
df.pfizer_mutation_edges$type <- "Pfizer"

figure5b1 <-ggplot(data=df.pfizer_mutation_edges, aes(x=Edge.Name, y=value, fill=muteType)) + 
  geom_bar(colour="black",stat="identity", size=.3) +  # Thinner lines
  xlab("") + ylab("") + # Set axis labels
  scale_fill_brewer(palette="Set1")+
  ggtitle("Pfizer") +     # Set title
  theme(text = element_text(lineheight=.8,size = 14),legend.position="none",
        axis.text.x = element_text(angle = 90,hjust = 1))+
  coord_flip() 

# option 2
figure5b1 <-ggplot(data=df.pfizer_mutation_edges, aes(x=Edge.Name, y=value, fill=muteType)) + 
  geom_bar(colour="black",stat="identity", size=.3) +  # Thinner lines
  xlab("") + ylab("") + # Set axis labels
  scale_fill_brewer(palette="Set1")+
  ggtitle("Pfizer") +     # Set title
  theme_linedraw()+
  theme(   legend.position="none",
           panel.grid.minor=element_blank(), # remove grid
           panel.grid.major=element_blank()
  )+
  coord_flip() 


# flu details
df.flu_edges<- read.csv('flu/FLU_Edges_mutation_num_per_edge.csv', header = T, stringsAsFactors = F)

df.flu_edges$Edge.Name <-gsub("G1", "G", df.flu_edges$Edge.Name)
df.flu_edges$Edge.Name <-gsub("G2", "G", df.flu_edges$Edge.Name)

# get the ratio beteewn each mutation type in flu
df.flu_mutation_edges <- ddply(df.flu_edges, c("Edge.Name"), summarise,
                               mutationNumber = sum(depth),
                               meanSynFR = sum(synonyms_FR)/mutationNumber,
                               meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                               meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                               meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

df.flu_mutation_edges<-df.flu_mutation_edges[-c(3,6,11),]



df.flu_mutation_edges <- subset(df.flu_mutation_edges,select = -c(mutationNumber))
df.flu_mutation_edges<- melt(df.flu_mutation_edges, id.vars = 'Edge.Name', variable.name='muteType')
df.flu_mutation_edges$type <- "Flu"


# option one
figure5b2 <-ggplot(data=df.flu_mutation_edges, aes(x=Edge.Name, y=value, fill=muteType)) + 
  geom_bar(colour="black",stat="identity", size=.3) +  # Thinner lines
  xlab("") + ylab("") + # Set axis labels
  scale_fill_brewer(palette="Set1")+
  coord_flip() +
  ggtitle("Flu") +     # Set title
  theme(text = element_text(lineheight=.8,size = 14),
        legend.position="none",axis.text.x = element_text(angle = 90,hjust = 1))

# option two
figure5b2 <-ggplot(data=df.flu_mutation_edges, aes(x=Edge.Name, y=value, fill=muteType)) + 
  geom_bar(colour="black",stat="identity", size=.3) +  # Thinner lines
  xlab("") + ylab("") + # Set axis labels
  scale_fill_brewer(palette="Set1")+
  #coord_flip() +
  ggtitle("Flu") +     # Set title
  theme_linedraw()+
  theme(   legend.position="none",
           panel.grid.minor=element_blank(), # remove grid
           panel.grid.major=element_blank()
      )


figure5b <- grid.arrange(figure5b1, figure5b2,ncol =1)
