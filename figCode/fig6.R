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

df.pfizer_mut_germ<- read.csv('pfizer/PFIZER_RS.csv', header = T, stringsAsFactors = F)

# get the ratio beteewn each mutation type in pfizer
df.pfizer_mutation_germ <- ddply(df.pfizer_mut_germ, c("isotype"), summarise,
                                 mutationNumber = sum(distance),
                                 meanSynFR = sum(synonyms_FR)/mutationNumber,
                                 meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                 meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                 meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

df.pfizer_mutation_germ <- subset(df.pfizer_mutation_germ,select = -c(mutationNumber))
df.pfizer_mutation_germ<- melt(df.pfizer_mutation_germ, id.vars = 'isotype', variable.name='muteType')
df.pfizer_mutation_germ$type <- "Pfizer"



figure3c1 <-ggplot(data=df.pfizer_mutation_germ, aes(x=isotype, y=value, fill=muteType)) + 
  geom_bar(colour="black",stat="identity", size=.3) +  # Thinner lines
  xlab("") + ylab("") + # Set axis labels
  scale_fill_brewer(palette="Set1")+
  ggtitle("Pfizer") +     # Set title
  theme(text = element_text(lineheight=.8,size = 14),
        axis.text.x = element_text(angle = 90,hjust = 1))+
  coord_flip() 

#############################################################################################

##### `~~~~~~ pfizer all data ~~~~~ #########
df.pfizer_mut_germ<- read.csv('pfizer/PFIZER_RS.csv', header = T, stringsAsFactors = F)

## naive
df.pfizer_germ_mutation_naive<- df.pfizer_mut_germ[df.pfizer_mut_germ$isotype=="naive_IGHM" ,]
df.pfizer_germ_mutation_naive <- ddply(df.pfizer_germ_mutation_naive, c("distance"), summarise,
                                       mutationNumber = sum(distance),
                                       meanSynFR = sum(synonyms_FR)/mutationNumber,
                                       meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                       meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                       meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.pfizer_germ_mutation_naive <- subset(df.pfizer_germ_mutation_naive,select = -c(mutationNumber))
df.pfizer_germ_mutation_naive<- melt(df.pfizer_germ_mutation_naive, id.vars = c('distance'), variable.name='muteType')
df.pfizer_germ_mutation_naive$isotype <- "naive IGHM"

# IGAM
df.pfizer_germ_mutation_IGHM<- df.pfizer_mut_germ[df.pfizer_mut_germ$isotype=="IGHM" ,]
df.pfizer_germ_mutation_IGHM <- ddply(df.pfizer_germ_mutation_IGHM, c("distance"), summarise,
                                      mutationNumber = sum(distance),
                                      meanSynFR = sum(synonyms_FR)/mutationNumber,
                                      meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                      meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                      meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.pfizer_germ_mutation_IGHM <- subset(df.pfizer_germ_mutation_IGHM,select = -c(mutationNumber))
df.pfizer_germ_mutation_IGHM<- melt(df.pfizer_germ_mutation_IGHM, id.vars = c('distance'), variable.name='muteType')
df.pfizer_germ_mutation_IGHM$isotype <- "IGHM"

# IGHA
df.pfizer_germ_mutation_IGHA<- df.pfizer_mut_germ[df.pfizer_mut_germ$isotype=="IGHM" ,]
df.pfizer_germ_mutation_IGHA <- ddply(df.pfizer_germ_mutation_IGHA, c("distance"), summarise,
                                      mutationNumber = sum(distance),
                                      meanSynFR = sum(synonyms_FR)/mutationNumber,
                                      meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                      meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                      meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.pfizer_germ_mutation_IGHA <- subset(df.pfizer_germ_mutation_IGHA,select = -c(mutationNumber))
df.pfizer_germ_mutation_IGHA<- melt(df.pfizer_germ_mutation_IGHA, id.vars = c('distance'), variable.name='muteType')
df.pfizer_germ_mutation_IGHA$isotype <- "IGHA"

# IGHG
df.pfizer_germ_mutation_IGHG<- df.pfizer_mut_germ[df.pfizer_mut_germ$isotype=="IGHM" ,]
df.pfizer_germ_mutation_IGHG <- ddply(df.pfizer_germ_mutation_IGHG, c("distance"), summarise,
                                      mutationNumber = sum(distance),
                                      meanSynFR = sum(synonyms_FR)/mutationNumber,
                                      meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                      meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                      meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.pfizer_germ_mutation_IGHG <- subset(df.pfizer_germ_mutation_IGHG,select = -c(mutationNumber))
df.pfizer_germ_mutation_IGHG<- melt(df.pfizer_germ_mutation_IGHG, id.vars = c('distance'), variable.name='muteType')
df.pfizer_germ_mutation_IGHG$isotype <- "IGHG"

df.pfizer_germ_mutation <- rbind(df.pfizer_germ_mutation_naive,
                                 df.pfizer_germ_mutation_IGHM,
                                 df.pfizer_germ_mutation_IGHA,
                                 df.pfizer_germ_mutation_IGHG)

df.pfizer_germ_mutation$type <- "Pfizer" 

figure6a <- ggplot(data=df.pfizer_germ_mutation, aes(x=distance, y=value, fill=muteType))  +
  facet_grid(isotype ~ type ,scale = "free_y") + 
  geom_bar(colour="black", stat="identity",size=.3) +   # Thinner lines
  scale_fill_brewer(palette="Set1")+
  ylab("") + # Set axis labels
  scale_y_continuous(breaks = c(0,0.5,1)) +
  scale_x_continuous(breaks = seq(from = 1,to = 20,by = 1),limits = c(0,21))+
  xlab("") +
  #scale_x_continuous(breaks = c(0,10,20,30) +
  ggtitle("Mutation distribution in each depth") +  # Set title
  background_grid(major = 'y', minor = "y")+
  scale_colour_discrete(name  ="Mutation type",
                        breaks=c("meanSynFR", "meanNonSynFR","meanSynCDR", "meanNonSynCDR"),
                        labels=c("Syn FR", "NonSyn FR","Syn CDR", "NonSyn CDR"))+
  theme_linedraw()+
  theme(  panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank()
       )



##### `~~~~~~ pfizer combine isotype data ~~~~~ #########
df.pfizer_mut_germ<- read.csv('pfizer/PFIZER_Combine_isotype_RS.csv', header = T, stringsAsFactors = F)

## naive
df.pfizer_germ_mutation_naive<- df.pfizer_mut_germ[df.pfizer_mut_germ$isotype=="naive_IGHM",]
df.pfizer_germ_mutation_naive <- ddply(df.pfizer_germ_mutation_naive, c("distance"), summarise,
                                       mutationNumber = sum(distance),
                                       meanSynFR = sum(synonyms_FR)/mutationNumber,
                                       meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                       meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                       meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.pfizer_germ_mutation_naive <- subset(df.pfizer_germ_mutation_naive,select = -c(mutationNumber))
df.pfizer_germ_mutation_naive<- melt(df.pfizer_germ_mutation_naive, id.vars = c('distance'), variable.name='muteType')
df.pfizer_germ_mutation_naive$isotype <- "naive IGHM"

# IGAM
df.pfizer_germ_mutation_IGHM<- df.pfizer_mut_germ[df.pfizer_mut_germ$isotype=="IGHM" ,]
df.pfizer_germ_mutation_IGHM <- ddply(df.pfizer_germ_mutation_IGHM, c("distance"), summarise,
                                      mutationNumber = sum(distance),
                                      meanSynFR = sum(synonyms_FR)/mutationNumber,
                                      meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                      meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                      meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.pfizer_germ_mutation_IGHM <- subset(df.pfizer_germ_mutation_IGHM,select = -c(mutationNumber))
df.pfizer_germ_mutation_IGHM<- melt(df.pfizer_germ_mutation_IGHM, id.vars = c('distance'), variable.name='muteType')
df.pfizer_germ_mutation_IGHM$isotype <- "IGHM"

# IGHA
df.pfizer_germ_mutation_IGHA<- df.pfizer_mut_germ[df.pfizer_mut_germ$isotype=="IGHM" ,]
df.pfizer_germ_mutation_IGHA <- ddply(df.pfizer_germ_mutation_IGHA, c("distance"), summarise,
                                      mutationNumber = sum(distance),
                                      meanSynFR = sum(synonyms_FR)/mutationNumber,
                                      meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                      meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                      meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.pfizer_germ_mutation_IGHA <- subset(df.pfizer_germ_mutation_IGHA,select = -c(mutationNumber))
df.pfizer_germ_mutation_IGHA<- melt(df.pfizer_germ_mutation_IGHA, id.vars = c('distance'), variable.name='muteType')
df.pfizer_germ_mutation_IGHA$isotype <- "IGHA"

# IGHG
df.pfizer_germ_mutation_IGHG<- df.pfizer_mut_germ[df.pfizer_mut_germ$isotype=="IGHM" ,]
df.pfizer_germ_mutation_IGHG <- ddply(df.pfizer_germ_mutation_IGHG, c("distance"), summarise,
                                      mutationNumber = sum(distance),
                                      meanSynFR = sum(synonyms_FR)/mutationNumber,
                                      meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                      meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                      meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.pfizer_germ_mutation_IGHG <- subset(df.pfizer_germ_mutation_IGHG,select = -c(mutationNumber))
df.pfizer_germ_mutation_IGHG<- melt(df.pfizer_germ_mutation_IGHG, id.vars = c('distance'), variable.name='muteType')
df.pfizer_germ_mutation_IGHG$isotype <- "IGHG"

df.pfizer_germ_mutation <- rbind(df.pfizer_germ_mutation_naive,
                                 df.pfizer_germ_mutation_IGHM,
                                 df.pfizer_germ_mutation_IGHA,
                                 df.pfizer_germ_mutation_IGHG)

df.pfizer_germ_mutation$type <- "Pfizer" 

figure6b <- ggplot(data=df.pfizer_germ_mutation, aes(x=distance, y=value, fill=muteType))  +
  facet_grid(isotype ~ type ,scale = "free_y") + 
  geom_bar(colour="black", stat="identity",size=.3) +   # Thinner lines
  scale_fill_brewer(palette="Set1")+
  ylab("") + # Set axis labels
  scale_y_continuous(breaks = c(0,0.5,1)) +
  scale_x_continuous(breaks = seq(from = 1,to = 18,by = 1),limits = c(0,19))+
  xlab("") +
  #scale_x_continuous(breaks = c(0,10,20,30) +
  ggtitle("Combined data - Mutation distribution in each depth") +  # Set title
  background_grid(major = 'y', minor = "y")+
  scale_colour_discrete(name  ="Mutation type",
                        breaks=c("meanSynFR", "meanNonSynFR","meanSynCDR", "meanNonSynCDR"),
                        labels=c("Syn FR", "NonSyn FR","Syn CDR", "NonSyn CDR"))+
  theme_linedraw()+
  theme(  panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank()
  )

#############################################################
###################3FLU 

##### `~~~~~~ flu all data ~~~~~ #########
df.flu_mut_germ<- read.csv('flu/FLU_RS.csv', header = T, stringsAsFactors = F)

## naive
df.flu_germ_mutation_IGHD<-df.flu_mut_germ[df.flu_mut_germ$isotype=="IGHD" ,]
df.flu_germ_mutation_IGHD <- ddply(df.flu_germ_mutation_IGHD, c("distance"), summarise,
                                   mutationNumber = sum(distance),
                                   meanSynFR = sum(synonyms_FR)/mutationNumber,
                                   meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                   meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                   meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.flu_germ_mutation_IGHD <- subset(df.flu_germ_mutation_IGHD,select = -c(mutationNumber))
df.flu_germ_mutation_IGHD<- melt(df.flu_germ_mutation_IGHD, id.vars = c('distance'), variable.name='muteType')
df.flu_germ_mutation_IGHD$isotype <- "IGHD"

# IGAM
df.flu_germ_mutation_IGHM<- df.flu_mut_germ[df.flu_mut_germ$isotype=="IGHM" ,]
df.flu_germ_mutation_IGHM <- ddply(df.flu_germ_mutation_IGHM, c("distance"), summarise,
                                   mutationNumber = sum(distance),
                                   meanSynFR = sum(synonyms_FR)/mutationNumber,
                                   meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                   meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                   meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.flu_germ_mutation_IGHM <- subset(df.flu_germ_mutation_IGHM,select = -c(mutationNumber))
df.flu_germ_mutation_IGHM<- melt(df.flu_germ_mutation_IGHM, id.vars = c('distance'), variable.name='muteType')
df.flu_germ_mutation_IGHM$isotype <- "IGHM"

# IGHA
df.flu_germ_mutation_IGHA<- df.flu_mut_germ[df.flu_mut_germ$isotype=="IGHM" ,]
df.flu_germ_mutation_IGHA <- ddply(df.flu_germ_mutation_IGHA, c("distance"), summarise,
                                   mutationNumber = sum(distance),
                                   meanSynFR = sum(synonyms_FR)/mutationNumber,
                                   meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                   meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                   meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.flu_germ_mutation_IGHA <- subset(df.flu_germ_mutation_IGHA,select = -c(mutationNumber))
df.flu_germ_mutation_IGHA<- melt(df.flu_germ_mutation_IGHA, id.vars = c('distance'), variable.name='muteType')
df.flu_germ_mutation_IGHA$isotype <- "IGHA"

# IGHG
df.flu_germ_mutation_IGHG<- rbind(df.flu_mut_germ[df.flu_mut_germ$isotype=="IGHG-1",],df.flu_mut_germ[df.flu_mut_germ$isotype=="IGHG-2",])
df.flu_germ_mutation_IGHG <- ddply(df.flu_germ_mutation_IGHG, c("distance"), summarise,
                                   mutationNumber = sum(distance),
                                   meanSynFR = sum(synonyms_FR)/mutationNumber,
                                   meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                   meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                   meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.flu_germ_mutation_IGHG <- subset(df.flu_germ_mutation_IGHG,select = -c(mutationNumber))
df.flu_germ_mutation_IGHG<- melt(df.flu_germ_mutation_IGHG, id.vars = c('distance'), variable.name='muteType')
df.flu_germ_mutation_IGHG$isotype <- "IGHG"

df.flu_germ_mutation <- rbind(df.flu_germ_mutation_IGHD,
                              df.flu_germ_mutation_IGHM,
                              df.flu_germ_mutation_IGHA,
                              df.flu_germ_mutation_IGHG)

df.flu_germ_mutation$type <- "Flu" 

figure6c <- ggplot(data=df.flu_germ_mutation, aes(x=distance, y=value, fill=muteType))  +
  facet_grid(isotype ~ type ,scale = "free_y") + 
  geom_bar(colour="black", stat="identity",size=.3) +   # Thinner lines
  scale_fill_brewer(palette="Set1")+
  ylab("") + # Set axis labels
  scale_y_continuous(breaks = c(0,0.5,1)) +
  xlab("") +
  #scale_x_continuous(breaks = c(0,10,20,30,40)) +
  scale_x_continuous(breaks = seq(from = 1,to = 20,by = 1),limits = c(0,21))+
  ggtitle("Mutation distribution in each depth") +  # Set title
  background_grid(major = 'y', minor = "y")+
  scale_colour_discrete(name  ="Mutation type",
                        breaks=c("meanSynFR", "meanNonSynFR","meanSynCDR", "meanNonSynCDR"),
                        labels=c("Syn FR", "NonSyn FR","Syn CDR", "NonSyn CDR"))+ 
  theme_linedraw()+
  theme(  panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank()
  )



##### `~~~~~~ flu combine isotype data ~~~~~ #########

df.flu_mut_germ<- read.csv('flu/FLU_Combine_isotype_RS.csv', header = T, stringsAsFactors = F)

## naive
df.flu_germ_mutation_IGHD<-df.flu_mut_germ[df.flu_mut_germ$isotype=="IGHD" ,]
df.flu_germ_mutation_IGHD <- ddply(df.flu_germ_mutation_IGHD, c("distance"), summarise,
                                   mutationNumber = sum(distance),
                                   meanSynFR = sum(synonyms_FR)/mutationNumber,
                                   meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                   meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                   meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.flu_germ_mutation_IGHD <- subset(df.flu_germ_mutation_IGHD,select = -c(mutationNumber))
df.flu_germ_mutation_IGHD<- melt(df.flu_germ_mutation_IGHD, id.vars = c('distance'), variable.name='muteType')
df.flu_germ_mutation_IGHD$isotype <- "IGHD"

# IGAM
df.flu_germ_mutation_IGHM<- df.flu_mut_germ[df.flu_mut_germ$isotype=="IGHM" ,]
df.flu_germ_mutation_IGHM <- ddply(df.flu_germ_mutation_IGHM, c("distance"), summarise,
                                   mutationNumber = sum(distance),
                                   meanSynFR = sum(synonyms_FR)/mutationNumber,
                                   meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                   meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                   meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.flu_germ_mutation_IGHM <- subset(df.flu_germ_mutation_IGHM,select = -c(mutationNumber))
df.flu_germ_mutation_IGHM<- melt(df.flu_germ_mutation_IGHM, id.vars = c('distance'), variable.name='muteType')
df.flu_germ_mutation_IGHM$isotype <- "IGHM"

# IGHA
df.flu_germ_mutation_IGHA<- df.flu_mut_germ[df.flu_mut_germ$isotype=="IGHM" ,]
df.flu_germ_mutation_IGHA <- ddply(df.flu_germ_mutation_IGHA, c("distance"), summarise,
                                   mutationNumber = sum(distance),
                                   meanSynFR = sum(synonyms_FR)/mutationNumber,
                                   meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                   meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                   meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.flu_germ_mutation_IGHA <- subset(df.flu_germ_mutation_IGHA,select = -c(mutationNumber))
df.flu_germ_mutation_IGHA<- melt(df.flu_germ_mutation_IGHA, id.vars = c('distance'), variable.name='muteType')
df.flu_germ_mutation_IGHA$isotype <- "IGHA"

# IGHG
df.flu_germ_mutation_IGHG<- rbind(df.flu_mut_germ[df.flu_mut_germ$isotype=="IGHG-1",],df.flu_mut_germ[df.flu_mut_germ$isotype=="IGHG-2",])
df.flu_germ_mutation_IGHG <- ddply(df.flu_germ_mutation_IGHG, c("distance"), summarise,
                                   mutationNumber = sum(distance),
                                   meanSynFR = sum(synonyms_FR)/mutationNumber,
                                   meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                   meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                   meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)


df.flu_germ_mutation_IGHG <- subset(df.flu_germ_mutation_IGHG,select = -c(mutationNumber))
df.flu_germ_mutation_IGHG<- melt(df.flu_germ_mutation_IGHG, id.vars = c('distance'), variable.name='muteType')
df.flu_germ_mutation_IGHG$isotype <- "IGHG"

df.flu_germ_mutation <- rbind(df.flu_germ_mutation_IGHD,
                              df.flu_germ_mutation_IGHM,
                              df.flu_germ_mutation_IGHA,
                              df.flu_germ_mutation_IGHG)

df.flu_germ_mutation$type <- "Flu" 

figure6d <- ggplot(data=df.flu_germ_mutation, aes(x=distance, y=value, fill=muteType))  +
  facet_grid(isotype ~ type ,scale = "free_y") + 
  geom_bar(colour="black", stat="identity",size=.3) +   # Thinner lines
  scale_fill_brewer(palette="Set1")+
  ylab("") + # Set axis labels
  scale_y_continuous(breaks = c(0,0.5,1)) +
  xlab("") +
  #scale_x_continuous(breaks = c(0,10,20,30,40)) +
  scale_x_continuous(breaks = seq(from = 1,to = 20,by = 1),limits = c(0,21))+
  ggtitle("Combined - Mutation distribution in each depth") +  # Set title
  background_grid(major = 'y', minor = "y")+
  scale_colour_discrete(name  ="Mutation type",
                        breaks=c("meanSynFR", "meanNonSynFR","meanSynCDR", "meanNonSynCDR"),
                        labels=c("Syn FR", "NonSyn FR","Syn CDR", "NonSyn CDR"))+
  theme_linedraw()+
  theme(  panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank()
  )

grid.arrange( figure6a ,figure6b ,figure6c,figure6d,ncol =1)  










