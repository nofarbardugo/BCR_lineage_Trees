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


###### fig1 #######

# ~~~~~~pfizer details~~~~~~

df.pfizer_node_diveded_son<- read.csv('pfizer/PFIZER_Node_son_divide.csv', header = T, stringsAsFactors = F)


df.pfizer_node_diveded_son$bins2Depth <- df.pfizer_node_diveded_son$Depth
df.pfizer_mutation_edges[df.pfizer_mutation_edges$bins2Depth%%2==1,"bins2Depth"] <- df.pfizer_mutation_edges[df.pfizer_mutation_edges$bins2Depth%%2==1,"bins2Depth"] +1 

# get the ratio beteewn each mutation type in pfizer
df.pfizer_mutation_edges <- ddply(df.pfizer_node_diveded_son, c("Isotype.Name","bins2Depth"), summarise,
                                  Ratio = mean(son.differnt/Son.Number), # for each ratio calc avg
                                  SD = sd(son.differnt/Son.Number), # for each ratio calc SD
                                  SE = SD/(sqrt(length(bins2Depth)))                         
)

# remove IGAH - not relvent becase do not have IS 
df.pfizer_mutation_edges <- df.pfizer_mutation_edges[df.pfizer_mutation_edges$Isotype.Name!="IGHA",] 
df.pfizer_mutation_edges <- df.pfizer_mutation_edges[df.pfizer_mutation_edges$Ratio> 0.001,] 
#df.pfizer_mutation_edges <- df.pfizer_mutation_edges[-36, ]

#df.pfizer_mutation_edges$bins2Depth <- df.pfizer_mutation_edges$Depth
#df.pfizer_mutation_edges[df.pfizer_mutation_edges$bins2Depth%%2==1,"bins2Depth"] <- df.pfizer_mutation_edges[df.pfizer_mutation_edges$bins2Depth%%2==1,"bins2Depth"] +1  
#df.pfizer_mutation_edges_bin2 <- ddply(df.pfizer_mutation_edges, c("Isotype.Name","bins2Depth"), summarise,
#                                  sonsBin = sum(sons),
#                                  son.differntBin = sum(son.differnt),
#                                  Ratio = son.differntBin/sonsBin                                 
#)
  
df.pfizer_mutation_edges <- df.pfizer_mutation_edges[df.pfizer_mutation_edges$bins2Depth <36,]
figure1a1 <-ggplot(df.pfizer_mutation_edges, aes(x=bins2Depth, y=Ratio, colour=Isotype.Name)) +
  geom_line() + geom_point()  +
  geom_errorbar(aes(ymin=Ratio-SE, ymax=Ratio+SE), width=.4) +
  scale_colour_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","naive_IGHM" = "green")) +
  scale_x_continuous(breaks=seq(0, 40, 2))+
  xlab("Depth: S+NS") +
  #xlim(c(0,40)) +
  ylim(c(0,0.5))+
  ylab("Fraction of switches") +
  ggtitle("Pfizer - Switch Fraction") +
  theme(plot.title = element_text(size = rel(1.5)))+ theme_gray(17)


# ~~~~~~flu details~~~~~~
df.flu_node_diveded_son<- read.csv('flu/FLU_Node_son_divide.csv', header = T, stringsAsFactors = F)

# merge IGHG1 with IGHG2
df.flu_node_diveded_son[df.flu_node_diveded_son$Isotype.Name=="IGHG-1","Isotype.Name"]<-"IGHG"
df.flu_node_diveded_son[df.flu_node_diveded_son$Isotype.Name=="IGHG-2","Isotype.Name"]<-"IGHG"

# get bins of 2 
df.flu_node_diveded_son$bins2Depth <- df.flu_node_diveded_son$Depth
df.flu_node_diveded_son[df.flu_node_diveded_son$bins2Depth%%2==1,"bins2Depth"] <- df.flu_node_diveded_son[df.flu_node_diveded_son$bins2Depth%%2==1,"bins2Depth"] +1 

# get the ratio beteewn each mutation type in pfizer
df.flu_mutation_edges <- ddply(df.flu_node_diveded_son, c("Isotype.Name","bins2Depth"), summarise,
                                  Ratio = mean(son.differnt/Son.Number), # for each ratio calc avg
                                  SD = sd(son.differnt/Son.Number), # for each ratio calc SD
                                  SE = SD/(sqrt(length(bins2Depth)))                         
)

# get the ratio beteewn each mutation type in flu
#df.flu_mutation_edges <- ddply(df.flu_node_diveded_son, c("Isotype.Name","Depth"), summarise,
#                               sons = sum(Son.Number),
#                               son.differnt = sum(son.differnt),
#                               Ratio = son.differnt/sons
                               
#)

# edit data for geraph
IGAD <- df.flu_mutation_edges[df.flu_mutation_edges$Isotype.Name=="IGHD",]
IGAD  <-  IGAD[IGAD$bins2Depth <28,]  

IGAM <- df.flu_mutation_edges[df.flu_mutation_edges$Isotype.Name=="IGHM",]
IGAM  <-  IGAM[IGAM$bins2Depth <37,]  

IGAG <-df.flu_mutation_edges[df.flu_mutation_edges$Isotype.Name=="IGHG",]
IGAG  <-  IGAG[IGAG$bins2Depth <37,]  

df.flu_mutation_edges_B <- rbind(IGAD,IGAG,IGAM)
df.flu_mutation_edges_B <- df.flu_mutation_edges_B[df.flu_mutation_edges_B$bins2Depth <36,]
figure1a2 <-ggplot(df.flu_mutation_edges_B, aes(x=bins2Depth, y=Ratio, colour=Isotype.Name)) +
  geom_line() + geom_point()  +
  geom_errorbar(aes(ymin=Ratio-SE, ymax=Ratio+SE), width=.4)+
  scale_x_continuous(breaks=seq(0, 40, 2))+
  #xlim(c(0,40)) +
  scale_colour_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","IGHD" = "#FF3399","IGHE" = "purple")) +
  xlab("Depth: S+NS") + 
  ylab("Fraction of switches") +
  ggtitle("Flu - Switch Fraction") +
  theme(plot.title = element_text(size = rel(1.5)))+ theme_gray(17) 

figure1a <-  grid.arrange( figure1a1,figure1a2,ncol =1)  


##### pre vacine flu 

df.flu_node_diveded_son<- read.csv('flu/FLU_PreVacine_Node_son_divide.csv', header = T, stringsAsFactors = F)


# merge IGHG1 with IGHG2
df.flu_node_diveded_son[df.flu_node_diveded_son$Isotype.Name=="IGHG-1","Isotype.Name"]<-"IGHG"
df.flu_node_diveded_son[df.flu_node_diveded_son$Isotype.Name=="IGHG-2","Isotype.Name"]<-"IGHG"

# get bins of 2 
df.flu_node_diveded_son$bins2Depth <- df.flu_node_diveded_son$Depth
df.flu_node_diveded_son[df.flu_node_diveded_son$bins2Depth%%2==1,"bins2Depth"] <- df.flu_node_diveded_son[df.flu_node_diveded_son$bins2Depth%%2==1,"bins2Depth"] +1 
d <-df.flu_node_diveded_son[df.flu_node_diveded_son$Depth==21 & df.flu_node_diveded_son$Isotype.Name =="IGHD",]
# get the ratio beteewn each mutation type in pfizer
df.flu_mutation_edges <- ddply(df.flu_node_diveded_son, c("Isotype.Name","bins2Depth"), summarise,
                               Ratio = mean(son.differnt/Son.Number), # for each ratio calc avg
                               SD = sd(son.differnt/Son.Number), # for each ratio calc SD
                               SE = SD/(sqrt(length(Depth)))                         
)

# get the ratio beteewn each mutation type in flu
#df.flu_mutation_edges <- ddply(df.flu_node_diveded_son, c("Isotype.Name","Depth"), summarise,
#                               sons = sum(Son.Number),
#                               son.differnt = sum(son.differnt),
#                               Ratio = son.differnt/sons

#)

# edit data for geraph
IGAD <- df.flu_mutation_edges[df.flu_mutation_edges$Isotype.Name=="IGHD",]
IGAD  <-  IGAD[IGAD$bins2Depth <28,]  

IGAM <- df.flu_mutation_edges[df.flu_mutation_edges$Isotype.Name=="IGHM",]
IGAM  <-  IGAM[IGAM$bins2Depth <37,]  

IGAG <-df.flu_mutation_edges[df.flu_mutation_edges$Isotype.Name=="IGHG",]
IGAG  <-  IGAG[IGAG$bins2Depth <37,]  

df.flu_mutation_edges_B <- rbind(IGAD,IGAG,IGAM)
df.flu_mutation_edges_B <- df.flu_mutation_edges_B[df.flu_mutation_edges_B$bins2Depth <36,]
figure1a2 <-ggplot(df.flu_mutation_edges_B, aes(x=bins2Depth, y=Ratio, colour=Isotype.Name)) +
  geom_line() + geom_point()  +
  geom_errorbar(aes(ymin=Ratio-SE, ymax=Ratio+SE), width=.4)+
  scale_x_continuous(breaks=seq(0, 40, 2))+
  #xlim(c(0,40)) +
  scale_colour_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","IGHD" = "#FF3399","IGHE" = "purple")) +
  xlab("Depth: S+NS") + 
  ylab("Fraction of switches") +
  ggtitle("Flu - pre vacine - Switch Fraction") +
  theme(plot.title = element_text(size = rel(1.5)))+ theme_gray(17) 


##### post vacine flu - as a function of a time point
df.flu_postVacineFrac<- read.csv('flu/FLU_isotype_Fraction_With_Point_Time.csv', header = T, stringsAsFactors = F)

df.Frac <- ddply(df.flu_postVacineFrac, c("Isotype.Name","Child.Name.","depth","Time.Point"), summarise,
                               number = length(Isotype.Name)                 
)

df.FracAll <- ddply(df.flu_postVacineFrac, c("Isotype.Name","depth","Time.Point"), summarise,
                 number = length(Isotype.Name)                 
)

df.Frac$totalSons <- 0
##  enter the fracion for each row in specific time point 
for(i in 1: nrow(df.Frac))
{
  df.Frac[i,"totalSons"] <- df.FracAll[df.FracAll$Isotype.Name==df.Frac[i,"Isotype.Name"] &
                                        df.FracAll$depth ==df.Frac[i,"depth"]&
                                        df.FracAll$Time.Point ==df.Frac[i,"Time.Point"],"number"]
}

df.Frac$fraction <- df.Frac$number/df.Frac$totalSon
df.Frac <- df.Frac[-grep('IGHA', df.Frac$Isotype.Name),]

ggplot(df.Frac, aes(x=depth, y=fraction, colour=Child.Name.)) +
  facet_grid(Isotype.Name ~ Time.Point) +
  geom_line() + geom_point()  +
  #scale_colour_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","naive_IGHM" = "green")) +
  xlab("Depth: S+NS") +
  xlim(c(0,40)) +
  ylim(c(0,0.5))+
  ylab("Fraction of switches") +
  ggtitle("Pfizer - Switch Fraction") +
  theme(plot.title = element_text(size = rel(1.5)))+ theme_gray(17)

####################~~~~~~~~~~~~~~~~~~~ fig2 ~~~~~~~~~~~~~~~~~~`#####################

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



# pfizer
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

figure2a <-ggplot(pfizer_ds_event, aes(x=depth,colour = isotype)) + geom_density() +
 xlab("Mutation number") + 
  xlim(c(0,30))+
  ylab("Density - Number of nodes") +
  ggtitle("Histogram of mutations per isotype") +
  scale_fill_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","Naive_IGHM" = "green")) +
  theme(text = element_text(lineheight=.8, face="bold",size = 16)) + theme_gray(17) 
  theme(plot.title = element_text(size = rel(1.5)))+ theme_gray(17) 


figure2a <-ggplot(pfizer_ds_event, aes(x=depth,fill = isotype)) + geom_histogram(aes(y=..density..),alpha=.7,adjust = 0.1,position ="dodge",binwidth = 1 ) +
  xlab("Mutation number") + 
  xlim(c(0,30))+
  ylab("Density - Number of nodes") +
  ggtitle("Histogram of mutations per isotype") +
  scale_fill_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","Naive_IGHM" = "green")) +
  theme(text = element_text(lineheight=.8, face="bold",size = 16)) + theme_gray(17) 
theme(plot.title = element_text(size = rel(1.5)))+ theme_gray(17) 

###########################################################################
pfizer_isotypeSize <- ddply(pfizer_ds_event, c("isotype"), summarise,
                                  repeatTimes = length(depth)
                                  
)

pfizer_ds_event_to_point <- ddply(pfizer_ds_event, c("isotype","depth"), summarise,
                            repeatTimes = length(depth)

)

pfizer_ds_event_to_point$frequency <- 
  pfizer_ds_event_to_point[pfizer_ds_event_to_point$isotype==pfizer_isotypeSize$isotype,]

figure21 <-ggplot(pfizer_ds_event_to_point, aes(x=depth, y=repeatTimes, colour=isotype)) +
  geom_line() + geom_point()  +
  xlab("Mutation number") + 
  xlim(c(0,30))+
  ylab("Density - Number of nodes") +
  ggtitle("Histogram of mutations per isotype") +
  scale_fill_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","Naive_IGHM" = "green")) +
  theme(text = element_text(lineheight=.8, face="bold",size = 16)) + theme_gray(17) 
theme(plot.title = element_text(size = rel(1.5)))+ theme_gray(17) 

ggplot(df.pfizer_mutation_edges,  +
 
geom_errorbar(aes(ymin=Ratio-SE, ymax=Ratio+SE), width=.4) +

####################################################################################################

# flu
flu_dfIGHG1 <- read.csv('flu/FLU_event_IGHG1.csv', header = T, stringsAsFactors = F)
flu_dfIGHG1$isotype <- "IGHG"

flu_dfIGHG2<- read.csv('flu/FLU_event_IGHG2.csv', header = T, stringsAsFactors = F)
flu_dfIGHG2$isotype <- "IGHG"

flu_dfIGHA <- read.csv('flu/FLU_event_IGHA.csv', header = T, stringsAsFactors = F)
flu_dfIGHA$isotype <- "IGHA"


#flu_dfIGHE <- read.csv('flu/FLU_event_IGHE.csv', header = T, stringsAsFactors = F)
#flu_dfIGHE$isotype <- "IGHE"

flu_dfIGHM <- read.csv('flu/FLU_event_IGHM.csv', header = T, stringsAsFactors = F)
flu_dfIGHM$isotype <- "IGHM"


flu_dfIGHD <- read.csv('flu/FLU_event_naive_IGHD.csv', header = T, stringsAsFactors = F)
flu_dfIGHD$isotype <- "IGHD"

# merge all into one dataFream
flu_ds_event <- rbind(flu_dfIGHG1,flu_dfIGHG2,flu_dfIGHA,flu_dfIGHM,flu_dfIGHD)#,flu_dfIGHE,)

figure2b <-ggplot(flu_ds_event, aes(x=depth,fill = isotype)) + geom_histogram(aes(y=..density..),alpha=.7,adjust = 0.1,position ="dodge",binwidth = 1 ) +
  xlab("Mutation number") + 
  xlim(c(0,30)) +
  ylab("Density - Number of nodes") +
  ggtitle("Histogram of mutations per isotype") +
  scale_fill_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHG1" = "#CC6666","IGHG2" = "brown","IGHM" ="blue","IGHD" = "#FF3399","IGHE" = "purple")) +
  theme(text = element_text(lineheight=.8, face="bold",size = 16)) + theme_gray(17) 

figure2 <- grid.arrange(figure2a, figure2b , ncol=2) 

############### ~~~~~~~~~~~fig 3 ~~~~~~~~~~~~~~~~ ################################

#~~~~~ fig3A

# pfizer details
df.pfizer_depthVSchildren<- read.csv('pfizer/PFIZER_mutation_num_per_isotype_with_certain_children_numberWithInternal.csv', header = T, stringsAsFactors = F)
df.pfizer_depthVSchildren$mutPerChild <- df.pfizer_depthVSchildren$depth/df.pfizer_depthVSchildren$Children.Number
df.pfizer_depthVSchildren$type <- "Pfizer"

# flu details
df.flu_depthVSchildren<- read.csv('flu/FLU_mutation_num_per_isotype_with_certain_children_numberWithInternal.csv', header = T, stringsAsFactors = F)
df.flu_depthVSchildren$mutPerChild <- df.flu_depthVSchildren$depth/df.flu_depthVSchildren$Children.Number
df.flu_depthVSchildren <- df.flu_depthVSchildren[df.flu_depthVSchildren$Children.Number >0 ,] 
df.flu_depthVSchildren$type <- "Flu"

# merge both data set
df.combine_depthVSchildren <- rbind(df.pfizer_depthVSchildren,df.flu_depthVSchildren)
# remove rows with ziro child (not supposed to be but append in igG2)
df.combine_depthVSchildren <- df.combine_depthVSchildren[df.combine_depthVSchildren$Children.Number >0 ,] 

# get avg of the ratio between mutation num to children num
#df.combine_depthVSchildrenwitoutZiroLength <- df.combine_depthVSchildren[df.combine_depthVSchildren$depth >0 ,] 
df.combine_depthVSchildrenMueserd <- ddply(df.combine_depthVSchildren, c("type","Isotype.Name"), summarise,
                                           isotypeNumber    = length(Isotype.Name),
                                           childrenNumber = sum(Children.Number),
                                           mutationNumber = sum(depth),
                                           meanDepthVSchildren = mutationNumber/childrenNumber,
                                           #meanDepthVSchildren = mean(mutPerChild),
                                           meanSynFR = sum(synonyms_FR)/mutationNumber,
                                           meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                                           meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                                           meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber,
                                           meanSynonymus = sum(Synonymus)/mutationNumber,
                                           meanNonSynonymus = sum(nonSynonymus)/mutationNumber
)

# get the ratio beteewn each mutation type in flu
df.flu_mutation <- ddply(df.flu_depthVSchildren, c("Isotype.Name"), summarise,
                         mutationNumber = sum(depth),
                         meanSynFR = sum(synonyms_FR)/mutationNumber,
                         meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                         meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                         meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

df.flu_mutation <- subset(df.flu_mutation,select = -c(mutationNumber))
df.flu_mutation<- melt(df.flu_mutation, id.vars = 'Isotype.Name', variable.name='muteType')
df.flu_mutation$type <- "Flu"

# get the ratio beteewn each mutation type in pfizer
df.pfizer_mutation <- ddply(df.pfizer_depthVSchildren, c("Isotype.Name"), summarise,
                            mutationNumber = sum(depth),
                            meanSynFR = sum(synonyms_FR)/mutationNumber,
                            meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                            meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                            meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

df.pfizer_mutation <- subset(df.pfizer_mutation,select = -c(mutationNumber))
df.pfizer_mutation<- melt(df.pfizer_mutation, id.vars = 'Isotype.Name', variable.name='muteType')
df.pfizer_mutation$type <- "Pfizer"

# combine flu and pfizer
df.combine_mutation <- rbind(df.pfizer_mutation,df.flu_mutation)



figure3A <- ggplot(data=df.combine_depthVSchildrenMueserd, aes(x=Isotype.Name, y=meanDepthVSchildren, fill=Isotype.Name)) + 
  facet_wrap(~type, ncol=2, scale="free_x") + #facet_grid(type ~ . ,scale = "free_y") +
  geom_bar(colour="black", stat="identity",position=position_dodge(),size=.3) +   # Thinner lines
  scale_fill_manual(values=c("IGHA" = "orange","IGHG-1" = "#CC6666","IGHG-2" = "brown",
                             "IGHM" ="blue","IGHD" = "#FF3399","IGHE" = "purple",
                             "naive_IGHM" = "green","IGHG" = "red")) +
  #   scale_fill_hue(name="Isotype name") +      # Set legend title
  ylab("Avg mutation per child") + # Set axis labels
  xlab("") + 
  ggtitle("Average mutations number per child in each isotype") +  # Set title
  background_grid(major = 'y', minor = "y") + # add thin horizontal
  theme(text = element_text(lineheight=.8,size = 14),
        legend.position="none",axis.text.x = element_text(angle = 90,hjust = 1))


#~~~~~ fig3B

figure3B1 <-ggplot(data=df.pfizer_mutation, aes(x=Isotype.Name, y=value, fill=muteType)) + 
  geom_bar(colour="black",stat="identity", size=.3) +  # Thinner lines
  xlab("") + ylab("") + # Set axis labels
  ggtitle("Pfizer") +     # Set title
  scale_fill_brewer(palette="Set1")+
  coord_flip() +# convert exes
  theme(text = element_text(lineheight=.8,size = 14),
        #axis.text.x = element_text(angle = 90,hjust = 1))+
        axis.text.x = element_text(hjust = 1))+
  scale_colour_discrete(name  ="Mutation type",
                        breaks=c("meanSynFR", "meanNonSynFR","meanSynCDR", "meanNonSynCDR"),
                        labels=c("Syn FR", "NonSyn FR","Syn CDR", "NonSyn CDR")) + coord_flip() 

figure3B2 <-ggplot(data=df.flu_mutation, aes(x=Isotype.Name, y=value, fill=muteType)) + 
  geom_bar(colour="black",stat="identity", size=.3) +  # Thinner lines
  scale_fill_brewer(palette="Set1")+
  xlab("") + ylab("") + # Set axis labels
  ggtitle("Flu") +     # Set title
  coord_flip() +# convert exes
  theme(text = element_text(lineheight=.8,size = 14),
        legend.position="none",axis.text.x = element_text(angle = 90,hjust = 1))

a <- ggdraw() +
  draw_plot(figure3A, 0, 0, 0.5, 1) +
  draw_plot(figure3B2, 0.5, 0, .2, 1) +
  draw_plot(figure3B1, 0.68, 0, .33, 1) +
  draw_plot_label(c("A", "B"),c(0, 0.55), c(0.98, 0.98), size = 15)


############### ~~~~~~~~~~~fig 4 ~~~~~~~~~~~~~~~~ ################################

#~~~~~ fig4A

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

df.flu_Igm_mutation <- df.flu_Ige_mutation[df.flu_Igm_mutation$depth<9,]
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
                                      df.flu_ige_mutationRatio_depth,
                                      df.flu_igm_mutationRatio_depth,
                                      df.flu_igd_mutationRatio_depth)

df.flu_combine_depth_mutaion$type <- "Flu" 

# combine pfizer with flu 
df.combine_depth_mutaion <- rbind(df.pfizer_combine_depth_mutaion,df.flu_combine_depth_mutaion) 

figure4a1 <- ggplot(data=df.pfizer_combine_depth_mutaion, aes(x=depth, y=value, fill=muteType))  +
  facet_grid(isotype ~ type ,scale = "free_y") + 
  geom_bar(colour="black", stat="identity",size=.3) +   # Thinner lines
  scale_fill_brewer(palette="Set1")+
  ylab("") + # Set axis labels
  scale_y_continuous(breaks = c(0,0.5,1)) +
  xlab("") +
  scale_x_continuous(breaks = c(0,10,20,30,40,50)) +
  
  ggtitle("Mutation distribution in each depth") +  # Set title
  background_grid(major = 'y', minor = "y")+
  scale_colour_discrete(name  ="Mutation type",
                        breaks=c("meanSynFR", "meanNonSynFR","meanSynCDR", "meanNonSynCDR"),
                        labels=c("Syn FR", "NonSyn FR","Syn CDR", "NonSyn CDR"))+ theme_gray(15)

figure4a2 <- ggplot(data=df.flu_combine_depth_mutaion, aes(x=depth, y=value, fill=muteType))  + 
  facet_grid(isotype ~ type ,scale = "free_y") + 
  geom_bar(colour="black", stat="identity",size=.3) +   # Thinner lines
  scale_fill_brewer(palette="Set1")+
  ylab("") + # Set axis labels
  scale_y_continuous(breaks = c(0,0.5,1)) +
  xlab("Depth") +
  scale_x_continuous(breaks = c(0,10,20,30,40,50)) +
  ggtitle("") +  # Set title
  background_grid(major = 'y', minor = "y")+
  scale_colour_discrete(name  ="Mutation type",
                        breaks=c("meanSynFR", "meanNonSynFR","meanSynCDR", "meanNonSynCDR"),
                        labels=c("Syn FR", "NonSyn FR","Syn CDR", "NonSyn CDR"))+ theme_gray(15)

figure4a <- grid_arrange_shared_legend(figure4a1 + xlim(c(0,40)), figure4a2 + xlim(c(0,40)))



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

figure4b1 <-ggplot(data=df.pfizer_mutation_edges, aes(x=Edge.Name, y=value, fill=muteType)) + 
  geom_bar(colour="black",stat="identity", size=.3) +  # Thinner lines
  xlab("") + ylab("") + # Set axis labels
  scale_fill_brewer(palette="Set1")+
  ggtitle("Pfizer") +     # Set title
  theme(text = element_text(lineheight=.8,size = 14),legend.position="none",
        axis.text.x = element_text(angle = 90,hjust = 1))+
  coord_flip() 


# flu details
df.flu_edges<- read.csv('flu/FLU_Edges_mutation_num_per_edge.csv', header = T, stringsAsFactors = F)


# get the ratio beteewn each mutation type in flu
df.flu_mutation_edges <- ddply(df.flu_edges, c("Edge.Name"), summarise,
                               mutationNumber = sum(depth),
                               meanSynFR = sum(synonyms_FR)/mutationNumber,
                               meanNonSynFR = sum(nonSynonyms_FR)/mutationNumber,
                               meanSynCDR =  sum(synonyms_CDR)/mutationNumber,
                               meanNonSynCDR  = sum(nonSynonyms_CDR)/mutationNumber
)

df.flu_mutation_edges <- subset(df.flu_mutation_edges,select = -c(mutationNumber))
df.flu_mutation_edges<- melt(df.flu_mutation_edges, id.vars = 'Edge.Name', variable.name='muteType')
df.flu_mutation_edges$type <- "Flu"

figure4b2 <-ggplot(data=df.flu_mutation_edges, aes(x=Edge.Name, y=value, fill=muteType)) + 
  geom_bar(colour="black",stat="identity", size=.3) +  # Thinner lines
  xlab("") + ylab("") + # Set axis labels
  scale_fill_brewer(palette="Set1")+
  coord_flip() +
  ggtitle("Flu") +     # Set title
  theme(text = element_text(lineheight=.8,size = 14),
        legend.position="none",axis.text.x = element_text(angle = 90,hjust = 1))


figure4b <- grid.arrange(figure4b1, figure4b2,ncol =2)







