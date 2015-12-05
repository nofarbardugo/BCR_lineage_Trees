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
############### ~~~~~~~~~~~fig 4 ~~~~~~~~~~~~~~~~ ################################

#~~~~~ fig4

# pfizer details
df.pfizer_depthVSchildren<- read.csv('pfizer/PFIZER_mutation_num_per_isotype_with_certain_children_numberWithInternal.csv', header = T, stringsAsFactors = F)
df.pfizer_depthVSchildren$mutPerChild <- df.pfizer_depthVSchildren$depth/df.pfizer_depthVSchildren$Children.Number
df.pfizer_depthVSchildren$type <- "Pfizer"

# flu details
df.flu_depthVSchildren<- read.csv('flu/FLU_mutation_num_per_isotype_with_certain_children_numberWithInternal.csv', header = T, stringsAsFactors = F)
df.flu_depthVSchildren$mutPerChild <- df.flu_depthVSchildren$depth/df.flu_depthVSchildren$Children.Number
df.flu_depthVSchildren <- df.flu_depthVSchildren[df.flu_depthVSchildren$Children.Number >0 ,] 
df.flu_depthVSchildren[df.flu_depthVSchildren$Isotype.Name=="IGHG-1" |df.flu_depthVSchildren$Isotype.Name=="IGHG-2","Isotype.Name"] <-"IGHG" 
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

figure4A <- ggplot(data=df.combine_depthVSchildrenMueserd, aes(x=Isotype.Name, y=meanDepthVSchildren, fill=Isotype.Name)) + 
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
  theme_linedraw()+
  theme(  legend.title=element_blank(),
          legend.position="none", # legend position
          legend.justification=c(1,1),
         panel.grid.minor=element_blank(),# remove grid
         panel.grid.major=element_blank()
          
  ) 



#~~~~~ fig4B

figure4B1 <-ggplot(data=df.pfizer_mutation, aes(x=Isotype.Name, y=value, fill=muteType)) + 
  geom_bar(colour="black",stat="identity", size=.3) +  # Thinner lines
  xlab("") + ylab("") + # Set axis labels
  ggtitle("Pfizer") +     # Set title
  scale_fill_brewer(palette="Set1")+
  scale_colour_discrete(name  ="Mutation type",
                        breaks=c("meanSynFR", "meanNonSynFR","meanSynCDR", "meanNonSynCDR"),
                        labels=c("Syn FR", "NonSyn FR","Syn CDR", "NonSyn CDR")) +
  theme_linedraw()+
  theme( 
         panel.grid.minor=element_blank(), # remove grid
         panel.grid.major=element_blank())



figure4B2 <-ggplot(data=df.flu_mutation, aes(x=Isotype.Name, y=value, fill=muteType)) + 
  geom_bar(colour="black",stat="identity", size=.3) +  # Thinner lines
  scale_fill_brewer(palette="Set1")+
  xlab("") + ylab("") + # Set axis labels
  ggtitle("Flu") +     # Set title
 # coord_flip() +# convert exes
  theme_linedraw()+
  theme( legend.position="none",
         panel.grid.minor=element_blank(), # remove grid
         panel.grid.major=element_blank())


grid_arrange_shared_legend(figure4B1,figure4B2)

#  theme(text = element_text(lineheight=.8,size = 14),
#        legend.position="none",axis.text.x = element_text(angle = 90,hjust = 1))

a1 <- grid.arrange(figure4B2, figure4B1,  ncol=1) 
a <- grid.arrange(figure4A, a1,  ncol=2) 
#grid.arrange(figure4B2, figure4B1,ncol =2)
a <- ggdraw() +
  draw_plot(figure4A, 0, 0, 0.5, 1) +
  draw_plot(figure4B2, 0.5, 0, .24, 1) +
  draw_plot(figure4B1, 0.74, 0, .25, 1) +
  draw_plot_label(c("A", "B"),c(0, 0.55), c(0.98, 0.98), size = 15)
