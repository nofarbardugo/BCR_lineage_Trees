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



# ~~~~~~pfizer details~~~~~~

df.pfizer_node_diveded_son<- read.csv('pfizer/PFIZER_Node_son_divide_per_isotype.csv', header = T, stringsAsFactors = F)

# get the ratio beteewn each mutation type in pfizer
df.pfizer_son_number <- ddply(df.pfizer_node_diveded_son, c("Isotype.Name","Depth"), summarise,
                              sons = length(Isotype.Name)                              
)

# get the ratio beteewn each mutation type in pfizer
df.pfizer_spec_son <- ddply(df.pfizer_node_diveded_son, c("Isotype.Name","Son.type","Depth"), summarise,
                            sons = length(Son.type)                              
)

# remove rowa with IGHA or '0' as father
df.pfizer_spec_son <- df.pfizer_spec_son[-grep('0', df.pfizer_spec_son$Isotype.Name),]
df.pfizer_spec_son <- df.pfizer_spec_son[-grep('IGHA', df.pfizer_spec_son$Isotype.Name),]

##  enter the fracion for each row 
for(i in 1: nrow(df.pfizer_spec_son))
{
  df.pfizer_spec_son[i,"totalSon"] <- df.pfizer_son_number[df.pfizer_son_number$Isotype.Name==df.pfizer_spec_son[i,"Isotype.Name"] &
                                                             df.pfizer_son_number$Depth ==df.pfizer_spec_son[i,"Depth"] ,"sons"]
}

# calc fracion (diffrentSon/total son)
df.pfizer_spec_son$fraction <- df.pfizer_spec_son$sons/df.pfizer_spec_son$totalSon

# create beans of 2:
df.pfizer_spec_son$bins2Depth <- df.pfizer_spec_son$Depth
df.pfizer_spec_son[df.pfizer_spec_son$bins2Depth%%2==1,"bins2Depth"] <- df.pfizer_spec_son[df.pfizer_spec_son$bins2Depth%%2==1,"bins2Depth"] +1 

df.pfizer_frac <- ddply(df.pfizer_spec_son, c("Isotype.Name","Son.type","bins2Depth"), summarise,
                                                                    Sons = sum(sons), #
                                                                    TotalSon =sum(totalSon)                              
                                  )                        

df.pfizer_frac<- df.pfizer_frac[df.pfizer_frac$Isotype.Name!=df.pfizer_frac$Son.type,]
df.pfizer_frac$fraction <- df.pfizer_frac$Sons/df.pfizer_frac$TotalSon
df.pfizer_frac$SE <- df.pfizer_frac$fraction* (1/sqrt(df.pfizer_frac$Sons))
  
figure2a1 <-ggplot(df.pfizer_frac, aes(x=bins2Depth, y=fraction, colour=Son.type)) +
  facet_wrap(~ Isotype.Name) +
  geom_line() + geom_point()  +
  geom_errorbar(aes(ymin=fraction-SE, ymax=fraction+SE), width=1) +
  scale_y_log10(limits = (c(0.001,1)), breaks=c(0.0001,0.001,0.01,0.1,1)) +
  scale_x_continuous( breaks=seq(0, 40, 2),limits = (c(0,40)))+
  scale_colour_manual(values=c("IGHA" = "orange","IGHG" = "red","IGHM" ="blue","naive_IGHM" = "green")) +
  xlab("Depth: S+NS") +
  ylab("Fraction of switches") +
  ggtitle("Pfizer - Switch Fraction") +
  theme_linedraw() +
  theme(panel.grid.major = element_line(colour = "blue",size = 0.1),
        panel.grid.minor = element_blank(),
        legend.position="none")
 # theme(  legend.title=element_blank(),
 #        legend.position=c(1,1), # legend position
 #         legend.justification=c(1,1),
  #        panel.grid.minor=element_blank(), # remove grid
  #       panel.grid.major=element_blank(),
  #        axis.text.x = element_text(angle=90, vjust=1),
          #axis.title.x=element_blank()
  #) +
  #guides(fill=guide_legend(nrow=2))


###########################################################################################################


# ~~~~~~flu details~~~~~~

# flu all time points
df.flu_node_diveded_son<- read.csv('flu/FLU_Node_son_divide_per_isotype.csv', header = T, stringsAsFactors = F)

df.flu_node_diveded_son[df.flu_node_diveded_son$Isotype.Name=="IGHG-1","Isotype.Name"]<-"IGHG"
df.flu_node_diveded_son[df.flu_node_diveded_son$Isotype.Name=="IGHG-2","Isotype.Name"]<-"IGHG"
df.flu_node_diveded_son[df.flu_node_diveded_son$Son.type=="IGHG-1","Son.type"]<-"IGHG"
df.flu_node_diveded_son[df.flu_node_diveded_son$Son.type=="IGHG-2","Son.type"]<-"IGHG"

# get the ratio beteewn each mutation type in flu
df.flu_son_number <- ddply(df.flu_node_diveded_son, c("Isotype.Name","Depth"), summarise,
                              sons = length(Isotype.Name)                              
)

# get the ratio beteewn each mutation type in flu
df.flu_spec_son <- ddply(df.flu_node_diveded_son, c("Isotype.Name","Son.type","Depth"), summarise,
                            sons = length(Son.type)                              
)

# remove rows with IGHA  as father | rempve rows with IGHE
df.flu_spec_son <- df.flu_spec_son[-grep('IGHA', df.flu_spec_son$Isotype.Name),]
df.flu_spec_son <- df.flu_spec_son[-grep('IGHE', df.flu_spec_son$Son.type),]
##  enter the fracion for each row 
for(i in 1: nrow(df.flu_spec_son))
{
  df.flu_spec_son[i,"totalSon"] <- df.flu_son_number[df.flu_son_number$Isotype.Name==df.flu_spec_son[i,"Isotype.Name"] &
                                                             df.flu_son_number$Depth ==df.flu_spec_son[i,"Depth"] ,"sons"]
}

# calc fracion (diffrentSon/total son)
df.flu_spec_son$fraction <- df.flu_spec_son$sons/df.flu_spec_son$totalSon

# create beans of 2:
df.flu_spec_son$bins2Depth <- df.flu_spec_son$Depth
df.flu_spec_son[df.flu_spec_son$bins2Depth%%2==1,"bins2Depth"] <- df.flu_spec_son[df.flu_spec_son$bins2Depth%%2==1,"bins2Depth"] +1 

df.flu_frac <- ddply(df.flu_spec_son, c("Isotype.Name","Son.type","bins2Depth"), summarise,
                        Sons = sum(sons), #
                        TotalSon =sum(totalSon)                              
)                        

df.flu_frac<- df.flu_frac[df.flu_frac$Isotype.Name!=df.flu_frac$Son.type,]
df.flu_frac$fraction <- df.flu_frac$Sons/df.flu_frac$TotalSon
df.flu_frac$SE <- df.flu_frac$fraction* (1/sqrt(df.flu_frac$Sons))

figure2a2 <-ggplot(df.flu_frac, aes(x=bins2Depth, y=fraction, colour=Son.type)) +
  facet_wrap(~ Isotype.Name) +
  geom_line() + geom_point()  +
  geom_errorbar(aes(ymin=fraction-SE, ymax=fraction+SE), width=1) +
  scale_y_log10(limits = (c(0.001,1)), breaks=c(0.0001,0.001,0.01,1)) +
  scale_x_continuous( breaks=seq(0, 40, 2),limits = (c(0,40)))+
  scale_colour_manual(values=c("IGHA" = "orange","IGHG-1" = "#CC6666","IGHG-2" = "brown",
                               "IGHM" ="blue","IGHD" = "#FF3399","IGHE" = "purple",
                               "naive_IGHM" = "green","IGHG" = "red")) +
  xlab("Depth: S+NS") +
  ylab("Fraction of switches") +
  ggtitle("flu - Switch Fraction") + theme_linedraw()+
  theme(  legend.title=element_blank(),
          legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  ) +
  guides(fill=guide_legend(nrow=2))


# ~~~~~~~flu pre vacine time points
df.flu_preVacine_node_diveded_son<- read.csv('flu/FLU_preVacine_Node_son_divide_per_isotype.csv', header = T, stringsAsFactors = F)

df.flu_preVacine_node_diveded_son[df.flu_preVacine_node_diveded_son$Isotype.Name=="IGHG-1","Isotype.Name"]<-"IGHG"
df.flu_preVacine_node_diveded_son[df.flu_preVacine_node_diveded_son$Isotype.Name=="IGHG-2","Isotype.Name"]<-"IGHG"
df.flu_preVacine_node_diveded_son[df.flu_preVacine_node_diveded_son$Son.type=="IGHG-1","Son.type"]<-"IGHG"
df.flu_preVacine_node_diveded_son[df.flu_preVacine_node_diveded_son$Son.type=="IGHG-2","Son.type"]<-"IGHG"

# get the ratio beteewn each mutation type in flu
df.flu_preVacine_son_number <- ddply(df.flu_preVacine_node_diveded_son, c("Isotype.Name","Depth"), summarise,
                           sons = length(Isotype.Name)                              
)

# get the ratio beteewn each mutation type in flu
df.flu_preVacine_spec_son <- ddply(df.flu_preVacine_node_diveded_son, c("Isotype.Name","Son.type","Depth"), summarise,
                         sons = length(Son.type)                              
)

# remove rows with IGHA  as father | rempve rows with IGHE
df.flu_preVacine_spec_son <- df.flu_preVacine_spec_son[-grep('IGHA', df.flu_preVacine_spec_son$Isotype.Name),]
df.flu_preVacine_spec_son <- df.flu_preVacine_spec_son[-grep('IGHE', df.flu_preVacine_spec_son$Son.type),]
##  enter the fracion for each row 
for(i in 1: nrow(df.flu_preVacine_spec_son))
{
  df.flu_preVacine_spec_son[i,"totalSon"] <- df.flu_preVacine_son_number[df.flu_preVacine_son_number$Isotype.Name==df.flu_preVacine_spec_son[i,"Isotype.Name"] &
                                                       df.flu_preVacine_son_number$Depth ==df.flu_preVacine_spec_son[i,"Depth"] ,"sons"]
}

# calc fracion (diffrentSon/total son)
df.flu_preVacine_spec_son$fraction <- df.flu_preVacine_spec_son$sons/df.flu_preVacine_spec_son$totalSon

# create beans of 2:
df.flu_preVacine_spec_son$bins2Depth <- df.flu_preVacine_spec_son$Depth
df.flu_preVacine_spec_son[df.flu_preVacine_spec_son$bins2Depth%%2==1,"bins2Depth"] <- df.flu_preVacine_spec_son[df.flu_preVacine_spec_son$bins2Depth%%2==1,"bins2Depth"] +1 

df.flu_preVacine_frac <- ddply(df.flu_preVacine_spec_son, c("Isotype.Name","Son.type","bins2Depth"), summarise,
                     Sons = sum(sons), #
                     TotalSon =sum(totalSon)                              
)                        

df.flu_preVacine_frac<- df.flu_preVacine_frac[df.flu_preVacine_frac$Isotype.Name!=df.flu_preVacine_frac$Son.type,]
df.flu_preVacine_frac$fraction <- df.flu_preVacine_frac$Sons/df.flu_preVacine_frac$TotalSon
df.flu_preVacine_frac$SE <- df.flu_preVacine_frac$fraction* (1/sqrt(df.flu_preVacine_frac$Sons))

figure2a3 <-ggplot(df.flu_preVacine_frac, aes(x=bins2Depth, y=fraction, colour=Son.type)) +
  facet_wrap(~ Isotype.Name) +
  geom_line() + geom_point()  +
  geom_errorbar(aes(ymin=fraction-SE, ymax=fraction+SE), width=1) +
  scale_y_log10(limits = (c(0.001,1)), breaks=c(0.0001,0.001,0.01,0.1,1)) +
  scale_x_continuous( breaks=seq(0, 40, 2),limits = (c(0,40)))+
  scale_colour_manual(values=c("IGHA" = "orange","IGHG-1" = "#CC6666","IGHG-2" = "brown",
                               "IGHM" ="blue","IGHD" = "#FF3399","IGHE" = "purple",
                               "naive_IGHM" = "green","IGHG" = "red")) +
  xlab("Depth: S+NS") +
  ylab("Fraction of switches") +
  ggtitle("flu - Switch Fraction") +
  theme_linedraw()+
  theme(panel.grid.major = element_line(colour = "blue",size = 0.1),
        panel.grid.minor = element_blank(),
        legend.position="none")
#  theme(  legend.title=element_blank(),
#          legend.position=c(1,1), # legend position
#          legend.justification=c(1,1),
#          panel.grid.minor=element_blank(), # remove grid
#          panel.grid.major=element_blank(),
#          axis.text.x = element_text(angle=90, vjust=1),
#          axis.title.x=element_blank()
#  ) +
#  guides(fill=guide_legend(nrow=2))


##### post vacine flu - as a function of a time point
df.flu_postVacineFrac<- read.csv('flu/FLU_isotype_Fraction_With_Point_Time.csv', header = T, stringsAsFactors = F)
df.flu_postVacineFrac <- df.flu_postVacineFrac[-grep('IGHE', df.flu_postVacineFrac$Child.Name.  ),]
df.flu_postVacineFrac[df.flu_postVacineFrac$Isotype.Name=="IGHG-1","Isotype.Name"]<-"IGHG"
df.flu_postVacineFrac[df.flu_postVacineFrac$Isotype.Name=="IGHG-2","Isotype.Name"]<-"IGHG"
df.flu_postVacineFrac[df.flu_postVacineFrac$Child.Name.=="IGHG-1","Child.Name."]<-"IGHG"
df.flu_postVacineFrac[df.flu_postVacineFrac$Child.Name.=="IGHG-2","Child.Name."]<-"IGHG"

df.Frac <- ddply(df.flu_postVacineFrac, c("Isotype.Name","Child.Name.","Time.Point"), summarise,
                 number = length(Isotype.Name)                 
)

df.FracAll <- ddply(df.flu_postVacineFrac, c("Isotype.Name","Time.Point"), summarise,
                    number = length(Isotype.Name)                 
)

df.Frac$totalSons <- 0
##  enter the fracion for each row in specific time point 
for(i in 1: nrow(df.Frac))
{
  df.Frac[i,"totalSons"] <- df.FracAll[df.FracAll$Isotype.Name==df.Frac[i,"Isotype.Name"] &
                                     
                                         df.FracAll$Time.Point ==df.Frac[i,"Time.Point"],"number"]
}

df.Frac$fraction <- df.Frac$number/df.Frac$totalSon
df.Frac <- df.Frac[-grep('IGHA', df.Frac$Isotype.Name),]
df.Frac<- df.Frac[df.Frac$Isotype.Name!=df.Frac$Child.Name.,]

figure2a4 <-ggplot(df.Frac, aes(x=Time.Point, y=fraction, fill=Child.Name.)) +
  facet_wrap( ~Isotype.Name) +
  geom_bar(colour="black", stat="identity",position=position_dodge(),size=.3)+
  scale_fill_manual(values=c("IGHA" = "orange","IGHG-1" = "#CC6666","IGHG-2" = "brown",
                             "IGHM" ="blue","IGHD" = "#FF3399","IGHE" = "purple",
                             "naive_IGHM" = "green","IGHG" = "red")) +
  xlab("Depth: S+NS") +
  ylab("Fraction of switches") +
  ggtitle("Flu post vacine - Switch Fraction per time point") +
  theme(plot.title = element_text(size = rel(1.5)))+ theme_gray(17)+
  theme(text = element_text(lineheight=.8,size = 14),
        legend.position="none",axis.text.x = element_text(angle = 90,hjust = 1))


figure3 <- grid.arrange(figure2a1,figure2a3, ncol=1) 

