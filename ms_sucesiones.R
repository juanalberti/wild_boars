library(ggplot2)
library(BiodiversityR)
library(tidyr)
library(cowplot)
library(car)
library(nlme)
library(multcomp)
library(dplyr)
library(rstatix)

riquezaMC <- read.csv("MarChiquita.csv",dec=",")
riquezaC1exp <- read.csv("Canal1exp.csv",dec=",")
riquezaC1mues <- read.csv("Canal1samp.csv",dec=",")
riquezaCTmues <- read.csv("CamposdelTuyu.csv",dec=",")
marchi <- read.csv("wild guinea pig experiment.csv",dec=",", fileEncoding = "latin1")

##### ALPHA DIVERSITY IN 6X6m2 #####

##### MAR CHIQUITA #####

# EXPERIMENTAL

riquezaMC
riquezaMC[is.na(riquezaMC)] <- 0
str(riquezaMC)
riquezaMC$Condition<-as.factor(riquezaMC$Condition)

#community data matrix
rspsmc<-riquezaMC[,6:52]

#calculating richness
riquezaMC$riqueza<-rowSums(rspsmc)

#t-test
mriqmc<-aov(riqueza~Treatment, riquezaMC)
shapiro.test(resid(mriqmc))
leveneTest(mriqmc)
t.test(riquezaMC$riqueza~riquezaMC$Treatment, riquezaMC, var.equal=T)

#figure
RMCE<-ggplot(riquezaMC, aes(Condition,riqueza,fill=Condition))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ 
  geom_boxplot(outlier.size = 0, 
               position=position_dodge(width=0.85))+ 
  labs(x="", y="", title = "Mar Chiquita")+ 
  scale_y_continuous(limits=c(0, 20))+
  scale_x_discrete( labels = c("with disturb"="Control", "without disturb" = "Exclosure"))+ 
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c("olivedrab3", "salmon2"))+
  theme_cowplot()+
  theme(plot.title = element_text(hjust =  0.5),
        strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "none",
        legend.text = element_text(size=rel(1)))

RMCE

##### CANAL 1 #####

#EXPERIMENTAL

riquezaC1exp
riquezaC1exp[is.na(riquezaC1exp)] <- 0
str(riquezaC1exp)

riquezaC1exp$Condition<-as.factor(riquezaC1exp$Condition)

#community data matrix
rspsc1<-riquezaC1exp[,6:30]

#calculating richness
riquezaC1exp$riqueza<-rowSums(rspsc1)

#t-test
mriqc1ex<-aov(riqueza~Treatment, riquezaC1exp)
shapiro.test(resid(mriqc1ex))
leveneTest(mriqc1ex)
t.test(riquezaC1exp$riqueza~riquezaC1exp$Treatment, riquezaC1exp, var.equal=T)

#figure
RC1E<-ggplot(riquezaC1exp, aes(Condition,riqueza,fill=Condition))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ 
  geom_boxplot(outlier.size = 0, 
               position=position_dodge(width=0.85))+ 
  labs(x="", y="", title = "Canal 1")+ 
  scale_y_continuous(limits=c(0, 20))+
  scale_x_discrete( labels = c("with disturb"="", "without disturb" = ""))+ 
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c("olivedrab3", "salmon2"))+
  theme_cowplot()+
  theme(plot.title = element_text(hjust =  0.5),
        strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "none",
        legend.text = element_text(size=rel(1)))

RC1E


#SAMPLING

riquezaC1mues
riquezaC1mues[is.na(riquezaC1mues)] <- 0
str(riquezaC1mues)

riquezaC1mues$Condition<-as.factor(riquezaC1mues$Condition)

#community data matrix
rspsc1m<-riquezaC1mues[,5:20]

#calculating richness
riquezaC1mues$riqueza<-rowSums(rspsc1m)

#t-test
mriqc1mue<-aov(riqueza~Condition, riquezaC1mues)
shapiro.test(resid(mriqc1mue))
leveneTest(mriqc1mue)
t.test(riquezaC1mues$riqueza~riquezaC1mues$Condition, riquezaC1mues, var.equal=T)

#figure
RC1M<-ggplot(riquezaC1mues, aes(Condition,riqueza,fill=Condition))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ 
  geom_boxplot(outlier.size = 0,
               position=position_dodge(width=0.85))+ 
  labs(x="", y="Richness", title = "Canal 1")+
  scale_x_discrete( labels = c("with disturb"="", "without disturb" = ""))+ 
  scale_y_continuous(limits=c(0, 20))+
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c("seagreen2", "firebrick1"))+
  theme_cowplot()+
  theme(plot.title = element_text(hjust =  0.5),
        strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "none",
        legend.text = element_text(size=rel(1)))


RC1M


##### CAMPOS DEL TUYU#####

#SAMPLING

riquezaCTmues
riquezaCTmues[is.na(riquezaCTmues)] <- 0
str(riquezaCTmues)

riquezaCTmues$Condition<-as.factor(riquezaCTmues$Condition)

#community datamatrix
rspsctm<-riquezaCTmues[,5:35]

#calculating richness
riquezaCTmues$riqueza<-rowSums(rspsctm)

#t-test
mriqct<-aov(log(riqueza)~Condition, riquezaCTmues)
shapiro.test(resid(mriqct))
leveneTest(mriqct)
t.test(riquezaCTmues$riqueza~riquezaCTmues$Condition, riquezaCTmues, var.equal=T)

#figure
RCTM<-ggplot(riquezaCTmues, aes(Condition,riqueza,fill=Condition))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ 
  geom_boxplot(outlier.size = 0, 
               position=position_dodge(width=0.85))+ 
  labs(x="", y="Richness", title = "Campos del Tuy?")+
  scale_x_discrete( labels = c("with disturb"="With disturbances", "without disturb" = "Without disturbances"))+ 
  guides(fill=guide_legend(title=""))+
  scale_y_continuous(limits=c(0, 20))+
  scale_fill_manual(values = c("seagreen2", "firebrick1"))+
  theme_cowplot()+
  theme(plot.title = element_text(hjust =  0.5),
        strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "none",
        legend.text = element_text(size=rel(1)))

RCTM


##### ALL RICHNESS FIGURES TOGETHER (FIG 2) #####


RICHNESS<-plot_grid(RC1M, RC1E, RCTM, RMCE)


RICHNESS

FIG2<-ggdraw(RICHNESS)+
  draw_plot_label(
    c("Sampling", "Experimental"),
    x= c(0.18, 0.68),
    y= c(0.06, 0.06),
    size = c(15, 15),
    colour = c("black", "black")
  )

FIG2


##############################
##### BETA DIVERSITY #####

##### MAR CHIQUITA #####

#community data matrix
rspsmc<-riquezaMC[,6:42]
#dissmilarity matrix
dissimilitudMCe <- vegdist(rspsmc, method="jaccard")
#treatment differences
beta.comuMCe<-betadisper(dissimilitudMCe,riquezaMC$Condition) 
(permu.comuMCe<-permutest(beta.comuMCe, pairwise=T, permutations = how(nperm = 9999)))
#new table with distances between groups
distMCexp<-data.frame(grupo=beta.comuMCe$group, dist=beta.comuMCe$distances)%>%
  separate(grupo,sep="\\.", into = c("Treatment"))

#figure
DMCE<-ggplot(distMCexp, aes(Treatment,dist,fill=Treatment))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ #primero se grafican las barras de error con los m?nimos y m?ximos
  geom_boxplot(outlier.size = 0, 
               position=position_dodge(width=0.85))+ 
  labs(x="Treatment", y="Distance to the centroid", title = "Mar Chiquita experimental")+ 
  scale_x_discrete( labels = c("with disturb"="Control", "without disturb" = "Exclosure"))+ # renombro y reordeno los r?tulos del eje x
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c("olivedrab3", "salmon2"))+
  theme_cowplot()+
  theme(plot.title = element_text(hjust =  0.5),
        strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "none",
        legend.text = element_text(size=rel(1)))
DMCE

##### CANAL 1 #####

### EXPERIMENTAL 

#community data matrix
rspsc1ex<-riquezaC1exp[,6:30]


dissimilitudC1e <- vegdist(rspsc1ex, method="jaccard")
beta.comuC1e<-betadisper(dissimilitudC1e,riquezaC1exp$Condition) 
(permu.comuC1e<-permutest(beta.comuC1e, pairwise=T, permutations = how(nperm = 999)))
distC1exp<-data.frame(grupo=beta.comuC1e$group, dist=beta.comuC1e$distances)%>%
  separate(grupo,sep="\\.", into = c( "Treatment"))

#figure
DC1E<-ggplot(distC1exp, aes(Treatment,dist,fill=Treatment))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ 
  geom_boxplot(outlier.size = 0, 
               position=position_dodge(width=0.85))+ 
  labs(x="Treatment", y="Distance to the centroid", title="Canal 1 experimental")+ 
  scale_x_discrete( labels = c("with disturb"="Control", "without disturb" = "Exclosure"))+ 
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c("olivedrab3", "salmon2"))+
  theme_cowplot()+
  theme(plot.title = element_text(hjust =  0.5),
        strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "none",
        legend.text = element_text(size=rel(1)))
DC1E

### SAMPLING

#community data matrix
rspsc1m<-riquezaC1mues[,5:20]

dissimilitudC1m <- vegdist(rspsc1m, method="jaccard")
beta.comuC1m<-betadisper(dissimilitudC1m,riquezaC1mues$Condition) #de las diferencias en varianzas por grupo (adem
(permu.comuC1m<-permutest(beta.comuC1m, pairwise=T, permutations = how(nperm = 999)))
distC1mue<-data.frame(grupo=beta.comuC1m$group, dist=beta.comuC1m$distances)%>%
  separate(grupo,sep="\\.", into = c( "Condition"))

#figure
DC1M<-ggplot(distC1mue, aes(Condition,dist,fill=Condition))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+
  geom_boxplot(outlier.size = 0, 
               position=position_dodge(width=0.85))+ 
  labs(x="Condition", y="Beta diversity", title = "Canal 1 sampling")+ 
  scale_x_discrete( labels = c("with disturb"="With disturbances", "without disturb" = "Without disturbances"))+ 
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c("seagreen2", "firebrick1"))+
  theme_cowplot()+
  theme(plot.title = element_text(hjust =  0.5),
        strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "none",
        legend.text = element_text(size=rel(1)))
DC1M

##### CAMPOS DEL TUYU #####

#community data matrix
rspsctm<-riquezaCTmues[,5:35]

dissimilitudCTm <- vegdist(rspsctm, method="jaccard")
beta.comuCTm<-betadisper(dissimilitudCTm,riquezaCTmues$Condition) #de las diferencias en varianzas por grupo (adem
(permu.comuCTm<-permutest(beta.comuCTm, pairwise=T, permutations = how(nperm = 999)))
distCTmue<-data.frame(grupo=beta.comuCTm$group, dist=beta.comuCTm$distances)%>%
  separate(grupo,sep="\\.", into = c( "Condition"))

#figure
DCTM<-ggplot(distCTmue, aes(Condition,dist,fill=Condition))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ 
  geom_boxplot(outlier.size = 0, 
               position=position_dodge(width=0.85))+ 
  labs(x="Condition", y="Beta diversity", title="Campos del Tuy? sampling")+ 
  scale_x_discrete( labels = c("with disturb"="With disturbances", "without disturb" = "Without disturbances"))+ 
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c("seagreen2", "firebrick1"))+
  theme_cowplot()+
  theme(plot.title = element_text(hjust =  0.5),
        strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "none",
        legend.text = element_text(size=rel(1)))
DCTM



##############################
##### WILD GUINEA PIG EXPERIMENT #####
marchi

marchi[is.na(marchi)] <- 0
str(marchi)
marchi$trat<-as.factor(marchi$trat)

##### ALPHA DIVERSITY #####

#community data matrix
spsMC<-marchi[,8:24]
#calculating richness
marchi$riq<-rowSums(decostand(spsMC, "pa"))

# anova
modelo1<-aov(riq~trat, data=marchi)
anova(modelo1)

shapiro.test(modelo1$residuals)
leveneTest(modelo1)

marchi %>% group_by(Tratamiento) %>%
  summarise_at(vars(riq), list(name= mean))

#figure
RIQ<-ggplot(marchi, aes(trat,riq,fill=trat))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ #primero se grafican las barras de error con los m?nimos y m?ximos
  geom_boxplot(outlier.size = 0, #sin outliers
               position=position_dodge(width=0.85))+ #para que no se superpongan las barras
  labs(x="", y="Richness")+ #defino la leyenda de los ejes
  scale_x_discrete(labels=c("a"= "Control", "b"= "Wild guinea pig \nexclosure", "c"= "Wild boar \nexclosure"))+
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c("peachpuff", "slateblue1", "lightgreen"))+
  theme_cowplot()+
  theme(strip.text.x = element_text(size = rel(1)),
        axis.title=element_text(size=rel(1)),
        axis.text=element_text(size=rel(1)),
        legend.position = "none",
        legend.text = element_text(size=rel(1)))


RIQ


##### BETA DIVERSITY #####

grupo<-interaction(marchi$trat, marchi$Tratamiento)

#dissimilarity matrix
dissimilitudSUCMC <- vegdist(spsMC, method="bray")
beta.comuSUCMC<-betadisper(dissimilitudSUCMC,grupo) 
(aa<-permutest(beta.comuSUCMC, pairwise=T, permutations = how(nperm = 9999)))

distSUCMC<-data.frame(grupo=beta.comuSUCMC$group, dist=beta.comuSUCMC$distances)%>%
  separate(grupo,sep="\\.", into = c("trat","Tratamiento"))

distSUCMC %>% group_by(Tratamiento) %>%
  summarise_at(vars(dist), list(names=mean))

#figure
DIST<-ggplot(distSUCMC, aes(trat,dist,fill=trat))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ #primero se grafican las barras de error con los m?nimos y m?ximos
  geom_boxplot(outlier.size = 0, #sin outliers
               position=position_dodge(width=0.85))+ #para que no se superpongan las barras
  labs(x="", y="Beta diversity")+ #defino la leyenda de los ejes
  scale_x_discrete(labels=c("a"= "Control", "b"= "Wild guinea pig \nexclosure", "c"= "Wild boar \nexclosure"))+
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c("peachpuff", "slateblue1", "lightgreen"))+
  theme_cowplot()+
  theme(strip.text.x = element_text(size = rel(1)),
        axis.title=element_text(size=rel(1)),
        axis.text=element_text(size=rel(1)),
        legend.position = "none",
        legend.text = element_text(size=rel(1)))

DIST



##### FIGURE 3 #####
FIG3<-plot_grid(RIQ ,DIST, labels = c("A", "B"))

FIG3

##### SUBORDINATED SPECIES COVER #####

marchiSIN<-marchi%>% dplyr::select(-Spartina, -Festuca.arundinacea)
spssin<-marchiSIN[,8:22]
str(marchiSIN)
#subordinated species cover
marchiSIN$covS<-rowSums(spssin)

#anova

modelos<-aov(sqrt(covS)~trat, data=marchiSIN)
anova(modelos)
TukeyHSD(modelos)
shapiro.test(modelos$residuals)
leveneTest(modelos)

#figure
COVSIN<-ggplot(marchiSIN, aes(trat,covS,fill=trat))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ #primero se grafican las barras de error con los m?nimos y m?ximos
  geom_boxplot(outlier.size = 0, #sin outliers
               position=position_dodge(width=0.85))+ #para que no se superpongan las barras
  labs(x="", y="Subordinate \nplant cover (%)")+ #defino la leyenda de los ejes
  guides(fill=guide_legend(title=""))+
  scale_x_discrete(labels=c("a"= "Control", "b"= "Wild guinea pig \nexclosure", "c"= "Wild boar \nexclosure"))+
  scale_fill_manual(values = c("peachpuff", "slateblue1", "lightgreen"))+
  theme_cowplot()+
  theme(strip.text.x = element_text(size = rel(1)),
        axis.title=element_text(size=rel(1)),
        axis.text=element_text(size=rel(1)),
        legend.position = "none",
        legend.text = element_text(size=rel(1)))


COVSIN


##### DOMINANT SPECIES COVER #####

marchiSF<-marchi[1:8]
marchiSF$festuca<-marchi$Festuca.arundinacea
spsSF<-marchiSF[,8:9]
marchiSF$covD<-rowSums(spsSF)

summary(marchiSF$covD ~marchiSF$Tratamiento)

#anova
modelod<-aov(covD~trat, data=marchiSF)
anova(modelod)
shapiro.test(modelod$residuals)
leveneTest(modelod)
TukeyHSD(modelod)

#figure
COVSF<-ggplot(marchiSF, aes(trat,covD,fill=trat))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ #primero se grafican las barras de error con los m?nimos y m?ximos
  geom_boxplot(outlier.size = 0, #sin outliers
               position=position_dodge(width=0.85))+ #para que no se superpongan las barras
  labs(x="", y="Dominant \nplant cover (%)")+ #defino la leyenda de los ejes
  guides(fill=guide_legend(title=""))+
  scale_x_discrete(labels=c("a"= "", "b"= "", "c"=""))+
  scale_fill_manual(values = c("peachpuff", "slateblue1", "lightgreen"))+
  theme_cowplot()+
  theme(strip.text.x = element_text(size = rel(1)),
        axis.title=element_text(size=rel(1)),
        axis.text=element_text(size=rel(1)),
        legend.position = "none",
        legend.text = element_text(size=rel(1)))


COVSF

##### BARE SOIL #####

summary(marchi$Suelo.desnudo ~marchi$Tratamiento)

#anova
modeloSD<-aov(Suelo.desnudo~trat, data=marchi)
anova(modeloSD)
TukeyHSD(modeloSD)
shapiro.test(modeloSD$residuals)
leveneTest(modeloSD)




SD<-ggplot(marchi, aes(trat,Suelo.desnudo,fill=trat))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ #primero se grafican las barras de error con los m?nimos y m?ximos
  geom_boxplot(outlier.size = 0, #sin outliers
               position=position_dodge(width=0.85))+ #para que no se superpongan las barras
  labs(x="", y="Bare soil (%)")+ #defino la leyenda de los ejes
  scale_x_discrete(labels=c("a"= "Control", "b"= "Wild guinea pig \nexclosure", "c"= "Wild boar \nexclosure"))+
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c("peachpuff", "slateblue1", "lightgreen"))+
  theme_cowplot()+
  theme(strip.text.x = element_text(size = rel(1)),
        axis.title=element_text(size=rel(1)),
        axis.text=element_text(size=rel(1)),
        legend.position = "none",
        legend.text = element_text(size=rel(1)))

SD

##### FIGURE 4 #####

cover<-plot_grid(COVSF, COVSIN, nrow=2, labels = c("A", "B")) 


FIG4<-plot_grid(cover, SD, ncol=2, labels= c ("", "C"))

FIG4
