library(ggplot2)
library(tidyr)
library(cowplot)
library(car)
library(nlme)
library(multcomp)
library(vegan)
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(readr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))

riquezaMC <- read_delim("MarChiquita.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
riquezaC1exp <- read_delim("Canal1exp.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
riquezaC1mues <- read_delim("Canal1samp.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
riquezaCTmues <- read_delim("CamposdelTuyu.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
marchi <- read_delim("wild guinea pig experiment.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

###############################
##### ALPHA DIVERSITY IN 6X6m #####

###############################
##### MAR CHIQUITA SITE #####

#EXPERIMENTAL


riquezaMC[is.na(riquezaMC)] <- 0
riquezaMC$Condition<-as.factor(riquezaMC$Condition)

#community data matrix
rspsmc<-riquezaMC[,6:41]

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

##### CANAL 1 SITE #####

#EXPERIMENTAL

riquezaC1exp[is.na(riquezaC1exp)] <- 0
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

RC1E


#SAMPLING

riquezaC1mues[is.na(riquezaC1mues)] <- 0
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
  scale_x_discrete( labels = c("with disturb"="With disturbances", "without disturb" = "Withput disturbances"))+ 
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


##### CAMPOS DEL TUYU SITE #####

#SAMPLING

riquezaCTmues[is.na(riquezaCTmues)] <- 0
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
  labs(x="", y="Richness", title = "Campos del Tuyú")+
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


RICHNESS<-plot_grid(RC1M, RCTM, 
                    RC1E, RMCE)


RICHNESS

FIG2<-ggdraw(RICHNESS)+
  draw_plot_label(
    c("Sampling", "Experiment"),
    x= c(1, 1),
    y= c(0.9, 0.5),
    size = c(15, 15),
    angle = c(270, 270),
    colour = c("black", "black")
  )

FIG2

###############################
##### NATIVE OR EXOTIC SPECIES?#####

###############################
##### MAR CHIQUITA SITE #####
natexoMC<-riquezaMC %>% select(-Disturbio) %>% 
  pivot_longer(cols = 6:41, names_to= "especies", values_to = "sino" ) %>% 
  reframe(especies=unique(especies)) 
unique(natexoMC$especies)

natexoMC$natexo<-c("si", "si", "si", "no", "si", "no", "si",
                   "no", "si", "si", "si", "si", "si", "si",
                   "si", "no", "si", "si", "si", "si", "si", "si",
                   "si", "si", "no", "no", "si", "no", "si",
                   "si", "no", "si", "no", "si", "si", "no") 

NATEXOMC<-riquezaMC %>% select(-Disturbio) %>% 
  pivot_longer(cols = 6:41, names_to= "especies", values_to = "sino" ) %>% 
  left_join(natexoMC) %>% 
  filter(sino != 0) %>% 
  group_by(Condition, Plot, natexo) %>% 
  summarize(total = sum(sino) ) %>% 
  pivot_wider( names_from = "natexo", values_from = "total") 

PMC<-NATEXOMC %>%   pivot_longer(cols = 3:4, names_to = "native", values_to = "richness")
c1<-aggregate(PMC$richness, list( PMC$Condition, PMC$Plot), sum) 

aggregate(PMC$richness, list(PMC$native), mean)


NEMC<-aov(richness~Condition*native, PMC)
shapiro.test(resid(NEMC))
leveneTest(NEMC)
anova(NEMC)
TukeyHSD(NEMC)

NEMCfig<-ggplot(PMC, aes(native,richness,fill=Condition))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ #primero se grafican las barras de error con los mínimos y máximos
  geom_boxplot(outlier.size = 0, #sin outliers
               position=position_dodge(width=0.85))+ #para que no se superpongan las barras
  labs(x="", y="")+ #defino la leyenda de los ejes
  scale_x_discrete( labels = c("no"="Exotic", "si" = "Native"))+ # renombro y reordeno los rótulos del eje x
  scale_y_continuous(limits=c(0,15))+
  labs(x="", y="", title = "Mar Chiquita")+
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c("olivedrab3", "salmon2"), 
                    labels=c("with disturb"="Control", "without disturb"="Exclosure"))+
  theme_cowplot()+
  theme(plot.title = element_text(hjust =  0.5),
        strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "top",
        legend.text = element_text(size=rel(1)))

NEMCfig

##### CANAL 1 SITE #####
#EXPERIMENTAL

natexoC1e<-riquezaC1exp %>% select(-SP1) %>% 
  pivot_longer(cols = 6:29, names_to= "especies", values_to = "sino" ) %>% 
  reframe(especies=unique(especies)) 
unique(natexoC1e$especies)

natexoC1e$natexo<-c("si", "si", "no", "si", "si",
                   "si", "si", "si", "si", "si",
                   "si", "si", "no", "no", "no", 
                   "si", "no", "si", "si", "no", 
                   "si", "si", "no", "si") 

NATEXOC1E<-riquezaC1exp %>% select(-SP1) %>% 
  pivot_longer(cols = 6:29, names_to= "especies", values_to = "sino" ) %>% 
  left_join(natexoC1e) %>% 
  group_by(Condition, Plot, natexo) %>% 
  summarize(total = sum(sino) ) %>% 
  pivot_wider( names_from = "natexo", values_from = "total") 

PC1E<-NATEXOC1E %>% 
  pivot_longer(cols = 3:4, names_to = "native", values_to = "richness")
c2<-aggregate(PC1E$richness, list( PC1E$Condition, PC1E$Plot), sum) 

aggregate(PC1E$richness, list(PC1E$native), mean)


NEC1E<-aov(richness~Condition*native, PC1E)
shapiro.test(resid(NEC1E))
leveneTest(NEC1E)
anova(NEC1E)
TukeyHSD(NEC1E)

NEC1Efig<-ggplot(PC1E, aes(native,richness,fill=Condition))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ #primero se grafican las barras de error con los mínimos y máximos
  geom_boxplot(outlier.size = 0, #sin outliers
               position=position_dodge(width=0.85))+ #para que no se superpongan las barras
  scale_x_discrete( labels = c("no"="Exotic", "si" = "Native"))+ # renombro y reordeno los rótulos del eje x
  scale_y_continuous(limits=c(0,15))+
  labs(x="", y="Richness", title = "Canal 1")+
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c("olivedrab3", "salmon2"), 
                    labels=c("with disturb"="Control", "without disturb"="Exclosure"))+
  theme_cowplot()+
  theme(plot.title = element_text(hjust =  0.5),
        strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "top",
        legend.text = element_text(size=rel(1)))

NEC1Efig

# SAMPLING
natexoC1m<-riquezaC1mues %>% select(-SP1, -SP2, -SP6, -SP7) %>% 
  pivot_longer(cols = 5:16, names_to= "especies", values_to = "sino" ) %>% 
  reframe(especies=unique(especies)) 
unique(natexoC1m$especies)

natexoC1m$natexo<-c("si", "si", "si", "si", "si", "si",
                    "si", "no", "si", "no", "no", "si") 

NATEXOC1M<-riquezaC1mues %>% select(-SP1, -SP2, -SP6, -SP7) %>% 
  pivot_longer(cols = 5:16, names_to= "especies", values_to = "sino" ) %>% 
  left_join(natexoC1m) %>% 
  group_by(Condition, Plot, natexo) %>% 
  summarize(total = sum(sino) ) %>% 
  pivot_wider( names_from = "natexo", values_from = "total")


PC1M<-NATEXOC1M %>% pivot_longer(cols = 3:4, names_to = "native", values_to = "richness")
c3<-aggregate(PC1M$richness, list( PC1M$Condition, PC1M$Plot), sum) 

aggregate(PC1M$richness, list(PC1M$native), mean)


NEC1M<-aov(richness~Condition*native, PC1M)
shapiro.test(resid(NEC1M))
leveneTest(NEC1M)
anova(NEC1M)
TukeyHSD(NEC1M)

NEC1Mfig<-ggplot(PC1M, aes(native,richness,fill=Condition))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ #primero se grafican las barras de error con los mínimos y máximos
  geom_boxplot(outlier.size = 0, #sin outliers
               position=position_dodge(width=0.85))+ #para que no se superpongan las barras
  labs(x="", y="Richness", title = "Canal 1")+ #defino la leyenda de los ejes
  scale_x_discrete( labels = c("no"="", "si" = ""))+ # renombro y reordeno los rótulos del eje x
  scale_y_continuous(limits=c(0,15))+
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c("seagreen2", "firebrick1"), 
                    labels=c("with disturb"="With disturbances", "without disturb"="Without disturbances"))+
  theme_cowplot()+
  theme(plot.title = element_text(hjust =  0.5),
        strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "top",
        legend.text = element_text(size=rel(1)))
NEC1Mfig

##### CAMPOS DEL TUYÚ #####
natexoCDT<-riquezaCTmues %>% select(-SP1, -SP2, -`Disturbio (%)`) %>% 
  pivot_longer(cols = 5:33, names_to= "especies", values_to = "sino" ) %>% 
  reframe(especies=unique(especies)) 
unique(natexoCDT$especies)

natexoCDT$natexo<-c("si", "si", "no", "si", "si", "no",
                    "si", "si", "si", "si", "si", "si",
                    "si", "si", "si", "si", "si", "no",
                    "si", "si", "no", "si", "si", "si",
                    "si", "no", "no", "si", "si") 

NATEXOCDT<-riquezaCTmues %>% select(-SP1, -SP2, -`Disturbio (%)`) %>% 
  pivot_longer(cols = 5:33, names_to= "especies", values_to = "sino" ) %>% 
  left_join(natexoCDT) %>% 
  group_by(Condition, Plot, natexo) %>% 
  summarize(total = sum(sino) ) %>% 
  pivot_wider( names_from = "natexo", values_from = "total")

PCDT<-NATEXOCDT %>%  pivot_longer(cols = 3:4, names_to = "native", values_to = "richness")
c4<-aggregate(PCDT$richness, list(PCDT$Condition, PCDT$Plot), sum) 

aggregate(PCDT$richness, list(PCDT$native), mean)

NECDT<-aov(richness~Condition*native, PCDT)
shapiro.test(resid(NECDT))
leveneTest(NECDT)
anova(NECDT)

#model structuring the variance
modzero<-glmmTMB(richness~0+Condition*native,  data=PCDT, dispformula = ~native) #dispformula para estructurar la varianza
modzero1<-glmmTMB(richness~native*Condition,  data=PCDT, dispformula = ~Condition) #dispformula para estructurar la varianza
modzero2<-glmmTMB(richness~native+Condition,  data=PCDT, dispformula = ~Condition) #dispformula para estructurar la varianza
modzero3<-glmmTMB(richness~native,  data=PCDT, dispformula = ~Condition) #dispformula para estructurar la varianza
modzero4<-glmmTMB(richness~Condition,  data=PCDT, dispformula = ~Condition) #dispformula para estructurar la varianza


anova(modzero1, modzero2)
anova(modzero2, modzero3)
anova(modzero2, modzero4)

#figure
NECDTMfig<-ggplot(PCDT, aes(native,richness,fill=Condition))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ #primero se grafican las barras de error con los mínimos y máximos
  geom_boxplot(outlier.size = 0, #sin outliers
               position=position_dodge(width=0.85))+ #para que no se superpongan las barras
  labs(x="", y="", title = "Campos del Tuyú")+ #defino la leyenda de los ejes
  scale_x_discrete( labels = c("no"="", "si" = ""))+ # renombro y reordeno los rótulos del eje x
  scale_y_continuous(limits=c(0,15))+
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c("seagreen2", "firebrick1"), 
                    labels=c("with disturb"="With disturbances", "without disturb"="Without disturbances"))+
  theme_cowplot()+
  theme(plot.title = element_text(hjust =  0.5),
        strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "top",
        legend.text = element_text(size=rel(1)))
NECDTMfig

##### ALL RICHNESS OF NATIVE OR EXOTIC SPECIES FIGURES TOGHETHER (SUP MAT 5) #####

PREVENANCE<-plot_grid(NEC1Mfig,NECDTMfig,  NEC1Efig, NEMCfig)


PREVENANCE

FIGSM5<-ggdraw(PREVENANCE)+
  draw_plot_label(
    c("Sampling", "Experiment"),
    x= c(1, 1),
    y= c(0.9, 0.5),
    size = c(15, 15),
    angle = c(270, 270),
    colour = c("black", "black")
  )

FIGSM5
##############################
##### BETA DIVERSITY #####

##############################
##### MAR CHIQUITA SITE #####

#community data matrix
rspsmc<-riquezaMC[,6:41]
#dissmilarity matrix
dissimilitudMCe <- vegdist(rspsmc, method="jaccard")
#treatment differences
beta.comuMCe<-betadisper(dissimilitudMCe,riquezaMC$Condition) 
permu.comuMCe<-permutest(beta.comuMCe, pairwise=T, permutations = how(nperm = 9999))

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
beta.comuC1e<-betadisper(dissimilitudC1e,riquezaC1exp$Condition) #de las diferencias en varianzas por grupo (adem
permu.comuC1e<-permutest(beta.comuC1e, pairwise=T, permutations = how(nperm = 999))
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
permu.comuC1m<-permutest(beta.comuC1m, pairwise=T, permutations = how(nperm = 999))
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
permu.comuCTm<-permutest(beta.comuCTm, pairwise=T, permutations = how(nperm = 999))
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

marchi[is.na(marchi)] <- 0

marchi$`Festuca arundinacea`<-as.numeric(marchi$`Festuca arundinacea`)
marchi$trat<-as.factor(marchi$trat)

##############################
##### ALPHA DIVERSITY #####

#community data matrix
spsMC<-marchi[,8:24]

#calculating richness
marchi$riq<-rowSums(decostand(spsMC, "pa"))

# anova
modelo1<-aov(riq~trat, data=marchi)
anova(modelo1)

plot(resid(modelo1))
qqnorm(residuals(modelo1))
qqline(residuals(modelo1))
shapiro.test(modelo1$residuals)
leveneTest(modelo1)
par(mfrow=c(2,2))
plot(modelo1)
par(mfrow=c(1,1))
TukeyHSD(modelo1)


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
aa<-permutest(beta.comuSUCMC, pairwise=T, permutations = how(nperm = 9999))
res<-adonis2(dissimilitudSUCMC~grupo, nperm=999)

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

marchiSIN<-marchi%>% dplyr::select(-Spartina, -`Festuca arundinacea`)

marchiSIN[is.na(marchiSIN)] <- 0
spssin<-marchiSIN[,8:22]
#subordinated species cover
marchiSIN$covS<-rowSums(spssin)

#anova

modelos<-aov(sqrt(covS)~trat, data=marchiSIN)
anova(modelos)
TukeyHSD(modelos)
plot(resid(modelos))
qqnorm(residuals(modelos))
qqline(residuals(modelos))
shapiro.test(modelos$residuals)
leveneTest(modelos)
par(mfrow=c(2,2))
plot(modelo1)
par(mfrow=c(1,1))

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
#Spartina densiflora and Festuca arundinacea#

marchiSF<-marchi[1:8]
marchiSF$festuca<-marchi$`Festuca arundinacea`
spsSF<-marchiSF[,8:9]
marchiSF$covD<-rowSums(spsSF)

#anova
modelod<-aov(covD~trat, data=marchiSF)
anova(modelod)
plot(resid(modelod))
qqnorm(residuals(modelod))
qqline(residuals(modelod))
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

#anova
modeloSD<-aov(`Suelo desnudo`~trat, data=marchi)
anova(modeloSD)
TukeyHSD(modeloSD)
plot(resid(modeloSD))
qqnorm(residuals(modeloSD))
qqline(residuals(modeloSD))
shapiro.test(modeloSD$residuals)
leveneTest(modeloSD)




SD<-ggplot(marchi, aes(trat, `Suelo desnudo`, fill=trat))+
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



##############################
##### SUPPLEMENTARY MATERIAL 4 #####
##############################
##### RICHNESS FOR ALL SITES TOGETHER #####

RIQexp<-riquezaC1exp  %>% select(Condition, Treatment,riqueza) %>% 
  mutate(Sitio="Canal 1 experimental", Tipo="Experimental") %>% 
  bind_rows(riquezaMC  %>% select(Condition, Treatment,riqueza)%>% 
              mutate(Sitio="Mar Chiquita",Tipo="Experimental"))


RIQexp$Sitio<-as.factor(RIQexp$Sitio)
RIQexp$Condition<-as.factor(RIQexp$Condition)
RIQexp$Tipo<-as.factor(RIQexp$Tipo)

m.riqexp<-aov(riqueza~Condition*Sitio, RIQexp)
shapiro.test(resid(m.riqexp))
leveneTest(m.riqexp)
anova(m.riqexp)

REXP<-ggplot(RIQexp, aes(Sitio,riqueza,fill=Condition))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ #primero se grafican las barras de error con los mínimos y máximos
  geom_boxplot(outlier.size = 0, #sin outliers
               position=position_dodge(width=0.85))+ #para que no se superpongan las barras
  labs(x="", y="Richness")+ #defino la leyenda de los ejes
  guides(fill=guide_legend(title="Treatment"))+
  scale_fill_manual(values = c("olivedrab3", "salmon2"), 
                    labels=c("with disturb"= "Control","without disturb"="Exclosure"))+
  scale_x_discrete(labels=c("Canal 1 experimental"="Canal 1", "Mar Chiquita"="Mar Chiquita"))+
  theme_cowplot()+
  theme(strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "top")
REXP


RIQsamp<-riquezaC1mues  %>% select(Condition,riqueza) %>% 
  mutate(Sitio="Canal 1 samplin", Tipo="Sampling") %>% 
  bind_rows(riquezaCTmues  %>% select(Condition,riqueza)%>% 
              mutate(Sitio="Campos del Tuyú",Tipo="Sampling"))


RIQsamp$Sitio<-as.factor(RIQsamp$Sitio)
RIQsamp$Condition<-as.factor(RIQsamp$Condition)
RIQsamp$Tipo<-as.factor(RIQsamp$Tipo)

m.riqsamp<-aov(sqrt(riqueza)~Condition*Sitio, RIQsamp)
shapiro.test(resid(m.riqsamp))
leveneTest(m.riqsamp)
anova(m.riqsamp)



RSAMP<-ggplot(RIQsamp, aes(Sitio,riqueza,fill=Condition))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ #primero se grafican las barras de error con los mínimos y máximos
  geom_boxplot(outlier.size = 0, #sin outliers
               position=position_dodge(width=0.85))+ #para que no se superpongan las barras
  labs(x="", y="Richness")+ #defino la leyenda de los ejes
  scale_fill_manual(values = c("seagreen2", "firebrick1"), 
                    labels=c("with disturb"="With disturbances", "without disturb"="Without disturbances"))+
  scale_x_discrete(labels=c("Canal 1 samplin"="Canal 1", "Campos del Tuyú"="Campos del Tuyú"))+
  guides(fill=guide_legend(title="Condition"))+
  theme_cowplot()+
  theme(strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "top")

RSAMP




##### BETA DIVERSITY FOR ALL SITES TOGETHER #####

DISTEXP<-distC1exp  %>% mutate(Sitio="Canal 1") %>% 
  bind_rows(distMCexp %>%  mutate(Sitio="Mar Chiquita")) 

m.distexp<-aov((dist)^2~Treatment*Sitio, DISTEXP)
shapiro.test(resid(m.distexp))
leveneTest(m.distexp)
anova(m.distexp)

DEXP<-ggplot(DISTEXP, aes(Sitio,dist,fill=Treatment))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ #primero se grafican las barras de error con los mínimos y máximos
  geom_boxplot(outlier.size = 0, #sin outliers
               position=position_dodge(width=0.85))+ #para que no se superpongan las barras
  labs(x="", y="Distance to the centroid")+ #defino la leyenda de los ejes
  guides(fill=guide_legend(title="Treatment"))+
  scale_fill_manual(values = c("olivedrab3", "salmon2"), 
                    labels=c("with disturb"= "Control","without disturb"="Exclosure"))+
#  scale_x_discrete(labels=c("Canal 1 experimental"="Canal 1", "Mar Chiquita"="Mar Chiquita"))+
  guides(fill=guide_legend(title="Treatment"))+
  theme_cowplot()+
  theme(strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        legend.position = "top")

DEXP


DISTSAMP<-distC1mue  %>% mutate(Sitio="Canal 1") %>% 
  bind_rows(distCTmue %>%  mutate(Sitio="Campos del Tuyú")) 

m.distsamp<-aov(dist~Condition*Sitio, DISTSAMP)
shapiro.test(resid(m.distsamp))
leveneTest(m.distsamp)
anova(m.distsamp)

DSAMP<-ggplot(DISTSAMP, aes(Sitio,dist,fill=Condition))+
  stat_boxplot(geom='errorbar',coef=100000,position=position_dodge(width=0.85))+ #primero se grafican las barras de error con los mínimos y máximos
  geom_boxplot(outlier.size = 0, #sin outliers
               position=position_dodge(width=0.85))+ #para que no se superpongan las barras
  labs(x="", y="Distance to the centroid")+ #defino la leyenda de los ejes
  guides(fill=guide_legend(title="Condition"))+
  scale_fill_manual(values = c("seagreen2", "firebrick1"), 
                    labels=c("with disturb"= "With disturbances","without disturb"="Without disturbances"))+
  theme_cowplot()+
  theme(strip.text.x = element_text(size = rel(1.2)),
        axis.title=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1)),
        #legend.position = "none",
        legend.text = element_text(size=rel(1)),
        legend.position = "top")

DSAMP



ayb<-plot_grid(RSAMP, DSAMP,
                REXP, DEXP)
ayb

SUP4<-ggdraw(ayb)+
  draw_plot_label(
    c("Sampling", "Experimental"),
    x= c(1, 1),
    y= c(0.9, 0.5),
    size = c(15, 15),
    angle = c(270, 270),
    colour = c("black", "black")
  )
SUP4

