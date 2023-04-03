#Analysis for interaction evenness, restricted connectance, and IGP
#Temporal constancy in the structure of a spider-focused food web with high rates of intraguild predation
#David H. Wise, Robin M. Mores, Jennifer M. Pajada-De La O, Matthew A. McCary 
#April 3 2023

#load packages
#library
library(ggplot2)
library(dplyr)
library(tidyverse)

#=========================================================
##Figures for restricted connectance

#=====import data======
#relative pathname
rc_web <- file.path(".", "Data", "restricted_connectance.csv")
print(rc_web)

#import data
rc<- read_csv(rc_web)

##=====================
##filter by the years
year1_rc<-
  rc%>%
  filter(Category == "Years")%>%
  rename( "Year" = Test)


#ggplot function
year1_rc%>%
  ggplot(aes(x = Year, y=Median, min=Lower, ymax=Upper)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=1.5, width=0.1)+
  geom_point(size = 5)+
  geom_pointrange()+
  ylim(0.1, 0.3)+
  xlab(NULL) +
  ylab("Restricted Connectance")+
  #geom_hline(aes(yintercept=0.5))+
  theme(axis.text = element_text(size=16, face="bold", colour = "black"),
        axis.title.y = element_text(size=16, face="bold", colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_rect(fill= "white"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))

##filter by seasons
seasons1_rc<-
  rc%>%
  filter(Category == "Seasons")%>%
  rename( "Season" = Test)

#ggplot function
seasons1_rc%>%
  mutate(Season = fct_relevel(Season, "Spring", "Summer", "Fall"))%>%
  ggplot(aes(x = Season, y=Median, min=Lower, ymax=Upper)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=1.5, width=0.1)+
  geom_point(size = 5)+
  geom_pointrange()+
  ylim(0.1, 0.3)+
  xlab(NULL) +
  ylab("Restricted Connectance")+
  theme(axis.text = element_text(size=16, face="bold", colour = "black"),
        axis.title.y = element_text(size=16, face="bold", colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_rect(fill= "white"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))

#============================================================
#Figures for interaction evenness

#=====import data======
#relative pathname
ie <- file.path(".", "Data", "interaction_evenness.csv")
print(ie)

#import data
int.e<- read_csv(ie)

##=====================
##filter by the years
year1_ie<-
  int.e%>%
  filter(Category == "Years")%>%
  rename( "Year" = Test)

#ggplot function
year1_ie%>%
  ggplot(aes(x = Year, y=Median, min=Lower, ymax=Upper)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=1.5, width=0.1)+
  geom_point(size = 5)+
  geom_pointrange()+
  ylim(0.8, 1.0)+
  xlab(NULL) +
  ylab("Interaction Evenness")+
  #geom_hline(aes(yintercept=0.5))+
  theme(axis.text = element_text(size=16, face="bold", colour = "black"),
        axis.title.y = element_text(size=16, face="bold", colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_rect(fill= "white"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))

##filter by seasons
seasons1_ie<-
  int.e%>%
  filter(Category == "Seasons")%>%
  rename( "Season" = Test)

#ggplot function
seasons1_ie%>%
  mutate(Season = fct_relevel(Season, "Spring", "Summer", "Fall"))%>%
  ggplot(aes(x = Season, y=Median, min=Lower, ymax=Upper)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=1.5, width=0.1)+
  geom_point(size = 5)+
  geom_pointrange()+
  ylim(0.8, 1.0)+
  xlab(NULL) +
  ylab("Interaction Evenness")+
  theme(axis.text = element_text(size=16, face="bold", colour = "black"),
        axis.title.y = element_text(size=16, face="bold", colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_rect(fill= "white"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))


#=================================================
##Intraguild predation pooled by year and season

#=====import data======
#relative pathname
ig <- file.path(".", "Data", "intraguild_predation.csv")
print(ig)

#import data
igp<- read_csv(ig)

##=====================
##filter by the years
year1_ig<-
  igp%>%
  filter(Category == "Years")%>%
  rename( "Year" = Test)

#ggplot function
year1_ig%>%
  ggplot(aes(x = Year, y=Percentage)) + 
  geom_point(size = 5)+
  geom_errorbar(aes(ymin=Percentage-Error, ymax=Percentage+Error), size=1.5, width=0.1)+
  ylim(0, 75)+
  xlab(NULL) +
  ylab("Percent IGP")+
  #geom_hline(aes(yintercept=0.5))+
  theme(axis.text = element_text(size=22, face="bold", colour = "black"),
        axis.title.y = element_text(size=22, face="bold", colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_rect(fill= "white"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))

##filter by seasons
seasons1_ig<-
  igp%>%
  filter(Category == "Seasons")%>%
  rename( "Season" = Test)

#ggplot function
seasons1_ig%>%
  mutate(Season = fct_relevel(Season, "Spring", "Summer", "Fall"))%>%
  ggplot(aes(x = Season, y=Percentage)) + 
  geom_point(size = 5)+
  geom_errorbar(aes(ymin=Percentage-Error, ymax=Percentage+Error), size=1.5, width=0.1)+
  #geom_pointrange()+
  ylim(0, 75)+
  xlab(NULL) +
  ylab("Percent IGP")+
  #geom_hline(aes(yintercept=0.5))+
  theme(axis.text = element_text(size=22, face="bold", colour = "black"),
        axis.title.y = element_text(size=22, face="bold", colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_rect(fill= "white"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))

#==================================================
##Figures for IGP for season and season not pooled

#=====import data======
#relative pathname
ig <- file.path(".", "Data", "IGP-season-by-year.csv")
print(ig)

#import data
igp<- read_csv(ig)

##=====================
#ggplot function
igp%>%
  mutate(Season = fct_relevel(Season, "Spring", "Summer", "Fall"))%>%
  ggplot(aes(x = Season, y=Percentage)) + 
  geom_point(size = 5)+
  geom_errorbar(aes(ymin=Percentage-Error, ymax=Percentage+Error), size=1.5, width=0.1)+
  ylim(0, 100)+
  xlab(NULL) +
  ylab("Percent IGP")+
  #geom_hline(aes(yintercept=0.5))+
  theme(axis.text = element_text(size=22, face="bold", colour = "black"),
        axis.title.y = element_text(size=22, face="bold", colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_rect(fill= "white"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 24, color = 'black', face = "bold"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))+
  facet_wrap(~Year, ncol = 1)