library(tidyverse)
library(readxl)
library(ggplot2)
library(glue)
library(ggtext)
library(ggpubr)
library(scales)

data<-read.csv("data_qPCR_ARVC_validation_clinical_all_boxplot.txt", sep = "\t", header = T)

data_group=data %>% 
 filter(type %in% c("Control", "BrS","ARVC"))

T0_n = data_group %>% 
  filter(type == "Control") %>% 
  nrow()
T1_n = data_group %>% 
  filter(type == "BrS") %>% 
  nrow()
T2_n = data_group %>% 
  filter(type == "ARVC") %>% 
  nrow()

data_group$type <- factor(data_group$type, 
                              levels = c("Control","BrS","ARVC"))
p1 <- data_group %>%
  ggplot(aes(x=type, hsamiR92a3p, fill=type))+
  labs(title="hsa-miR92a-3p", x=NULL, y="Relative miRNA Expression")+
    geom_boxplot(show.legend = F, outlier.shape = NA, 
               alpha=0.3, width=0.6, coef=0)+

  scale_x_discrete(breaks = c("Control","BrS","ARVC"), 
                   labels = c(glue("Control 
                                   (N={T0_n})"), 
                              glue("BrS 
                                   (N={T1_n})"),
                              glue("ARVC 
                                   (N={T2_n})")))+
  
  scale_fill_manual(values=c("seagreen3","coral","#144272"))+

  theme(axis.text = element_text(size = 19, color = "black"),
      axis.title = element_text(size = 20, color = "black"),
      panel.background = element_rect(fill="white"),
      axis.line = element_line(color = "black"))
  
compare_means(hsamiR92a3p ~ type, data = data_group, method = "wilcox.test")
compare_means(hsamiR92a3p ~ type, data = data_group, method = "kruskal.test")

my_comparisons <- list( c("Control","BrS"), c("Control","ARVC"), c("BrS","ARVC"))
p1 = p1 + stat_compare_means(comparisons = my_comparisons, size=5)
          #stat_compare_means(label.y = 13, size=5) 
p1

p1=p1+geom_point(aes(colour=factor(type), 
                 fill = factor(type)), shape=21, size = 4, alpha = .7,
                 show.legend = F,position=position_jitter(width=0.15)) + 
  scale_fill_manual(values=c("seagreen3", "coral","#144272")) + 
  scale_colour_manual(values=c("black", "black","black"))
  
p1
  


data<-read.csv("data_qPCR_ARVC_validation_clinical_all_boxplot.txt", sep = "\t", header = T)
data_score=data %>% 
  filter(scoreARVC %in% c("Low risk", "Interm. risk","High risk"))

T0_n = data_score %>% 
  filter(scoreARVC == "Low risk") %>% 
  nrow()
T1_n = data_score %>% 
  filter(scoreARVC == "Interm. risk") %>% 
  nrow()
T2_n = data_score %>% 
  filter(scoreARVC == "High risk") %>% 
  nrow()

data_score$scoreARVC=factor(data_score$scoreARVC, 
                            levels = c("Low risk","Interm. risk","High risk"))
p2 <- data_score %>%
 
  ggplot(aes(x=scoreARVC, hsamiR92a3p, fill = scoreARVC))+
  geom_boxplot(show.legend = F, outlier.shape = NA, alpha=0.3, width=0.6, coef=0)+
  geom_jitter(aes(color = scoreARVC), alpha = .7, size=4, 
              position=position_jitter(width=0.15),show.legend = F)+
  guides(fill = "none")+
  
  stat_summary(fun = median, show.legend=F, geom='crossbar', width=0.6, 
               size=0.15, colour="black")+
  labs(title="hsa-miR92a-3p", x=NULL, y="Relative miRNA Expression")+ 
  scale_x_discrete(breaks = c("Low risk","Interm. risk","High risk"), 
                   labels = c(glue("Low risk
                                   (N={T0_n})"), 
                              glue("Interm. risk
                                   (N={T1_n})"),
                              glue("High risk
                                   (N={T2_n})")))+
  scale_fill_manual(values=c("#144272","#144272","#144272"))+
  scale_color_manual(values = c("#144272","#144272","#144272"))+
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16, face = "bold",hjust = 0.5))

compare_means(hsamiR92a3p ~ scoreARVC, data = data_score, method = "wilcox.test")
#compare_means(hsamiR92a3p ~ scoreARVC, data = data_score, method = "kruskal.test")

my_comparisons <- list( c("Low risk","Interm. risk"), c("Low risk","High risk"),
                        c("Interm. risk","High risk"))
p2 = p2 + stat_compare_means(comparisons = my_comparisons, size=5)
#stat_compare_means(label.y = 13, size=5) 
p2

