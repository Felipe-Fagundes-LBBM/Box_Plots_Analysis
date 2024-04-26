library(tidyverse)
library(readxl)
library(ggplot2)
library(glue)
library(ggtext)
library(ggpubr)

################# Zuric and Polish by nonARVC ################

setwd('c:/Users/User/Downloads')

data_arvc<-read.csv("data_qPCR_ARVC_validation_clinical_all.txt", sep = "\t", header = T)
tail(data_arvc)

data_arvc=data_arvc %>% 
  filter(cohort %in% c("Polish", "non-ARVC", "Zurich"))

T1_n = data_arvc %>% 
  filter(cohort == "non-ARVC") %>% 
  nrow()
T2_n = data_arvc %>% 
  filter(cohort == "Polish") %>% 
  nrow()
T3_n = data_arvc %>% 
  filter(cohort == "Zurich") %>% 
  nrow()


data_arvc$cohort <- factor(data_arvc$cohort, levels = c("non-ARVC", "Polish","Zurich"))

####################### functions for abbreviations ##################################

Tema_1 <- theme(axis.text = element_text(size = 19, color = "black"),
                axis.title = element_text(size = 20, color = "black"),
                panel.background = element_rect(fill="white"),
                axis.line = element_line(color = "black"),legend.position="top",
                legend.text = element_text(size = 19, color = "black"),
                legend.title = element_blank())

geom_caixa <- geom_boxplot(show.legend = F, outlier.shape = NA, alpha=0.5, width=0.7, coef=0)


geom_j <- geom_jitter(aes(shape = cohort, color = cohort), size = 2.5, show.legend = T, width=0.15)

guides <- guides(fill = "none")

stat_s <-  stat_summary(fun = median, show.legend=F, geom='crossbar', width=0.6, size=0.5, color='black')

scale_d <- scale_x_discrete(breaks = c("non-ARVC","Polish","Zurich"), 
                           labels = c(glue("non-ARVC (N={T1_n})"),
                                      glue("Polish (N={T2_n})"), 
                                      glue("Zurich (N={T3_n})")))

scale.fill <- scale_fill_grey(start = 1, end = 0)

scale.color <- scale_color_manual(values = c("gray", "steelblue3","coral3"))

compare_wilcox <- compare_means(hsa.miR.16.5p ~ cohort, data = data_arvc, method = "wilcox.test")


compare_kruskal <- compare_means(hsa.miR.16.5p ~ cohort, data = data_arvc, method = "kruskal.test")


my_compa <- my_comparisons <- list( c("non-ARVC", "Polish"), c("non-ARVC", "Zurich"),c("Polish", "Zurich"))

########################################################################################################



p1 <- data_arvc %>%
  ggplot(aes(x=cohort, hsa.miR.145.5p, fill = cohort))+
  geom_caixa +
  geom_j +
  guides(fill = "none")+
  stat_s +
  labs(x=NULL, y="hsa-miR-145-5p")+
  scale_d +
  scale_fill_grey(start = 1, end = 0)+
  scale.color +  
  Tema_1 

compare_wilcox
compare_kruskal
my_compa
p1 = p1 + stat_compare_means(comparisons = my_comparisons, size=7)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.53, size=7)     # Add global p-value
p1

p2 <- data_arvc %>%
  ggplot(aes(x=cohort, hsa.miR.16.5p, fill= cohort))+
  geom_caixa+
  geom_j +
  guides(fill = "none")+
  stat_s +
  labs(x=NULL, y="hsa-miR-16-5p")+
  scale_d +
  scale_fill_grey(start = 1, end = 0)+
  scale.color +  
  Tema_1


compare_wilcox
compare_kruskal

my_compa
p2 = p2 + stat_compare_means(comparisons = my_comparisons,size=7)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 7.8,size=7)     # Add global p-value
p2


p3 <- data_arvc %>%
  ggplot(aes(x=cohort, hsa.miR.15a.5p, fill= cohort))+
  geom_caixa +
  geom_j +
  guides(fill = "none")+
  stat_s +
  labs(x=NULL, y="hsa-miR-15a-5p")+
  scale_d +
  scale_fill_grey(start = 1, end = 0)+
  scale.color +  
  Tema_1

compare_wilcox
compare_kruskal

my_compa
p3 = p3 + stat_compare_means(comparisons = my_comparisons, size=7)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 5.3, size=7)     # Add global p-value
p3


p4 <- data_arvc %>%
  ggplot(aes(x=cohort, hsa.miR.92a.3p, fill= cohort))+
  geom_caixa +
  geom_j +
  guides(fill = "none")+
  stat_s +
  labs(x=NULL, y="hsa-miR-92a-3p")+
  scale_d +
  scale_fill_grey(start = 1, end = 0)+
  scale.color +  
  Tema_1

compare_wilcox
compare_kruskal

my_compa
p4 = p4 + stat_compare_means(comparisons = my_comparisons,size=7)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 12.5,size=7)     # Add global p-value
p4

p5 <- data_arvc %>%
  ggplot(aes(x=cohort, hsa.miR.19a.3p, fill= cohort))+
  geom_caixa +
  geom_j +
  guides(fill = "none")+
  stat_s +
  labs(x=NULL, y="hsa-miR-19a-3p")+
  scale_d +
  scale_fill_grey(start = 1, end = 0)+
  scale.color +  
  Tema_1

compare_wilcox
compare_kruskal

my_compa
p5 = p5 + stat_compare_means(comparisons = my_comparisons,size=7)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1.35,size=7)     # Add global p-value
p5

p6 <- data_arvc %>%
  ggplot(aes(x=cohort, hsa.miR.29a.3p, fill= cohort))+
  geom_caixa +
  geom_j +
  guides(fill = "none")+
  stat_s +
  labs(x=NULL, y="hsa-miR-29a-3p")+
  scale_d +
  scale_fill_grey(start = 1, end = 0)+
  scale.color +  
  Tema_1

compare_wilcox
compare_kruskal

my_compa
p6 = p6 + stat_compare_means(comparisons = my_comparisons, size=7)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.33, size=7)     # Add global p-value
p6

ggarrange(p1, p2, p3, p4, p5, p6, 
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)

