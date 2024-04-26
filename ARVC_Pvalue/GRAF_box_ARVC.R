library(tidyverse)
library(readxl)
library(ggplot2)
library(glue)
library(ggtext)
library(ggpubr)
library(hrbrthemes)
library(rstatix)
library(scales)
library(ggbreak)

install.packages("ggtext")
install.packages("hrbrthemes")
install.packages("ggbreak")

################# Zuric and Polish ARVC for diagnosis and prognosis ################

setwd("C:/Users/User/Documents/Prof andré")
list.files()

data_arvc <- read.csv("MASTER_data_antonia(alt).csv", sep = ";", header = T, dec = ",")

tail(data_arvc)

dim(data_arvc)


data_arvc_filted <- na.omit(data_arvc[ , c(3,7,9,11,13,15,17)]) 

dplyr::glimpse(data_arvc_filted)


data_arvc_to_box <- data_arvc_filted %>% tidyr::gather("hsamiR" ,"value",  -TERCEIL )
data_arvc_to_box_3 <- subset(data_arvc_to_box, hsamiR == "hsamiR1455p" | hsamiR == "hsamiR15a5p"| hsamiR == "hsamiR165p")

estat.test_3 <- data_arvc_to_box_3 %>% group_by(hsamiR) %>% dunn_test(value ~TERCEIL ) %>% 
  adjust_pvalue(method = "fdr") %>% add_significance() %>% add_xy_position(x="TERCEIL", dodge = .8)

View(estat.test_filted)

estat.test_filted_3 <- estat.test_3 %>% filter(p.adj.signif != "ns")


p1 <- data_arvc_to_box_3 %>%
  ggplot(aes(x = TERCEIL, y = value )) + 
  geom_boxplot(aes( fill = TERCEIL),  show.legend = F, outlier.shape = NA, alpha=0.5, width=0.7, coef=0) +
  geom_jitter(aes(fill = TERCEIL), size = 1.5, show.legend = F, width=0.15) +
  guides(fill = "none")+
  facet_wrap(~hsamiR, scales = "fixed", ncol = 3) + 
  #coord_cartesian (ylim = c(0,5)) +
 # scale_y_continuous(limits = c(0,11), breaks = pretty(c(0,2,11), n=0.5)) +
 # stat_compare_means(method = "kruskal.test", size=5, label.y = 11) +
  # stat_summary(fun = median, show.legend=F, geom='crossbar', width=0.6, size=0.5, color='black')+
  # labs(x=NULL, y="hsa-miR-92a3p", title ="Kruskal-Wallis 6.22e-14") + #trocar o miR
  # scale_x_discrete(breaks = c("HC","BrS","Low", "Interm", "High"), 
  #                  labels = c(glue("HC (N={Low_n})"),
  #                             glue("BrS (N={Interm_n})"), 
  #                             glue("Low (N={High_n})"), 
  #                             glue("Interm (N={T4_n})"),
  #                       glue("High (N={T5_n})")))+
  scale_fill_manual(values = c("HC"="#0000CD", "BrS"= "#32CD32", "Low"= "#DC143C", "Interm"= "#DC143C", "High"= "#DC143C")) +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black")) +
  stat_pvalue_manual(estat.test_filted_3, lable = "p.adj.signif", tip.length = 0.005, y.position = "y.position")


p1 + scale_y_break(c(0.5), scales = "free") 


compare_means(hsamir92a3p ~ cohort, data = data_arvc, method = "wilcox.test") #trocar o miR
compare_means(hsamir92a3p ~ cohort, data = data_arvc, method = "kruskal.test") #trocar o miR
my_comparisons <- list( c("HC", "BrS"),c("HC", "Low"),c("HC", "Interm"),c("HC", "High"))
p1 = p1 + stat_compare_means(comparisons = my_comparisons, size=4)+ # ajustar o size
  stat_compare_means(label.y = 39.5, size=7)     # ajustar o size e o label.y se necessario

p1

?ggbreak

