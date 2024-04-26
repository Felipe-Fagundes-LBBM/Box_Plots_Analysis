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
load("/Users/dennysonfonseca/Desktop/Enviroment(Graf_box_vivan).RData")
setwd("C:/Users/User/Documents/Prof andré")
list.files()

data_arvc <- read.csv("MASTER_data_antonia(alt).csv", sep = ";", header = T, dec = ",")
tail(data_arvc)
dim(data_arvc)
data_arvc_filted <- na.omit(data_arvc[ , c(3,7,9,11,13,15,17)]) #trocar o miR que vc estiver usando2
data_arvc_to_box <- data_arvc_filted |> tidyr::gather("hsamiR", "value", -TERCEIL)

#### 1
data_arvc_to_box_hsamiR1455p <- subset(data_arvc_to_box, hsamiR == "hsamiR1455p")
estat.test_hsamiR1455p <- data_arvc_to_box_hsamiR1455p |> group_by(hsamiR) |> dunn_test(value ~TERCEIL ) |> 
  adjust_pvalue(method = "fdr") |> add_significance() |> add_xy_position(x="TERCEIL", dodge = .8)
estat.test_filted_hsamiR1455p <- estat.test_hsamiR1455p|> filter(p.adj.signif != "ns")

#### 2
data_arvc_to_box_hsamiR92a3p <- subset(data_arvc_to_box, hsamiR == "hsamiR92a3p")
estat.test_hsamiR92a3p <- data_arvc_to_box_hsamiR92a3p |> group_by(hsamiR) |> dunn_test(value ~TERCEIL ) |> 
  adjust_pvalue(method = "fdr") |> add_significance() |> add_xy_position(x="TERCEIL", dodge = .8)
estat.test_filted_hsamiR92a3p <- estat.test_hsamiR92a3p |> filter(p.adj.signif != "ns")

#### 3
data_arvc_to_box_hsamir19a3p <- subset(data_arvc_to_box, hsamiR == "hsamir19a3p")
estat.test_hsamir19a3p <- data_arvc_to_box_hsamir19a3p |> group_by(hsamiR) |> dunn_test(value ~TERCEIL ) |> 
  adjust_pvalue(method = "fdr") |> add_significance() |> add_xy_position(x="TERCEIL", dodge = .8)
estat.test_filted_hsamir19a3p <- estat.test_hsamir19a3p |> filter(p.adj.signif != "ns")

#### 4
data_arvc_to_box_hsamiR15a5p <- subset(data_arvc_to_box, hsamiR == "hsamiR15a5p")
estat.test_hsamiR15a5p <- data_arvc_to_box_hsamiR15a5p |> group_by(hsamiR) |> dunn_test(value ~TERCEIL ) |> 
  adjust_pvalue(method = "fdr") |> add_significance() |> add_xy_position(x="TERCEIL", dodge = .8)
estat.test_filted_hsamiR15a5p <- estat.test_hsamiR15a5p |> filter(p.adj.signif != "ns")

#### 5
data_arvc_to_box_hsamiR29a3p <- subset(data_arvc_to_box, hsamiR == "hsamiR29a3p")
estat.test_hsamiR29a3p <- data_arvc_to_box_hsamiR29a3p |> group_by(hsamiR) |> dunn_test(value ~TERCEIL ) |> 
  adjust_pvalue(method = "fdr") |> add_significance() |> add_xy_position(x="TERCEIL", dodge = .8)
estat.test_filted_hsamiR29a3p <- estat.test_hsamiR29a3p |> filter(p.adj.signif != "ns")

#### 6
data_arvc_to_box_hsamiR165p <- subset(data_arvc_to_box, hsamiR == "hsamiR165p")
estat.test_hsamiR165p <- data_arvc_to_box_hsamiR165p |> group_by(hsamiR) |> dunn_test(value ~TERCEIL ) |> 
  adjust_pvalue(method = "fdr") |> add_significance() |> add_xy_position(x="TERCEIL", dodge = .8)
estat.test_filted_hsamiR165p <- estat.test_hsamiR165p |> filter(p.adj.signif != "ns")

stat_all <- rbind(estat.test_filted_hsamiR1455p, estat.test_filted_hsamiR92a3p, estat.test_filted_hsamir19a3p,
                 estat.test_filted_hsamiR15a5p, estat.test_filted_hsamiR29a3p, estat.test_filted_hsamiR165p)

data_arvc_to_box |>
  ggplot(aes(x = TERCEIL, y = value)) + 
  geom_boxplot(aes( fill = TERCEIL),  show.legend = T, outlier.shape = NA, alpha=0.5, width=0.7, coef=0,  lwd = 0.4) +
  geom_jitter(aes(color = TERCEIL), size = 1.5, show.legend = F, width=0.15, alpha = .5) +
  guides(fill = "none")+
  facet_wrap(~hsamiR, scales = "free", ncol = 3) + 
  labs(x=NULL, y="hsa-miR") + #trocar o miR
  scale_fill_manual(values = c("HC"="#4C0033", "BrS"= "#D8D8D8", "Low"= "#144272", "Interm"= "#205295", "High"= "#2C74B3")) +
  scale_color_manual(values = c("HC"="#4C0033", "BrS"= "#D8D8D8", "Low"= "#144272", "Interm"= "#205295", "High"= "#2C74B3")) +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 14, family = "Tahoma", face = "bold")) +
  stat_pvalue_manual(stat_all, lable = "p.adj.signif", tip.length = 0.005, y.position = "y.position")



