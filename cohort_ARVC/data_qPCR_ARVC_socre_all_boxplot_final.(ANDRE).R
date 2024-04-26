data<-read.csv("data_qPCR_ARVC_validation_clinical_all_boxplot.txt", sep = "\t", header = T)
names(data)
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
  geom_point(aes(colour=factor(type), 
                 fill = factor(type)), shape=21, size = 3, alpha = .7,
             show.legend = F,position=position_jitter(width=0.15))+
  scale_x_discrete(breaks = c("Control","BrS","ARVC"), 
                   labels = c(glue("Control 
                                   (N={T0_n})"), glue("BrS 
                                   (N={T1_n})"), glue("ARVC 
                                   (N={T2_n})")))+
  scale_fill_manual(values=c("seagreen3","coral","#144272"))+
  scale_colour_manual(values=c("black", "black","black"))+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16, face = "bold",hjust = 0.5))

compare_means(hsamiR92a3p ~ type, data = data_group, method = "wilcox.test")
#compare_means(hsamiR92a3p ~ type, data = data_group, method = "kruskal.test")
my_comparisons <- list( c("Control","BrS"), c("Control","ARVC"), c("BrS","ARVC"))
p1 = p1 + stat_compare_means(comparisons = my_comparisons, size=5)
#stat_compare_means(label.y = 13, size=5) 
p1

p2 <- data_group %>%
  ggplot(aes(x=type, hsamiR165p, fill=type))+
  labs(title="hsa-miR16-5p", x=NULL, y="Relative miRNA Expression")+
  geom_boxplot(show.legend = F, outlier.shape = NA, 
               alpha=0.3, width=0.6, coef=0)+
  geom_point(aes(colour=factor(type), 
                 fill = factor(type)), shape=21, size = 3, alpha = .7,
             show.legend = F,position=position_jitter(width=0.15))+
  scale_x_discrete(breaks = c("Control","BrS","ARVC"), 
                   labels = c(glue("Control 
                                   (N={T0_n})"), glue("BrS 
                                   (N={T1_n})"), glue("ARVC 
                                   (N={T2_n})")))+
  scale_fill_manual(values=c("seagreen3","coral","#144272"))+
  scale_colour_manual(values=c("black", "black","black"))+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16, face = "bold",hjust = 0.5))

compare_means(hsamiR165p ~ type, data = data_group, method = "wilcox.test")
#compare_means(hsamiR92a3p ~ type, data = data_group, method = "kruskal.test")
my_comparisons <- list( c("Control","BrS"), c("Control","ARVC"), c("BrS","ARVC"))
p2 = p2 + stat_compare_means(comparisons = my_comparisons, size=5)
#stat_compare_means(label.y = 13, size=5) 
p2

p3 <- data_group %>%
  ggplot(aes(x=type, hsamiR15a5p, fill=type))+
  labs(title="hsa-miR15a-5p", x=NULL, y="Relative miRNA Expression")+
  geom_boxplot(show.legend = F, outlier.shape = NA, 
               alpha=0.3, width=0.6, coef=0)+
  geom_point(aes(colour=factor(type), 
                 fill = factor(type)), shape=21, size = 3, alpha = .7,
             show.legend = F,position=position_jitter(width=0.15))+
  scale_x_discrete(breaks = c("Control","BrS","ARVC"), 
                   labels = c(glue("Control 
                                   (N={T0_n})"), glue("BrS 
                                   (N={T1_n})"), glue("ARVC 
                                   (N={T2_n})")))+
  scale_fill_manual(values=c("seagreen3","coral","#144272"))+
  scale_colour_manual(values=c("black", "black","black"))+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16, face = "bold",hjust = 0.5))

compare_means(hsamiR15a5p ~ type, data = data_group, method = "wilcox.test")
#compare_means(hsamiR15a5p ~ type, data = data_group, method = "kruskal.test")
my_comparisons <- list( c("Control","BrS"), c("Control","ARVC"), c("BrS","ARVC"))
p3 = p3 + stat_compare_means(comparisons = my_comparisons, size=5)
#stat_compare_means(label.y = 13, size=5) 
p3

p4 <- data_group %>%
  ggplot(aes(x=type, hsamiR1455p, fill=type))+
  labs(title="hsa-miR145-5p", x=NULL, y="Relative miRNA Expression")+
  geom_boxplot(show.legend = F, outlier.shape = NA, 
               alpha=0.3, width=0.6, coef=0)+
  geom_point(aes(colour=factor(type), 
                 fill = factor(type)), shape=21, size = 3, alpha = .7,
             show.legend = F,position=position_jitter(width=0.15))+
  scale_x_discrete(breaks = c("Control","BrS","ARVC"), 
                   labels = c(glue("Control 
                                   (N={T0_n})"), glue("BrS 
                                   (N={T1_n})"), glue("ARVC 
                                   (N={T2_n})")))+
  scale_fill_manual(values=c("seagreen3","coral","#144272"))+
  scale_colour_manual(values=c("black", "black","black"))+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16, face = "bold",hjust = 0.5))

compare_means(hsamiR1455p ~ type, data = data_group, method = "wilcox.test")
#compare_means(hsamiR1455p ~ type, data = data_group, method = "kruskal.test")
my_comparisons <- list( c("Control","BrS"), c("Control","ARVC"), c("BrS","ARVC"))
p4 = p4 + stat_compare_means(comparisons = my_comparisons, size=5)
#stat_compare_means(label.y = 13, size=5) 
p4

p5 <- data_group %>%
  ggplot(aes(x=type, hsamiR19a3p, fill=type))+
  labs(title="hsa-miR19a-3p", x=NULL, y="Relative miRNA Expression")+
  geom_boxplot(show.legend = F, outlier.shape = NA, 
               alpha=0.3, width=0.6, coef=0)+
  geom_point(aes(colour=factor(type), 
                 fill = factor(type)), shape=21, size = 3, alpha = .7,
             show.legend = F,position=position_jitter(width=0.15))+
  scale_x_discrete(breaks = c("Control","BrS","ARVC"), 
                   labels = c(glue("Control 
                                   (N={T0_n})"), glue("BrS 
                                   (N={T1_n})"), glue("ARVC 
                                   (N={T2_n})")))+
  scale_fill_manual(values=c("seagreen3","coral","#144272"))+
  scale_colour_manual(values=c("black", "black","black"))+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16, face = "bold",hjust = 0.5))

compare_means(hsamiR19a3p ~ type, data = data_group, method = "wilcox.test")
#compare_means(hsamiR19a3p ~ type, data = data_group, method = "kruskal.test")
my_comparisons <- list( c("Control","BrS"), c("Control","ARVC"), c("BrS","ARVC"))
p5 = p5 + stat_compare_means(comparisons = my_comparisons, size=5)
#stat_compare_means(label.y = 13, size=5) 
p5

p6 <- data_group %>%
  ggplot(aes(x=type, hsamiR29a3p, fill=type))+
  labs(title="hsa-miR29a-3p", x=NULL, y="Relative miRNA Expression")+
  geom_boxplot(show.legend = F, outlier.shape = NA, 
               alpha=0.3, width=0.6, coef=0)+
  geom_point(aes(colour=factor(type), 
                 fill = factor(type)), shape=21, size = 3, alpha = .7,
             show.legend = F,position=position_jitter(width=0.15))+
  scale_x_discrete(breaks = c("Control","BrS","ARVC"), 
                   labels = c(glue("Control 
                                   (N={T0_n})"), glue("BrS 
                                   (N={T1_n})"), glue("ARVC 
                                   (N={T2_n})")))+
  scale_fill_manual(values=c("seagreen3","coral","#144272"))+
  scale_colour_manual(values=c("black", "black","black"))+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16, face = "bold",hjust = 0.5))

compare_means(hsamiR29a3p ~ type, data = data_group, method = "wilcox.test")
#compare_means(hsamiR29a3p ~ type, data = data_group, method = "kruskal.test")
my_comparisons <- list( c("Control","BrS"), c("Control","ARVC"), c("BrS","ARVC"))
p6 = p6 + stat_compare_means(comparisons = my_comparisons, size=5)
#stat_compare_means(label.y = 13, size=5) 
p6
###################

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
p7 <- data_score %>%
  
  ggplot(aes(x=scoreARVC, hsamiR92a3p, fill = scoreARVC))+
  geom_boxplot(show.legend = F, outlier.shape = NA, alpha=0.3, width=0.6, coef=0)+
  geom_point(aes(colour=factor(scoreARVC), 
                 fill = factor(scoreARVC)), shape=21, size = 3, alpha = .7,
             show.legend = F,position=position_jitter(width=0.15))+
  guides(fill = "none")+
  stat_summary(fun = median, show.legend=F, geom='crossbar', width=0.6, 
               size=0.15, colour="black")+
  labs(title="hsa-miR92a-3p", x=NULL, y="Relative miRNA Expression")+ 
  scale_x_discrete(breaks = c("Low risk","Interm. risk","High risk"), 
                   labels = c(glue("Low risk
                                   (N={T0_n})"), glue("Interm. risk
                                   (N={T1_n})"), glue("High risk
                                   (N={T2_n})")))+
  scale_fill_manual(values=c("#144272","#144272","#144272"))+
  scale_color_manual(values = c("black","black","black"))+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

compare_means(hsamiR92a3p ~ scoreARVC, data = data_score, method = "wilcox.test")
#compare_means(hsamiR92a3p ~ scoreARVC, data = data_score, method = "kruskal.test")

my_comparisons <- list( c("Low risk","Interm. risk"), c("Low risk","High risk"),
                        c("Interm. risk","High risk"))
p7 = p7 + stat_compare_means(comparisons = my_comparisons, size=5)
#stat_compare_means(label.y = 13, size=5) 
p7

p8 <- data_score %>%
  
  ggplot(aes(x=scoreARVC, hsamiR165p, fill = scoreARVC))+
  geom_boxplot(show.legend = F, outlier.shape = NA, alpha=0.3, width=0.6, coef=0)+
  geom_point(aes(colour=factor(scoreARVC), 
                 fill = factor(scoreARVC)), shape=21, size = 3, alpha = .7,
             show.legend = F,position=position_jitter(width=0.15))+
  guides(fill = "none")+
  stat_summary(fun = median, show.legend=F, geom='crossbar', width=0.6, 
               size=0.15, colour="black")+
  labs(title="hsa-miR16-5p", x=NULL, y="Relative miRNA Expression")+ 
  scale_x_discrete(breaks = c("Low risk","Interm. risk","High risk"), 
                   labels = c(glue("Low risk
                                   (N={T0_n})"), glue("Interm. risk
                                   (N={T1_n})"), glue("High risk
                                   (N={T2_n})")))+
  scale_fill_manual(values=c("#144272","#144272","#144272"))+
  scale_color_manual(values = c("black","black","black"))+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

compare_means(hsamiR165p ~ scoreARVC, data = data_score, method = "wilcox.test")
#compare_means(hsamiR165p ~ scoreARVC, data = data_score, method = "kruskal.test")

my_comparisons <- list( c("Low risk","Interm. risk"), c("Low risk","High risk"),
                        c("Interm. risk","High risk"))
p8 = p8 + stat_compare_means(comparisons = my_comparisons, size=5)
#stat_compare_means(label.y = 13, size=5) 
p8

p9 <- data_score %>%
  
  ggplot(aes(x=scoreARVC, hsamiR15a5p, fill = scoreARVC))+
  geom_boxplot(show.legend = F, outlier.shape = NA, alpha=0.3, width=0.6, coef=0)+
  geom_point(aes(colour=factor(scoreARVC), 
                 fill = factor(scoreARVC)), shape=21, size = 3, alpha = .7,
             show.legend = F,position=position_jitter(width=0.15))+
  guides(fill = "none")+
  stat_summary(fun = median, show.legend=F, geom='crossbar', width=0.6, 
               size=0.15, colour="black")+
  labs(title="hsa-miR15a-5p", x=NULL, y="Relative miRNA Expression")+ 
  scale_x_discrete(breaks = c("Low risk","Interm. risk","High risk"), 
                   labels = c(glue("Low risk
                                   (N={T0_n})"), glue("Interm. risk
                                   (N={T1_n})"), glue("High risk
                                   (N={T2_n})")))+
  scale_fill_manual(values=c("#144272","#144272","#144272"))+
  scale_color_manual(values = c("black","black","black"))+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

compare_means(hsamiR15a5p ~ scoreARVC, data = data_score, method = "wilcox.test")
#compare_means(hsamiR15a5p ~ scoreARVC, data = data_score, method = "kruskal.test")

my_comparisons <- list( c("Low risk","Interm. risk"), c("Low risk","High risk"),
                        c("Interm. risk","High risk"))
p9 = p9 + stat_compare_means(comparisons = my_comparisons, size=5)
#stat_compare_means(label.y = 13, size=5) 
p9

p10 <- data_score %>%
  
  ggplot(aes(x=scoreARVC, hsamiR1455p, fill = scoreARVC))+
  geom_boxplot(show.legend = F, outlier.shape = NA, alpha=0.3, width=0.6, coef=0)+
  geom_point(aes(colour=factor(scoreARVC), 
                 fill = factor(scoreARVC)), shape=21, size = 3, alpha = .7,
             show.legend = F,position=position_jitter(width=0.15))+
  guides(fill = "none")+
  stat_summary(fun = median, show.legend=F, geom='crossbar', width=0.6, 
               size=0.15, colour="black")+
  labs(title="hsa-miR145-5p", x=NULL, y="Relative miRNA Expression")+ 
  scale_x_discrete(breaks = c("Low risk","Interm. risk","High risk"), 
                   labels = c(glue("Low risk
                                   (N={T0_n})"), glue("Interm. risk
                                   (N={T1_n})"), glue("High risk
                                   (N={T2_n})")))+
  scale_fill_manual(values=c("#144272","#144272","#144272"))+
  scale_color_manual(values = c("black","black","black"))+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

compare_means(hsamiR1455p ~ scoreARVC, data = data_score, method = "wilcox.test")
#compare_means(hsamiR1455p ~ scoreARVC, data = data_score, method = "kruskal.test")

my_comparisons <- list( c("Low risk","Interm. risk"), c("Low risk","High risk"),
                        c("Interm. risk","High risk"))
p10 = p10 + stat_compare_means(comparisons = my_comparisons, size=5)
#stat_compare_means(label.y = 13, size=5) 
p10

p11 <- data_score %>%
  
  ggplot(aes(x=scoreARVC, hsamiR19a3p, fill = scoreARVC))+
  geom_boxplot(show.legend = F, outlier.shape = NA, alpha=0.3, width=0.6, coef=0)+
  geom_point(aes(colour=factor(scoreARVC), 
                 fill = factor(scoreARVC)), shape=21, size = 3, alpha = .7,
             show.legend = F,position=position_jitter(width=0.15))+
  guides(fill = "none")+
  stat_summary(fun = median, show.legend=F, geom='crossbar', width=0.6, 
               size=0.15, colour="black")+
  labs(title="hsa-miR19a-3p", x=NULL, y="Relative miRNA Expression")+ 
  scale_x_discrete(breaks = c("Low risk","Interm. risk","High risk"), 
                   labels = c(glue("Low risk
                                   (N={T0_n})"), glue("Interm. risk
                                   (N={T1_n})"), glue("High risk
                                   (N={T2_n})")))+
  scale_fill_manual(values=c("#144272","#144272","#144272"))+
  scale_color_manual(values = c("black","black","black"))+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

compare_means(hsamiR19a3p ~ scoreARVC, data = data_score, method = "wilcox.test")
#compare_means(hsamiR19a3p ~ scoreARVC, data = data_score, method = "kruskal.test")

my_comparisons <- list( c("Low risk","Interm. risk"), c("Low risk","High risk"),
                        c("Interm. risk","High risk"))
p11 = p11 + stat_compare_means(comparisons = my_comparisons, size=5)
#stat_compare_means(label.y = 13, size=5) 
p11

p12 <- data_score %>%
  
  ggplot(aes(x=scoreARVC, hsamiR29a3p, fill = scoreARVC))+
  geom_boxplot(show.legend = F, outlier.shape = NA, alpha=0.3, width=0.6, coef=0)+
  geom_point(aes(colour=factor(scoreARVC), 
                 fill = factor(scoreARVC)), shape=21, size = 3, alpha = .7,
             show.legend = F,position=position_jitter(width=0.15))+
  guides(fill = "none")+
  stat_summary(fun = median, show.legend=F, geom='crossbar', width=0.6, 
               size=0.15, colour="black")+
  labs(title="hsa-miR29a-3p", x=NULL, y="Relative miRNA Expression")+ 
  scale_x_discrete(breaks = c("Low risk","Interm. risk","High risk"), 
                   labels = c(glue("Low risk
                                   (N={T0_n})"), glue("Interm. risk
                                   (N={T1_n})"), glue("High risk
                                   (N={T2_n})")))+
  scale_fill_manual(values=c("#144272","#144272","#144272"))+
  scale_color_manual(values = c("black","black","black"))+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

compare_means(hsamiR29a3p ~ scoreARVC, data = data_score, method = "wilcox.test")
#compare_means(hsamiR29a3p ~ scoreARVC, data = data_score, method = "kruskal.test")

my_comparisons <- list( c("Low risk","Interm. risk"), c("Low risk","High risk"),
                        c("Interm. risk","High risk"))
p12 = p12 + stat_compare_means(comparisons = my_comparisons, size=5)
#stat_compare_means(label.y = 13, size=5) 
p12


