## Load libraries----
library(Seurat)
library(ggplot2)
library(dplyr)

## Load Seurat object----
load("results/objects/obj.Rdata")

## Cluster composition KO v WT----

# pull data out of meta.data
composition <- (prop.table(table(obj$annotation,
                                 obj$genotype),
                           margin = 2))*100
composition <- as.data.frame(composition)
colnames(composition)<- c("Cluster",
                          "Genotype",
                          "Percent_total")
composition <- composition %>% 
  mutate_if(is.numeric, round, digits = 2)

# plot
p <- ggplot(composition, aes(x = Cluster, y = Percent_total)) +
  geom_bar(
    aes(fill = Genotype),
    stat = "identity",
    position = position_dodge(0.8),
    width = 0.7) +
  ylab("% of total cells") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = c(0.85, 0.75),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title.x = element_blank())+
  ylim(0,30) +
  RotatedAxis() +
  scale_fill_manual("legend",
                    values = c("WT" = "#4878CD",
                               "KO" = "#D75438"))

# export
pdf(file="results/composition/cluster_composition.pdf", 
    width=3, 
    height=4,
    pointsize = 12,
    useDingbats = F)
p 
dev.off()

