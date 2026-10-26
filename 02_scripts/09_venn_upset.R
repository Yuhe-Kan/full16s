# Fig1D
# ROOT DIRECTORY (to modify on your computer)
rm(list = ls())
getwd()
path.root <- "E:/b_myproject/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

ps = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))
ntaxa(ps)#查看otu总数

library(VennDiagram)
library(UpSetR)
library(MicrobiotaProcess)
library(ggplot2)

vennlist <- get_vennlist(obj=ps, factorNames="Group")
vennp <- venn.diagram(vennlist,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("#FDDED7", "#A28CC2","#CCD376","#5CB0C3"),
                      cat.col=c("#FDDED7", "#A28CC2","#CCD376","#5CB0C3"),
                      alpha = 0.85, 
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 1.4,
                      fontcolor = "black",
                      cat.cex = 1.7,
                      cat.just = list(c(0.5, -5), c(0.5, -5), c(0.5, 0.6), c(0.5, 1)),
                      cat.default.pos = "outer",
                      cat.dist=0.1,
                      cat.size = 5,
                      margin = 0, 
                      lwd = 4,
                      lty ='dotted',
                      imagetype = "tiff")
grid::grid.draw(vennp)
ggsave(file.path(path.plot,"05_venn_upset/05_venn.jpeg"), plot = vennp, width = 10, height = 8)

upsetda <- get_upset(obj=ps, factorNames="Group")
upset_plot <- upset(upsetda, sets=unique(as.vector(sample_data(ps)$Group)), 
      sets.bar.color = "#56B4E9",
      order.by = "freq", 
      empty.intersections = "on",
      #set_size.show = T,
      main.bar.color = c("black","black","black","black","#A28CC2","#A28CC2","#A28CC2","black","#A28CC2",
                         "black","black","black","black","black","black"))
jpeg(file.path(path.plot,"05_venn_upset/05_upset.jpeg"), width = 8, height = 4, units = "in", res = 300)
print(upset_plot)
dev.off()
#______________________________________________________________________________________________________________________________

