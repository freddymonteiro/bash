library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(viridis)

dat <- read.csv(file = "~/Desktop/PLOTS/assembly_sizes_NLRS_2.csv",
                header=TRUE,
                sep=",",
                colClasses= c('character','numeric','numeric','numeric','numeric',
                              'numeric','numeric','numeric','numeric','numeric',
                              'numeric','character','character','numeric'))
dat$Accessions <- factor(dat$Accessions, levels= dat$Accessions) #FIX dat1$Accessions category order

dat3 <- read.csv(file='~/Desktop/PLOTS/architectures-summary.csv', #uses processed matrix: Architectures counts, Arch ID + Pfam domains, Count/TOTAL NLR fraction
                 header=TRUE,
                 sep=',',
                 colClasses=c('character','character','numeric','character','numeric'))
dat3$Arch_name <- factor(dat3$Arch_name, levels= dat3$Arch_name) #FIX dat3$Archname category order
meltdat3 <- melt(dat3, 
                 id.vars="TOTAL",
                 measure.vars = c("Arch_count"))

meltdat3$TOTAL <- factor(meltdat3$TOTAL, levels= meltdat3$TOTAL) #FIX dat3$Archname category order

#PLOT6 Total number of domain architectures per accession 
plot6 <- ggplot(dat, aes(x=dat$Accessions, y=dat$uniarch, fill = dat$Collection)) +
  geom_bar(stat="identity") +
  geom_text(aes(y=dat$uniarch, label=dat$uniarch), nudge_x=0, nudge_y=2, color="black", size=3, angle=0) +
  scale_fill_manual(values=c("red2","mediumpurple3","gold","orange")) +
  theme(panel.background = element_rect(fill = "white", color = "white", size = 0.25)) +
  theme(panel.grid.major=element_blank()) +
  theme(panel.grid.minor=element_blank()) +
  theme(axis.line = element_line(size = 0.25, colour = "black")) +
  theme(axis.ticks = element_line(size = 0.25)) +
  theme(axis.ticks.length = unit(.2, "cm")) +
  theme(plot.title=element_text(color="black", size=24, vjust=0, hjust=0.025)) +
  theme(axis.text.x=element_text(size=14,color="black",angle=90, hjust =0.8, vjust=0.5)) +
  theme(axis.text.y=element_text(size=14,color="black")) +
  theme(axis.title.x=element_text(size=18,color="black", vjust=1.25)) +
  theme(axis.title.y=element_text(size=18,color="black", vjust=1.25)) +
  theme(legend.position = "none") +
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(size = 14, colour = "black")) +
  theme(legend.key = element_rect(fill = "white", colour = "black",size = 0.5)) +
  labs(title="Number of domain arrangements", x="Accessions", y="Total number of arrangements") +
  scale_y_continuous(expand = c(0,0), limits=c(0,55))
plot6

#PLOT7 Number of Genes Vs Different Architectures
plot7 <- ggplot(dat3, aes(x=dat3$Arch_name, y=dat3$TOTAL, fill=dat3$Fraction)) +
  geom_bar(stat="identity") +
  theme(panel.background = element_rect(fill = "white", color = "white", size = 0.25)) +
  theme(panel.grid.major=element_blank()) +
  theme(panel.grid.minor=element_blank()) +
  theme(axis.line = element_line(size = 0.25, colour = "black")) +
  theme(axis.ticks = element_line(size = 0.25)) +
  theme(axis.ticks.length = unit(.2, "cm")) +
  theme(plot.title=element_text(color="black", size=24, vjust=0, hjust=0.025)) +
  theme(axis.text.x=element_text(size=7,color="black",angle=90, hjust =0.8, vjust=0.5)) +
  theme(axis.text.y=element_text(size=14,color="black")) +
  theme(axis.title.x=element_text(size=18,color="black", vjust=1.25)) +
  theme(axis.title.y=element_text(size=18,color="black", vjust=1,25)) +
  theme(legend.position = "none") +
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(size = 14, colour = "black")) +
  theme(legend.key = element_rect(fill = "white", colour = "black",size = 0.5)) +
  labs(title="PanNLRs domain architectures (semilog)", x="Domain Arrangements/Architectures", y="Log10 number of genes") +
  scale_y_continuous(trans='log10') +
  scale_fill_hue(h = c(10, 10000))
#  scale_fill_viridis(discrete=TRUE, option ="magma", direction = -1, begin = 0,1, end = 0.92,0)
plot7

#PLOT8 Number of Domain architectures sharing same total number of genes
plot8 <- ggplot(meltdat3, aes(x=meltdat3$TOTAL, y=meltdat3$value, fill=meltdat3$variable)) +
  geom_bar(stat="identity") +
  theme(panel.background = element_rect(fill = "white", color = "white", size = 2)) +
  theme(panel.background = element_rect(fill = "white", color = "white", size = 0.25)) +
  theme(panel.grid.major=element_blank()) +
  theme(panel.grid.minor=element_blank()) +
  theme(axis.line = element_line(size = 0.25, colour = "black")) +
  theme(axis.ticks = element_line(size = 0.25)) +
  theme(axis.ticks.length = unit(.2, "cm")) +
  theme(plot.title=element_text(color="black", size=24, vjust=0, hjust=0.025)) +
  theme(axis.text.x=element_text(size=14,color="black",angle=90, hjust =0.8, vjust=0.5)) +
  theme(axis.text.y=element_text(size=14,color="black")) +
  theme(axis.title.x=element_text(size=18,color="black", vjust=1.25)) +
  theme(axis.title.y=element_text(size=18,color="black", vjust=1.25)) +
  theme(legend.position = "none") +
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(size = 14, colour = "black")) +
  theme(legend.key = element_rect(fill = "white", colour = "black",size = 0.5)) +
  labs(title="Number of Architectures sharing number of genes", x="Number of genes", y="Number of domain architectures") +
  scale_y_continuous(expand = c(0,0), limits=c(0,90))
plot8

grid.newpage()
grid.draw(rbind(ggplotGrob(plot7), ggplotGrob(plot6), ggplotGrob(plot8), size = "first"))