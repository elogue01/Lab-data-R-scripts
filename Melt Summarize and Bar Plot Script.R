EM <- read.delim("~/Desktop/Extracellular matrix remodeling cluster_AM RNAseq.txt", stringsAsFactors=FALSE) 
EM$Cohort = as.factor(EM$Cohort)
library(reshape2)
data_long <- melt(data=EM, id.var=c("Sample", "Cohort"),
                  measure.vars=c("ITGA6", "A2M", "SDC4", "MMP2", "MMP9", "MMP12", "COL1A1", 
                                 "COL23A1", "COL6A1", "COL6A2", "COL8A2", "CYP1B1", 
                                 "ITGAE", "ITGAX", "ITGA11", "TGFBR1"),
                  variable.name="Genes")
names(data_long)[names(data_long)=="value"] <- "Expression"

library(plyr)
cdata <- ddply(data_long, c("Cohort", "Genes"), summarise,
               N    = length(Expression),
               mean = mean(Expression),
               sd   = sd(Expression),
               se   = sd / sqrt(N)
)
cdata

library(ggplot2)
ggplot(cdata, aes(x=Genes, y=mean, fill=Cohort)) + 
        geom_bar(position=position_dodge(), stat="identity",
                 colour="black", # Use black outlines,
                 size=.3) +      # Thinner lines
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                      size=.3,    # Thinner lines
                      width=.2,
                      position=position_dodge(.9)) +
        coord_cartesian(ylim=c(0,28)) +
        xlab("Genes") +
        ylab("Relative Expression") +
        ggtitle("Extracellular Matrix Reorganizing Genes") +
        scale_y_continuous(breaks=0:20*4) +
        theme_classic()