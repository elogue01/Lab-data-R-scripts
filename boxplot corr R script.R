panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
        usr <- par("usr"); on.exit(par(usr)) 
        par(usr = c(0, 1, 0, 1)) 
        r <- abs(cor(x, y)) 
        txt <- format(c(r, 0.123456789), digits=digits)[1] 
        txt <- paste(prefix, txt, sep="") 
        if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
        
        test <- cor.test(x,y) 
        # borrowed from printCoefmat
        Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                         symbols = c("***", "**", "*", ".", " ")) 
        
        text(0.5, 0.5, txt, cex = cex * r) 
        text(.8, .8, Signif, cex=cex, col=2) 
}

library(ggplot2)

conc <- read.delim("~/Desktop/CHIT1 conc merge pos SM.txt", stringsAsFactors=FALSE) 
conc$Cohort = as.factor(conc$Cohort)

p <- ggplot(conc, aes(x=Cohort, y= Calc..CHIT1.Conc..BAL)) + 
        geom_boxplot() 
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
        #scale_y_log10()  +
        ## ylim(0,80) +
        labs(x="Cohort", y = "Relative CYP1B1 expression") +
        theme_classic()

chit1.conc = aov(formula = Calc..CHIT1.Conc..BAL ~ Cohort, data=conc)
summary(chit1.conc)
TukeyHSD(chit1.conc)


Cyp1b1.expression <- read.delim("~/Desktop/AM downreg bacterial defense genes.txt", stringsAsFactors=FALSE) 
Cyp1b1.expression$Cohort = as.factor(Cyp1b1.expression$Cohort)

# subset the HIV-NS data
neg_ns = Cyp1b1.expression[Cyp1b1.expression$Cohort  == "HIV-NS",]

# subset the HIV-SM data
neg_sm = Cyp1b1.expression[Cyp1b1.expression$Cohort  == "HIV-SM",]

# subset the HIV+NS data
pos_ns = Cyp1b1.expression[Cyp1b1.expression$Cohort  == "HIV+NS",]

# subset the HIV+SM data
pos_sm = Cyp1b1.expression[Cyp1b1.expression$Cohort  == "HIV+SM",]

# subset the HIV+SM ART data
pos_sm_art = Cyp1b1.expression[Cyp1b1.expression$Cohort  == "HIV+ SM ART",]

# subset the HIV+NS ART data
pos_ns_art = Cyp1b1.expression[Cyp1b1.expression$Cohort  == "HIV+ NS ART",]

# subset the HIV- COPD data
neg_copd = Cyp1b1.expression[Cyp1b1.expression$Cohort  == "HIV- COPD",]

# subset the HIV- COPD data
pos_sm_copd = Cyp1b1.expression[Cyp1b1.expression$Cohort  == "HIV+SM ART COPD",]

#Using ANOVA to do significance testing
chit1 = aov(formula = relative.CHIT1 ~ Cohort, data=Cyp1b1.expression)
summary(chit1)
TukeyHSD(chit1)

pla2g7 = aov(formula = relative.PLA2G7 ~ Cohort, data=Cyp1b1.expression)
summary(pla2g7)
TukeyHSD(pla2g7)

cyp1b1 = aov(formula = relative.CYP1B1 ~ Cohort, data=Cyp1b1.expression)
summary(cyp1b1)
TukeyHSD(cyp1b1)

mmp12 = aov(formula = relative.MMP12 ~ Cohort, data=Cyp1b1.expression)
summary(mmp12)
TukeyHSD(mmp12)

col6a2 = aov(formula = relative.COL6A2 ~ Cohort, data=Cyp1b1.expression)
summary(col6a2)
TukeyHSD(col6a2)

tgfbr1 = aov(formula = relative.TGFBR1 ~ Cohort, data=Cyp1b1.expression)
summary(tgfbr1)
TukeyHSD(tgfbr1)

#Using T tests to do significance testing
var.test(neg_ns$Calc..CHIT1.Conc..BAL,pos_sm$Calc..CHIT1.Conc..BAL)

t.test(neg_ns$Calc..CHIT1.Conc..BAL,pos_sm$Calc..CHIT1.Conc..BAL, alternative = "less", var.equal=FALSE, paired=FALSE)



var.test(pos_ns$Calc..CHIT1.Conc..BAL,pos_sm$Calc..CHIT1.Conc..BAL)

wilcox.test(pos_ns$Calc..CHIT1.Conc..BAL,pos_sm$Calc..CHIT1.Conc..BAL, alternative = "less", var.equal=TRUE, paired=FALSE)


var.test(neg_sm$Calc..CHIT1.Conc..BAL,pos_sm$Calc..CHIT1.Conc..BAL)

wilcox.test(neg_sm$Calc..CHIT1.Conc..BAL,pos_sm$Calc..CHIT1.Conc..BAL, alternative = "l", var.equal=TRUE, paired=FALSE)


var.test(neg_sm$Calc..CHIT1.Conc..BAL,pos_sm_art$Calc..CHIT1.Conc..BAL)

wilcox.test(neg_sm$Calc..CHIT1.Conc..BAL,pos_sm_art$Calc..CHIT1.Conc..BAL, alternative = "less", var.equal=FALSE, paired=FALSE)

pairs(Cyp1b1.expression[,3:16], lower.panel=panel.smooth, upper.panel=panel.cor)

cor.test(Cyp1b1.expression$Relative.Quantity , Cyp1b1.expression$MCP.1)

library(ggplot2)

conc <- read.delim("~/Desktop/CHIT1 conc 1st v 2nd bronch.txt", stringsAsFactors=FALSE) 
conc$Cohort = as.factor(conc$Cohort)
conc$Bronch = as.factor(conc$Bronch)

ggplot(data=conc, aes(x=Bronch, y=Calc..CHIT1.Conc..BAL, colour=Cohort, group=Sample)) + 
        geom_line() +
        geom_point() +
        labs(x="Bronchoscopy", y = "CHIT1 Conc. (ng/mL)") +
        theme_classic()

p <- ggplot(Cyp1b1.expression, aes(x=Cohort, y= relative.ITGA11)) + 
        geom_boxplot() 
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
       # scale_y_log10()  +
        ## ylim(0,80) +
        labs(title="ITGA11 relative expression",x="Cohort", y = "ITGA11 expression") +
        theme_classic()

#Plot linear regression lines

plot(Cyp1b1.expression$relative.CHIT1, Cyp1b1.expression$relative.PLA2G7,  type = "p", main = "CHIT1 vs PLA2G7 expression", xlab = "relative CHIT1 expression", ylab ="relative PLA2G7 expression")
abline(fit <- lm(relative.PLA2G7 ~ relative.CHIT1, data=Cyp1b1.expression), col="red") # regression line (y~x) 
legend("topleft", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)

plot(Cyp1b1.expression$relative.CHIT1, Cyp1b1.expression$relative.CYP1B1,  type = "p", main = "CHIT1 vs CYP1B1 expression", xlab = "relative CHIT1 expression", ylab ="relative CYP1B1 expression")
abline(fit <- lm(relative.CYP1B1 ~ relative.CHIT1, data=Cyp1b1.expression), col="red") # regression line (y~x) 
legend("topleft", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)

plot(Cyp1b1.expression$relative.CHIT1, Cyp1b1.expression$relative.MMP12,  type = "p", main = "CHIT1 vs MMP12 expression", xlab = "relative CHIT1 expression", ylab ="relative MMP12 expression")
abline(fit <- lm(relative.MMP12 ~ relative.CHIT1, data=Cyp1b1.expression), col="red") # regression line (y~x) 
legend("topleft", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)

plot(Cyp1b1.expression$relative.CHIT1, Cyp1b1.expression$relative.COL6A2,  type = "p", main = "CHIT1 vs COL6A2 expression", xlab = "relative CHIT1 expression", ylab ="relative COL6A2 expression")
abline(fit <- lm(relative.COL6A2 ~ relative.CHIT1, data=Cyp1b1.expression), col="red") # regression line (y~x) 
legend("topleft", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)

plot(Cyp1b1.expression$relative.CHIT1, Cyp1b1.expression$relative.ITGA11,  type = "p", main = "CHIT1 vs ITGA11 expression", xlab = "relative CHIT1 expression", ylab ="relative ITGA11 expression")
abline(fit <- lm(relative.ITGA11 ~ relative.CHIT1, data=Cyp1b1.expression), col="red") # regression line (y~x) 
legend("topleft", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)

plot(Cyp1b1.expression$relative.CHIT1, Cyp1b1.expression$relative.TGFBR1,  type = "p", main = "CHIT1 vs TGFBR1 expression", xlab = "relative CHIT1 expression", ylab ="relative TGFBR1 expression")
abline(fit <- lm(relative.TGFBR1 ~ relative.CHIT1, data=Cyp1b1.expression), col="red") # regression line (y~x) 
legend("topleft", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)


plot(Cyp1b1.expression$relative.CYP1B1, Cyp1b1.expression$relative.MMP12,  type = "p", main = "CYP1B1 vs MMP12 expression", xlab = "relative CYP1B1 expression", ylab ="relative MMP12 expression")
abline(fit <- lm(relative.MMP12 ~ relative.CYP1B1, data=Cyp1b1.expression), col="red") # regression line (y~x) 
legend("topleft", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)

plot(Cyp1b1.expression$relative.CYP1B1, Cyp1b1.expression$relative.COL6A2,  type = "p", main = "CYP1B1 vs COL6A2 expression", xlab = "relative CYP1B1 expression", ylab ="relative COL6A2 expression")
abline(fit <- lm(relative.COL6A2 ~ relative.CYP1B1, data=Cyp1b1.expression), col="red") # regression line (y~x) 
legend("topleft", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)

plot(Cyp1b1.expression$relative.CYP1B1, Cyp1b1.expression$relative.ITGA11,  type = "p", main = "CYP1B1 vs ITGA11 expression", xlab = "relative CYP1B1 expression", ylab ="relative ITGA11 expression")
abline(fit <- lm(relative.ITGA11 ~ relative.CYP1B1, data=Cyp1b1.expression), col="red") # regression line (y~x) 
legend("topleft", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)

plot(Cyp1b1.expression$relative.CYP1B1, Cyp1b1.expression$relative.TGFBR1,  type = "p", main = "CYP1B1 vs TGFBR1 expression", xlab = "relative CYP1B1 expression", ylab ="relative TGFBR1 expression")
abline(fit <- lm(relative.TGFBR1 ~ relative.CYP1B1, data=Cyp1b1.expression), col="red") # regression line (y~x) 
legend("topleft", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)


par(mfrow=c(2,2))
plot(fit)

seqvqPCR <- read.delim("~/Desktop/CHIT1 RNAseq v RTqPCR.txt", stringsAsFactors=FALSE) 
seqvqPCR $Cohort = as.factor(seqvqPCR $Cohort)

plot(seqvqPCR$RNAseq, seqvqPCR$RTqPCR,  type = "p", main = "CHIT1 expression: RNA-seq vs RT-qPCR", xlab = "CHIT1 RNA-seq expression", ylab ="CHIT1 RT-qPCR expression")
abline(fit <- lm(RTqPCR ~ RNAseq, data=seqvqPCR), col="red") # regression line (y~x) 
legend("topleft", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)

par(mfrow=c(2,2))
plot(fit)

seqvqPCR <- read.delim("~/Desktop/CYP1B1 RNAseq v RTqPCR.txt", stringsAsFactors=FALSE) 
seqvqPCR $Cohort = as.factor(seqvqPCR $Cohort)

plot(seqvqPCR$RNAseq, seqvqPCR$RTqPCR,  type = "p", main = "CYP1B1 expression: RNA-seq vs RT-qPCR", xlab = "CYP1B1 RNA-seq expression", ylab ="CYP1B1 RT-qPCR expression")
abline(fit <- lm(RTqPCR ~ RNAseq, data=seqvqPCR), col="red") # regression line (y~x) 
legend("topleft", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)

ggplot(Cyp1b1.expression, aes(x=bronch, y=Calc..CHIT1.Conc..BAL, group=sex, shape=sex)) + 
        geom_line(size=1.5) + 
        geom_point(size=3, fill="white") +
        scale_shape_manual(values=c(22,21))

Cyp1b1.expression <- read.delim("~/Desktop/CYP1B1 expression w ART.txt", stringsAsFactors=FALSE) 
Cyp1b1.expression$Cohort = as.factor(Cyp1b1.expression$Cohort)


p <- ggplot(Cyp1b1.expression, aes(x=Cohort, y= relative.qRT.PCR)) + 
        geom_boxplot() 
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
        scale_y_log10()  +
        ## ylim(0,80) +
        labs(title="qRT-PCR relative CYP1B1 expression",x="Cohort", y = "Relative CYP1B1 expression") +
        theme_classic()

p <- ggplot(Cyp1b1.expression, aes(x=Cohort, y= Calc..CHIT1.Conc..BAL)) + 
        geom_boxplot() 
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
        ##scale_y_log10()  +
        ## ylim(0,80) +
        labs(title="CHIT1 BAL Conc.",x="Cohort", y = "CHIT1 Conc. (ng/mL)") +
        theme_classic()

p <- ggplot(Cyp1b1.expression, aes(x=Cohort, y= Relative.CHIT1.expression)) + 
        geom_boxplot() 
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
        #scale_y_log10()  +
        ## ylim(0,80) +
        labs(title="Relative CHIT1 expression",x="Cohort", y = "Relative CHIT1 expression") +
        theme_classic()

p <- ggplot(Cyp1b1.expression, aes(x=Cohort, y= Normalized.CHIT1.Conc..BAL)) + 
        geom_boxplot() 
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
        #scale_y_log10()  +
        ## ylim(0,80) +
        labs(title="Normalized CHIT1 BAL Conc.",x="Cohort", y = "Normalized CHIT1 Conc. (ng/mL)") +
        theme_classic()

p <- ggplot(Cyp1b1.expression, aes(x=Cohort, y= Serum.conc.)) + 
        geom_boxplot() 
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
        scale_y_log10()  +
        ## ylim(0,80) +
        labs(title="Serum CHIT1 Conc.",x="Cohort", y = "CHIT1 Conc. (ng/mL)") +
        theme_classic()

