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

Cyp1b1.expression <- read.delim("~/Desktop/AM relative expression w smoking info.txt", stringsAsFactors=FALSE) 
Cyp1b1.expression$Cohort = as.factor(Cyp1b1.expression$Cohort)

pairs(Cyp1b1.expression[c(3,5,6,15,16,17,18,20,19,21,22,23,26,27)], lower.panel=panel.smooth, upper.panel=panel.cor)
cor(Cyp1b1.expression[c(3,5,6,15,16,17,18,20,19,21,22,23,26,27)])

il37 <- read.delim("~/Desktop/IL-37 corr.txt", stringsAsFactors=FALSE) 
il37$Cohort = as.factor(il37$Cohort)

pairs(il37[,3:16], lower.panel=panel.smooth, upper.panel=panel.cor)
cor(il37[,3:16])

cor.test(il37$relative.CAMP , il37$relative.IL37)

plot(il37$relative.CAMP, il37$relative.IL37,  type = "p", main = "CAMP vs IL37 expression", xlab = "relative CAMP expression", ylab ="relative IL37 expression")
abline(fit <- lm(relative.IL37 ~ relative.CAMP, data=il37), col="red") # regression line (y~x) 
legend("toprigh", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)


spir <- read.delim("~/Desktop/CHIT1 expression v spirometry.txt", stringsAsFactors=FALSE) 
spir$Cohort = as.factor(spir$Cohort)

pairs(spir[,3:16], lower.panel=panel.smooth, upper.panel=panel.cor)
cor(spir[,3:16])

plot(spir$Calc..CHIT1.Conc..BAL, spir$FEV1 ,  type = "p", main = "CHIT1 Conc IN BALF vs FEV1", xlab = "CHIT1 Conc. (ng/mL)", ylab ="FEV1")
abline(fit <- lm(FEV1  ~ Calc..CHIT1.Conc..BAL, data=spir), col="red") # regression line (y~x) 
legend("toprigh", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)

plot(spir$Calc..CHIT1.Conc..BAL, spir$X..FEV1 ,  type = "p", main = "CHIT1 Conc IN BALF vs %FEV1", xlab = "CHIT1 Conc. (ng/mL)", ylab ="%FEV1")
abline(fit <- lm(X..FEV1  ~ Calc..CHIT1.Conc..BAL, data=spir), col="red") # regression line (y~x) 
legend("toprigh", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)

plot(spir$Calc..CHIT1.Conc..BAL, spir$FVC,  type = "p", main = "CHIT1 Conc IN BALF vs FVC", xlab = "CHIT1 Conc. (ng/mL)", ylab ="FVC")
abline(fit <- lm(FVC  ~ Calc..CHIT1.Conc..BAL, data=spir), col="red") # regression line (y~x) 
legend("toprigh", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)