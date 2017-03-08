
# Function for correlation in pairs plots
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

# Load data
msd_12 <- read.delim("~/Desktop/CHIT1 RNA v Protein v cytokine.txt", stringsAsFactors=TRUE) 
# msd_12$Relative.Conc. = as.numeric(msd_12$Relative.Conc.)
# remove relative conc,std dev and CV
# msd = msd_12[,1:8]

pairs(msd_12[,3:16], lower.panel=panel.smooth, upper.panel=panel.cor, main = "Cytokine levels for unstimulated AM")

cor(msd_12[,3:16])

cor.test(msd_12[,3], msd_12[,12])

cor.test(msd_12[,4], msd_12[,12])

cor.test(msd_12[,5], msd_12[,12])

# subset the unstimulated data
msd_unstim = msd_10wide[msd_10wide$Stimulation  == "Unstimulated",]

# pairwise plots and correlations for unstimulated data
pairs(msd_unstim[,4:15], lower.panel=panel.smooth, upper.panel=panel.cor, main = "Cytokine levels for unstimulated AM")

cor(msd_unstim[,4:15])

# make the table wide to prepare for pca 
library(reshape2)
msd_wide <- dcast(msd_12, Sample + Stimulation + Cohort + Pack.years + Viral.load ~ Assay, value.var="Calc..Conc..Mean")

# remove data for Eotaxin-3 and IL-13 
msd_10wide = msd_wide[, c(1:5, 7:10, 12:17)]

# remove data beta-glu and stool stimulations
msd_10wide = msd_10wide[msd_10wide$Stimulation %in% c("LPS", "LTA",  "PGN",  "PolyI:C",  "R848",  "recFLA",  "Unstimulated"),]

msd_10wide$Viral.load = as.numeric(msd_10wide$Viral.load)

# pairwise plots and correlations for all cohorts
pairs(msd_10wide[,4:15], lower.panel=panel.smooth, upper.panel=panel.cor, main = "Cytokine levels for stimulated AM")

cor(msd_10wide[,4:15])

library(Ez)
ezCor(msd_10wide[,4:15])

# subset the unstimulated data
msd_unstim = msd_10wide[msd_10wide$Stimulation  == "Unstimulated",]

# pairwise plots and correlations for unstimulated data
pairs(msd_unstim[,4:15], lower.panel=panel.smooth, upper.panel=panel.cor, main = "Cytokine levels for unstimulated AM")

cor(msd_unstim[,4:15])

cor.test(msd_unstim[,4], msd_unstim[,9])

# subset the LPS data
msd_lps = msd_10wide[msd_10wide$Stimulation  == "LPS",]

# pairwise plots and correlations for LPS data
pairs(msd_lps[,4:15], lower.panel=panel.smooth, upper.panel=panel.cor, main = "Cytokine levels for LPS stimulated AM")

cor(msd_lps[,4:15])

# subset the LTA data
msd_lta = msd_10wide[msd_10wide$Stimulation  == "LTA",]

# pairwise plots and correlations for LTA data
pairs(msd_lta[,4:15], lower.panel=panel.smooth, upper.panel=panel.cor, main = "Cytokine levels for LTA stimulated AM")

cor(msd_lps[,4:15])

# subset the R848 data
msd_r848 = msd_10wide[msd_10wide$Stimulation  == "R848",]

# pairwise plots and correlations for R848 data
pairs(msd_r848[,4:15], lower.panel=panel.smooth, upper.panel=panel.cor, main = "Cytokine levels for R848 stimulated AM")

cor(msd_r848[,4:15])

# subset the PGN data
msd_pgn = msd_10wide[msd_10wide$Stimulation  == "PGN",]

# pairwise plots and correlations for PGN data
pairs(msd_pgn[,4:15], lower.panel=panel.smooth, upper.panel=panel.cor, main = "Cytokine levels for PGN stimulated AM")

cor(msd_pgn[,4:15])

# subset the PolyI:C data
msd_polyic = msd_10wide[msd_10wide$Stimulation  == "PolyI:C",]

# pairwise plots and correlations for PolyI:C data
pairs(msd_polyic[,4:15], lower.panel=panel.smooth, upper.panel=panel.cor, main = "Cytokine levels for PolyI:C stimulated AM")

cor(msd_polyic[,4:15])

# subset the recFLA data
msd_recfla = msd_10wide[msd_10wide$Stimulation  == "recFLA",]

# pairwise plots and correlations for recFLA data
pairs(msd_recfla[,4:15], lower.panel=panel.smooth, upper.panel=panel.cor, main = "Cytokine levels for recFLA stimulated AM")

cor(msd_recfla[,4:15])

# log transform 
log.unstim <- log(msd_unstim[, 6:15])
log.unstim$Pack.years <- msd_unstim[, 4]
msd.cohort_unstim <- msd_unstim[, 3]


# apply PCA  scale = TRUE is highly 
# advisable, but default is FALSE. 
msd_unstim.pca <- prcomp(log.unstim,
                         center = TRUE,
                         scale. = TRUE) 

# print method
print(msd_unstim.pca)

# plot method
plot(msd_unstim.pca, type = "l")

# summary method
summary(msd_unstim.pca)

library(ggbiplot)
g <- ggbiplot(msd_unstim.pca, obs.scale = 1, var.scale = 1, 
              groups = msd.cohort_unstim, ellipse = TRUE, 
              circle = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'vertical', 
               legend.position = 'right') +
        labs(title="PolyI:C AM stimulation PCA")
print(g)



# subset the HIV-NS cohort
msd_norm = msd_10wide[msd_10wide$Cohort == "HIV-NS",]

# pairwise plots and correlations for HIV-NS cohorts
pairs(msd_norm[,4:13], lower.panel=panel.smooth, upper.panel=panel.cor, main = "Cytokine levels for stimulated HIV-NS AM")

cor(msd_norm[,4:13])


# log transform 
log.msd <- log(msd_norm[, 6:13])
msd.stim <- msd_norm[, 2]


# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
msd.pca <- prcomp(log.msd,
                 center = TRUE,
                 scale. = TRUE) 

# print method
print(msd.pca)

# plot method
plot(msd.pca, type = "l")

# summary method
summary(msd.pca)

library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(msd.pca, obs.scale = 1, var.scale = 1, 
              groups = msd.stim, ellipse = TRUE, 
              circle = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top') +
        labs(title="HIV-NS AM stimulation PCA")
print(g)

# subset the HIV-SM cohort
msd_negsm = msd_10wide[msd_10wide$Cohort == "HIV-SM",]

# pairwise plots and correlations for HIV-SM cohorts
pairs(msd_negsm[,4:13], upper.panel = NULL, main = "Cytokine levels for stimulated HIV-SM AM")

cor(msd_negsm[,4:13])

# subset the HIV+NS cohort
msd_posns = msd_10wide[msd_10wide$Cohort == "HIV+NS",]

# pairwise plots and correlations for HIV+NS cohorts
pairs(msd_posns[,4:13], upper.panel = NULL, main = "Cytokine levels for stimulated HIV+NS AM")

cor(msd_posns[,4:13])

# subset the HIV+SM cohort
msd_possm = msd_10wide[msd_10wide$Cohort == "HIV+SM",]

# pairwise plots and correlations for HIV+SM cohorts
pairs(msd_possm[,4:13], upper.panel = NULL, main = "Cytokine levels for stimulated HIV+SM AM")

cor(msd_possm[,4:13])

# log transform HIV+SM
log.msd_possm <- log(msd_possm[, 4:13])
msd.stim_possm <- msd_possm[, 2]


# apply PCA to HIV+SM- scale. = TRUE is highly 
# advisable, but default is FALSE. 
msd_possm.pca <- prcomp(log.msd_possm,
                  center = TRUE,
                  scale. = TRUE) 

# print method
print(msd_possm.pca)

# plot method
plot(msd_possm.pca, type = "l")

# summary method
summary(msd_possm.pca)

library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(msd_possm.pca, obs.scale = 1, var.scale = 1, 
              groups = msd.stim_possm, ellipse = TRUE, 
              circle = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'vertical', 
               legend.position = 'right') +
        labs(title="HIV+SM AM stimulation PCA")
print(g)

# subset the PolyI:C stimulated data
msd_polyic = msd_10wide[msd_10wide$Stimulation  == "PolyI:C",]

# pairwise plots and correlations for PolyI:C stimulated data
pairs(msd_polyic[,4:13], upper.panel = NULL, main = "Cytokine levels for PolyI:C stimulated AM")

cor(msd_polyic[,4:13])

cor.test(msd_polyic[,8], msd_polyic[,11])

# log transform HIV+SM
log.msd_polyic <- log(msd_polyic[, 4:13])
msd.cohort_polyic <- msd_polyic[, 3]


# apply PCA to HIV+SM- scale. = TRUE is highly 
# advisable, but default is FALSE. 
msd_polyic.pca <- prcomp(log.msd_polyic,
                        center = TRUE,
                        scale. = TRUE) 

# print method
print(msd_polyic.pca)

# plot method
plot(msd_polyic.pca, type = "l")

# summary method
summary(msd_polyic.pca)

library(ggbiplot)
g <- ggbiplot(msd_polyic.pca, obs.scale = 1, var.scale = 1, 
              groups = msd.cohort_polyic, ellipse = TRUE, 
              circle = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'vertical', 
               legend.position = 'right') +
        labs(title="PolyI:C AM stimulation PCA")
print(g)
