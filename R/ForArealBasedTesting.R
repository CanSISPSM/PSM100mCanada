# This R script is used to conduct Pearson correlation analysis and graphing
# after effort on skewness testing and normalization on soil survey polygon-based
# soil property data from three of the independent statistical validation sites.
# Input: CSV files stored in the data folder of this R porject.
# Output: various scatter-stats-plots for further use.
# Author: Xiaoyuan Geng
# Copyright (c) 2025, Agriculture and Agri-Food Canada.
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sub-license,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# Disclaimer:
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
# 
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggstatsplot)
library(correlation)
library(epiR)
library(outliers)
# Define a function for maximum and mininum based normalization
min_max_normalization <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
# Define a function for Z score normalization
z_score_normalization <- function(x) {
  return ((x - mean(x)) / sd(x))
}
getwd()
# data folder is relative of the R project workspace here
subSLSOC<-read.csv("data/SLSOC_v2.csv",header=TRUE)
subSLBD<-read.csv("data/SLBD_v2.csv",header=TRUE)
subSLSand<-read.csv("data/SLSand_v2.csv",header=TRUE)
subSLClay<-read.csv("data/SLClay_v2.csv",header=TRUE)
str(subSLSOC)
subSLSOC <- subSLSOC %>% select(X0_30cm,SOC10m,SOC100m, SOC250m,SOC250mv2)
subSLSOC[subSLSOC <= 0]<-NA
nrow(subSLSOC)
subSLSOC<-subSLSOC[complete.cases(subSLSOC), ]
nrow(subSLSOC)
subSLSOCSWTest<-subSLSOC %>% shapiro_test(X0_30cm,SOC10m,SOC100m, SOC250m,SOC250mv2)
subSLSOCSWTest
outlier(subSLSOC$X0_30cm)
outlier(subSLSOC$SOC10m)
ggdensity(subSLSOC$X0_30cm, fill = "lightgray")
subSLSOC$X0_30cm<-log(subSLSOC$X0_30cm)
ggdensity(subSLSOC$X0_30cm, fill = "lightgray")
shapiro_test(subSLSOC$X0_30cm)
ggdensity(subSLSOC$SOC10m, fill = "lightgray")
subSLSOC$SOC10m<-log(subSLSOC$SOC10m)
ggdensity(subSLSOC$SOC10m, fill = "lightgray")
ggdensity(subSLSOC$SOC100m, fill = "lightgray")
subSLSOC$SOC100m<-log(subSLSOC$SOC100m)
ggdensity(subSLSOC$SOC250m, fill = "lightgray")
subSLSOC$SOC250m<-log(subSLSOC$SOC250m)
ggdensity(subSLSOC$SOC250m, fill = "lightgray")
subSLSOC$SOC250mv2<-log(subSLSOC$SOC250mv2)
ggdensity(subSLSOC$SOC250mv2, fill = "lightgray")
# Generate scatter stats plots based on transformed data
ggscatterstats(
  data = subSLSOC,
  x = SOC100m,
  y = X0_30cm,
  xlab = "Log(SOC%) 100m grids",
  ylab = "Log(SOC%) soil survey",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = subSLSOC,
  x = SOC100m,
  y = SOC10m,
  xlab = "Log(SOC%) 100m grids",
  ylab = "Log(SOC%) 10m grids",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = subSLSOC,
  x = X0_30cm,
  y = SOC250mv2,
  xlab = "Log(SOC%) 250m grids",
  ylab = "Log(SOC%) soil survey",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = subSLSOC,
  x = SOC100m,
  y = SOC250mv2,
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
shapiro_test(subSLSOC$SOC10m)
subSLSOCcor <- subSLSOC %>% select(X0_30cm,SOC10m,SOC100m, SOC250m, SOC250mv2)
# correlation for all variables
round(cor(subSLSOCcor), digits = 2)
subSLSOCccc <- epiR::epi.ccc(subSLSOC$SOC10m, subSLSOC$SOC100m, ci="z-transform", conf.level=0.95)
subSLSOCccc.result <- subSLSOCccc$rho.c
subSLSOCccc.result
subSLSOCcorall <- subSLSOC %>% select(X0_30cm,SOC10m,SOC100m, SOC250mv2)
subSLSOCCorTab<-correlation::correlation(subSLSOCcorall,
                          include_factors = TRUE, method = "auto"
 )
subSLSOCCorTab
subSLSOCCorTab[1,1]
str(subSLBD)
subSLBD <- subSLBD %>% select(X0_30cm,BD10m30,BD100m30, BD250m30,BD250mv2)
subSLBD[subSLBD<= 0]<-NA
nrow(subSLBD)
outlier(subSLBD$BD10m30)
outlier(subSLBD$BD100m30)
subSLBD<-subSLBD[complete.cases(subSLBD), ]
nrow(subSLBD)
subSLBDSWTest<-subSLBD %>% shapiro_test(X0_30cm,BD10m30,BD100m30, BD250m30)
subSLBDSWTest
subSLBD$X0_30cm<-log(subSLBD$X0_30cm)
shapiro_test(subSLBD$X0_30cm)
ggdensity(subSLBD$BD10m30, fill = "lightgray")
subSLBD$BD10m30<-log(subSLBD$BD10m30)
ggdensity(subSLBD$BD10m30, fill = "lightgray")
ggdensity(subSLBD$BD100m30, fill = "lightgray")
subSLBD$BD100m30<-log(subSLBD$BD100m30)
ggdensity(subSLBD$BD250mv2, fill = "lightgray")
subSLBD$BD250mv2<-log(subSLBD$BD250mv2)
ggdensity(subSLBD$BD250mv2, fill = "lightgray")
ggscatterstats(
  data = subSLBD,
  x = BD100m30,
  y = X0_30cm,
  xlab = "Log(BD g/cm3) 100m grids",
  ylab = "Log(BD g/cm3) soil survey",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = subSLBD,
  x = BD100m30,
  y = BD10m30,
  xlab = "Log(BD g/cm3) 100m grids",
  ylab = "Log(BD g/cm3) 10m grids",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)

ggscatterstats(
  data = subSLBD,
  x = BD100m30,
  y = BD250mv2,
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
subSLBDcor <- subSLBD %>% select(X0_30cm,BD10m30,BD100m30, BD250m30,BD250mv2)
# correlation for all variables
round(cor(subSLBDcor), digits = 2)
subSLBDccc <- epiR::epi.ccc(subSLBD$BD10m30, subSLBD$BD100m30, ci="z-transform", conf.level=0.95)
subSLBDccc.result <- subSLBDccc$rho.c
subSLBDccc.result
subSLBDcorall <- subSLBD %>% select(X0_30cm,BD10m30,BD100m30, BD250m30,BD250mv2)
correlation::correlation(subSLBDcorall,
                         include_factors = TRUE, method = "auto"
)
str(subSLSand)
subSLSand <- subSLSand %>% select(X0_30_cm,Sand10m30,Sand100m30, Sand250m30,Sand250mv2)
subSLSand[subSLSand<= 0]<-NA
nrow(subSLSand)
subSLSand<-subSLSand[complete.cases(subSLSand), ]
nrow(subSLSand)
subSLSandSWTest <- subSLSand %>% shapiro_test(X0_30_cm,Sand10m30,Sand100m30, Sand250m30, Sand250mv2)
subSLSandSWTest
subSLSand$X0_30_cm<-log(subSLSand$X0_30_cm)
ggdensity(subSLSand$X0_30_cm, fill = "lightgray")
shapiro_test(subSLSand$X0_30_cm)
ggdensity(subSLSand$Sand10m30, fill = "lightgray")
subSLSand$Sand10m30<-log(subSLSand$Sand10m30)
ggdensity(subSLSand$Sand10m30, fill = "lightgray")
ggdensity(subSLSand$Sand100m30, fill = "lightgray")
subSLSand$Sand100m30<-log(subSLSand$Sand100m30)
ggdensity(subSLSand$Sand100m30, fill = "lightgray")
ggdensity(subSLSand$Sand250mv2, fill = "lightgray")
subSLSand$Sand250mv2<-log(subSLSand$Sand250mv2)
ggdensity(subSLSand$Sand250mv2, fill = "lightgray")

ggscatterstats(
  data = subSLSand,
  x = Sand100m30,
  y = X0_30_cm,
  xlab = "Log(Sand%) 100m grids",
  ylab = "Log(Sand%) soil survey",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = subSLSand,
  x = Sand100m30,
  y = Sand10m30,
  xlab = "Log(Sand%) 100m grids",
  ylab = "Log(Sand%) 10m grids",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)

ggscatterstats(
  data = subSLSand,
  x = Sand100m30,
  y = Sand250mv2,
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
subSLSandcor <- subSLSand %>% select(X0_30_cm,Sand10m30,Sand100m30, Sand250m30, Sand250mv2)
# correlation for all variables
round(cor(subSLSandcor), digits = 2)
subSLSandccc <- epiR::epi.ccc(subSLSand$Sand10m30, subSLSand$Sand100m30, ci="z-transform", conf.level=0.95)
subSLSandccc.result <- subSLSandccc$rho.c
subSLSandccc.result
subSLSandcorall <- subSLSand %>% select(X0_30_cm,Sand10m30,Sand100m30, Sand250m30, Sand250mv2)
correlation::correlation(subSLSandcorall,
                         include_factors = TRUE, method = "auto"
)
str(subSLClay)
subSLClay <- subSLClay %>% select(X0_30_cm,Clay10m30,Clay100m30, Clay250m30,Clay250mv2)

subSLClay[subSLClay<= 0]<-NA
nrow(subSLClay)
outlier(subSLClay$Clay10m30)
outlier(subSLClay$Clay100m30)
subSLClay<-subSLClay[complete.cases(subSLClay), ]
nrow(subSLClay)
subSLClaySWTest <- subSLClay %>% shapiro_test(X0_30_cm,Clay10m30,Clay100m30, Clay250m30)
subSLClaySWTest
subSLClay$X0_30_cm<-log(subSLClay$X0_30_cm)
ggdensity(subSLClay$X0_30_cm, fill = "lightgray")
shapiro_test(subSLClay$X0_30_cm)
ggdensity(subSLClay$Clay10m30, fill = "lightgray")
subSLClay$Clay10m30<-log(subSLClay$Clay10m30)
ggdensity(subSLClay$Clay10m30, fill = "lightgray")
ggdensity(subSLClay$Clay100m30, fill = "lightgray")
subSLClay$Clay100m30<-log(subSLClay$Clay100m30)
ggdensity(subSLClay$Clay100m30, fill = "lightgray")
ggdensity(subSLClay$Clay250mv2, fill = "lightgray")
subSLClay$Clay250mv2<-log(subSLClay$Clay250mv2)
ggdensity(subSLClay$Clay250mv2, fill = "lightgray")
ggscatterstats(
  data = subSLClay,
  x = Clay100m30,
  y = X0_30_cm,
  xlab = "Log(Clay%) 100m grids",
  ylab = "Log(Clay%) soil survey",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = subSLClay,
  x = Clay100m30,
  y = Clay10m30,
  xlab = "Log(Clay%) 100m grids",
  ylab = "Log(Clay%) 10m grids",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
SLClayPlot <- ggscatterstats(
  data = subSLClay,
  x = Clay250mv2,
  y = X0_30_cm,
  xlab = "Log(Clay%) 250m grids",
  ylab = "Log(Clay%) soil survey",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
SLClayPlot
extract_stats(SLClayPlot)

ggscatterstats(
  data = subSLClay,
  x = Clay100m30,
  y = Clay250mv2,
  xlab = "Log(Clay%) 100m grids",
  ylab = "Log(Clay%) 250m grids",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
subSLClaycor <- subSLClay %>% select(X0_30_cm,Clay10m30,Clay100m30, Clay250m30, Clay250mv2)
# correlation for all variables
round(cor(subSLClaycor), digits = 2)
subSLClayccc <- epiR::epi.ccc(subSLClay$Clay10m30, subSLClay$Clay100m30, ci="z-transform", conf.level=0.95)
subSLClayccc.result <- subSLClayccc$rho.c
subSLClayccc.result
subSLClaycorall <- subSLClay %>% select(X0_30_cm,Clay10m30,Clay100m30, Clay250m30, Clay250mv2)
correlation::correlation(subSLClaycorall,
                         include_factors = TRUE, method = "auto"
)
# For the site of Bread Albane
BASOC<-read.csv("data/BASOC_v2.csv",header=TRUE)
BABD<-read.csv("data/BABD_v2.csv",header=TRUE)
BASand<-read.csv("data/BASand_v2.csv",header=TRUE)
BAClay<-read.csv("data/BAClay_v2.csv",header=TRUE)
str(BASOC)
BASOC<- BASOC %>% select(X0_30cm,SOC10m30,SOC100m30, SOC250m30,SOC250mv2)

BASOC[BASOC<= 0]<-NA
nrow(BASOC)
BASOC<-BASOC[complete.cases(BASOC), ]
nrow(BASOC)
BASOCSWTest<-BASOC %>% shapiro_test(X0_30cm,SOC10m30,SOC100m30, SOC250m30, SOC250mv2)
BASOCSWTest
ggdensity(BASOC$X0_30cm, fill = "lightgray")
BASOC$X0_30cm<-log(BASOC$X0_30cm)
ggdensity(BASOC$X0_30cm, fill = "lightgray")
shapiro_test(BASOC$X0_30cm)
ggdensity(BASOC$SOC10m30, fill = "lightgray")
BASOC$SOC10m30<-log(BASOC$SOC10m30)
ggdensity(BASOC$SOC10m30, fill = "lightgray")
shapiro_test(BASOC$SOC10m30)
ggdensity(BASOC$SOC100m30, fill = "lightgray")
BASOC$SOC100m30<-log(BASOC$SOC100m30)
ggdensity(BASOC$SOC100m30, fill = "lightgray")
ggdensity(BASOC$SOC250mv2, fill = "lightgray")
BASOC$SOC250mv2<-log(BASOC$SOC250mv2)
ggdensity(BASOC$SOC250mv2, fill = "lightgray")

ggscatterstats(
  data = BASOC,
  x = SOC100m30,
  y = X0_30cm,
  xlab = "Log(SOC%) 100m grids",
  ylab = "Log(SOC%) soil survey",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = BASOC,
  x = SOC100m30,
  y = SOC250mv2,
  xlab = "Log(SOC%) 100m grids",
  ylab = "Log(SOC%) 250m grids",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
BASOCcor <- BASOC %>% select(X0_30cm,SOC10m30,SOC100m30, SOC250m30,SOC250mv2)
# correlation for all variables
round(cor(BASOCcor), digits = 2)
BASOCccc <- epiR::epi.ccc(BASOC$SOC10m30, BASOC$SOC100m30, ci="z-transform", conf.level=0.95)
BASOCccc.result <- BASOCccc$rho.c
BASOCccc.result
BASOCcorall <- BASOC %>% select(X0_30cm,SOC10m30,SOC100m30, SOC250m30,SOC250mv2)
correlation::correlation(BASOCcorall,
                         include_factors = TRUE, method = "auto"
)
## Bulk density case
str(BABD)
BABD <- BABD %>% select(X0_30cm,BD10m30,BD100m30, BD250m30,BD250mv2)

BABD[BABD<= 0]<-NA
nrow(BABD)
BABD<-BABD[complete.cases(BABD), ]
nrow(BABD)
BABDSWTest<-BABD %>% shapiro_test(X0_30cm,BD10m30,BD100m30, BD250m30, BD250mv2)
BABDSWTest
BABD$X0_30cm<-log(BABD$X0_30cm)
ggdensity(BABD$X0_30cm, fill = "lightgray")
shapiro_test(BABD$X0_30cm)
ggdensity(BABD$BD10m30, fill = "lightgray")
BABD$BD10m30<-log(BABD$BD10m30)
ggdensity(BABD$BD10m30, fill = "lightgray")
shapiro_test(BABD$BD10m30)
ggdensity(BABD$BD100m30, fill = "lightgray")
BABD$BD100m30<-log(BABD$BD100m30)
ggdensity(BABD$BD100m30, fill = "lightgray")
shapiro_test(BABD$BD100m30)
ggdensity(BABD$BD250mv2, fill = "lightgray")
BABD$BD250mv2<-log(BABD$BD250mv2)
ggdensity(BABD$BD250mv2, fill = "lightgray")

ggscatterstats(
  data = BABD,
  x = BD100m30,
  y = X0_30cm,
  xlab = "Log(BD g/cm3) 100m grids",
  ylab = "Log(BD g/cm3) soil survey",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = BABD,
  x = BD100m30,
  y = BD250mv2,
  xlab = "Log(BD g/cm3) 100m grids",
  ylab = "Log(BD g/cm3) 250m grids",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
BABDcor <- BABD %>% select(X0_30cm,BD10m30,BD100m30, BD250m30, BD250mv2)
# correlation for all variables
round(cor(BABDcor), digits = 2)
BABDccc <- epiR::epi.ccc(BABD$BD10m30, BABD$BD100m30, ci="z-transform", conf.level=0.95)
BABDccc.result <- BABDccc$rho.c
BABDccc.result
BABDcorall <- BABD %>% select(X0_30cm,BD10m30,BD100m30, BD250m30, BD250mv2)
correlation::correlation(BABDcorall,
                         include_factors = TRUE, method = "auto"
)
str(BASand)
BASand <- BASand %>% select(X0_30cm,Sand10m30,Sand100m30, Sand250m30,Sand250mv2)
BASand[BASand<= 0]<-NA
nrow(BASand)
BASand<-BASand[complete.cases(BASand), ]
nrow(BASand)
BASandSWTest <- BASand %>% shapiro_test(X0_30cm,Sand10m30,Sand100m30, Sand250m30, Sand250mv2)
BASandSWTest
BASand$X0_30_cm<-log(BASand$X0_30cm)
shapiro_test(BASand$X0_30cm)
ggdensity(BASand$X0_30cm, fill = "lightgray")
ggdensity(BASand$Sand10m30, fill = "lightgray")
BASand$Sand10m30<-log(BASand$Sand10m30)
ggdensity(BASand$Sand10m30, fill = "lightgray")
shapiro_test(BASand$Sand10m30)
ggdensity(BASand$Sand100m30, fill = "lightgray")
BASand$Sand100m30<-log(BASand$Sand100m30)
ggdensity(BASand$Sand100m30, fill = "lightgray")
ggdensity(BASand$Sand250mv2, fill = "lightgray")
BASand$Sand250mv2<-log(BASand$Sand250mv2)
ggdensity(BASand$Sand250mv2, fill = "lightgray")
ggscatterstats(
  data = BASand,
  x = Sand100m30,
  y = X0_30cm,
  xlab = "Log(Sand %) 100m grids",
  ylab = "Log(Sand %) soil survey",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = BASand,
  x = Sand100m30,
  y = Sand250mv2,
  xlab = "Log(Sand %) 100m grids",
  ylab = "Log(Sand %) 250m grids",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
BASandcor <- BASand %>% select(X0_30cm,Sand10m30,Sand100m30, Sand250m30, Sand250mv2)
# correlation for all variables
round(cor(BASandcor), digits = 2)
BASandccc <- epiR::epi.ccc(BASand$Sand10m30, BASand$Sand100m30, ci="z-transform", conf.level=0.95)
BASandccc.result <- BASandccc$rho.c
BASandccc.result
BASandcorall <- BASand %>% select(X0_30_cm,Sand10m30,Sand100m30, Sand250m30, Sand250mv2)
correlation::correlation(BASandcorall,
                         include_factors = TRUE, method = "auto"
)
str(BAClay)
BAClay <- BAClay %>% select(X0_30cm,Clay10m30,Clay100m30, Clay250m30,Clay250mv2)
BAClay[BAClay<= 0]<-NA
nrow(BAClay)
BAClay<-BAClay[complete.cases(BAClay), ]
nrow(BAClay)
BAClaySWTest<-BAClay %>% shapiro_test(X0_30cm,Clay10m30,Clay100m30, Clay250m30)
BAClaySWTest
BAClay$X0_30_cm<-log(BAClay$X0_30cm)
BAClay$Clay10m30<-log(BAClay$Clay10m30)
shapiro_test(BAClay$X0_30cm)
ggdensity(BAClay$X0_30cm, fill = "lightgray")
shapiro_test(BAClay$Clay10m30)
ggdensity(BAClay$Clay10m30, fill = "lightgray")
BAClay$Clay10m30<-log(BAClay$Clay10m30)
ggdensity(BAClay$Clay10m30, fill = "lightgray")
shapiro_test(BAClay$Clay10m30)
ggdensity(BAClay$Clay100m30, fill = "lightgray")
BAClay$Clay100m30<-log(BAClay$Clay100m30)
ggdensity(BAClay$Clay100m30, fill = "lightgray")
shapiro_test(BAClay$Clay100m30)
ggscatterstats(
  data = BAClay,
  x = Clay100m30,
  y = X0_30_cm,
  xlab = "Log(Clay %) 100m grids",
  ylab = "Log(Clay %) soil survey",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = BAClay,
  x = Clay100m30,
  y = Clay250mv2,
  xlab = "Log(Sand %) 100m grids",
  ylab = "Log(Sand %) 250m grids",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
BAClaycor <- BAClay %>% select(X0_30cm,Clay10m30,Clay100m30, Clay250m30, Clay250mv2)
# correlation for all variables
round(cor(BAClaycor), digits = 2)
BAClayccc <- epiR::epi.ccc(BAClay$Clay10m30, BAClay$Clay100m30, ci="z-transform", conf.level=0.95)
BAClayccc.result <- BAClayccc$rho.c
BAClayccc.result
BAClaycorall <- BAClay %>% select(X0_30cm,Clay10m30,Clay100m30, Clay250m30, Clay250mv2)
correlation::correlation(BAClaycorall,
                         include_factors = TRUE, method = "auto"
)
# For the site of Grassland National Park
WBGNPSOC<-read.csv("data/GNPSOC_v2.csv",header=TRUE)
WBGNPBD<-read.csv("data/GNPBD_v2.csv",header=TRUE)
WBGNPSand<-read.csv("data/GNPSand_v2.csv",header=TRUE)
WBGNPClay<-read.csv("data/GNPClay_v2.csv",header=TRUE)

str(WBGNPSOC)
WBGNPSOC <- WBGNPSOC %>% select(X0_30cm,SOC50m30,SOC100m30, SOC250m30,SOC250mv2)
WBGNPSOC[WBGNPSOC<= 0]<-NA
nrow(WBGNPSOC)
WBGNPSOC<-WBGNPSOC[complete.cases(WBGNPSOC), ]
nrow(WBGNPSOC)
WBGNPSOCSWTest<-WBGNPSOC %>% shapiro_test(X0_30cm,SOC50m30,SOC100m30, SOC250m30)
WBGNPSOCSWTest
# outlier(WBGNPSOC$X0_30cm)
# outlier(WBGNPSOC$SOC10m)
ggdensity(WBGNPSOC$X0_30cm, fill = "lightgray")
WBGNPSOC$X0_30_cm<-log(WBGNPSOC$X0_30cm)
ggdensity(WBGNPSOC$X0_30cm, fill = "lightgray")
shapiro_test(WBGNPSOC$X0_30cm)
ggdensity(WBGNPSOC$SOC50m30, fill = "lightgray")
WBGNPSOC$SOC50m30<-log(WBGNPSOC$SOC50m30)
ggdensity(WBGNPSOC$SOC50m30, fill = "lightgray")
ggdensity(WBGNPSOC$SOC100m30, fill = "lightgray")
WBGNPBD$BD100m30<-log(WBGNPBD$BD100m30)
ggdensity(WBGNPSOC$SOC100m30, fill = "lightgray")
# WBGNPSOC$SOC50m<-log(WBGNPSOC$SOC50m30)
# ggdensity(WBGNPSOC$SOC50m30, fill = "lightgray")
shapiro_test(WBGNPSOC$SOC50m30)
ggscatterstats(
  data = WBGNPSOC,
  x = SOC100m30,
  y = X0_30cm,
  xlab = "Log(SOC %) 100m grids",
  ylab = "Log(SOC %) soil survey",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = WBGNPSOC,
  x = SOC100m30,
  y = SOC250mv2,
  xlab = "Log(SOC %) 100m grids",
  ylab = "Log(SOC %) 250m grids",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
WBGNPSOCcor <- WBGNPSOC %>% select(X0_30cm,SOC50m30,SOC100m30, SOC250m30, SOC250mv2)
# correlation for all variables
round(cor(WBGNPSOCcor), digits = 2)
WBGNPSOCccc <- epiR::epi.ccc(WBGNPSOC$SOC50m30, WBGNPSOC$SOC100m30, ci="z-transform", conf.level=0.95)
WBGNPSOCccc.result <- WBGNPSOCccc$rho.c
WBGNPSOCccc.result
WBGNPSOCcorall <- WBGNPSOC %>% select(X0_30_cm,SOC50m30,SOC100m30, SOC250m30, SOC250mv2)
correlation::correlation(WBGNPSOCcorall,
                         include_factors = TRUE, method = "auto"
)
str(WBGNPBD)
WBGNPBD <- WBGNPBD %>% select(X0_30cm,BD50m30,BD100m30, BD250m30,BD250mv2)
WBGNPBD[WBGNPBD<= 0]<-NA
nrow(WBGNPBD)
WBGNPBD<-WBGNPBD[complete.cases(WBGNPBD), ]
nrow(WBGNPBD)
WBGNPBDSWTest<-WBGNPBD %>% shapiro_test(X0_30cm,BD50m30,BD100m30, BD250m30)
WBGNPBDSWTest
WBGNPBD$X0_30cm<-log(WBGNPBD$X0_30cm)
ggdensity(WBGNPBD$X0_30cm, fill = "lightgray")
shapiro_test(WBGNPBD$X0_30cm)
ggdensity(WBGNPBD$BD50m30, fill = "lightgray")
WBGNPBD$BD50m30<-log(WBGNPBD$BD50m30)
ggdensity(WBGNPBD$BD50m30, fill = "lightgray")
shapiro_test(WBGNPBD$BD50m30)
ggdensity(WBGNPBD$BD100m30, fill = "lightgray")
WBGNPBD$BD100m30<-log(WBGNPBD$BD100m30)
ggdensity(WBGNPBD$BD100m30, fill = "lightgray")
ggscatterstats(
  data = WBGNPBD,
  x = BD100m30,
  y = X0_30cm,
  xlab = "Log(BD g/cm3) 100m grids",
  ylab = "Log(BD g/cm3) soil survey",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = WBGNPBD,
  x = BD100m30,
  y = BD250mv2,
  xlab = "Log(BD g/cm3) 100m grids",
  ylab = "Log(BD g/cm3) 250m grids",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
WBGNPBDcor <- WBGNPBD %>% select(X0_30cm,BD50m30,BD100m30, BD250m30, BD250mv2)
# correlation for all variables
round(cor(WBGNPBDcor), digits = 2)
WBGNPBDccc <- epiR::epi.ccc(WBGNPBD$BD50m30, WBGNPBD$BD100m30, ci="z-transform", conf.level=0.95)
WBGNPBDccc.result <- WBGNPBDccc$rho.c
WBGNPBDccc.result
WBGNPBDcorall <- WBGNPBD %>% select(X0_30cm,BD50m30,BD100m30, BD250m30, BD250mv2)
correlation::correlation(WBGNPBDcorall,
                         include_factors = TRUE, method = "auto"
)
str(WBGNPSand)
WBGNPSand <- WBGNPSand %>% select(X0_30cm,Sand50m30,Sand100m30, Sand250m30,Sand250mv2)
WBGNPSand[WBGNPSand<= 0]<-NA
nrow(WBGNPSand)
WBGNPSand<-WBGNPSand[complete.cases(WBGNPSand), ]
nrow(WBGNPSand)
WBGNPSandSWTest <- WBGNPSand %>% shapiro_test(X0_30cm,Sand50m30,Sand100m30, Sand250m30, Sand250mv2)
WBGNPSandSWTest
WBGNPSand$X0_30_cm<-log(WBGNPSand$X0_30cm)
shapiro_test(WBGNPSand$X0_30cm)
ggdensity(WBGNPSand$X0_30cm, fill = "lightgray")
ggdensity(WBGNPSand$Sand50m30, fill = "lightgray")
WBGNPSand$Sand50m30<-log(WBGNPSand$Sand50m30)
ggdensity(WBGNPSand$Sand50m30, fill = "lightgray")
shapiro_test(WBGNPSand$Sand50m30)
ggdensity(WBGNPSand$Sand100m30, fill = "lightgray")
WBGNPSand$Sand100m30<-log(WBGNPSand$Sand100m30)
ggdensity(WBGNPSand$Sand100m30, fill = "lightgray")
ggscatterstats(
  data = WBGNPSand,
  x = Sand100m30,
  y = X0_30_cm,
  xlab = "Log(Sand %) 100m grids",
  ylab = "Log(Sand %) soil survey",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = WBGNPSand,
  x = Sand100m30,
  y = Sand250mv2,
  xlab = "Log(Sand %) 100m grids",
  ylab = "Log(Sand %) 250m grids",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
WBGNPSandcor <- WBGNPSand %>% select(X0_30cm,Sand50m30,Sand100m30, Sand250m30, Sand250mv2)
# correlation for all variables
round(cor(WBGNPSandcor), digits = 2)
WBGNPSandccc <- epiR::epi.ccc(WBGNPSand$Sand50m30, WBGNPSand$Sand100m30, ci="z-transform", conf.level=0.95)
WBGNPSandccc.result <- WBGNPSandccc$rho.c
WBGNPSandccc.result
WBGNPSandcorall <- WBGNPSand %>% select(X0_30cm,Sand50m30,Sand100m30, Sand250m30, Sand250mv2)
correlation::correlation(WBGNPSandcorall,
                         include_factors = TRUE, method = "auto"
)
str(WBGNPClay)
WBGNPClay <- WBGNPClay %>% select(X0_30cm,Clay50m30,Clay100m30, Clay250m30,Clay250mv2)
WBGNPClay[WBGNPClay<= 0]<-NA
nrow(WBGNPClay)
WBGNPClay<-WBGNPClay[complete.cases(WBGNPClay), ]
nrow(WBGNPClay)
WBGNPClaySWTest<-WBGNPClay %>% shapiro_test(X0_30cm,Clay50m30,Clay100m30, Clay250m30)
WBGNPClaySWTest
WBGNPClay$X0_30cm<-log(WBGNPClay$X0_30cm)
ggdensity(WBGNPClay$X0_30cm, fill = "lightgray")
WBGNPClay$Clay50m30<-log(WBGNPClay$Clay50m30)
shapiro_test(WBGNPClay$X0_30cm)
shapiro_test(WBGNPClay$Clay50m30)
ggdensity(WBGNPClay$Clay50m30, fill = "lightgray")
WBGNPClay$Clay50m30<-log(WBGNPClay$Clay50m30)
ggdensity(WBGNPClay$Clay50m30, fill = "lightgray")
shapiro_test(WBGNPClay$Clay50m30)
ggdensity(WBGNPClay$Clay100m30, fill = "lightgray")
WBGNPClay$Clay100m30<-log(WBGNPClay$Clay100m30)
ggdensity(WBGNPClay$Clay100m30, fill = "lightgray")
ggscatterstats(
  data = WBGNPClay,
  x = Clay100m30,
  y = X0_30cm,
  xlab = "Log(Clay %) 100m grids",
  ylab = "Log(Clay %) soil survey",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = WBGNPClay,
  x = Clay50m30,
  y = Clay100m30,
  xlab = "Log(Clay %) 50m grids",
  ylab = "Log(Clay %) 100m grids",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
ggscatterstats(
  data = WBGNPClay,
  x = Clay100m30,
  y = Clay250mv2,
  xlab = "Log(Clay %) 100m grids",
  ylab = "Log(Clay %) 100m grids",
  bf.message = FALSE,
  marginal = FALSE # remove histograms
)
WBGNPClaycor <- WBGNPClay %>% select(X0_30cm,Clay50m30,Clay100m30, Clay250m30, Clay250mv2)
# correlation for all variables
round(cor(WBGNPClaycor), digits = 2)
WBGNPClayccc <- epiR::epi.ccc(WBGNPClay$Clay50m30, WBGNPClay$Clay100m30, ci="z-transform", conf.level=0.95)
WBGNPClayccc.result <- WBGNPClayccc$rho.c
WBGNPClayccc.result
WBGNPClaycorall <- WBGNPClay %>% select(X0_30_cm,Clay50m30,Clay100m30, Clay250m30, Clay250mv2)
correlation::correlation(WBGNPClaycorall,
                         include_factors = TRUE, method = "auto"
)
