# This R script is to compare dominate components (i.e. soil great groups) for 
# SLC polygons to dominate soil great groups from CanSIS 100m grids
# Input:  Dominate soil great group per SLC polygon from SLC3.2 polygons.
#         Dominate soil great group per SLC polygon from CanSIS 100m grids.
# Output: overall agreement, overall disagreement, overall quantiy disagreement,
#         overall allocation disagreement, and etc.
# Author: Juanxia He
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
library(qs)
library(yardstick)
library(diffeR)
library(dplyr)
# Make sure that there is R sub-folder within the R project
source("R/FunctionsOfClassificationAssessment.R")
### input -------------
slc_cmp_dominate <- qread("data/soil_type_assessment_result/dominate_cmp_per_poly.qs") # 4530
psm_ext_slc_dominate <- qread("data/soil_type_assessment_result/psm_dominate_values_in_slc_poly_all.qs") # 4456
### preprocess -------------
str(slc_cmp_dominate) # GG_ID, int
str(psm_ext_slc_dominate) # GG, character
psm_ext_slc_dominate <- psm_ext_slc_dominate %>% mutate(GG = as.integer(GG))

slc_GG_PSM_GG_ID_dominate <- slc_cmp_dominate %>% left_join(psm_ext_slc_dominate, by = "POLY_ID") # 4530
slc_GG_PSM_GG_ID_dominate <- slc_GG_PSM_GG_ID_dominate %>% filter(!is.na(GG)) # 4456

names(slc_GG_PSM_GG_ID_dominate) # "POLY_ID"        "GG_ID"          "PERCENT_merged" "GG" 
names(slc_GG_PSM_GG_ID_dominate) <- c("POLY_ID", "SLC_GG", "PERCENT_SLC", "PSM_GG")
str(slc_GG_PSM_GG_ID_dominate)

length(unique(slc_GG_PSM_GG_ID_dominate$SLC_GG))
length(unique(slc_cmp_dominate$GG_ID))
length(unique(slc_GG_PSM_GG_ID_dominate$PSM_GG))
length(unique(psm_ext_slc_dominate$GG))

levels_value <- as.character(seq(1, 31)) # values of GG_ID, characters
levels_value <- factor(levels_value)
levels_value

slc_GG_PSM_GG_ID_dominate <- slc_GG_PSM_GG_ID_dominate %>%
  mutate(SLC_GG_3 = factor(as.character(SLC_GG), levels = levels_value),
         PSM_GG_3 = factor(as.character(PSM_GG), levels = levels_value))
str(slc_GG_PSM_GG_ID_dominate)
summary(slc_GG_PSM_GG_ID_dominate$SLC_GG - as.integer(as.character(slc_GG_PSM_GG_ID_dominate$SLC_GG_3)))
summary(slc_GG_PSM_GG_ID_dominate$PSM_GG - as.integer(as.character(slc_GG_PSM_GG_ID_dominate$PSM_GG_3)))
table(slc_GG_PSM_GG_ID_dominate$SLC_GG)
table(slc_GG_PSM_GG_ID_dominate$SLC_GG_3)
table(slc_GG_PSM_GG_ID_dominate$PSM_GG)
table(slc_GG_PSM_GG_ID_dominate$PSM_GG_3)

str(slc_GG_PSM_GG_ID_dominate)
summary(levels(slc_GG_PSM_GG_ID_dominate$SLC_GG_3) == levels(slc_GG_PSM_GG_ID_dominate$PSM_GG_3))

### Comparison -------------
conf_matrix <- yardstick::conf_mat(slc_GG_PSM_GG_ID_dominate, "SLC_GG_3", "PSM_GG_3") # data, truth, estimate

yardstick::accuracy(slc_GG_PSM_GG_ID_dominate, "SLC_GG_3", "PSM_GG_3") # same as pontius's Overall Agreement

str(conf_matrix)

accu_all <- pontius(conf_matrix$table)
accu_all

qsave(accu_all, file = "data/soil_type_assessment_result/pontius_result_SLC_all_vs_PSM.qs")
qsave(slc_GG_PSM_GG_ID_dominate, file = "data/soil_type_assessment_result/dominate_GG_from_PSM_and_SLC_CMP.qs")
qsave(conf_matrix, file = "data/soil_type_assessment_result/confusion_matrix_dominate_GG_from_PSM_and_SLC_CMP.qs")

# join to slc for showing
library(tidyterra)

slc32 <- vect("data/soil_type_assessment_input/slc32_projected.shp") 
slc32_2 <- slc32 %>% left_join(slc_GG_PSM_GG_ID_dominate, by = "POLY_ID")
unique(slc32_2$PSM_GG_3)
slc32_2 <- slc32_2 %>% filter(!is.na(PSM_GG_3))
dim(slc32_2)
writeVector(slc32_2, file = "data/soil_type_assessment_result/slc_GG_PSM_GG_ID_dominate.shp")

