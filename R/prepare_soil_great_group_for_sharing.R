
### This script is used to extract dominate soil great group types from SLC3.2 polygons and from CanSIS 100m great group grids. 
## Input:
#   SLC polygon shapefile: the projection of this shapefile is LCC-WGS84
#   SLC component table
#   SLC soil name table
# the lookup table that relates soil great group types in SLC3.2 to pixel values in CanSIS 100m great group grids 
# 
## output
# dominate_cmp_per_poly.qs # POLY_ID, GG_ID, PERCENT_merged # dominate soil great group type for each SLC polygon
# psm_dominate_values_in_slc_poly_all.qs # POLY_ID, GG # dominate soil great group type for PSM pixels within each SLC polygon
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

library(terra)
library(foreign) # read.dbf
library(tidyterra) 
library(qs)
library(dplyr)


### input ----------------
slc32 <- vect("data/soil_type_assessment_input/slc32_projected.shp") # LCC WGS84, projected from F:/ca_all_slc_v3r2/ca_all_slc_v3r2.shp
cmp <- read.dbf("data/soil_type_assessment_input/ca_all_slc_v3r2_cmp.dbf", as.is = T)
soilname_table <- read.dbf("data/soil_type_assessment_input/soil_name_canada_v2r20231107.dbf", as.is = T)
lookup_table <- read.csv("data/soil_type_assessment_input/GG_ID_GTGRP_lookup_table.csv") # GG_code == G_GROUP4

psm <- rast("data/soil_type_assessment_input/GreatGroupSoil100mv1.tif")

str(cmp) # POLY_ID, CMP, PERCENT, SOIL_ID, CMP_ID
str(soilname_table) # SOIL_ID, G_GROUP3
str(lookup_table) # GTGRP, GG_CODE, GG_ID

### preprocess of SLC -----------------
# join lookup_table and soilname_table to cmp in order to get GG_ID for each cmp for each SLC polygon

unique(soilname_table$G_GROUP3)
soilname_table <- soilname_table %>% filter(G_GROUP3 != "-") # 14323
# str(soilname_table)

# join GG_ID (in lookup table) to soilname_table 
str(lookup_table) # GG_code: character, e.g. MB, EB, SB
soilname_table <- soilname_table %>% left_join(lookup_table, by = join_by(G_GROUP3 == GG_code)) #  "SM"  "RG"  "SBR" not in lookup table
soilname_table <- soilname_table %>% dplyr::select(SOIL_ID, SOIL_CODE, ORDER3, G_GROUP3, GTGRP, GG_ID) # use GG_ID
unique(soilname_table$GG_ID)
soilname_table <- soilname_table %>% filter(!is.na(GG_ID))
str(soilname_table) # 14316

# join soilname table to cmp table
str(cmp)
str(soilname_table)

cmp_join <- cmp %>% left_join(soilname_table, by = "SOIL_ID") # good to join based on SOIL_ID
unique(cmp_join$GG_ID)
cmp_join <- cmp_join %>% filter(!is.na(GG_ID))
str(cmp_join)
cmp_join <- cmp_join %>% select(POLY_ID, CMP, PERCENT, SOIL_CODE = SOIL_CODE.x, SOIL_ID, CMP_ID, ORDER3, G_GROUP3, GTGRP, GG_ID)
str(cmp_join) # component information for each polygon, including percent, GG_ID, need to join it to SLC
write.csv(cmp_join, file = "soil_type_assessment_result/cmp_join_GG_ID.csv", row.names = F) # one polygon usually have more than 1 row

# now cmp_join has POLY_ID, CMP, GG_ID, PERCENT

# calculate dominate soil great group for each SLC polygon. 
unique_polyID <- unique(cmp_join$POLY_ID) # 4530 polygon

dominate_cmp_per_poly <- data.frame() # dominate GG_ID per polygon

for (i in 1:length(unique_polyID)){
  t1 <- cmp_join %>% filter(POLY_ID == unique_polyID[i])
  t2 <- t1 %>% group_by(GG_ID) %>% mutate(PERCENT_merged = sum(PERCENT)) %>% data.frame
  t2 <- t2 %>% mutate(flag = paste0(POLY_ID,"-", GG_ID)) %>% group_by(flag) %>% slice(1) %>% data.frame() %>% select(POLY_ID, GG_ID, PERCENT_merged)
  t2 <- t2 %>% arrange(desc(PERCENT_merged)) 
  
  dominate_cmp_per_poly <- rbind(dominate_cmp_per_poly, t2[1,])
}

length(unique(dominate_cmp_per_poly$POLY_ID)) # 4530
qsave(dominate_cmp_per_poly, file = "data/soil_type_assessment_result/dominate_cmp_per_poly.qs")
write.csv(dominate_cmp_per_poly, file = "data/soil_type_assessment_result/dominate_cmp_per_poly.csv", row.names = F)

### preprocess of CanSIS PSM100m -----------------

cmp_join_GG_ID <- read.csv("data/soil_type_assessment_result/cmp_join_GG_ID.csv")
slc32_all <-  slc32 %>% filter(POLY_ID %in% cmp_join_GG_ID$POLY_ID) 
dim(slc32_all) # 4530

psm_ext_slc_dominate <- data.frame()

for (i in 1: nrow(slc32_all)){
  poly <- slc32_all[i, ]
  tmp <- extract(psm, poly) # extract() is faster than rasterize, mask, crop
  tmp <- na.omit(tmp)
 
  if (nrow(tmp) != 0){
    tmp0 <- data.frame(table(tmp$pred)) %>% mutate(pred = as.character(Var1))
    
    tmp1 <- tmp0 %>% arrange(desc(Freq)) %>% slice(1) %>% dplyr::select(pred)
    psm_ext_slc_dominate <- rbind(psm_ext_slc_dominate, data.frame(POLY_ID = slc32_all[i, ]$POLY_ID, GG = tmp1$pred))
  }
}
 

str(psm_ext_slc_dominate) # 4456 obj, POLY_ID, GG
psm_ext_slc_dominate[duplicated(psm_ext_slc_dominate$POLY_ID), ]
length(unique(psm_ext_slc_dominate$POLY_ID)) # no duplicated POLY_ID
qsave(psm_ext_slc_dominate, file = "data/soil_type_assessment_result/psm_dominate_values_in_slc_poly_all.qs")
write.csv(psm_ext_slc_dominate, file = "data/soil_type_assessment_result/psm_dominate_values_in_slc_poly_all.csv", row.names = F)
