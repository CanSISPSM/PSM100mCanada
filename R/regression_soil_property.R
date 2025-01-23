
## This script generates RF predictions for soil properties
## 
## The steps include: creating RF models, generating predictions per tile, merging tile prediction results to form the final result 
## 
## The final results include predictions and quantile intervals for soil properties
##

library(ranger)
library(terra)
library(tictoc)
library(future)
library(furrr)
library(stringr)
library(qs)
library(dplyr)

getwd()

source("R/PSM_functions.R")

#========================================
#  Input
#========================================

# name of soil properties
soil_prop_name <- c("SOC", # 1
                    "BD",  # 2
                    "CEC", # 3
                    "pHCACL2", # 4
                    "Sand", # 5
                    "Clay", # 6
                    "Phosphorus_sat_index",
                    "CN_Ratio")
depth_interval <- c(0, 5, 15, 30, 60, 100) # soil depths                  

depth_field <- paste0("X", depth_interval[-length(depth_interval)], "_", depth_interval[-1], "_cm") # depth fields, e.g "X0_5_cm"    "X5_15_cm"   "X15_30_cm"  "X30_60_cm"  "X60_100_cm"

n_prop = 1 # index of soil property

output_folder <- soil_prop_name[n_prop] # set the output folder name for all the results for the soil property
if(!dir.exists(output_folder)){dir.create(output_folder)}
output_folder


# folder for prediction to new data
path_2_newdata = "newdata_tiles/newdata4prediction" # new data, used for generating per tile prediction
path_2_base_raster_folder = "newdata_tiles/base_tile_for_20_20" # used for converting dataframe to tif

out_path_root = paste0(output_folder,"/tmp") # save all the temporary (intermediate) data. The folder can be deleted once the final result is created.
if(!dir.exists(out_path_root)){dir.create(out_path_root)}
out_path_root

final_result_out_path = paste0(output_folder, "/final_result") # folder to save all the final results (i.e. tif files)
if(!dir.exists(final_result_out_path)){dir.create(final_result_out_path)}
final_result_out_path

lower = 0.05 # upper and lower quantile values 
upper = 0.95

tile_newdata_list <- list.files(path_2_newdata, pattern =  "*.qs", full.names = T) # list of file names for new data, results from regression_soil_property_step_1_create_tiles_for_prediction.R

#========================================
#  RF prediction
#========================================

# for (n_layer in 1:length(depth_field)){

# input: read in regression matrix, factor levels ------------

n_layer = 1 # set n_layer to 1,2,3,4,5 to generate predictions for all the layers

# load regression matrix for each layer
paste0(output_folder, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_regression_matrix.qs")
regre_matrix <- qread(paste0(output_folder, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_regression_matrix.qs")) # regression matrix
summary(regre_matrix[,1])


# extract levels for all factor covariates 
factor_cov_index <- which(sapply(regre_matrix, is.factor))
names(regre_matrix)[factor_cov_index] # e.g surficial geology, land cover 
levels_all <- sapply(regre_matrix[, factor_cov_index], levels)
levels_all

# step 1:  generate RF models for BC-RF stage 1 (i.e. predictions to mean), BC-RF stage 2(i.e. predictions to residuals), and quantile predictions for uncertainty------

formula_string <- as.formula(paste(depth_field[n_layer], ' ~ ', paste(names(regre_matrix %>% dplyr::select(-depth_field[n_layer])), collapse="+")))
formula_string # e.g. X5_15_cm ~ LandSatNDVI + TWI + RHSP + ...

#  BC-RF stage1 model
rf_model_mean <- ranger(formula_string, 
                        data = regre_matrix, 
                        # num.trees = rf_model_tuning$trees,
                        # mtry = rf_model_tuning$mtry,
                        # min.node.size = rf_model_tuning$min_n, 
                        importance = "impurity", # GINI for classification, variance of responses for regression
                        respect.unordered.factors = T,
                        write.forest = T
)

rf_model_mean
rf_model_mean$variable.importance

# BC-RF stage 1 OOB prediction residuals 
resi <- regre_matrix[ , 1] - rf_model_mean$predictions # resi = observation - prediction, thus the final should be prediction (stage 1 RF) + resi (prediction on resi, i.e. stage 2 RF)
# summary(resi)  
# summary(regre_matrix[,1]) 

# BC-RF stage2 model based on residuals
resi_regre_matrix <- cbind(resi, regre_matrix)
formula_string_resi <- as.formula(paste("resi", ' ~ ', paste(names(resi_regre_matrix)[-1], collapse="+")))
formula_string_resi # e.g. resi ~ X5_15_cm + LandSatNDVI + TWI + RHSP + ...

rf_model_mean_resi <- ranger(formula_string_resi, 
                             data = resi_regre_matrix, 
                             importance = "impurity",
                             respect.unordered.factors = T,
                             write.forest = T
)
rf_model_mean_resi

qsave(rf_model_mean, paste0(out_path_root, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_rf_model_stage_1_mean.qs"))
qsave(rf_model_mean_resi, paste0(out_path_root, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_rf_model_stage_2_mean_resi.qs"))

# RF model for quantile predictions for uncertainty
rf_model_quantile <- ranger(formula_string, 
                            data = regre_matrix, 
                            importance = "impurity",
                            respect.unordered.factors = T,
                            write.forest = T,
                            quantreg = T
)
rf_model_quantile
qsave(rf_model_quantile, paste0(out_path_root, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_rf_model_quantile_uncertainty.qs"))

# step 2:  generate BC-RF predictions and quantile predictions (dataframe in qs format) for each tile ------

# predict RF mean based on stage 1 RF model
out_tmp <- paste0(out_path_root, "/layer_", depth_field[n_layer]) # folder to save intermediate results for each layer, e.g. SOC/tmp/layer_X60_100_cm
if(!dir.exists(out_tmp)){dir.create(out_tmp)}
out_tmp

tic()
plan(multisession, workers = 4) 
tile_newdata_list %>% future_map(rfRegressionModelPredict2TilesWithBCmodel,
                                 rf_model_stage1 = rf_model_mean,
                                 rf_model_stage2 = rf_model_mean_resi,
                                 Y_name = names(regre_matrix)[1],
                                 out_path = out_tmp,
                                 all_factor_levels = levels_all, 
                                 .options = furrr_options(seed = 123))
plan(sequential)
toc()
cat("done") 

# predict uncertainty (quantiles) for each tile 
tic()
plan(multisession, workers = 2) 
tile_newdata_list %>% future_map(rfRregressionQuantileModelPredict2Tiles4Uncertainty,
                                 rf_quantile_model = rf_model_quantile,
                                 quantile_lower = lower,
                                 quantile_upper = upper,
                                 out_path = out_tmp,
                                 all_factor_levels = levels_all,
                                 .options = furrr_options(seed = 123))
plan(sequential) 
toc()
cat("done") 

# step 3:  convert tile predictions in dataframe format (saved as .qs files) to tif files, and merge tile results to one ------

RF_prediction_final_files <- list.files(out_tmp, pattern = "RF_prediction_final_", full.names = T) # list of .qs files containing BC-RF predictions
RF_prediction_stage1_mean_files <- list.files(out_tmp, pattern = "RF_stage1_mean_prediction_", full.names = T) # list of .qs files containing RF stage1 prediction (mean)

# BC-RF results in tif format but saved as .qs
tic() 
plan(multisession, workers = 20) 
RF_prediction_final_files %>% future_map(qsDF2tifFun,  
                                         path_newdata_with_NA = path_2_newdata,
                                         path_base_raster = path_2_base_raster_folder, 
                                         result_file_name_inti = "RF_final", 
                                         all_factor_levels = levels_all,
                                         prediction_data_type = "numeric",
                                         .options = furrr_options(seed = 123))
plan(sequential)
toc()
cat("done")  

# RF stage1 predictions in tif format but saved as .qs
tic()
plan(multisession, workers = 20) 
RF_prediction_stage1_mean_files %>% future_map(qsDF2tifFun,  
                                               path_newdata_with_NA=path_2_newdata,
                                               path_base_raster = path_2_base_raster_folder, 
                                               result_file_name_inti = "RF_stage1_mean",
                                               all_factor_levels = levels_all,
                                               prediction_data_type = "numeric",
                                               .options = furrr_options(seed = 123))
plan(sequential)
toc()
cat("done")

# merge BC-RF per-tile results to one
tile_pred_final_tif <- list.files(out_tmp,  pattern = "RF_final*", full.names = T)
tile_pred_final_tif_sprc <- lapply(tile_pred_final_tif, function(x) {unwrap(qread(x))})
tile_pred_final_tif_sprc <- sprc(tile_pred_final_tif_sprc)
final_pred <- merge(tile_pred_final_tif_sprc)
final_pred[final_pred < 0] = 0
if (n_prop %in% c(1, 5, 6)){ # for SOC, Sand,Clay, adjustment needed
  writeRaster(final_pred, paste0(final_result_out_path, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_prediction_RF_BC_0.tif"), datatype = "FLT4S", NAflag = -999)
} else { # for BD, CEC, pH, no adjustment needed
  writeRaster(final_pred, paste0(final_result_out_path, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_prediction_RF_BC.tif"), datatype = "FLT4S", NAflag = -999) 
}


if(n_prop == 1) { # for SOC
  final_pred[final_pred > 60] = 60
  writeRaster(final_pred, paste0(final_result_out_path, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_prediction_RF_BC.tif") , datatype = "FLT4S", NAflag = -999)#, datatype = "INT2S", NAflag = -999)
}
if(n_prop %in% c(5, 6)){ # for sand, clay
  final_pred[final_pred > 100] = 100
  writeRaster(final_pred, paste0(final_result_out_path, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_prediction_RF_BC.tif"), datatype = "FLT4S", NAflag = -999) #, datatype = "INT2S", NAflag = -999)
}

# RF stage1 prediction mean
tile_pred_mean_tif <- list.files(out_tmp,  pattern = "RF_stage1_mean_tile*", full.names = T)
tile_pred_mean_tif_sprc <- lapply(tile_pred_mean_tif, function(x) {unwrap(qread(x))})
tile_pred_mean_tif_sprc <- sprc(tile_pred_mean_tif_sprc)
mean_pred <- merge(tile_pred_mean_tif_sprc)
mean_pred[mean_pred < 0] = 0
if(n_prop == 1) { # for SOC
  mean_pred[mean_pred > 60] = 60
}
writeRaster(mean_pred, paste0(final_result_out_path, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_prediction_stage1_mean.tif"), datatype = "FLT4S", NAflag = -999) #, datatype = "INT2S", NAflag = -999)

# for prediction of uncertainty (quantiles)
quantile_upper_pred_files <- list.files(out_tmp, pattern = "quantile_prediction_upper*", full.names = T)
tic()
plan(multisession, workers = 20) # workers = 20 for 100m PSM
quantile_upper_pred_files %>% future_map(qsDF2tifFun,  
                                         path_newdata_with_NA = path_2_newdata,
                                         path_base_raster = path_2_base_raster_folder, 
                                         result_file_name_inti = "quantile_upper_prediction",
                                         all_factor_levels = levels_all,
                                         prediction_data_type = "numeric",
                                         .options = furrr_options(seed = 123))


plan(sequential)
toc()
cat("done")

quantile_lower_pred_files <- list.files(out_tmp, pattern = "quantile_prediction_lower*", full.names = T)
tic()
plan(multisession, workers = 20) 
quantile_lower_pred_files %>% future_map(qsDF2tifFun,  
                                         path_newdata_with_NA = path_2_newdata,
                                         path_base_raster = path_2_base_raster_folder, 
                                         result_file_name_inti = "quantile_lower_prediction",
                                         all_factor_levels = levels_all,
                                         prediction_data_type = "numeric",
                                         .options = furrr_options(seed = 123))

plan(sequential)
toc()
cat("done")

quantile_median_pred_files <- list.files(out_tmp, pattern = "quantile_prediction_median*", full.names = T)
tic()
plan(multisession, workers = 20) 
quantile_median_pred_files %>% future_map(qsDF2tifFun,  
                                          path_newdata_with_NA = path_2_newdata,
                                          path_base_raster = path_2_base_raster_folder, 
                                          result_file_name_inti = "quantile_median_prediction",
                                          all_factor_levels = levels_all,
                                          prediction_data_type = "numeric",
                                          .options = furrr_options(seed = 123))

plan(sequential)
toc()
cat("done")

# merge lower, upper, median quantile per-tile results to one separately
tile_pred_quantile_upper_tif <- list.files(out_tmp,  pattern = "quantile_upper_prediction*", full.names = T)
tile_pred_quantile_upper_final_tif_sprc <- lapply(tile_pred_quantile_upper_tif, function(x) {unwrap(qread(x))})
tile_pred_quantile_upper_final_tif_sprc <- sprc(tile_pred_quantile_upper_final_tif_sprc)
final_quantile_upper_pred <- merge(tile_pred_quantile_upper_final_tif_sprc)
if(n_prop == 1) { # for SOC
  final_quantile_upper_pred[final_quantile_upper_pred > 60] = 60
}
if(n_prop %in% c(5, 6)){ # for sand
  final_quantile_upper_pred[final_quantile_upper_pred > 100] = 100
}
writeRaster(final_quantile_upper_pred, paste0(final_result_out_path, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_quantile_upper_prediction.tif"), datatype = "FLT4S", NAflag = -999) #, datatype = "INT2S", NAflag = -999)

tile_pred_quantile_lower_tif <- list.files(out_tmp,  pattern = "quantile_lower_prediction*", full.names = T)
tile_pred_quantile_lower_final_tif_sprc <- lapply(tile_pred_quantile_lower_tif, function(x) {unwrap(qread(x))})
tile_pred_quantile_lower_final_tif_sprc <- sprc(tile_pred_quantile_lower_final_tif_sprc)
final_quantile_lower_pred <- merge(tile_pred_quantile_lower_final_tif_sprc)
final_quantile_lower_pred[final_quantile_lower_pred < 0] = 0
writeRaster(final_quantile_lower_pred, paste0(final_result_out_path, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_quantile_lower_prediction.tif") , datatype = "FLT4S", NAflag = -999) #, datatype = "INT2S", NAflag = -999)

tile_pred_quantile_median_tif <- list.files(out_tmp,  pattern = "quantile_median_prediction*", full.names = T)
tile_pred_quantile_median_final_tif_sprc <- lapply(tile_pred_quantile_median_tif, function(x) {unwrap(qread(x))})
tile_pred_quantile_median_final_tif_sprc <- sprc(tile_pred_quantile_median_final_tif_sprc)
final_quantile_median_pred <- merge(tile_pred_quantile_median_final_tif_sprc)
writeRaster(final_quantile_median_pred, paste0(final_result_out_path, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_quantile_median_prediction.tif") , datatype = "FLT4S", NAflag = -999)#, datatype = "INT2S", NAflag = -999)

# calculate quantile intervals: upper - lower
quantile_interval <- final_quantile_upper_pred - final_quantile_lower_pred 
writeRaster(quantile_interval, paste0(final_result_out_path, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_quantile_interval.tif"), datatype = 'FLT4S', NAflag = -999)

# calculate relative quantile intervals based on the quantile interval value of training data
train_date_quantile_interval <- quantile(regre_matrix[, 1], probs = upper) - quantile(regre_matrix[, 1], probs = lower)  # the quantile interval value of training data
relative_quantile_interval <- (quantile_interval)/train_date_quantile_interval
writeRaster(relative_quantile_interval, paste0(final_result_out_path, "/", soil_prop_name[n_prop], "_", depth_field[n_layer], "_relative_quantile_interval_based_on_training_data_quantile_interval.tif"), datatype='FLT4S', NAflag = -999)


# unlink(out_tmp, recursive = TRUE) 
# cat('done')
# } # for n_layer loop



