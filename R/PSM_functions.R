
#' @description function removeFileExtension() removes file extension(no extension, no path)
#' @param fname file name (file basename) without path
#' @param extensionName file extension to be removed, e.g ".csv", ".tif
#' @return a file name without file extension and path
#' 
#' @example removeFileExtension(name1.tif, ".tif")
#' 
removeFileExtension <- function(fname, extensionName){ # fname = name1.tif, extensionName = ".tif"
  return(substr(fname,1,nchar(fname)-nchar(extensionName)))
}



#' @description function rfOneModelPredict2Tiles() generates RF classification predictions for each tile based on a RF model----
#' @param x full path for newdata for each tile. 
#' @param rf_model: a ranger random forest classification model
#' @param all_factor_levels: a list storing factor levels of all factor covariates
#' @param out_path: folder to save results
#' 
#' @export classification predictions (in qs format)
#' 
rfOneModelPredict2Tiles <- function(x, # name of new data without NA (full path)
                                    rf_model,  # RF model
                                    all_factor_levels, # includes levels for  all factor covariates, like sg,landcover, landform type(IP), geomor
                                    out_path){
 # require(stringr)
  require(qs)
  require(ranger)

  newdata <- qread(x)  # data.frame, x is df saved in qs
  for (i in 1:length(all_factor_levels)){
    newdata[, c(names(all_factor_levels)[i])] = factor(newdata[, c(names(all_factor_levels)[i])], levels =  all_factor_levels[[i]])
  }

 # tile_name = str_split(x, "_")[[1]][length(str_split(x, "_")[[1]])]
 # tile_name  <- paste0("tile_", str_replace(basename(x), ".qs", ""))

 newdata <- newdata[complete.cases(newdata), ]

 tmp <- predictions(predict(rf_model, newdata, predict.all = T)) # one row for one pixel, one column for one tree
 tmp_vote_majority <- predictions(predict(rf_model, newdata))
 qsave(tmp, paste0(out_path, "/predict_all_value_", basename(x)))
 qsave(tmp_vote_majority, paste0(out_path, "/predict_voted_majority_value_", basename(x)))

}
#' 

#' @description function modeFrequenceyCalFunc() calculates frequencies of the predicted majority classes
#' @param x is the predicted all values from predict.all = T, we also need the majority qs files for all tiles in the same folder
#' @export frequency values of majority-voted classes(mode)
#' 

#' @description function qsDF2tifFun() converts predicted results in a dataframe (saved as .qs files) to tif.
#' 
#' @param x a file name with full path. It saves predicted results for a tile in a .qs file
#' @param path_newdata_with_NA a file(.qs) name with full path. It contains new data (in dataframe format) to be predicted.
#' @param path_base_raster a file (.tif) name with full path. It contains a terra rast object used as the base (template) to convert predicted values to tiff files. 
#' @param result_file_name_inti a string that will be added tot the output file name  
#' @param all_factor_levels a list contains all levels of all factoro covariates
#' 
#' @export  a tif file but it is saved as a .qs file (using wrap()). It contains the prediction result(.tif file) per tile.
#'  
qsDF2tifFun <- function(x, # predictions for a tile generated from predcit(rf model, newdata, predict.all = T)
                        path_newdata_with_NA,
                        path_base_raster,
                        result_file_name_inti,
                        all_factor_levels,
                        prediction_data_type){
  # require(terra)
  # require(stringr)

  out_prediction = qread(x)

  tile_name = str_split(x, "_")[[1]][length(str_split(x, "_")[[1]])]
  tile_name  <- paste0("tile_", str_replace(tile_name, ".qs", ""))

  newdata <- qread(paste0(path_newdata_with_NA, "/", tile_name, ".qs")) # with NA in the data frame

  if(!missing(all_factor_levels)){
    for (i in 1:length(all_factor_levels)){
      newdata[, c(names(all_factor_levels)[i])] = factor(newdata[, c(names(all_factor_levels)[i])], levels =  all_factor_levels[[i]])
    }
  }

  tmp_v <- newdata[, 1] * NA # note: number of NA values of tmp may be different from that of base.raster.
  
  if(prediction_data_type == "factor"){ tmp_v[complete.cases(newdata)] <- as.integer(as.character(out_prediction))}
  if(prediction_data_type == "numeric"){ tmp_v[complete.cases(newdata)] <- out_prediction}
    
  base_raster <- rast(paste0(path_base_raster, "/", tile_name, ".tif")) # no need to set to NA
  values(base_raster) = tmp_v # length of tmp_v same as number of pixel of base_raster
  names(base_raster) = "pred"

  out_qs = paste0(dirname(x), "/", result_file_name_inti, "_", tile_name,  ".qs")
  qsave(wrap(base_raster), file = out_qs)
}


#' @description function rfRegressionModelPredict2TilesWithBCmodel() generates predictions per tile based on BC-RF models.
#' @param x a file name with full path for tiled new data. It is generated by regression_soil_property_step_1_create_tiles_for_prediction.R
#' @param rf_model_stage1 stage 1 RF model generated from ranger()
#' @param rf_model_stage2 stage 2 RF model (based on stage 1 OOB prediction error) generated from ranger()
#' @param Y_name the response variable name, e.g. X0_5_cm, X5_15-cm
#' @param out_path a folder name with path to save the function output (i.e. prediction results)
#' @param all_factor_levels a list contains all levels of all factor covariates
#' @return
#' @export  two qs files for each tile are saved. One qs file contains stage 1 RF predictions (dataframe),
#'  and the other qs file contains BC-RF final predictions (dataframe)
#'  
 rfRegressionModelPredict2TilesWithBCmodel <- function(x, # name of new data without NA (full path)
                                                       rf_model_stage1,  # path to RF model
                                                       rf_model_stage2,
                                                       Y_name,
                                                       out_path,
                                                       all_factor_levels # includes levels for  all factor covariates, like sg,landcover, landform type(IP), geomor
                                                       ){
   require(stringr)
   require(qs)
   require(ranger)

   newdata <- qread(x)  # data.frame, x is tile.newdata.list[i]

   if(!missing(all_factor_levels)){
     for (i in 1:length(all_factor_levels)){
       newdata[, c(names(all_factor_levels)[i])] = factor(newdata[, c(names(all_factor_levels)[i])], levels =  all_factor_levels[[i]])
     }
   }

   tile_name = str_split(x, "_")[[1]][length(str_split(x, "_")[[1]])]
   tile_name  <- paste0("tile_", str_replace(tile_name, ".qs", ""))

   newdata <- newdata[complete.cases(newdata), ]

   rf_stage1_pred_new_data <- data.frame(predictions(predict(rf_model_stage1, newdata))) 
   names(rf_stage1_pred_new_data) = Y_name

   rf_stage2_pred_new_data  <- data.frame(predictions(predict(rf_model_stage2, data.frame(rf_stage1_pred_new_data, newdata)))) # the second RF (resi) to new data
   names(rf_stage2_pred_new_data) = "pred_resi"

   rf_pred_new_data_final <- rf_stage1_pred_new_data + rf_stage2_pred_new_data # final = pred + resi

   # pred_all = list(rf_stage1_pred = rf_stage1_pred_new_data, rf_pred_final = rf_pred_new_data_final)
   # qsave(pred_all, paste0(out_path, "/predict_stage1_stage2_", tile_name, ".qs"))
   qsave(rf_stage1_pred_new_data[,1], paste0(out_path, "/RF_stage1_mean_prediction_", tile_name, ".qs"))
   qsave(rf_pred_new_data_final[,1], paste0(out_path, "/RF_prediction_final_", tile_name, ".qs"))
 }

 #' @description function rfRregressionQuantileModelPredict2Tiles4Uncertainty() generates quantile predictions per tile based on a RF model.
 #' @param x a file name with full path for tiled new data. It is generated by regression_soil_property_step_1_create_tiles_for_prediction.R
 #' @param rf_quantile_model a RF quantile model generated from ranger()
 #' @param quantile_lower lower quantile value, e.g. 0.05
 #' @param quantile_upper upper quantile value, e.g. 0.95
 #' @param out_path a folder name with path to save the function output (i.e. prediction results)
 #' @param all_factor_levels a list contains all levels of all factoro covariates
 #' @return
 #' @export  three .qs files for each tile are saved for upper, lower, and median quantile predictions
 #' 
rfRregressionQuantileModelPredict2Tiles4Uncertainty <- function(x,
                                                                rf_quantile_model,
                                                                quantile_lower,
                                                                quantile_upper,
                                                                out_path,
                                                                all_factor_levels){
  require(stringr)
  require(qs)
  require(ranger)

  newdata <- qread(x)  # data.frame, x is tile.newdata.list[i], levels of factor covariates not set

  if(!missing(all_factor_levels)){
    for (i in 1:length(all_factor_levels)){
      newdata[, c(names(all_factor_levels)[i])] = factor(newdata[, c(names(all_factor_levels)[i])], levels =  all_factor_levels[[i]])
    }
  }


  tile_name = str_split(x, "_")[[1]][length(str_split(x, "_")[[1]])]
  tile_name  <- paste0("tile_", str_replace(tile_name, ".qs", ""))

  newdata <- newdata[complete.cases(newdata), ]
  # colnumes(rf_model_quantile_pred): the first column is for lower, the second for median, the third for upper
  rf_model_quantile_pred <- predictions(predict(rf_quantile_model,
                                    newdata,
                                    type = "quantiles",
                                    quantiles = c(quantile_lower, 0.5, quantile_upper)))


  qsave(rf_model_quantile_pred[,1], paste0(out_path, "/quantile_prediction_lower_", tile_name, ".qs"))
  qsave(rf_model_quantile_pred[,2], paste0(out_path, "/quantile_prediction_median_", tile_name, ".qs"))
  qsave(rf_model_quantile_pred[,3], paste0(out_path, "/quantile_prediction_upper_", tile_name, ".qs"))
}




