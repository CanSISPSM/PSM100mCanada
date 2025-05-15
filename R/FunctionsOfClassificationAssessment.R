
# Some R functions for agreement and disagreement evaluations
# https://github.com/EcoDyn/rsacc/tree/master/R
#' Calculates quantity, allocation, exchange and shift disagreement indexes.
#'
#' function to calculate the disagreement measures from:
#' Pontius Jr, Robert Gilmore, and Millones, Marco. "Death to Kappa: birth of quantity disagreement and allocation disagreement for accuracy assessment." International Journal of Remote Sensing 32.15 (2011): 4407-4429 AND Pontius Jr, Robert Gilmore, and Santacruz, Alí. "Quantity, exchange, and shift components of difference in a square contingency table." International Journal of Remote Sensing, 35.21 (2014): 7543–7554.
#'
#' @param confmat a confusion matrix created by conf_mat(), or matrix with named columns and rows.
#' @return Overall agreement, disagreement, Quantity Disagreement, Allocation Disagreement, Exchange Disagreement, Shift Disagreement, and metrics at class level.
#' @importFrom diffeR overallExchangeD
#' @importFrom diffeR overallShiftD
#' @importFrom diffeR diffTablej
#' @export
pontius <- function(confmat){
  confsumms <- conf_summs(confmat)
  
  # Class quantity disagreement vector for J classes as in Pontius et al 2011
  Qvec <- abs(confsumms$Pig - confsumms$Pgj)
  
  # Overall Quantity Disagreement as in Pontius et al 2011
  Q <- sum(Qvec/2)/100
  
  # Allocation disagreement as in Pontius et al 2011
  Amat <- cbind(confsumms$Pig-confsumms$Pgg,confsumms$Pgj-confsumms$Pgg)
  Avec <- unname(2 * apply(Amat,1,min))
  A <- sum(Avec)/2/100
  
  # Overall Agreement (Percent Correctly Classified - PCC) as in Pontius et al 2011
  PCC <- sum(confsumms$Pgg)/100
  
  # Overall Disagreement as in Pontius et al 2011
  D = 1 - PCC
  
  # Overall Exchange Disagreement as in Pontius and Santacruz 2014
  E <- overallExchangeD(confsumms$confpmat)/100
  
  # Overall Shift Disagreement as in Pontius and Santacruz 2014
  S <- overallShiftD(confsumms$confpmat)/100
  
  # Decomponsing metrics per class as in Pontius and Santacruz 2014
  class_metrics <- diffTablej(unclass(confsumms$confpmat))
  
  # Droping E and S from diffTablej, alreade computed to be shown in our 'output style'
  class_metrics <- class_metrics[-nrow(class_metrics), ]
  
  # Converting number display to absolute instead of relative (default of diffTablej)
  # and round to 4 decimal digits (droping meaningful factor column for this calculation)
  temp <- round(class_metrics[,2:ncol(class_metrics)]/100,4)
  
  # Joining classes again
  class_metrics <- cbind(class_metrics$Category, temp)
  names(class_metrics) <- c("Classes", names(temp))
  rownames(class_metrics) <- NULL
  
  
  # Present results as a nice table:
  pontius2011_df <- data.frame("Value" = round(c(PCC,D,Q,A),4),
                               row.names = c("Overall Agreement",
                                             "Overall Disagreement",
                                             "Overall Quantity Disagreement",
                                             "Overall Allocation Disagreement"))
  
  pontius2014_df <- data.frame("Value" = round(c(E,S),4),
                               row.names = c("Overall Exchange Disagreement",
                                             "Overall Shift Disagreement"))
  
  
  
  return(c(list("Pontius et al. 2011 Disagreement Metrics:" = pontius2011_df),
           list("Decomposing your Allocation Disagreement as in Pontius and Santacruz 2014:" = pontius2014_df),
           list("Disagreement Metrics at class level:" = class_metrics)))
  
  
}
#' Check the validity of remote sensing classification data
#'
#' This function takes classification and validation inputs, given as *raster* or *Spatial* objects, and tests if projection are the same.
#'
#' @param map classification results as a *raster* or *SpatialPolygons* object.
#' @param val reference data for validating the classification. Accepts *raster*, *SpatialPolygons* or *SpatialPoints* objects.
#' @param reproj logical flag indicating if classification data should be reprojected to match reference data.
#' @return The sum of \code{x} and \code{y}
#'
#' @importFrom raster projection
#'
check_inputs <- function(map,val,reproj=FALSE){
  # Returns object classes so user knows what she is working with
  message(paste("Map data is a", class(map)))
  message(paste("Reference data is a", class(val)))
  
  # Inform user of map and val projections and test if they
  # are the same. Throw error if reproj = FALSE
  
  if (raster::projection(map) != raster::projection(val)){
    message(paste("Map projection: ",raster::projection(map)))
    message(paste("Reference projection: ",raster::projection(val)))
    if (reproj == TRUE){
      return(FALSE)
    } else {
      return(FALSE)
      stop("Error! Map projections are not equal. Use reproj=TRUE.")
    }
  } else {
    message(paste("Projections are the same: ",raster::projection(map)))
    return(TRUE)
  }
}
#' Build a confusion matrix from classified *raster* or *Spatial* objects.
#'
#' This function builds confusion matrices by overlaying *raster* or *Spatial* objects, automatically choosing the best method depending on input class.
#'
#' @param map classification results as a *raster* or *SpatialPolygons* object.
#' @param val reference data for validating the classification. Accepts *raster*, *SpatialPolygons* or *SpatialPoints* objects.
#' @param map_field Which column of the classified *SpatialPolygonsDataFrame* has class labels? Only used if map input is a *SpatialPolygon* object
#' @param val_field Which column of the validation *SpatialPolygonsDataFrame* or *SpatialPointsDataFrame* has class labels? Only used if validation input is a *Spatial* object.
#' @param reproj logical flag indicating if classification data should be reprojected to match reference data.
#' @param na_val are there data values that should be considered as NODATA? Specified value will be replaced by NA.
#' @param use_extract use *extract* method instead of *resample* method, See notes.
#'
#' @details When the classification results are a *raster* object and the reference data is a *SpatialPolygonsDataFrame*, the default method is to rasterize the reference data to match the classification grid. This can be slow and memory-heavy for large rasters. Setting 'use_extract = TRUE' will extract pixel values from the classification raster, based on vector layer instead. Thus method may be slower than the 'rasterize" method if polygons cover a large portion of the classification raster.
#'
#' @return The sum of \code{x} and \code{y}
#'
#' @importFrom sp spTransform
#' @importFrom sp over
#' @importFrom raster projection
#' @importFrom raster resample
#' @importFrom raster projectRaster
#' @importFrom raster extract
#' @importFrom stats na.omit



#' @export
conf_mat <- function(map, val, map_field=NA, val_field=NA, na_val=NA, reproj=FALSE, use_extract = FALSE){
  
  # Check if all required parameters are given
  if (inherits(map,"Spatial") & is.na(map_field)){
    stop("Classification dataset is a vector. Please specify the field name holding class names using 'val_field'.")
  }
  if (inherits(val,"Spatial") & is.na(val_field)){
    stop("Reference dataset is a vector. Please specify the field name holding class names using 'val_field'.")
  }
  
  # Reprojects val data to match map
  check <- check_inputs(map,val,reproj = reproj)
  if (check == FALSE){
    if (reproj == FALSE){
      stop("Projections are not the same!")} else {
        message("Reprojecting validation to match map.")
        
        # If val is a Spatial object
        if (inherits(val,'Spatial')){
          val <- sp::spTransform(val,CRSobj = raster::projection(map))
        }
        
        # If val is a Raster object
        if (inherits(val,'Raster')){
          val <- raster::projectRaster(val,crs = raster::projection(map), method = "ngb")
        }
      }
  }
  # Build confusion matrix for Raster vs Raster
  if (class(map) == "RasterLayer" & class(val) == "RasterLayer"){
    map_val <- raster::resample(map,val)
    
    # Replace specified na_val with NA
    if (!is.na(na_val)){
      val[val == na_val] <- NA
      map_val[map_val == na_val] <- NA
    }
    
    # Compute confusion matrix
    cmat <- na.omit(raster::crosstab(map_val,val))
    cdf <- matrix(cmat$Freq,
                  nrow=length(unique(cmat$Var1)),
                  ncol=length(unique(cmat$Var2)),
                  byrow=F)
    valnames <- na.omit(as.character(unique(cmat$Var2)))
    mapnames <- na.omit(as.character(unique(cmat$Var1)))
    colnames(cdf) <- valnames
    rownames(cdf) <- mapnames
    return(cdf)
  }
  
  # Build confusion matrix for Raster vs SpatialPolygons
  if (class(map) == "RasterLayer" & inherits(val,"SpatialPolygons")) {
    if (use_extract == FALSE){
      
      ### rasterize method
      message("Using 'rasterize' method. Might be slow and/or result in a memory error if classification raster is too large. Consider using 'use_extract = TRUE'.")
      rasval <- raster::rasterize(val, map, field = as.numeric(val@data[,val_field]))
      
      # Replace specified na_val with NA
      if (!is.na(na_val)){
        rasval[rasval == na_val] <- NA
        map[map == na_val] <- NA
      }
      
      cmat <- raster::crosstab(map,rasval)
      if (inherits(cmat, 'data.frame')){
        names(cmat) <- c('map','val','freq')
        cdf <- matrix(cmat$freq,
                      nrow=length(unique(cmat$map)),
                      ncol=length(unique(cmat$val)),
                      byrow=F)
        valnames <- na.omit(as.character(unique(cmat$val)))
        mapnames <- na.omit(as.character(unique(cmat$map)))
        colnames(cdf) <- valnames
        rownames(cdf) <- mapnames
        return(cdf)
      } else {
        cmatx <- as.matrix(cmat)
        cmatx[,as.character(sort(as.integer(colnames(cmatx))))]
        return(cmatx)
      }
      
    }
    if (use_extract == TRUE) {
      
      ### extract method
      message("Using 'extract' method. Might be slow if validation polygons cover a large portion of the classified raster. Consider using 'use_extract = FALSE'.")
      classdf <- raster::extract(map,val,df=T)
      valdf <- data.frame(class_val = val@data[,c(val_field)],
                          ID = sapply(slot(val, "polygons"), function(x) slot(x, "ID")))
      cvdf <- merge(classdf,valdf,by='ID')
      cmat <-  table(cvdf[,-1])
      return(cmat)
    }
  }
  
  # Build confusion matrix for Raster vs SpatialPoints
  if (class(map) == "RasterLayer" & inherits(val,"SpatialPoints")) {
    
    # Replace specified na_val with NA
    if (!is.na(na_val)){
      val@data[val@data[,val_field] == na_val,val_field] <- NA
      map[map == na_val] <- NA
    }
    
    mapvals <- raster::extract(map,val)
    cmat <-  table(val@data[,val_field],mapvals)
  }
  
  # Build confusion matrix for SpatialPolygons vs SpatialPolygons
  if (inherits(map, "SpatialPolygons") & inherits(val,"SpatialPolygons")) {
    stop("Not supported yet! :-(")
  }
  
  # Build confusion matrix for SpatialPolygons vs SpatialPoints
  if (inherits(map, "SpatialPolygons") & inherits(val,"SpatialPoints")) {
    mapvals <- sp::over(val, map[,map_field])
    cmat <- table(val@data[,val_field],mapvals[,map_field])
  }
  return(cmat)
}
#' Calculates summary statistics and totals from a confusion matrix
#'
#' This function calculates diagonal sums, row sums, and similar summaries from the confusion matrix, which are needed to calculate accuracy metrics.
#'
#' @param confmat a confusion matrix created by conf_mat(), or matrix with named columns and rows.
#' @return Several useful statistics.

conf_summs <- function(confmat){
  rn <- nrow(confmat) # number of rows
  cn <- ncol(confmat) # number of columns
  
  # Compute J - number of classes
  J <- if(rn==cn) {nrow(confmat)} else
  {stop("Unequal number of map and reference classes!")}
  
  # Class Names
  classnames <- colnames(confmat)
  
  # Total number of observations
  N <- sum(confmat)
  
  # Convert matrix to proportions to avoid integer overflowing
  confpmat <- confmat/N*100
  
  # Pig (colsums) x+j
  Pig <- colSums(confpmat)
  
  # Pgj (rowsums) xi+
  Pgj <- rowSums(confpmat)
  
  # Pgg - diagonal
  Pgg <- diag(confpmat)
  
  # Get all results into a list
  
  summs_list <- list(
    "confmat" =  confmat,
    "confpmat" = confpmat,
    "class_names" = classnames,
    "J" =  J,
    "N" =  N,
    "Pig" = Pig,
    "Pgj" = Pgj,
    "Pgg" = Pgg)
  
  return(summs_list)
}
#' Per class accuracies and *kappa* index of agreement
#'
#' Calculate per-class Omission and Comission errors and the Kappa index of agreement.
#'
#' @param confmat a confusion matrix created by conf_mat(), or matrix with named columns and rows.
#' @return Overall Accuracy, Kappa Index of Agreement, Per-class Omission and Comission Errors
#' @export
kia <- function(confmat){
  # Calculate necessary confusion matrix summaries
  confsumms <- conf_summs(confmat)
  
  # Expected agreement
  Eg <- sum(confsumms$Pig*confsumms$Pgj)
  
  # N * Diagonal
  E <- 100 * sum(confsumms$Pgg)
  
  # Denominator N² * Eg
  D <- 100^2 - Eg
  
  # Kappa Index of Agreement
  khat <- (E-Eg)/D
  
  # Overall Accuracy (Percent Correctly Classified - PCC)
  PCC <- sum(confsumms$Pgg)/100
  
  # Per class accuracies: producer's accuracy / omission error
  pro_omm <- 1 - (confsumms$Pgg/confsumms$Pig)
  
  # Per class accuracies: user's accuracy / comission error
  use_com <- 1 - (confsumms$Pgg/confsumms$Pgj)
  
  # Format output as a data frame
  khat_df <- data.frame(Accuracy_results = round(c(PCC,1-PCC,khat),4),
                        row.names = c("Overall Accuracy",
                                      "Overall Error",
                                      "Kappa Index of Agreement"))
  
  perclass_df <- data.frame("Omission Error"= round(pro_omm,4),
                            "Comission Error" = round(use_com,4),
                            row.names=confsumms$class_names)
  
  
  return(list("Overall Accuracy"= khat_df,
              "Class Accuracy" = perclass_df))
}

