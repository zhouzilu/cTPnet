#' Imputation of 12 surface proteins from scRNA-seq data.
#'
#' @param data data file: can be Seurat object, a matrix or a dataframe.
#' @param data_type specify type of data.Need to be one of the four
#' options: Seurat2, Seurat3, matrix, dataframe. You can check your
#' Seurat version by sessionInfo()'
#' @param model_file_path str: indicate the path to the trained cTPnet pytorch
#' model for the prediction
#' @return If input is Seurat object, the returned value will also be a Seurat
#' object; if input is matrix or dataframe, the returned value will be a
#' dataframe with rows as proteins and columns as cells.
cTPnet=function(data,data_type='Seurat2',model_file_path,dprotein=24){
  cat('Start data preprocessing...\n')
	if(data_type=='Seurat2'){
		X=preprocess_seurat2(data, dprotein)
	}else if(data_type=='Seurat3'){
		X=preprocess_seurat3(data, dprotein)
	}else if (data_type=='matrix'|data_type=='dataframe'){
		X=preprocess_matrix(data, dprotein)
	}else{
		stop('Error: unrecognizable data_type argument. Need to be one of the four\n
		      options: Seurat2, Seurat3, matrix, dataframe. You can check your \n
		      Seurat version by sessionInfo()\n')
	}
  cat('Start imputation. Running python ...\n')
	ctpnet <- reticulate::import("ctpnet", convert = F)
	y_pred=reticulate::py_to_r(ctpnet$predict$predict(X,model_file_path,dprotein))
	cat('Postprocess...\n')
	if(data_type=='Seurat2'){
	  data=postprocess_seurat2(data,y_pred)
	}else if(data_type=='Seurat3'){
	  data=postprocess_seurat3(data,y_pred)
	}else{
	  data=postprocess_matrix(y_pred)
	}
	cat('Done!\n')
	return(data)
}



