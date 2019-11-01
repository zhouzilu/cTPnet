postprocess_seurat2=function(seurat_data,y_pred){
  rownames(y_pred)=paste0('ctpnet_',rownames(y_pred))
  y_pred=y_pred-apply(y_pred,1,min)
  seurat_data <- SetAssayData(seurat_data, assay.type = "cTPnet", slot = "data",
                              new.data = y_pred)
  return(seurat_data)
}

postprocess_seurat3=function(seurat_data,y_pred){
  y_pred=y_pred-apply(y_pred,1,min)
  seurat_data[["cTPnet"]] <- CreateAssayObject(data = as.matrix(y_pred))
  return(seurat_data)
}

postprocess_matrix=function(y_pred){
  rownames(y_pred)=paste0('ctpnet_',rownames(y_pred))
  y_pred=y_pred-apply(y_pred,1,min)
  return(y_pred)
}
