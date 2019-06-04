postprocess_seurat2=function(seurat_data,y_pred){
  rownames(y_pred)=paste0('cTPnet-',rownames(y_pred))
  seurat_data <- SetAssayData(seurat_data, assay.type = "cTPnet", slot = "data",
                              new.data = y_pred)
  return(seurat_data)
}

postprocess_seurat3=function(seurat_data,y_pred){
  seurat_data[["cTPnet"]] <- CreateAssayObject(data = seurat_data)
  return(seurat_data)
}

postprocess_matrix=function(y_pred){
  rownames(y_pred)=paste0('cTPnet-',rownames(y_pred))
  return(y_pred)
}
