preprocess_seurat2=function(seurat_data){
  data(gene_names)
	data=as.data.frame(as.matrix(seurat_data@raw.data))
	data=data[shared_gene,]
	data[is.na(data)]=0.0
	rownames(data)=shared_gene
	return(data)
}

preprocess_seurat3=function(seurat_data){
	data(gene_names_seurat3)
	data=as.data.frame(as.matrix(seurat_data$RNA@counts))
	data=data[shared_gene,]
	data[is.na(data)]=0.0
	rownames(data)=shared_gene
	return(data)
}

preprocess_matrix=function(matrix_data){
  data(gene_names)
	data=as.data.frame(as.matrix(matrix_data))
	data=data[shared_gene,]
	data[is.na(data)]=0.0
	rownames(data)=shared_gene
	return(data)
}
