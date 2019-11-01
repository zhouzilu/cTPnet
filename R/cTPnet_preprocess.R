preprocess_seurat2=function(seurat_data,dprotein){
    if (dprotein==12){
    	data(gene_names_protein12)
    }else if (dprotein==24){
    	data(gene_names_protein24)
    }else{
    	stop('Error: unrecognizable number of protein. Need to be 12 or 24\n')
    }
	data=as.data.frame(as.matrix(seurat_data@raw.data))
	data=data[shared_gene,]
	data[is.na(data)]=0.0
	rownames(data)=shared_gene
	return(data)
}

preprocess_seurat3=function(seurat_data,dprotein){
    if (dprotein==12){
    	data(gene_names_seurat3_protein12)
    }else if (dprotein==24){
    	data(gene_names_protein24)
    }else{
    	stop('Error: unrecognizable number of protein. Need to be 12 or 24\n')
    }
	data=as.data.frame(as.matrix(seurat_data$RNA@counts))
	data=data[shared_gene,]
	data[is.na(data)]=0.0
	rownames(data)=shared_gene
	return(data)
}

preprocess_matrix=function(matrix_data,dprotein){
    if (dprotein==12){
    	data(gene_names_protein12)
    }else if (dprotein==24){
    	data(gene_names_protein24)
    }else{
    	stop('Error: unrecognizable number of protein. Need to be 12 or 24\n')
    }
	data=as.data.frame(as.matrix(matrix_data))
	data=data[shared_gene,]
	data[is.na(data)]=0.0
	rownames(data)=shared_gene
	return(data)
}
