# cTPnet
single cell Transcriptome to Protein prediction with deep neural network (cTP-net) is a transfer learning framework to impute surface protein abundances from scRNA-seq data by learning from existing single-cell multi-omic resources. While single cell RNA sequencing (scRNA-seq) is invaluable for studying cell populations, cell-surface proteins are often integral markers of cellular function and serve as primary targets for therapeutic intervention.

## Manuscript

[Surface protein imputation from single cell transcriptomes by deep neural networks]()


## Questions & Problems

If you have any questions or problems when using cTPnet or ctpnetpy, please feel free to open a new issue [here](https://github.com/zhouzilu/cTPnet/issues). You can also email the maintainers of the corresponding packages -- the contact information is shown under Developers & Maintainers.

## Pipeline overview

This cTPnet package includes two analysis tools: (1) **SAVERX**, we stronly recommend denoise the scRNA-seq data before impute the surface protein abundance, and (2) **cTP-net**, which impute the surface protein abundance based on previously trained model. 

<p align="center">
  <img src="https://raw.githubusercontent.com/zhouzilu/cTPnet/master/figure/FIG_pkg.jpg" width='800'>
  </p>

  **Figure 1.**  *(a)* Overview of cTP-net analysis pipeline, which learns a mapping from the denoised scRNA-seq data to the relative abundance of surface proteins, capturing multi-gene features that reflect the cellular environment and related processes. *(b)* For three example proteins, cross-cell scatter and correlation of CITE-seq measured abundances vs. (1) raw RNA count, (2) SAVER-X denoised RNA level, and (3) cTP-net predicted protein abundance.

## Running cTP-net
For different Seurat version, we developed separate vignette see below:

### Seurat v2
  **cTP-net R notebook** with step-by-step demonstration and rich display is available [***here***](http://rawgit.com/zhouzilu/cTPnet/master/vignette/cTPnet_vignette_v2.html). Corresponding **Rmd script** is available [***here***](https://github.com/zhouzilu/master/blob/master/vignette/cTPnet_vignette_v2.Rmd).

### Seurat v3
  **cTP-net R notebook** with step-by-step demonstration and rich display is available [***here***](http://rawgit.com/zhouzilu/cTPnet/master/vignette/cTPnet_vignette_v3.html). Corresponding **Rmd script** is available [***here***](https://github.com/zhouzilu/cTPnet/blob/master/vignette/cTPnet_vignette_v3.Rmd).


## Citation

Please cite cTP-net.

* **cTP-net**: [link]()
<br>

## Developers & Maintainers

* [Zilu Zhou](https://statistics.wharton.upenn.edu/profile/zhouzilu/) (zhouzilu at pennmedicine dot upenn dot edu)
  <br>
  Genomics and Computational Biology Graduate Group, University of Pennsylvania

* [Nancy R. Zhang](https://statistics.wharton.upenn.edu/profile/nzh/) (nzh at wharton dot upenn dot edu)
  <br>
  Department of Statistics, University of Pennsylvania
