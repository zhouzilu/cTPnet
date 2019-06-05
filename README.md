# cTPnet
single cell Transcriptome to Protein prediction with deep neural network (cTP-net) is a transfer learning framework to impute surface protein abundances from scRNA-seq data by learning from existing single-cell multi-omic resources. While single cell RNA sequencing (scRNA-seq) is invaluable for studying cell populations, cell-surface proteins are often integral markers of cellular function and serve as primary targets for therapeutic intervention.

## Manuscript

([link]())


## Questions & Problems

If you have any questions or problems when using cTPnet or ctpnetpy, please feel free to open a new issue [here](https://github.com/zhouzilu/cTPnet/issues). You can also email the maintainers of the corresponding packages -- the contact information is shown under Developers & Maintainers.


## Installation

Install to R/RStudio
Install all packages in the latest version of [R](https://www.r-project.org/).
First, install the supporting Python package ctpnetpy. See the source code of the package [here](http://github.com/zhouzilu/ctpnetpy)

```python
pip install ctpnetpy
```

Next, open R and install the R package cTPnet
```r
devtools::install_github("zhouzilu/cTPnet")
```


## Pipeline overview

This cTPnet package includes two analysis tools: (1) **SAVERX**, we stronly recommend denoise the scRNA-seq data before impute the surface protein abundance, and (2) **cTP-net**, which impute the surface protein abundance based on previously trained model. 

### SAVERX pipeline

See Jingshu's github for more details

### cTP-net pipeline

<p align="center">
  <img src='https://raw.githubusercontent.com/zhouzilu/DENDRO/master/figure/Pkg_FIG-01.jpg' width='1000' height='600'>
  </p>

  **Figure 1.** A flowchart outlining the procedures of cTPnet. **Need update after submit manuscripts**.

### Running cTP-net

  **cTP-net R notebook** with step-by-step demonstration and rich display is available [***here***](http://rawgit.com/zhouzilu/DENDRO/master/vignette/DENDRO_vignette.html). Corresponding **Rmd script** is available [***here***](https://github.com/zhouzilu/DENDRO/blob/master/vignette/DENDRO_vignette.Rmd).


## Citation

Please cite cTP-net.

* **cTP-net**: [link](https://www.biorxiv.org/content/early/2018/10/30/457622)
<br>
  Genetic Heterogeneity Profiling by Single Cell RNA Sequencing ([GitHub](https://github.com/zhouzilu/DENDRO))

## Developers & Maintainers

* [Zilu Zhou](https://statistics.wharton.upenn.edu/profile/zhouzilu/) (zhouzilu at pennmedicine dot upenn dot edu)
  <br>
  Genomics and Computational Biology Graduate Group, University of Pennsylvania

* [Nancy R. Zhang](https://statistics.wharton.upenn.edu/profile/nzh/) (nzh at wharton dot upenn dot edu)
  <br>
  Department of Statistics, University of Pennsylvania