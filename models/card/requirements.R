install.packages('devtools')
install.packages('reticulate')
reticulate::install_miniconda(update = TRUE, force = FALSE)
reticulate::conda_install(packages = "numpy==1.26.1", pip = TRUE)
reticulate::conda_install(packages = "hdf5plugin==4.2.0", pip = TRUE)
reticulate::conda_install(packages = "h5py==3.10.0", pip = TRUE)
reticulate::conda_install(packages = "anndata==0.10.2", pip = TRUE)
install.packages('BiocManager')
BiocManager::install('Biobase', ask=F)
BiocManager::install('rhdf5', ask=F)
BiocManager::install('anndata', ask=F)
BiocManager::install('SummarizedExperiment', ask=F)
BiocManager::install('SingleCellExperiment', ask=F)
BiocManager::install('Seurat', ask=F)
BiocManager::install('TOAST', ask=F)
devtools::install_github('YingMa0107/CARD')
devtools::install_github('xuranw/MuSiC')
install.packages('R.utils')
install.packages('yaml')
install.packages('plumber')
install.packages('plyr')
install.packages('future')