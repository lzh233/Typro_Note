# 分析环境配置

```R
install.packages("BiocManager")
install.packages("devtools")
install.packages("tidyverse")
install.packages("Seurat")
install.packages("vegan")
install.packages("cowplot")
install.packages("ggnewscale")
install.packages("pheatmap")
install.packages("ape")
install.packages("car")
install.packages("nycflights13")
install.packages("patchwork")
install.packages("harmony")

BiocManager::install("SingleCellExperiment")
BiocManager::install("scater")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("ggsci")
BiocManager::install("ggtree")
BiocManager::install("ggtreeExtra")
BiocManager::install("phyloseq")
BiocManager::install("treeio")
BiocManager::install("glmGamPoi")
BiocManager::install("MAST")

#SingleR数据库
MonacoImmuneData()
HumanPrimaryCellAtlasData() 
BlueprintEncodeData() 
MouseRNAseqData() 
ImmGenData() 
DatabaseImmuneCellExpressionData()
#安装SeuratData, 整合了很多单细胞数据集, 使用AvailableData()可以查看数据集
devtools::install_github('satijalab/seurat-data')
#library(SeuratData)
#查看数存在的据集
#AvailableData()
#安装某一数据集
#InstallData("ifnb")
```

- **`BiocManager::install("celldex")`**

报错处理

```shell
vendor/boost/math/tools/config.hpp:408:13: fatal error: boost/detail/fenv.hpp: No such file or direc
#解决办法
#原因 是因为RSQLite报错
conda install -c conda-forge boost-cpp
```

