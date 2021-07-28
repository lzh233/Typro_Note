# Seurat-多组学数据整合

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-5-6
R version 4.0.3 (2020-10-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)
```

## 多组学数据

是指对同一个细胞同时测量多种指标， 如[CITE-seq](http://www.nature.com/nmeth/journal/v14/n9/full/nmeth.4380.html) ：是一种可以同时检测细胞表面蛋白和RNA的单细胞测序。一方面用带DNA标签的抗体来检测细胞表面蛋白，同时用10x的RNA单细胞测序来检测全部的mRNA。scRNA-seq+scATAC-seq，本例子主要是演示了如何在单细胞转录组cluster后的基础上去分析每个cell上的表面蛋白的含量

## 数据下载

***Here, we analyze a dataset of 8,617 cord blood mononuclear cells (CBMCs), where transcriptomic measurements are paired with abundance estimates for 11 surface proteins,***

```shell
#转录组数据
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz
#表面蛋白数据
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz
```

## 载入相关包并数据读入

**由于数据量太大，所以上服务器**

```r
library(Seurat)
library(tidyverse)
library(cowplot)
setwd("~/multi_data/")
#查看文件列表，两个文件
#'GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz'：表面蛋白的数据
#'GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz'：转录组的数据
#读入转录组的数据
cbmc.rna <- as.sparse(read.csv("GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",header = T,row.names = 1))
#仅仅提取人类的数据，CollapseSpeciesExpressionMatrix，当数据是由多物种组成的时候，可以提取某物种的数据，本数据分为人和老鼠，人的数据默认为HUMAN_老鼠尾MOUSE_，提取人类的前100的表达量基因
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna,prefix = "HUMAN_",controls = "MOUSE_",ncontrols = 100)
#读入表面蛋白的数据
cmbc.adt <- as.sparse(read.csv("GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz",header=T,row.names = 1))
#检查二者的cell是否一致
>all.equal(colnames(cbmc.rna),colnames(cmbc.adt))
TRUE
```

## 创建Seurat对象并将表面蛋白数据加入

```shell
#创建seurat对象
cbmc <- CreateSeuratObject(counts = cbmc.rna,project = "multi_test")
#将表面蛋白数据创建为seurat对象CreateAssayObject,并将其写入cbmc
adt_assay <- CreateAssayObject(counts = cmbc.adt)
cbmc[["ADT"]] <- adt_assay

###查看cbmc内的数据
>Assays(cbmc)
'RNA''ADT'
>cbmc[["ADT"]]
Assay data with 13 features for 8617 cells
First 10 features:
 CD3, CD4, CD8, CD45RA, CD56, CD16, CD10, CD11c, CD14, CD19 
>cbmc[["RNA"]]
Assay data with 36240 features for 8617 cells
First 10 features:
 ERCC-ERCC-00104, HUMAN-A1BG, HUMAN-A1BG-AS1, HUMAN-A1CF, HUMAN-A2M,
HUMAN-A2M-AS1, HUMAN-A2ML1, HUMAN-A4GALT, HUMAN-A4GNT, HUMAN-AAAS 
```

## 对转录组数据分析

```r
#设定分析数据
DefaultAssay(cbmc) <- "RNA"
#对数据标准分析，质控 归一化 高变基因筛选 标准化 降维 聚类.....
cbmc <- NormalizeData(cbmc,normalization.method = "LogNormalize",scale.factor = 10000,verbose = T)
cbmc <- FindVariableFeatures(cbmc,selection.method = "vst",nfeatures = 2000,verbose = T)
cbmc <- ScaleData(cbmc,features = VariableFeatures(cbmc),verbose = T)
cbmc <- RunPCA(cbmc,verbose = T)
cbmc <- FindNeighbors(cbmc,dims = 1:30,verbose = T)
cbmc <- FindClusters(cbmc,resolution = 0.8,verbose = T)
cbmc <- RunUMAP(cbmc,dims = 1:30,verbose = T)
p1 <- DimPlot(cbmc,label = T)
p1
```

![image-20210721093206893](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210721093206893.png)

## 表面蛋白数据的可视化

现在已经得到了转录组的结果，因此下一步可以通过对表面蛋白的数据进行分析，首先需要对表面蛋白的数据进行归一化处理

```R
#"CD19" %in% rownames(cbmc) 
#查看CD19 RNA在cell中的表达
DefaultAssay(cbmc) <- "RNA"
cd19.rna <- FeaturePlot(cbmc,reduction = "umap",features = "CD19",cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 RNA")
#查看CD19 蛋白在cell中的表达
DefaultAssay(cbmc) <- "ADT"
cd19.pro <- FeaturePlot(cbmc,reduction = "umap",features = "CD19",cols = c("gray", "red")) + ggtitle("CD19 pro")
cd19.rna + cd19.pror
```

![image-20210721094828924](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210721094828924.png)

```r
#可以使用Key为ADT和RNA数据添加前缀
Key(cbmc[["ADT"]])
Key(cbmc[["RNA"]])
##返回结果为前缀
'adt_'
'rna_'
##无需切换数据即可查看CD19的表达情况
p3 <- FeaturePlot(cbmc,features = "adt_CD19",cols = c("lightgrey", "darkgreen"))+ ggtitle("CD19 ADT")
p4 <- FeaturePlot(cbmc,features = "adt_CD19")+ ggtitle("CD19 RNA")
p3 + p4
```

![image-20210721100139495](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210721100139495.png)

## 查看表面蛋白在不同的cluster上的表达情况

```r
#CD19 是B cell的常见marker
p5 <- VlnPlot(cbmc,features = "adt_CD19")
p5
```

![image-20210721100921065](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210721100921065.png)



## 基因/蛋白在各个cluster的 表达情况

```R
#查看各个表面蛋白在不同cell中的表达情况
adt.markers <- FindMarkers(cbmc,ident.1 = 6,assay = "ADT")

>head(adt.markers)
A data.frame: 7 × 5
p_val	avg_log2FC	pct.1	pct.2	p_val_adj
<dbl>	<dbl>	<dbl>	<dbl>	<dbl>
CD19	2.067533e-215	1.2787751	1	1	2.687793e-214
CD45RA	8.106076e-109	0.4117172	1	1	1.053790e-107
CD4	1.123162e-107	-0.7255977	1	1	1.460110e-106
CD14	7.212876e-106	-0.5060496	1	1	9.376739e-105
CD3	1.639633e-87	-0.6565471	1	1	2.131523e-86
CD8	1.042859e-17	-0.3001131	1	1	1.355716e-16
CD11c	8.957964e-11	-0.4382277	1	1	1.164535e-0
#查看各个基因在不同cell中的表达情况
rna.marker <- FindMarkers(cbmc,ident.1 = 5,assay = "RNA")

>head(rna.marker)
A data.frame: 6 × 5
p_val	avg_log2FC	pct.1	pct.2	p_val_adj
<dbl>	<dbl>	<dbl>	<dbl>	<dbl>
AC109351.1	0	0.3203893	0.265	0.005	0
CTD-2090I13.1	0	2.0024376	0.972	0.062	0
DCAF5	0	0.6637418	0.619	0.055	0
DYNLL2	0	2.0387603	0.984	0.094	0
FAM186B	0	0.3000479	0.244	0.002	0
HIST2H2AB	0	1.3104432	0.812	0.013	0
```

## 查看各个基因/蛋白表达的相关性

```r
p5 <- FeatureScatter(cbmc,feature1 = "adt_CD19", feature2 = "adt_CD3")
p6 <- FeatureScatter(cbmc, feature1 = "adt_CD3", feature2 = "rna_CD3E")
p7 <- FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8")
#选择原始数据进行可视化"adt_CD4"和"adt_CD8"表达量之间的关系
p8 <- FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8", slot = "counts")
plot_grid(p5,p6,p7,p8,labels = c("adt_CD19 vs adt_CD3","adt_CD3 vs rna_CD3E","adt_CD4 vs adt_CD8","no-Normalize"))
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210721104209565.png" alt="image-20210721104209565" style="zoom:80%;" />

