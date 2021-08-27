# Seurat-多样本整合(rPCA)

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-7-21
sessionInfo("Seurat")
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)
version: Seurat_4.0.1
```

```
Tim S, Andrew Butler, Paul Hoffman , et al. Comprehensive integration of single cell data[J].Cell,2019.
```

### 数据整合原理

**见Seurat-多样本整合(CCA)**

https://www.jianshu.com/p/32ca61450fe9

通常一下情况可以选择使用***rPCA***的方法对数据进行整合

- 一个数据集中有大部分的cell无法与另一个数据集进行匹配( ***substantial fraction of cells in one dataset have no matching type in the other***)
- 数据集均来自同一平台，如10X( ***Datasets originate from the same platform***)
- 数据集非常的大(***There are a large number of datasets or cells to integrate***)

### 与基于CCA的数据整合的区别

- 在确定两个数据集共享的变异来源时，CCA更适合细胞类型比较保守，但是基因表达差异很大的数据，但是CCA可能会造成过度矫正，尤其是在非常大的数据集已经共享细胞较少的数据集

```
By identifying shared sources of variation between datasets, CCA is well-suited for identifying anchors when cell types are conserved, but there are very substantial differences in gene expression across experiments. CCA-based integration therefore enables integrative analysis when experimental conditions or disease states introduce very strong expression shifts, or when integrating datasets across modalities and species. However, CCA-based integration may also lead to overcorrection, especially when a large proportion of cells are non-overlapping across datasets.
```

- rPCA是一种运行速度更快的方法，同时也更为保守，因此整合后的数据可能并没有那么“整齐”，**同时在执行rPCA寻找锚点时，需要先对每个需要整合的数据运行PCA分析，**因此rPCA的方法更适合**上述三种情况**

```
RPCA-based integration runs significantly faster, and also represents a more conservative approach where cells in different biological states are less likely to ‘align’ after integration. 
```

## 使用rPCA的方法整合数据

### 载入数据与相关包

```r
rm(list = ls())
library(SeuratData)
library(Seurat)
library(tidyverse)
# install dataset
#InstallData("ifnb")
#载入ifnb数据
LoadData("ifnb")
```

### 分割对象并对每组数据进行处理

```r
#分割对象
ifnb.list <- SplitObject(ifnb,split.by = "stim")
#———查看数据----
> ifnb.list
$CTRL
An object of class Seurat 
14053 features across 6548 samples within 1 assay 
Active assay: RNA (14053 features, 0 variable features)

$STIM
An object of class Seurat 
14053 features across 7451 samples within 1 assay 
Active assay: RNA (14053 features, 0 variable features)
#————————----
#归一化，寻找高变基因
ifnb.list <- lapply(ifnb.list, function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x,selection.method = "vst",nfeatures = 2000)
})
#确定用于数据整合的features
features <- SelectIntegrationFeatures(ifnb.list)
#数据标准化和PCA分析
ifnb.list <- lapply(ifnb.list,function(x){
  x <- ScaleData(x,features = features)
  x <- RunPCA(x,features = features)
})
```

### 使用rPCA确定锚点并整合数据

```r
#确定锚点
ifnb.imm <- FindIntegrationAnchors(object.list = ifnb.list,
                                   anchor.features = features,
                                   reduction = "rpca")
#整合数据
ifnb.immdata <- IntegrateData(ifnb.imm)
```

### 整合数据标准分析

```r
#设定分析用的数据集
DefaultAssay(ifnb.immdata) <- "integrated"
#标准分析
ifnb.immdata <- ScaleData(ifnb.immdata, verbose = FALSE)
    ifnb.immdata <- RunPCA(ifnb.immdata, npcs = 30, verbose = FALSE)
ifnb.immdata <- RunUMAP(ifnb.immdata, reduction = "pca", dims = 1:30)
ifnb.immdata <- FindNeighbors(ifnb.immdata, reduction = "pca", dims = 1:30)
ifnb.immdata <- FindClusters(ifnb.immdata, resolution = 0.5)

p1 <- DimPlot(ifnb.immdata, reduction = "umap", group.by = "stim")
p2 <- DimPlot(ifnb.immdata, reduction = "umap", group.by = "seurat_annotations", label = TRUE,
              repel = TRUE)
p1 + p2
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210721151929453.png" alt="image-20210721151929453" style="zoom:80%;" />

从上图可以看出，与CCA的整合相比，rPCA的整合结果要略差于CCA，但是可以通过调整`FindIntegrationAnchors()`的参数进一步优化整合结果

### 优化rPCA寻找锚点的结果

通过调整`k.anchor`参数, 默认为5

**`k.anchor`**： How many neighbors (k) to use when picking anchors, 5

```r
#调整k.anchor
ifnb.imm <- FindIntegrationAnchors(object.list = ifnb.list,
                                   anchor.features = features,
                                   reduction = "rpca",
                                   k.anchor = 20)
ifnb.immdata <- IntegrateData(ifnb.imm)


DefaultAssay(ifnb.immdata) <- "integrated"

# Run the standard workflow for visualization and clustering
ifnb.immdata <- ScaleData(ifnb.immdata, verbose = FALSE)
ifnb.immdata <- RunPCA(ifnb.immdata, npcs = 30, verbose = FALSE)
ifnb.immdata <- RunUMAP(ifnb.immdata, reduction = "pca", dims = 1:30)
ifnb.immdata <- FindNeighbors(ifnb.immdata, reduction = "pca", dims = 1:30)
ifnb.immdata <- FindClusters(ifnb.immdata, resolution = 0.5)

p1 <- DimPlot(ifnb.immdata, reduction = "umap", group.by = "stim")
p2 <- DimPlot(ifnb.immdata, reduction = "umap", group.by = "seurat_annotations", label = TRUE,
              repel = TRUE)
p1 + p2
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210721152947634.png" alt="image-20210721152947634" style="zoom: 67%;" />

**可以看出整合的效果要好于默认参数**

## 使用SCTransform对数据标准化

**上服务器**

`NormalizeData()`、`FindVariableFeatures()`和`ScaleData()`被`SCTransform`替代

```r
library(glmGamPoi)
LoadData("ifnb")
ifnb.list <- SplitObject(ifnb, split.by = "stim")
#执行SCTransform标准化
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform, method = "glmGamPoi")
#选择用于整合的features
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
#Prepare an object list normalized with sctransform for integration.
?PrepSCTIntegration
#PCA分析
ifnb.list <- lapply(X = ifnb.list, FUN = RunPCA, features = features)
```

### 数据整合

- 在常规分析中，使用少量的PC既能关注到关键的生物学差异，又能够不引入更多的技术差异，相当于一种保守性的做法。是的，它会失去一些生物差异信息，但是同时又在常规手段中比较安全。

- 这里使用的`sctransform`，显然更“自信”一些，能提取更多的生物差异，并且兼顾不引入技术误差

常规分析中的`FindVariableFeatures`默认得到2000个高变异基因（HVGs），而**这里的`sctransform`因为使用了更多的PCs，算法也更优化，所以默认会得到3000个HVGs**。sctransform认为：新增加的这1000个基因就包含了之前没有检测到的微弱的生物学差异。而且，即使使用全部的全部的基因去做下游分析，得到的结果也是和`sctransform`这3000个基因的结果相似

```R
#确定锚点
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, 
                                         normalization.method = "SCT",
                                         anchor.features = features, 
                                         dims = 1:30, 
                                         reduction = "rpca", 
                                         k.anchor = 5)
#数据整合
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)
```

### 整合后的数据UMAP分析

```r
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_annotations", label = TRUE,
    repel = TRUE)
cowplot::plot_grid(p1,p2,nrow = 2)
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210721163603858.png" alt="image-20210721163603858" style="zoom:80%;" />
