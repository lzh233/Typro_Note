# Seurat-多样本整合

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-5-6
sessionInfo("Seurat")
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)
version: Seurat_4.0.1
```

## 资料来源

https://mp.weixin.qq.com/s/dz0dJNf2slhRNd-9bG188g

## 工作目录设置与查看数据

**目的**：将会使用3种方法对数据进行合并, 并比较区别

- 直接将数据读入构建10个数据的10x对象
- 使用R基础函数`merge()`进行直接合并
- 使用`Seurat`提供的**CCA+MNN**进行数据整合

```r
> setwd("D:\\Desktop\\s_note\\data\\singel_cell\\IntegrateData")
#共有10组数据，来自Immune Landscape of Viral- and Carcinogen-Driven Head and Neck Cancer，数据集GEO编号：GSE139324
> dir()
 [1] "GSM4138110" "GSM4138111" "GSM4138128" "GSM4138129" "GSM4138148" "GSM4138149" "GSM4138162" "GSM4138163"
 [9] "GSM4138168" "GSM4138169"
#每组数据均为10x下机数据标准格式
> dir(".\\GSM4138110")
[1] "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"
```

## 载入包并读入数据

```R
library(Seurat)
library(tidyverse)
library(patchwork)
library(celldex)
library(SingleR)
#由于配置不太行，所以就只读入三组吧......
#构建文件列表
data.list <- dir() %>% str_c(".\\",.)
data.list <- data.list[2:4]
names(data.list) = c('HNC01TI', 'HNC10PBMC', 'HNC10TIL')
#names(data.list) = c('HNC01PBMC', 'HNC01TIL', 'HNC10PBMC', 'HNC10TIL', 'HNC20PBMC', 
#               'HNC20TIL', 'PBMC1', 'PBMC2', 'Tonsil1', 'Tonsil2')
```

使用`Read10X()`读入数据

```r
counts <- Read10X(data.dir = data.list)
```

## 构建seurat对象

```r
#构建对象，并进行简单质控
scRNA<- CreateSeuratObject(counts = counts,project = "TEST2",min.cells = 3,min.features = 200)
> scRNA
An object of class Seurat 
17975 features across 4431 samples within 1 assay 
Active assay: RNA (17975 features, 0 variable features)
```

## 直接组合数据进行下游分析

### 数据质控

```r
#进行下游分析
#质控，
#计算每个cell中线粒体基因的含量
scRNA[["mt.percent"]] <- PercentageFeatureSet(scRNA,pattern = "^MT-")
#计算每个cell中红细胞基因的含量
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
#使用match函数匹配上述HB.gene是否在基因列表中，有则返回基因索引值，没有则返回NA，后利用str_subset()剔除NA将得到的index转换为数值，提取基因名
#str_subset("str",".")表示输出除了NA外所有字符
HB_m <- rownames(scRNA@assays$RNA)[str_subset(match(HB.genes, rownames(scRNA@assays$RNA)),".") %>% as.numeric()]
#HB_m <- <- rownames(scRNA@assays$RNA)[HB.genes %in% rownames(scRNA@assays$RNA)]
#每个cell中的红细胞数量
scRNA[["HB.gene"]] <- PercentageFeatureSet(scRNA,features = HB_m)
#对meta增加project列
scRNA[["project"]] <- rep("Test1",nrow(scRNA@meta.data))
#进行质控，设定质控参数
minGene=500
maxGene=3000
pctMT=10
pcHB=0
scRNA <- subset(scRNA,subset = 
                   nFeature_RNA > minGene & 
                   nFeature_RNA < maxGene & 
                   percent.mt <= pctMT & 
                   percent.HB <= 0)
```

### 数据归一化以及高变基因筛选

```r
#数据归一化与标准化
scRNA <- NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
#高变基因筛选
scRNA <- FindVariableFeatures(scRNA,selection.method = "vst",nfeatures = 2000)
#查看高变基因
plot1 <- VariableFeaturePlot(scRNA)
top10 <- head(VariableFeatures(scRNA),10)
LabelPoints(plot = plot1,points = top10,repel = T)
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210517095516349.png" alt="image-20210517095516349" style="zoom:60%;" />

### PCA分析

```r
#数据标准化
scRNA <- ScaleData(scRNA,features = VariableFeatures(immune.combine))
#PCA，细胞并没有很好地融合在一起
scRNA <- RunPCA(scRNA)
DimPlot(scRNA,reduction = "pca",group.by = "orig.ident")
#确定后续分析选用的主成分数量
scRNA <- JackStraw(scRNA,num.replicate = 100)
scRNA <- ScoreJackStraw(scRNA,dims = 1:20)
JackStrawPlot(scRNA,dims = 1:20)
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210517101043438.png" alt="image-20210517101043438" style="zoom:50%;" />

### 聚类分析

```r
scRNA <- FindNeighbors(scRNA,dims = 1:30)
scRNA <- FindClusters(scRNA,resolution = 0.5)
table(Idents(scRNA))
#查看聚类情况
> table(Idents(scRNA))
  0   1   2   3   4   5   6   7   8   9  10  11  12  13 
506 485 483 428 416 376 273 268 246 170 167 107  76  40 
```

### 细胞类型注释

```r
#SingleR的安装
BiocManager::install("SingleR")
#安装数据库
BiocManager::install("celldex")
#下载数据库
refdata <- MonacoImmuneData()
ref <- HumanPrimaryCellAtlasData() 
ref <- BlueprintEncodeData() 
ref <- MouseRNAseqData() 
ref <- ImmGenData() 
ref <- DatabaseImmuneCellExpressionData() 
ref <- NovershternHematopoieticData() 
```

```R
#载入需要注释的数据
testdata <- GetAssayData(scRNA, slot="data")
#载入数据库
refdata <- MonacoImmuneData()
#提取cluster
clusters <- scRNA@meta.data$seurat_clusters
#细胞类型鉴定
celltype <- SingleR(test = testdata,
                     ref= refdata,
                     label =refdata$label.main,
                     clusters = clusters1,
                     assay.type.test = "logcounts", 
                     assay.type.ref = "logcounts")
#提取注释信息
cell.name <- data.frame(clusters=rownames(celltype1),cellname=celltype1$labels)
#对cell类型进行注释
#细胞注释
new.cluster.ids <- cell.name
names(new.cluster.ids) <- levels(testdata)
testdata <- RenameIdents(testdata,new.cluster.ids)
#写入metadata
testdata[["cell.type"]] <- Idents(testdata)
#可视化
p1 <- DimPlot(scRNA,reduction = "umap",group.by = "cell.type",label = T)
p2 <- DimPlot(scRNA,reduction = "umap",group.by = "orig.ident",label = T)
p1+p2
#可以看出相同或的cell并没有很好的融合在一起
```

![image-20210517104410463](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210517104410463.png)



```r
#查看t-sne结果，同样批次效应明显
scRNA <- RunTSNE(scRNA,dims = 1:30)
p3 <- DimPlot(scRNA,reduction = "tsne",group.by = "cell.type",label = T)
p4 <- DimPlot(scRNA,reduction = "tsne",group.by = "orig.ident",label = T)
p3 + p4
```

![image-20210517105140500](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210517105140500.png)

## 使用seurat提供的数据整合算法

### 分割数据并整合

使用`Read10x()`读入数据后构建seurat对象后使用`SplitObject()`根据指定的因子分割对象

```r
counts <- Read10X(data.dir = data.list)
scRNA2 <- CreateSeuratObject(counts = counts,project = "TEST2",min.cells = 3,min.features = 200)
scRNA2 <- SplitObject(scRNA2,split.by = "orig.ident")
#返回包含3个seurat对象的列表
> scRNA2
$HNC01TI
An object of class Seurat 
17975 features across 1298 samples within 1 assay 
Active assay: RNA (17975 features, 0 variable features)

$HNC10PBMC
An object of class Seurat 
17975 features across 1750 samples within 1 assay 
Active assay: RNA (17975 features, 0 variable features)

$HNC10TIL
An object of class Seurat 
17975 features across 1383 samples within 1 assay 
Active assay: RNA (17975 features, 0 variable features)
#整合3组数据
#分别对每一组数据标准化，并筛选高变基因
scRNA2 <- lapply(scRNA2, function(obj){
  #归一化并分别选取高变基因
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj,selection.method = "vst",nfeatures = 2000)
})
#整合数据(详见多样本整合)
features <- SelectIntegrationFeatures(object.list = scRNA2)
immune.anchors <- FindIntegrationAnchors(object.list = scRNA2,anchor.features = features)
immune.combine <- IntegrateData(anchorset = immune.anchors)

#将下游分析改为整合后的原始数据，integrated为整合后的保留2000高变基因的数据
DefaultAssay(immune.combine) <- "RNA"
#查看数据情况
> immune.combine
An object of class Seurat 
19975 features across 4431 samples within 2 assays 
Active assay: RNA (17975 features, 0 variable features)
 1 other assay present: integrated
```

### 数据质控

```r
counts <- Read10X(data.dir = data.list)
scRNA2 <- CreateSeuratObject(counts = counts,project = "TEST2",min.cells = 3,min.features = 200)

scRNA2[["percent.mt"]] <- PercentageFeatureSet(scRNA2,pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- rownames(scRNA2@assays$RNA)[str_subset(match(HB.genes, rownames(scRNA2@assays$RNA)),".") %>% as.numeric()]
scRNA2[["percent.HB"]] <- PercentageFeatureSet(scRNA2,features = HB_m)

minGene=500
maxGene=3000
pctMT=10
pcHB=0

immune.combine <- subset(scRNA2,subset = 
                           nFeature_RNA > minGene & 
                           nFeature_RNA < maxGene & 
                           percent.mt <= pctMT & 
                           percent.HB <= 0)
```

### 归一化与高变基因筛选

```R
scRNA2 <- SplitObject(scRNA2,split.by = "orig.ident")
scRNA2 <- lapply(scRNA2, function(obj){
  #归一化并分别选取高变基因
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj,selection.method = "vst",nfeatures = 2000)
})
```

### 数据整合

```R
features <- SelectIntegrationFeatures(object.list = scRNA2)
immune.anchors <- FindIntegrationAnchors(object.list = scRNA2,anchor.features = features)
immune.combine <- IntegrateData(anchorset = immune.anchors)
```

### 数据分析

```R
immune.combine <- ScaleData(immune.combine, verbose = T)
immune.combine <- RunPCA(immune.combine, npcs = 30, verbose = T)
immune.combine <- RunUMAP(immune.combine, reduction = "pca", dims = 1:30)
immune.combine <- FindNeighbors(immune.combine, reduction = "pca", dims = 1:30)
immune.combine <- FindClusters(immune.combine, resolution = 0.5)
immune.combine <- RunTSNE(immune.combine,dims = 1:30)
#注释
testdata <- GetAssayData(immune.combine, slot="data")
#载入数据库
refdata <- MonacoImmuneData()
#提取cluster
clusters <- immune.combine@meta.data$seurat_clusters
#细胞类型鉴定
celltype <- SingleR(test = testdata,
                    ref= refdata,
                    label =refdata$label.main,
                    clusters = clusters,
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")
#提取注释信息
cell.name <- data.frame(clusters=rownames(celltype),cellname=celltype$labels)
#细胞注释
new.cluster.ids <- cell.name[,2]
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.cluster.ids)
#注释信息写入meta.data
pbmc$cell.type <- Idents(pbmc)
#将信息写入meta.data
immune.combine[["cell.type"]] <- clusters
#可视化
p1 <- DimPlot(immune.combine, reduction = "umap", group.by = "orig.ident",label = T)
p2 <- DimPlot(immune.combine, reduction = "umap", label = TRUE, repel = TRUE,group.by = "cell.type")
p1 + p2
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210517142546190.png" alt="image-20210517142546190" style="zoom:67%;" />

```R
p3 <- DimPlot(immune.combine, reduction = "tsne", group.by = "orig.ident",label = T)
p4 <- DimPlot(immune.combine, reduction = "tsne", label = TRUE, repel = TRUE,group.by = "cell.type")
p3 + p4
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210517144030574.png" alt="image-20210517144030574" style="zoom:67%;" />
