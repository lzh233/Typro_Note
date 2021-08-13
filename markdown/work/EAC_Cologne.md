# EAC_Cologne

## 基本分析

```r
#packages
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(SingleR)
library(celldex)
#读入原始数据
OE33_1 <- LoadH5Seurat("convert_h5ad/OE33_1_new.h5seurat",verbose = F)
OE33_2  <- LoadH5Seurat("convert_h5ad/OE33_2_new.h5seurat",verbose = F)
#提取表达矩阵
mtx_OE33_1 <- OE33_1@assays$RNA@counts
mtx_OE33_2 <- OE33_2@assays$RNA@counts
#分别构建10X对象
OE33_1_s <- CreateSeuratObject(counts = mtx_OE33_1,,project = "OE33_1",min.cells = 3,min.features = 200)
OE33_2_s <- CreateSeuratObject(counts = mtx_OE33_2,,project = "OE33_2",min.cells = 3,min.features = 200)
#合并数据
OE33_merge <- merge(OE33_1_s,OE33_2_s)
#质控指标计算
#红细胞含量
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.gene.ex <- rownames(OE33_merge)[rownames(OE33_merge) %in% HB.genes]
OE33_merge[["HB.percent"]] <- PercentageFeatureSet(OE33_merge,pattern = HB.gene.ex)
#计算线粒体含量
OE33_merge[["mt.percent"]] <- PercentageFeatureSet(OE33_merge,pattern = "^MT-")
#查看基本指标
VlnPlot(OE33_merge,features = c('nCount_RNA','nFeature_RNA','HB.percent','mt.percent'))
##查看count和feature的关系
plot1 <- FeatureScatter(OE33_merge,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
#查看count和线粒体含量的关系
plot2 <- FeatureScatter(OE33_merge,feature1 = "nCount_RNA",feature2 = "mt.percent")
#可视化
plot1 + plot2
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803141732161.png" alt="image-20210803141732161" style="zoom:50%;" />

```r
#质控
#进行质控，设定质控参数
minGene=200
maxGene=4000
pctMT=30
pcHB=0
OE33_merge <- subset(OE33_merge,subset = 
                   nFeature_RNA > minGene & 
                   nFeature_RNA < maxGene & 
                   mt.percent <= pctMT & 
                   HB.percent == 0)
OE33_merge
'An object of class Seurat 
16160 features across 3485 samples within 1 assay 
Active assay: RNA (16160 features, 0 variable features)''
```

```r
#数据归一化
OE33_merge <- NormalizeData(OE33_merge,normalization.method = "LogNormalize",scale.factor = 10000)
#高变基因
OE33_merge <- FindVariableFeatures(OE33_merge,selection.method = "vst",nfeatures = 2000)
#标准化
OE33_merge <- ScaleData(OE33_merge,verbose = F)
#PCA
OE33_merge <- RunPCA(OE33_merge,verbose = F)
DimPlot(OE33_merge,reduction = "pca")
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803141854622.png" alt="image-20210803141854622" style="zoom:80%;" />

```r
#确定主成分轴
ElbowPlot(OE33_merge)
#聚类
OE33_merge <- FindNeighbors(OE33_merge,dims = 1:20,verbose = F)
OE33_merge <- FindClusters(OE33_merge,resolution = 0.6)
#tsne/umap
OE33_merge <- RunTSNE(OE33_merge)
OE33_merge <- RunUMAP(OE33_merge,verbose = F,dims = 1:20)
```

```r
DimPlot(OE33_merge,reduction = "tsne")
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803141948520.png" alt="image-20210803141948520" style="zoom:80%;" />

```r
DimPlot(OE33_merge,reduction = "umap")
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803142009965.png" alt="image-20210803142009965" style="zoom:80%;" />

```r
#注释
   'BlueprintEncodeData, DatabaseImmuneCellExpressionData,
    HumanPrimaryCellAtlasData, ImmGenData, MonacoImmuneData,
    MouseRNAseqData, NovershternHematopoieticData'
ref <- HumanPrimaryCellAtlasData()
Ana_data <- GetAssayData(OE33_merge,slot = "data")
clusters <- OE33_merge@meta.data$seurat_clusters
cell_type <- SingleR(test = Ana_data, 
                     ref = ref, 
                    label = ref$label.main,
                    clusters = clusters)
cell.name <- data.frame(clusters=rownames(cell_type),cellname=cell_type$labels)
new.cluster.ids <- cell.name[,2]
names(new.cluster.ids) <- levels(OE33_merge)
OE33_merge <- RenameIdents(OE33_merge,new.cluster.ids)
OE33_merge[["cell.type"]] <- Idents(OE33_merge)

p1 <- DimPlot(OE33_merge,reduction = "umap",group.by = "orig.ident")
p2 <- DimPlot(OE33_merge,reduction = "umap")
cowplot::plot_grid(p1,p2,nrow = 2)
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803142130236.png" alt="image-20210803142130236" style="zoom:80%;" />



## therapy resistance

```
Single-cell Transcriptome Analyses Reveal Molecular Signals to Intrinsic and Acquired Paclitaxel Resistance in Esophageal Squamous Cancer Cells
```

### 数据分析部分

### 差异基因和therapy resistance相关通路的分析

首先使用三种算法，`Deseq`、`Monocle`、`SCDE`得到三者overlap的DEGs，后通过IPA(Ingenuity Pathway Analysis)([跟着Cell学作图| 11.Ingenuity Pathway Analysis(IPA) - 知乎 (zhihu.com)](https://zhuanlan.zhihu.com/p/374326403))

```
IPA（Ingenuity Pathway Analysis，通路分析软件）是一款基于云计算的图形化界面生物信息学软件，能够从生物学通路角度将组学数据进行分析、整合和理解，适用于转录组学、蛋白质组学、代谢组学等大数据分析，也适用于一些产生基因、化学物质列表的小规模实验。组学数据分析结果主要有**经典通路**、**上游转录调控**、**下游调控子效应**、**疾病与功能**以及**分子间相互作用网络**这5个结果。
```

![image-20210803143127859](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803143127859.png)

得到在`canonical pathways`中，有如下几个经典通路在紫杉醇抗药性细胞里被抑制

![image-20210803143830127](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803143830127.png)

在上游转录调控分析里，结果如下

![image-20210803144133150](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803144133150.png)

![image-20210803144154024](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803144154024.png)

在细胞互做网络分析中，TP53调控细胞死亡和周期，Histone H3 和CLDN7.....



![image-20210803144352184](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803144352184.png)

### Gene-network modules underlying diverse paclitaxel resistance identified by WGCNA

首先整体比较了单细胞数据的转录组差异KYSE-30 cells (30_S1 and 30_S2) and Taxol-R cells (Taxol-R_S1 and Taxol-R_S2). 基因聚类成了七个模块

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803145235692.png" alt="image-20210803145235692" style="zoom:80%;" />

比较了各个亚群与不同基因模块之间的相关性

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803145548473.png" alt="image-20210803145548473" style="zoom:80%;" />

 `DAVID-KEGG `分析表明 proteasome通路在Turquoise模块中受到显著调控



<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803145657002.png" alt="image-20210803145657002" style="zoom:80%;" />

关键基因网络分析得到一批关键基因

![image-20210803150037864](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803150037864.png)

![image-20210803161200845](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803161200845.png)



Turquoise模块中存在的许多基因被报道过

![image-20210803160012687](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210803160012687.png)







