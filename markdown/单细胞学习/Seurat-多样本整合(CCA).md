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

```
Tim S, Andrew Butler, Paul Hoffman , et al. Comprehensive integration of single cell data[J].Cell,2019.
```

### 数据整合原理

https://mp.weixin.qq.com/s/dz0dJNf2slhRNd-9bG188g

1、使用CCA分析将两个数据集降维到同一个低维空间，因为CCA降维之后的空间距离不是相似性而是相关性，所以相同类型与状态的细胞可以克服技术偏倚重叠在一起。CCA降维之后细胞在低维空间有了可以度量的“距离”，MNN(mutual nearest neighbor)算法以此找到两个数据集之间互相“距离”最近的细胞，**Seurat将这些相互最近邻细胞称为“锚点细胞”**。假设：

- A样本中的细胞A3与B样本中距离最近的细胞有3个（B1,B2,B3）
- B样本中的细胞B1与A样本中距离最近的细胞有4个（A1,A2,A3,A4）
- B样本中的细胞B2与A样本中距离最近的细胞有2个（A5,A6）
- B样本中的细胞B3与A样本中距离最近的细胞有3个（A1,A2,A7）

那么A3与B1是相互最近邻细胞，A3与B2、B3不是相互最近邻细胞，A3+B1就是A、B两个数据集中的锚点之一。实际数据中，两个数据集之间的锚点可能有几百上千个



<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210507103020730.png" alt="image-20210507103020730" style="zoom: 67%;" />

理想情况下相同类型和状态的细胞才能构成配对锚点细胞，但是异常的情况也会出现，如上图中query数据集中黑色的细胞团。它在reference数据集没有相同类型的细胞，但是它也找到了锚点配对细胞（红色连线）。Seurat会通过两步过滤这些不正确的锚点（DE）：

1. 在CCA低维空间找到的锚点，返回到基因表达数据构建的高维空间中验证，如果它们的转录特征相似性高则保留，否则过滤此锚点。

2. 检查锚点细胞所在数据集最邻近的30个细胞，查看它们重叠的锚点配对细胞的数量，重叠越多分值越高，代表锚点可靠性更高。

   经过层层过滤剩下的锚点细胞对，可以认为它们是相同类型和状态的细胞，它们之间的基因表达差异是技术偏倚引起的。Seurat计算它们的差异向量，然后用此向量校正这个锚点锚定的细胞子集的基因表达值。校正后的基因表达值即消除了技术偏倚，实现了两个单细胞数据集的整合。

## 数据读入

**本例子的目的是确定在刺激和不刺激两种状态下，在不同cell中表达相似的基因，并可视化**

**载入相关包**

```R
library(tidyverse)
library(cowplot)
library(Seurat)
#安装SeuratData, 整合了很多单细胞数据集, 使用AvailableData()可以查看数据集
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(patchwork)
#查看数存在的据集
#AvailableData()
#安装某一数据集
#InstallData("ifnb")
```

**载入数据集，`ifnb`**

```r
#数据集基本信息，数据分为两个处理，分别为刺激(Stim)和对照(Control)
ifnb.SeuratData     ifnb   3.1.0   IFNB-Stimulated and Control PBMCs   <NA>      human

#加载数据集
LoadData("ifnb")
> LoadData("ifnb")
An object of class Seurat 
14053 features across 13999 samples within 1 assay 
Active assay: RNA (14053 features, 0 variable features)
```

使用**`SplitObject()`**函数将数据按`metadata`中的数据进行分组, 本例是按照刺激`(Stim)`和对照`(Control)`将数据分成两组, 该函数的返回值为一个列表，包含了分出来的所有Seurat对象

```R
ifnb.list <- SplitObject(ifnb,split.by = "stim")
#查看结果
> class(ifnb.list)
[1] "list"
#数据被分为两组数据集，查看一下两组数据情况
> ifnb.list$CTRL@assays$RNA
Assay data with 14053 features for 6548 cells
Top 10 variable features:
 HBB, HBA2, HBA1, CCL3, CCL4, CXCL10, TXN, CCL7, CCL2, GNLY 
> ifnb.list$STIM@assays$RNA
Assay data with 14053 features for 7451 cells
Top 10 variable features:
 HBB, HBA2, HBA1, CCL4, APOBEC3B, CCL7, GNLY, CCL3, TXN, PPBP 
```

## 数据归一化并分别选取高变基因

使用`lapply()` `NormalizeData() `  `FindVariableFeatures()` ，进行批量归一化并分别选取高变基因

```R
ifnb.list <- lapply(ifnb.list, function(obj){
  #归一化并分别选取高变基因
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj,selection.method = "vst",nfeatures = 2000)
})
```

使用`SelectIntegrationFeatures()`，选择用于合并的features，**# select features that are repeatedly variable across datasets for integration？**

```r
features <- SelectIntegrationFeatures(object.list = ifnb.list)
```

### 寻找并确定锚点(anchors)

使用`FindIntegrationAnchors()`寻找`anchors`, 选项与参数：

**`object.list`**: 存放用于合并的seurat对象的列表

**`assay`** 指定使用列表中哪个seurat对象

**`reference`**: 指定数据整合过程中，哪个对象用于参考，默认`NULL`，当值为`NULL`时，所有的`anchors `都会被确定，如果指定了`reference`, 会首先从每个			  `reference`和`query`中确定`anchors`，后所有的`reference`会进行整合，然后`query`再与整合后的`reference`整合。

**`anchor.features`**: 指定使用的`features`

**`scale`**: 默认对数据执行标准化，如果过提前执行了则可设置为FALSE

**`normalization.method`**: 指定归一化方法，`LogNormalize` 或 `SCT`, 若是已经执行过则无需设置

**`sct.clip.range`**: 划定`anchors`相关性阈值

**`reduction`**： 指定降维方法，`cca`或`rpca`

**`l2.norm`**: CCA分析后对数据L2归一化？,默认 ture

**`dims`**: 选择的维数，默认 1:30

**`k.anchor`**： How many neighbors (k) to use when picking anchors, 5

**`k.filter`**: How many neighbors (k) to use when filtering anchors, 200

**`k.score`**: How many neighbors (k) to use when scoring anchors, 30

.........

```r
?FindIntegrationAnchors()
#选取anchors
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list,anchor.features = features)
#输出信息
Scaling features for provided objects
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
Finding all pairwise anchors
  |                                                  | 0 % ~calculating  Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 16393 anchors
Filtering anchors
	Retained 6756 anchors
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=04m 23s
```

### 整合数据

使用`IntegrateData()`对数据进行整合

```R
?IntegrateData()
immune.combine <- IntegrateData(anchorset = immune.anchors)
#输出信息
Merging dataset 1 into 2
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
 
#数据整合后，查看数据集
 > immune.combine@assays 
$RNA
Assay data with 14053 features for 13999 cells
First 10 features:
 AL627309.1, RP11-206L10.2, LINC00115, NOC2L, KLHL17, PLEKHN1, HES4, ISG15, AGRN, C1orf159 

$integrated
Assay data with 2000 features for 13999 cells
Top 10 variable features:
 HBB, HBA2, HBA1, CCL4, CCL3, CCL7, TXN, GNLY, PPBP, APOBEC3B 
```

使用`DefaultAssay()`函数设定默认分析使用的数据集

```r
DefaultAssay(immune.combine) <- "integrated"
```

### 标准流程分析

数据整合后，即可进行seurt标准分析，如，`umap` `clust`....., 详见基本分析流程，简单举例umap

```
immune.combine <- ScaleData(immune.combine, verbose = T)
immune.combine <- RunPCA(immune.combine, npcs = 30, verbose = T)
immune.combine <- RunUMAP(immune.combine, reduction = "pca", dims = 1:30)
immune.combine <- FindNeighbors(immune.combine, reduction = "pca", dims = 1:30)
immune.combine <- FindClusters(immune.combine, resolution = 0.5)

p1 <- DimPlot(immune.combine, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combine, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210507163346844.png" alt="image-20210507163346844" style="zoom:50%;" />

```r
DimPlot(immune.combine, reduction = "umap", split.by = "stim")
```



<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210507163659030.png" alt="image-20210507163659030" style="zoom:50%;" />

### 保守基因筛选

通过**`FindConservedMarkers()`**函数确定每个cluster上的保守基因（其表达不受外界的刺激影响，本函数对两个数据集/组执行差异基因表达测试并给出p-value与FC（以为cluster6：NKcell 为例子）

```r
#为什么切换回RNA，可能是由于整合后的数据仅仅保留2000个高变基因，忽略了那些保守基因，因此要切换为未筛选的数据进行保守基因的筛选？存疑
#We can explore these marker genes for each cluster and use them to annotate our clusters as specific cell types.(官方)
DefaultAssay(immune.combined) <- "RNA"
#创建一个list用于存放每个cluster的输出结果
marker.list <- list()
#创建一个list用于存放每个cluster的top10的输出结果
top10 <- list()
cluster_end <- dim(table(immune.combine@meta.data$seurat_clusters)) - 1
for (cluster in c(0:cluster_end)){
  #指定cluster名称，用于解决list索引和cluster编号冲突问题（list索引1开始，cluster则从0开始）
  cluster_num <- str_c("cluster_",cluster)
  marker.list[[cluster_num]] <- FindConservedMarkers(object = immune.combine,ident.1 = cluster,grouping.var = "stim")
  top10[[cluster_num]] <- head(marker.list[[cluster_num]])
}
#从差异倍数来看，表达量差异不大，为保守基因
> head(top10[[cluster_6]])
       CTRL_p_val CTRL_avg_log2FC CTRL_pct.1 CTRL_pct.2 CTRL_p_val_adj    STIM_p_val STIM_avg_log2FC STIM_pct.1
GNLY            0        6.006173      0.944      0.045              0  0.000000e+00        5.856524      0.957
FGFBP2          0        3.243588      0.505      0.020              0 4.224915e-162        2.187640      0.259
CLIC3           0        3.461957      0.597      0.024              0  0.000000e+00        3.540011      0.623
PRF1            0        2.650548      0.422      0.017              0  0.000000e+00        4.100285      0.864
CTSW            0        2.987507      0.531      0.029              0  0.000000e+00        3.136218      0.596
KLRD1           0        2.777231      0.495      0.019              0  0.000000e+00        2.865059      0.552
       STIM_pct.2 STIM_p_val_adj      max_pval minimump_p_val
GNLY        0.059   0.000000e+00  0.000000e+00              0
FGFBP2      0.015  5.937273e-158 4.224915e-162              0
CLIC3       0.031   0.000000e+00  0.000000e+00              0
PRF1        0.057   0.000000e+00  0.000000e+00              0
CTSW        0.035   0.000000e+00  0.000000e+00              0
KLRD1       0.027   0.000000e+00  0.000000e+00              0
```

```R
#对每个cluster进行注释
#meta.data中每个cell类型都注释过
immune.combine <- RenameIdents(immune.combine, `0` = "CD14 Mono",
                                `1` = "CD4 Naive T", 
                                `2` = "CD4 Memory T", 
                                `3` = "CD16 Mono", 
                                `4` = "B", 
                                `5` = "CD8 T", 
                                `6` = "NK", 
                                `7` = "T activated", 
                                `8` = "DC", 
                                `9` = "B Activated", 
                                `10` = "Mk", 
                                `11` = "pDC", 
                                `12` = "Eryth", 
                                `13` = "Mono/Mk Doublets", 
                                `14` = "HSPC")
DimPlot(immune.combined, label = TRUE)
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210511144430105.png" alt="image-20210511144430105" style="zoom:67%;" />

使用**`Dotplot()`**, 查看保守基因在不同cell中的表达情况

```r
#将cell的注释结果转换为因子
Idents(immune.combine) <- factor(Idents(immune.combine), 
                                  levels = c("HSPC", 
                                             "Mono/Mk Doublets",
                                             "pDC", 
                                             "Eryth", 
                                             "Mk", 
                                             "DC", 
                                             "CD14 Mono", 
                                             "CD16 Mono", 
                                             "B Activated", 
                                             "B", 
                                             "CD8 T", 
                                             "NK", 
                                             "T activated", 
                                             "CD4 Naive T", 
                                             "CD4 Memory T"))
#每个cluster上选择2-3个marker基因
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", 
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", 
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
#可视化
DotPlot(object = immune.combine,features = markers.to.plot,cols = c("blue","red"),dot.scale = 8,split.by = "stim")+RotatedAxis()
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210511153037696.png" alt="image-20210511153037696" style="zoom:67%;" />

### 不同刺激条件下的显著变化的基因筛选

我们取受刺激和受控制的原始T细胞（CD4 Naive T）和CD14单核细胞（CD14 Mono）群的平均表达量，并生成散点图，突出显示对干扰素刺激有显著反应的基因。并在散点图上寻找视觉异常值的基因

```r
#筛选出CD4 Naive T的基因表达数据
t.cell <- subset(immune.combine,idents = "CD4 Naive T")
#查看原来的因子
Idents(t.cell)
......
 [ reached getOption("max.print") -- omitted 1509 entries ]
Levels: CD4 Naive T
#将因子转换为stim(是否刺激)
Idents(t.cell) <- "stim"
...
 [ reached getOption("max.print") -- omitted 1509 entries ]
Levels: CTRL STIM
#使用AverageExpression()计算基因的平均表达量(log1p即log10(x + 1))
avg.t.cells <- as.data.frame(log1p((AverageExpression(t.cell)$RNA)))
#增加一列基因信息
avg.t.cells$gene <- rownames(avg.t.cells)
#查看表达情况
> head(avg.t.cells)
                     CTRL        STIM          gene
AL627309.1    0.000000000 0.000000000    AL627309.1
RP11-206L10.2 0.007066129 0.000000000 RP11-206L10.2
LINC00115     0.015736246 0.000000000     LINC00115
NOC2L         0.566377860 0.621794899         NOC2L
.....
#CD14 Mono细胞的基因表达水平筛选与上述一致
cd14.cell <- subset(immune.combine,idents = "CD14 Mono")
Idents(cd14.cell) <- "stim"
avg.cd14.cell <- as.data.frame(log1p(AverageExpression(cd14.cell)$RNA))
avg.cd14.cell$gene <- rownames(avg.cd14.cell) 
#可视化展示与分析
#指定基因的名称，可以自行通过绘图数据划定阈值筛选，这里为了展示直接指定
genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p.t.cell <- ggplot(avg.t.cells,aes(CTRL,STIM)) + geom_point(color = "black" ,size = 1) + theme_bw()
p.t.cell <- LabelPoints(plot = p.t.cell, points = genes.to.label, repel = TRUE) + ggtitle("T Cells")

p.cd14.cell <- ggplot(avg.cd14.cell,aes(CTRL,STIM)) + geom_point(color = "black" ,size = 1) + theme_bw()
p.cd14.cell <- LabelPoints(plot = p.cd14.cell, points = genes.to.label, repel = TRUE) + ggtitle("CD14 Cells")
#p.t.cell + p.cd14.cell 
plot_grid(p.t.cell,p.cd14.cell)
```

![image-20210511160152415](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210511160152415.png)

查看上述筛选出的基因是否是高变基因, 使用**`FindMarkers()`**函数筛选CTRL和STIM之间的高变基因

```r
#首先构建分组，即每个cell在两种状态下的分组，并将其写入meta.data
immune.combine@meta.data$celltype.stim <- str_c(immune.combine@meta.data$seurat_annotations,"_",immune.combine@meta.data$stim)
#查看新的分组,celltype.stim即为构建的新分组
> head(immune.combine@meta.data,3)
                  orig.ident nCount_RNA nFeature_RNA stim seurat_annotations
AAACATACATTTCC.1 IMMUNE_CTRL       3017          877 CTRL          CD14 Mono
AAACATACCAGAAA.1 IMMUNE_CTRL       2481          713 CTRL          CD14 Mono
AAACATACCTCGCT.1 IMMUNE_CTRL       3420          850 CTRL          CD14 Mono

                 integrated_snn_res.0.5 seurat_clusters     celltype.stim
AAACATACATTTCC.1                      0               0    CD14 Mono_CTRL
AAACATACCAGAAA.1                      0               0    CD14 Mono_CTRL
AAACATACCTCGCT.1                      0               0    CD14 Mono_CTRL
......
#将分组变量更改为新构建的celltype.stim
Idents(immune.combine) <- "celltype.stim"
#确定指定cell的高变基因ident.1和ident.2指定cell
tcell.interferon.response <- FindMarkers(immune.combine,ident.1 = "CD4 Naive T_CTRL",ident.2 = "CD4 Naive T_STIM")
#查看marker基因，可以看出上述筛选出的基因均在以下列表
> head(b.interferon.response,20)
                p_val avg_log2FC pct.1 pct.2     p_val_adj
ISG15    0.000000e+00 -4.2335570 0.167 0.993  0.000000e+00
IFI6     0.000000e+00 -4.0328939 0.065 0.938  0.000000e+00
IFIT3    0.000000e+00 -4.1628278 0.006 0.893  0.000000e+00
ISG20    0.000000e+00 -2.6453200 0.459 0.976  0.000000e+00
IFIT1   6.396139e-307 -3.8751542 0.015 0.844 8.988494e-303
LY6E    2.776011e-289 -3.0165895 0.148 0.893 3.901129e-285
MX1     1.918617e-266 -3.3277341 0.069 0.812 2.696232e-262
B2M     3.947113e-228 -0.5654633 1.000 1.000 5.546878e-224
IFIT2   2.194799e-193 -3.1122608 0.011 0.622 3.084351e-189
OAS1    1.255347e-183 -2.8898711 0.018 0.609 1.764139e-179
IFI44L  2.667595e-176 -2.7416684 0.025 0.603 3.748771e-172
.....
```

使用**`FeaturePlot()`**展示不同基因在不同刺激条件下的表达情况

```r
#ISG15 IFI6 均为刺激条件下上调的基因；CD3D为保守基因
FeaturePlot(immune.combine, features = c("ISG15", "CD3D", "IFI6"), split.by = "stim", max.cutoff = 3, 
            cols = c("grey", "red"))
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210511163129381.png" alt="image-20210511163129381" style="zoom:67%;" />

使用**`VlNplot()`**进行展示

```r
plots <- VlnPlot(immune.combine,features = c("ISG15", "CD3D", "IFI6"),split.by = "stim",group.by = "celltype.stim",pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210511165213954.png" alt="image-20210511165213954" style="zoom:67%;" />