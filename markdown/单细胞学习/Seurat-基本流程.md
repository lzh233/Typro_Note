# Seurat-基本流程

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-3-18
sessionInfo("Seurat")
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)
version: Seurat_4.0.1
```

```R
#查看输入文件列表
setwd("D:\\Desktop\\s_note\\data\\singel_cell\\seurat")
dir()
rm(list = ls())
> dir()
[1] "barcodes.tsv" "genes.tsv"    "matrix.mtx"  
```

## 单细胞数据分析的基本流程

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715115940415.png" alt="image-20210715115940415" style="zoom:70%;" />

## 数据读入

**载入相关包以及设置工作目录,使用官网的提供的数据集，10x下机处理好的数据**

```R
library("Seurat")
library("tidyverse")
library("patchwork")
gc()
```

直接读入10x数据，并创建Seurat对象，同时对基因和cell进行初步的QC

```R
#读入10X数据
pbmc.data <- Read10X(data.dir = "D:\\Desktop\\s_note\\data\\singel_cell\\seurat")
#创建Seurat对象，初步QC
#min.cell：所有cell中检测到数量小于3的基因剔除
#min.features：某个cell中检测到的基因数量小于200的细胞剔除
pbmc <- CreateSeuratObject(counts = pbmc.data,project = "pbmc3k",min.cells = 3,min.features = 200)
#查看导入后的数据是否正确
pbmc

#（输出结果）查看数据情况，cell数和基因数目
> pbmc
An object of class Seurat 
13714 features across 2700 samples within 1 assay 
Active assay: RNA (13714 features, 0 variable features)
```

**表达矩阵的直接读入(注意列为cell，行为基因)**

```R
setwd("setwd("D:\\Desktop\\s_note\\data\\singel_cell\\")
count.data <- read.delim("GSE111229_Mammary_Tumor_fibroblasts_768samples_rpkmNormalized.txt",header = T,row.names = 1,sep = "\t")
count.row <- CreateSeuratObject(counts = count.data,project = "cell_try",min.cells = 0,min.features = 0)
```

## **数据质控**

**质控原则：**

1. 标准数据QC过程，3个标准
   QC标准
   根据每个cell中发现的基因数量，进行过滤
   			a. 低质量的cell或空液滴通常有很低的gene数
   			b. doublets or multiplets cell 有很高的基因表达
2. 每个cell中检测到的基因总数
3.  线粒体,或其他指定需要过滤的基因的比例
      			a. 低质量与濒临死亡的cell通常有很高的线粒体含量
      			b. 使用内置函数**`PercentageFeatureSet()`**进行线cell线粒体gene比例计算，可以计算线粒体比例，本例以MT开头的gene为线粒体基因

```R
#计算各个cell中线粒体的含量
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc,pattern = "^MT-")
#还会计算红细胞基因的含量，具体质控方法间多样本整合数据分析-1
#计算结果查看
#计算储存在CreateSeuratObject对象的meta.data中，meta.data中已有的featurecount和rnacount是生成对象时计算的
head(pbmc@meta.data,6)

（输出结果）> head(pbmc@meta.data,6)
                 orig.ident nCount_RNA nFeature_RNA percent.mt
AAACATACAACCAC-1     pbmc3k       2419          779  3.0177759
AAACATTGAGCTAC-1     pbmc3k       4903         1352  3.7935958
AAACATTGATCAGC-1     pbmc3k       3147         1129  0.8897363

#补充：对其他类型的基因进行过滤（from周运来就是我，简书）

```

使用**`VlnPlot()`**函数绘制小提琴图，查看`nCount_RNA`、`nFeature_RNA`、`percent.mt`在每个cell中的数量分布情况

```R
#可视化统计结,可自行根据meta.data指定展示的对象
VlnPlot(pbmc,features = c("nCount_RNA","nFeature_RNA","percent.mt"))
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715094222158.png" alt="image-20210715094222158" style="zoom:80%;" />

使用**`FeatureScatter()`**查看各个feature之间的相关性，可以得出`nCount_RNA`与线粒体含量的关系，确定一个质控的范围

```R
#查看count和feature的关系
plot1 <- FeatureScatter(pbmc,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
#查看count和线粒体含量的关系
plot2 <- FeatureScatter(pbmc,feature1 = "nCount_RNA",feature2 = "percent.mt")
#可视化
plot1 + plot2
```

![image-20210715095343088](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715095343088.png)

使用**`subset()`**函数进行数据质控，可以比较一下质控前后数据的差异

```R
#根据上述结果删除了线粒体含量较高的cell(>5%)同时nfeatures的数量在200和2500之间
pbmc <- subset(pbmc,subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc

#质控结果比较
#质控后（64个cell被剔除）
> pbmc
An object of class Seurat 
13714 features across 2638 samples within 1 assay 
Active assay: RNA (13714 features, 0 variable features)
#质控前
> pbmc
An object of class Seurat 
13714 features across 2700 samples within 1 assay 
Active assay: RNA (13714 features, 0 variable features)
```

## 数据归一化

使用**`NormalizeData()`**函数归一化数据，Seurat提供了3种归一化方法，分别为：

1. `LogNormalize`:每个基因的count数除以该基因在所有cell中的count总和，在乘以缩放因子（默认10000）后取`log(x+1)`（本例所使用的方法）
2. `CLR`:Applies a centered log ratio transformation ??不明白
3. `RC`:即CPM值，count数除以该基因在所有cell中的count总和后乘以1000000

```R
pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",scale.factor = 10000)
```

## 高变基因的筛选

使用**`FindVariableFeatures()`**函数筛选高变基因，**`VariableFeatures()`**查看高变基因，Seurat提供3种方法筛选高变基因

1. `vst`：首先根据基因在所有细胞表达量的`log(variance) `和 `log(mean)`使用loess进行拟合，然后对每个基因的表达量进行标准化，即（基因表达量-平均表达量）/回归分析所得的每个基因的方差，后得到`平均表达量`和`标准方差`的散点图，通过划定特定的方差值来得到高变基因`clip.max`选项可以设定阈值，同时`loess.span`选项可以对拟合情况进行修改，`nfeatures`选项代表选择的基因数目。
2. `mean.var.plot (mvp)`：首先利用mean.function和 dispersion.function分别计算每个基因的平均表达量和离散情况；然后根据平均表达量将基因们分散到一定数量（默认是20个）的小区间（bin）中，并且计算每个bin中z-score。？
3. `dispersion (disp)`：选取离散值最高的基因

```R
#确认高表达的前2000个features
?FindVariableFeatures()
#vst使用均值和方差确定，具体看帮助文档
pbmc <- FindVariableFeatures(pbmc,selection.method = "vst",nfeatures = 2000)
#使用VariableFeatures创建对象，并获取差异基因前10
top10 <- head(VariableFeatures(pbmc),10)
plot3 <- VariableFeaturePlot(pbmc)
plot4 <- LabelPoints(plot = plot3,points = top10,repel = T)
plot3 + plot4
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715100355234.png" alt="image-20210715100355234" style="zoom:80%;" />

## 数据标准化

使用**`ScaleData()`**函数对数据进行标准化，默认使用只对前一步骤确认的前2000个高变基因进行标准化，后续分析同样一致，默认对数据同时进行中心化与z-score标准化`do.center=T`和`do.scale=T`，`features`指定对那些基因进行标准化，比较对全部基因标准化与前两个基因标准化对后续分析影响，使用**`vars.to.regress`**选项可以校正线粒体基因比例的影响，如

```R
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
```

- **全部数据标准化**

  ```R
  all.gene <- rownames(pbmc)
  pbmc <- ScaleData(pbmc,features = all.gene)
  #查看标准化后的数据
  > head(pbmc[["RNA"]]@scale.data[,1:3],3)
                AAACATACAACCAC-1 AAACATTGAGCTAC-1 AAACATTGATCAGC-1
  AL627309.1         -0.05812316      -0.05812316      -0.05812316
  AP006222.2         -0.03357571      -0.03357571      -0.03357571
  RP11-206L10.2      -0.04166819      -0.04166819      -0.04166819
  ```

- **前2000个高变基因标准化**

  ```R
  pbmc <- ScaleData(pbmc,features = VariableFeatures(pbmc))
  ```

## PCA分析

- **标准化后的全部数据**

使用**`RunPCA()`**函数进行PCA分析，数据需要提前进行标准化，首先使用全部基因进行PCAF分析

```R
pbmc <- RunPCA(pbmc,features = all.gene)
#每个cell的score值储存在以下位置
> head(pbmc@reductions$pca@cell.embeddings[,1:3],4)
                      PC_1       PC_2       PC_3
AAACATACAACCAC-1 -8.064775 -0.5646571 -0.5684279
AAACATTGAGCTAC-1 -0.118376  1.3762525 -3.6784549
AAACATTGATCAGC-1 -5.126320 -8.8074650 -6.7152749
AAACCGTGCTTCCG-1 18.846209  2.3536153 -2.2374983
#每个基因的score值（loading值储存在以下位置）
head(pbmc@reductions$pca@feature.loadings[,1:3],4)
                      PC_1          PC_2          PC_3
AL627309.1    0.0028691834  0.0008739089 -0.0017311397
AP006222.2    0.0003312499 -0.0008990262 -0.0002707583
RP11-206L10.2 0.0002819107  0.0006482639 -0.0023268238
RP11-206L10.9 0.0020072933 -0.0002714407 -0.0016980831
```

使用**`DimPlot()`**对PCA结果进行可视化，**`VizDimLoadings()`**和**`DimHeatmap()`**对不同基因在不同主成分轴上的loading进行可视化

```R
DimPlot(pbmc,reduction = "pca")
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715101341597.png" alt="image-20210715101341597" style="zoom:80%;" />

```R
#不同主成分轴上的loading，查看前2个主成分轴上的前22个features的loading值
VizDimLoadings(object = pbmc,reduction = "pca",nfeatures = 22,dims = 1:2)
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715101846227.png" alt="image-20210715101846227" style="zoom:80%;" />

```R
#前500个cell中的基因在PC1上的loading分布---探索异质性来源
#dim指定主成分轴数,cell指定细胞数目
DimHeatmap(pbmc,dims = 1:2,nfeatures = 22,reduction = "pca",cells = 500)
```

![image-20210715102119508](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715102119508.png)

- **使用前2000个基因PCA分析,结果基因一致**，***后续分析使用前<u>2000个高变基因</u>继续分析***

```R
pbmc <- RunPCA(pbmc)
DimPlot(pbmc,reduction = "pca")
```

![image-20210715102333512](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715102333512.png)

## 主成分轴数的确定

```R
#我们对矩阵进行了细胞水平质控，但是如何对基因进行质控，当然我们创建seurat 对象时候，已经设定一个cutoff. 基因必须早多少个细胞中表达. 但是如果直接用这些基因进行聚类和降维（T-sne/UMAP 可视化 可能引起偏差，所以需要先进行一下PCA 初步选取主要的主成分，对应着权重不同的基因，再用这些基因进行聚类，降维可视化效果更好。
```

1. 使用了**`JackStraw()`**函数确定主成分的轴数，内置算法为jackstraw法，JackStraw法，是随机地从所有的feature中选出1%，去计算PCA值，再将这部分基因的PCA值与实际观察到的PCA值进行比较，从而得到统计显著性，最终的结果是将每个基因与不同主成分的相关性的显著性水平计算出来。那么一个PC越显著，对结果的影响越大，那么它就应当更多的富集那些显著的feature（基因）。更多的原理参考文献[Macosko *et al*](http://www.cell.com/abstract/S0092-8674(15)00549-8)。`num.replicate`选项指定重复次数，使用**`ScoreJackStraw()`**函数进行分析，`dim`选项指定使用前几轴

```R
pbmc <- JackStraw(pbmc,num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc,dims = 1:20)
```

2. 使用**`JackStrawPlot()`**进行可视化操作，通过p-value选择合适的主成分轴，PC的显著性越强，则包含的显著性高的feature就越多，因此如果一个PC显著，则它在这张图上将显著地左偏。因此我们认为前10-12个PC较为显著，可以作为我们后续聚类的依据。

```R
JackStrawPlot(pbmc,dims = 1:20)
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715103735478.png" alt="image-20210715103735478" style="zoom:80%;" />

3. 使用**`ElbowPlot()`**函数可视化每个排序轴贡献度的分布，寻找贡献度开始趋于平稳的轴，**碎石图？**

```R
ElbowPlot(pbmc)
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715103829627.png" alt="image-20210715103829627" style="zoom:80%;" />

**关于PCA轴的选取：**

```r
#一是最显著的PC不一定能帮助我们区分一些很小的类，这些类由于本身数量较少，因此很容易淹没在噪声中；二是如果对取多少PC或哪些PC没有把握，则尽量往多了取，噪声虽然会更多，但总比信息丢失要好。
```

**本例子中结合上述两个分析，将前13个主成分作为后续分析选取的主成分**

## 聚类分析

聚类分析一般分为层次聚类(hclust)、k-mean划分、dbscan等等，Seurat提供的聚类算法为`graph-based`算法，结果保存在`object@ident`，Seurat的聚类分析步骤如下，

```R
step1：在PhenoGraph中，首先根据PCA的欧氏距离结果，构建KNN图，找到了各个近邻组成的小社区；然后根据近邻中两两住户（细胞）之间的交往次数（shared overlap）计算每条线的权重（术语叫Jaccard similarity） 。这些计算都是用包装好的函数FindNeighbors() 得到的，它的输入就是前面降维最终确定的主成分。

Step-2：得到权重后，为了对细胞进行聚类，使用了计算模块的算法（默认使用Louvain，另外还有包括SLM的其他三种），使用FindClusters() 进行聚类。其中包含了一个resolution的选项，它会设置一个”间隔“值，该值越大，间隔越大，得到的cluster数量越多。一般来说，这个值在细胞数量为3000左右时设为0.4-1.2 会有比较好的结果
链接：https://www.jianshu.com/p/4a2a2acfc10b
####关于knn的一个解释###
KNN算法，即K近邻算法是一种监督学习算法，本质上是要在给定的训练样本中找到与某一个测试样本A最近的K个实例，然后统计k个实例中所属类别计数最多的那个类，就是A的类别。
```

使用**`FindNeighbors()`**得到cell之间的关系，后使用**`FindClusters()`**函数进行聚类，`resolution`选项在细胞数量为**3000左右时设为0.4-1.2** 会有比较好的结果，**值越大得到的cluster就越多**

```R
pbmc <- FindNeighbors(pbmc,dims = 1:13)
pbmc <- FindClusters(pbmc,resolution = 0.5)

#查看聚类的情况
> head(Idents(pbmc),3)
AAACATACAACCAC-1 AAACATTGAGCTAC-1 AAACATTGATCAGC-1 
               4                3                2 
Levels: 0 1 2 3 4 5 6 7 8

#对聚类情况进行统计，看一下每个cluster中的cell数目
> table(Idents(pbmc))
  0   1   2   3   4   5   6   7   8 
696 477 446 343 307 165 158  32  14 
```

## 降维分析：t-sne

使用**`RunTSNE()`**函数可以进行t-sne分析，`dim`选项指定使用**PCA**分析的主成分数量，`feature`选项可以对指定数据降维，`reduction`可以指定降维方法输出的结果用于t-sne分析。

```R
pbmc <- RunTSNE(pbmc,dims = 1:13)
#label是否显示标签，pt.size点的大小
DimPlot(pbmc,reduction = "tsne",pt.size = 1,label = T)
#查看t-sne结果
> head(pbmc@reductions$tsne@cell.embeddings)
                       tSNE_1     tSNE_2
AAACATACAACCAC-1  -9.78163034  12.794337
AAACATTGAGCTAC-1  -6.99146947 -29.740188
AAACATTGATCAGC-1  -0.02430397  20.579422
AAACCGTGCTTCCG-1  27.96447311  -5.042012
AAACCGTGTATGCG-1 -27.70559808  19.731951
AAACGCACTGGTAC-1   1.03852106   9.407889
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715104552348.png" alt="image-20210715104552348" style="zoom:80%;" />

## 降维分析：UMAP

使用**`RunUMAP()`**函数进行UMAP降维，`dim`选项指定使用**PCA**分析的主成分数量

```R
pbmc <- RunUMAP(pbmc,dims = 1:13)
DimPlot(pbmc,reduction = "umap",pt.size = 1,label = T)

#查看UMAP结果
> head(pbmc@reductions$umap@cell.embeddings)
                    UMAP_1       UMAP_2
AAACATACAACCAC-1 -3.348948   2.84308863
AAACATTGAGCTAC-1 -4.933196 -11.56426882
AAACATTGATCAGC-1 -3.554357   5.64948679
AAACCGTGCTTCCG-1  9.946793   0.15786392
AAACCGTGTATGCG-1 -8.150798   0.04117789
AAACGCACTGGTAC-1 -2.074983   4.63934446
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715104828906.png" alt="image-20210715104828906" style="zoom:80%;" />

## marker基因的确定

使用**`FindMarkers()`**函数寻找Marker基因，

`ident.1`选项来选则要确定的Marker基因所在的轴，

`min.pct`选项表示基因的含量要大于指定的百分数，

`test.use`选项可以指定不同的差异基因确定的算法，

`logfc.threshold`选项确认foldchange的阈值。

`test.use`提供的方法分别为，`wilcox (默认)` 、`roc`、`bimod`、`t`、`negbinom`、`poisson`、`LR`、`MAST`、`DESeq2`,**可通过不同算法来确定差异基因**

不同方法的原理如下（待补充）

```R
?FindAllMarkers()
'wilcox：即Wilcoxon符号秩检验
roc：使用ROC曲线分析差异基因，得出的AUC值是该基因对分类模型的重要度，AUC=1或0时，说明该基因可以完美的区分两个cluster（分类器的AUC等于分类器随机选择的正实例的概率高于随机选择的负实例的概率），AUC=0.5时说明该基因对分类没有影响，使用‘predictive power即， (abs(AUC-0.5) * 2)来预测分类能力，越大说明该基因对模型的影响越大
t:使用t-test
bimod:	Likelihood-ratio test for single cell gene expression, (McDavid et al., Bioinformatics, 2013)'
```

```R
#使用wilcox法筛选cluster2的差异基因
clus1.marker <- FindMarkers(pbmc,ident.1 = 2,min.pct = 0.25,test.use = "wilcox")
#查看差异基因
#pct.1 为指定的cluster上的占比，pct.2为其余所有cluster上的占比
> head(clus1.marker,5)
            p_val avg_log2FC pct.1 pct.2    p_val_adj
LTB  1.068818e-82  1.2627653 0.982 0.647 1.465777e-78
IL32 8.430751e-81  1.0969979 0.946 0.472 1.156193e-76
LDHB 2.567403e-66  0.9531558 0.964 0.616 3.520936e-62
IL7R 2.097231e-65  1.1923612 0.756 0.331 2.876143e-61
CD3D 5.747103e-64  0.8837803 0.913 0.440 7.881577e-60

#使用roc法筛选cluster2的差异基因
clus1.marker <- FindMarkers(pbmc,ident.1 = 2,min.pct = 0.25,test.use = "roc",logfc.threshold = 0.25)
#查看差异基因， myAUC：
> head(clus2.marker,5)
     myAUC  avg_diff power avg_log2FC pct.1 pct.2
LTB  0.785 0.8752822 0.570  1.2627653 0.982 0.647
IL32 0.772 0.7603810 0.544  1.0969979 0.946 0.472
LDHB 0.754 0.6606772 0.508  0.9531558 0.964 0.616
CD3D 0.739 0.6125898 0.478  0.8837803 0.913 0.440
IL7R 0.727 0.8264818 0.454  1.1923612 0.756 0.331

#确定指定cluster间的差异基因
clus5.marker <- FindMarkers(pbmc,ident.1 = 5,ident.2 = c(0:2),min.pct = 0.25,test.use = "wilcox")
#查看差异基因，pct.1 为指定的cluster上的占比，pct.2为其余所有cluster上的占比
> head(clus5.marker,5)
               p_val avg_log2FC pct.1 pct.2     p_val_adj
FCGR3A 7.807758e-233   4.075321 0.976 0.069 1.070756e-228
CDKN1C 6.708315e-151   1.595292 0.509 0.011 9.199783e-147
RHOC   2.126300e-150   2.598861 0.867 0.101 2.916008e-146
HES4   4.075406e-144   1.725508 0.582 0.025 5.589011e-140
MS4A7  4.575958e-138   2.596997 0.800 0.090 6.275469e-134
```

使用**`FindAllMarkers()`**函数可以直接比较每个cluster上的差异基因,`only.pos=T`选项指定仅筛选上调的基因

```R
#确定所有cluster上的差异基因
cluster.marker.all <- FindAllMarkers(pbmc,logfc.threshold = 0.25,min.pct = 0.25,test.use = "wilcox",only.pos = T)
#根据log2foc值筛选出每个cluster上的top10差异基因
cluster.top <- 
	cluster.marker.all %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)

#查看所有差异基因
> head(clust.top)
# A tibble: 6 x 7
# Groups:   cluster [3]
      p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene  
      <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr> 
1 2.01e-116       1.06 0.914 0.59  2.75e-112 0       LDHB  
2 1.99e- 84       1.36 0.438 0.11  2.72e- 80 0       CCR7  
3 0.              5.58 0.996 0.217 0.        1       S100A9
4 0.              5.49 0.975 0.123 0.        1       S100A8
5 1.07e- 82       1.26 0.982 0.647 1.47e- 78 2       LTB   
6 3.39e- 55       1.28 0.422 0.115 4.64e- 51 2       AQP3  
```

## Marker基因表达水平的可视化

使用**`FeaturePlot()`**和**`VlnPlot()`**来可视化基因的表达情况

```R
#featuresz指定基因名称，一个或多个均可，默认降维方式为UMAP，通过reduction可以指定降维结果
FeaturePlot(pbmc,features = "CD3E")
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715111114358.png" alt="image-20210715111114358" style="zoom:80%;" />

```R
#指定展示多个基因
FeaturePlot(pbmc,features = c("CD3E","PF4","LDHB"))
```

![image-20210715111247853](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715111247853.png)

```R
#指定降维方式为t-sne
FeaturePlot(pbmc,features = "CD3E",reduction = "tsne",label = T)
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715111438810.png" alt="image-20210715111438810" style="zoom:67%;" />

```R
#查看某个基金在各个簇中的表达情况
VlnPlot(pbmc,features = "PF4")
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715111521686.png" alt="image-20210715111521686" style="zoom:67%;" />

```R
#同时展示多个基因
VlnPlot(pbmc,features = c("CD3E","PF4","LDHB","CCR7"))
```

![image-20210715111627396](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715111627396.png)

```R
#使用count值展示,slot指定使用的数据，log是否执行log转换
VlnPlot(pbmc,features = "LDHB",slot = "count",log = T)
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715111737550.png" alt="image-20210715111737550" style="zoom:67%;" />

## cell类型注释（具体见数据整合分析-2）

根据Marker基因确定了各个cluster的cell类型，**`RenameIdents()`**函数重新对cluster赋值

```R
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono","NK", "DC", "Platelet")
#将cell名赋值个各个cluster
new.cluster.ids <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.cluster.ids)

> levels(pbmc)
[1] "Naive CD4 T"  "CD14+ Mono"   "Memory CD4 T" "B"            "CD8 T"        "FCGR3A+ Mono" "NK"           "DC"          
[9] "Platelet"    
```

```R
#可视化重新注释后的umap
DimPlot(pbmc,reduction = "umap",label = T) + NoLegend()
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715112004445.png" alt="image-20210715112004445" style="zoom:67%;" />

```R
#其他的也会更改
VlnPlot(pbmc,features = c("CD3E","PF4","LDHB","CCR7"))
```

![image-20210715112023409](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210715112023409.png)