# Seurat-使用参考数据进行注释和mapping

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-7-21
R version 4.0.3 (2020-10-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)
```

## 构建参考数据对新数据进行整合

***In this vignette, we first build an integrated reference and then demonstrate how to leverage this reference to annotate new query datasets. Generating an integrated reference follows the same workflow described in more detail in the integration introduction [vignette](https://satijalab.org/seurat/articles/integration_introduction.html). Once generated, this reference can be used to analyze additional query datasets through tasks like cell type label transfer and projecting query cells onto reference UMAPs. Notably, this does not require correction of the underlying raw query data and can therefore be an efficient strategy if a high quality reference is available.***

## 数据载入

**本数据集是来自四种不同的测序方法的单细胞结果**（人胰岛细胞）

***CelSeq (GSE81076)***

***CelSeq2 (GSE85241),*** 

***Fluidigm C1 (GSE86469),*** 

***SMART-Seq2 (E-MTAB-5061)***

目的是使用***celseq, celseq2, smartseq2***构建参考数据集，来注释和整合***Fluidigm C1 (GSE86469),*** 

- 将多个不同的scna -seq数据集组装到一个参考数据集中
- 将细胞类型标签从参考数据集转移到新的查询数据集

```r
#载入包与数据
library(Seurat)
library(SeuratData)
#InstallData("panc8")
data("panc8")
#分割Seurat对象
pancreas.list <- SplitObject(panc8, split.by = "tech")
#取"celseq", "celseq2", "fluidigmc1", "smartseq2"作为后续分析对象
pancreas.list <- pancreas.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")]
>pancreas.list
$celseq
An object of class Seurat 
34363 features across 1004 samples within 1 assay 
Active assay: RNA (34363 features, 0 variable features)

$celseq2
An object of class Seurat 
34363 features across 2285 samples within 1 assay 
Active assay: RNA (34363 features, 0 variable features)

$fluidigmc1
An object of class Seurat 
34363 features across 638 samples within 1 assay 
Active assay: RNA (34363 features, 0 variable features)

$smartseq2
An object of class Seurat 
34363 features across 2394 samples within 1 assay 
Active assay: RNA (34363 features, 0 variable features)
```

## 数据标准化与高变基因的确定

```R
pancreas.list <- lapply(pancreas.list,function(x){
    x <- NormalizeData(x,verbose = FALSE )
    x <- FindVariableFeatures(x,selection.method = "vst",nfeatures = 2000,verbose = FALSE)
})
```

## 构建参考数据集

```r
#取出三个数据集构建参考数据集
reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]
#首先对这三个数据进行整合,构建参考数据集
pan.ant <- FindIntegrationAnchors(reference.list,dims = 1:30,verbose = FALSE)
pancreas.reference.Integra <- IntegrateData(anchorset = pan.ant, dims = 1:30,verbose = FALSE)
>pancreas.reference.Integra
An object of class Seurat 
36363 features across 5683 samples within 2 assays 
Active assay: integrated (2000 features, 2000 variable features)
 1 other assay present: RNA
```

## 参考数据集的基本分析

```r
DefaultAssay(pancreas.reference.Integra) <- "integrated"
#标准化
pancreas.reference.Integra <- ScaleData(pancreas.reference.Integra, verbose = FALSE)
#PCA
pancreas.reference.Integra <- RunPCA(pancreas.reference.Integra, npcs = 30, verbose = FALSE)
#UMAP
pancreas.reference.Integra.Integra <- RunUMAP(pancreas.reference.Integra, reduction = "pca", dims = 1:30, verbose = FALSE)
p1 <- DimPlot(pancreas.reference.Integra, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.reference.Integra, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
    NoLegend()
plot_grid(p1,p2,nrow = 2)
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210721120038593.png" alt="image-20210721120038593" style="zoom:90%;" />

## 使用整合后的参考数据集对其他数据集行cell注释

Seurat支持将参考数据集上的数据（或meta.data）投影到qurey数据集上，

这些方法可以对来自不同的个体、实验条件、测序技术甚至物种中收集来的数据进行整合，旨在识别出不同数据集之间的共享细胞状态(shared cell states),

**`data transfer` 和 `integration`有两个重要的区别**

- 前者不会对表达矩阵进行校正或调整***(Data transfer, Seurat does not correct or modify the query expression data.)***
- 前者进行PCA分析，后者为CCA***(data transfer, Seurat has an option (set by default) to project the PCA structure of a reference onto the query, instead of learning a joint structure with CCA. We generally suggest using this option when projecting data between scRNA-seq datasets)***

```R
#查看一下参考数据集的注释
pancreas.reference.Integra@meta.data %>% head()
>
A data.frame: 6 × 8
orig.ident	nCount_RNA	nFeature_RNA	tech	replicate	assigned_cluster	celltype	dataset
<chr>	<dbl>	<int>	<chr>	<chr>	<chr>	<chr>	<chr>
D101_5	D101	4615.810	1986	celseq	celseq	NA	gamma	celseq
D101_7	D101	29001.563	4209	celseq	celseq	NA	acinar	celseq
D101_10	D101	6707.857	2408	celseq	celseq	NA	alpha	celseq
D101_13	D101	8797.224	2964	celseq	celseq	NA	delta	celseq
D101_14	D101	5032.558	2264	celseq	celseq	NA	beta	celseq
D101_17	D101	13474.866	3982	celseq	celseq	NA	ductal	celseq
```
**使用`TransferData()`对分类数据进行转移（use the TransferData() function to classify the query cells based on reference data (a vector of reference cell type labels）,其返回值是一个矩阵，包括`pretict id`,`score`.....**

```r
#载入需要注释的数据集
pancreas.query <- pancreas.list[["fluidigmc1"]]
#寻找TransferAnchors
DefaultAssay(pancreas.reference.Integra) <- "integrated"
pancreas.anchors <- FindTransferAnchors(reference = pancreas.reference.Integra,
                                        query = pancreas.query,
                                        dims = 1:30,
                                        reference.reduction = "pca")
#将分类信息映射到query数据集
predict <- TransferData(anchorset = pancreas.anchors,
                        refdata = pancreas.reference.Integra$celltype,、
                        dims = 1:30)
#predict为预测的结果，每一列为得分，取得分结果最高的对应cell为该cell的cell类型
#将预测信息写入meta.data
pancreas.query <- AddMetaData(pancreas.query,metadata = predict)
#测试有多少cell是注释正确的 
c(pancreas.query@meta.data$predicted.id == pancreas.query@meta.data$celltype) %>% table()
#有617个cell是预测正确的，正确率在0.96左右
>
FALSE  TRUE 
   21   617 
```

## 查看指定基因在cell各个cell中的表达

```r
#预测结果
VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id") 
```

![image-20210721135646694](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210721135646694.png)

```r
#实际结果
VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "celltype")
```

![image-20210721135734809](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210721135734809.png)

## 将UMAP的结果映射到query数据集

参考数据集做`UMAP`分析, 后使用**`MapQuery()`**函数将UMAP的结果进行映射

***What is MapQuery doing?***
**`MapQuery()` is a wrapper around three functions: `TransferData()`, `IntegrateEmbeddings()`, and `ProjectUMAP()`. `TransferData()` is used to transfer cell type labels and impute the ADT values; `IntegrateEmbeddings() `is used to integrate reference with query by correcting the query’s projected low-dimensional embeddings; and finally `ProjectUMAP()`is used to project the query data onto the UMAP structure of the reference. The equivalent code for doing this with the intermediate functions is below:**

```r
pancreas.reference.Integra <- RunUMAP(pancreas.reference.Integra, di
                                      s = 1:30, 
                                      reduction = "pca", 
                                      return.model = TRUE,
                                      verbose = F)
#将query数据集映射到参考数据集
#refdata: 需要被transfer的数据
pancreas.query <- MapQuery(anchorset = pancreas.anchors, 
                           reference = pancreas.reference.Integra, 
                           query = pancreas.query,
                           refdata = list(celltype = "celltype"), 
                           reference.reduction = "pca", 
                           reduction.model = "umap")
```

查看映射结果

```r
p6 <- DimPlot(pancreas.reference.Integra, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p7 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
plot_grid(p6,p7,nrow = 2)
```

![image-20210721140034597](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210721140034597.png)

