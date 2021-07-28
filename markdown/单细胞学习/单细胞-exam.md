# Exam

### Q1 How is sequencing saturation calculated in CellRanger? (10)

```
Sequencing Saturation = 1 - (n_deduped_reads / n_reads)
n_deduped_reads = Number of unique (valid cell-barcode, valid UMI, gene)
n_reads = Total number of confidently mapped, valid cell-barcode, valid UMI reads.
```

### Q2 Explain the causes of high level of mitochondrial gene expression in scRNA-Seq data. (10)

```R
#细胞的质量较低，有较多的细胞发生了凋亡或溶解，也有可能是细胞本身原因所造成的，如肿瘤细胞的代谢和其坏死都会增加线粒体的比例，也有可能是细胞样品被死掉或溶解的细胞污染，导致线粒体的含量增加
```

### Q3 How to modify the STAR alignment parameters in CellRanger? How to modify them in CeleScope? (10)

```
通过修改cellranger的原代码对STAR的参数进行修改
```

### Q4 What is the difference between Seurat data integration and merging? When will you use `IntegrateData` instead of `merge`? (10)

```r
What is the difference between Seurat data integration and merging?
#merge为r中的基础函数，只是将两个数据集简单的合并到一起，而合并后的数据可能会存在批次效应（通过PCA tsne umap...均可看出），影响后续的分析
#Seurat包中的IntegrateData是利用CCA分析将两个数据集降维到同一个低维空间，进而找到两个数据集之间互相“距离”最近的细胞，Seurat将这些相互最近邻细胞称为“锚点细胞”，一般来说只用细胞类型一致和状态一致的细胞才会成为锚点细胞，然后Seurat会利用这些锚点细胞实现单细胞数据的整合
When will you use `IntegrateData` instead of `merge`?
#当需要对不同的单细胞数据进行整合时，如来自不同的实验，数据间可能存在批次效应时应该使用IntegrateData进行数据的整合分析
```

### Q5 Install SingleR through Conda and use SingleR to annotate [pbmc3k](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) data. (10)

```r
#安装SingleR和注释数据库celldex
BiocManager::install("SingleR")
BiocManager::install("celldex")
#下载数据库
MonacoImmuneData()
HumanPrimaryCellAtlasData() 
BlueprintEncodeData() 
MouseRNAseqData() 
ImmGenData() 
DatabaseImmuneCellExpressionData()

#细胞类群注释
library("SingleR")
library("celldex")
library("Seurat")
library("tidyverse")
library("patchwork")
#读入数据，质控降维聚类......
setwd("D:\\Desktop\\s_note\\data\\singel_cell\\seurat")
pbmc.data <- Read10X(data.dir = "D:\\Desktop\\s_note\\data\\singel_cell\\seurat")
pbmc <- CreateSeuratObject(counts = pbmc.data,project = "pbmc.seurt",min.cells = 3,min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc,pattern = "^MT-")
pbmc <- subset(pbmc,subset = nFeature_RNA > 200 & nFeature_RNA <2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc,selection.method = "vst",nfeatures = 2000)
pbmc <- ScaleData(pbmc,vars.to.regress = "percent.mt")
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc,dims = 1:13)
pbmc <- FindClusters(pbmc,resolution = 0.5)
#细胞注释
#加载参考数据库
ref.data <- MonacoImmuneData()
#载入需要被注释的数据
note.pbmc <- GetAssayData(pbmc, slot="data")
#提取clusters信息
clusters <- pbmc@meta.data$seurat_clusters
#细胞类型注释
celltype <- SingleR(test = note.pbmc,
                    ref = ref.data,
                    labels = ref.data$label.main,
                    clusters = clusters,
                    assay.type.test = "logcounts",
                    assay.type.ref = "logcounts")
type_note <- data.frame(cluster = rownames(celltype),cell_type=celltype$labels)
#细胞注释
new.cluster.ids <- type_note[,2]
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.cluster.ids)
#注释信息写入meta.data
pbmc$cell.type <- Idents(pbmc)
####查看注释结果####
> pbmc@meta.data %>% head()
                 orig.ident nCount_RNA nFeature_RNA percent.mt RNA_snn_res.0.5 seurat_clusters          cell_type
AAACATACAACCAC-1 pbmc.seurt       2419          779  3.0177759               0               0 CDT cells+ T cells
AAACATTGAGCTAC-1 pbmc.seurt       4903         1352  3.7935958               3               3            B cells
AAACATTGATCAGC-1 pbmc.seurt       3147         1129  0.8897363               2               2 CDT cells+ T cells
AAACCGTGCTTCCG-1 pbmc.seurt       2639          960  1.7430845               1               1          Monocytes
AAACCGTGTATGCG-1 pbmc.seurt        980          521  1.2244898               6               6           NK cells
AAACGCACTGGTAC-1 pbmc.seurt       2163          781  1.6643551               2               2 CDT cells+ T cells
> table(Idents(pbmc))
   CD4+ T cells       Monocytes         B cells         T cells        NK cells Dendritic cells 
           1140             657             343             311             156              31 
```

### Q6 Using pbmc3k data, write a script to:

1. Plot the average UMI per gene per cell on t-SNE plot. (10)
2. Generate an expression heatmap for the top 10 markers(10 each, total 30 genes) of Naive CD4+ T, CD8+ T and NK cells. (10)

```r
#1、Plot the average UMI per gene per cell on t-SNE plot. (10)
#细胞类群注释
library("Seurat")
library("tidyverse")
library("patchwork")
#读入数据，质控降维聚类......
setwd("D:\\Desktop\\s_note\\data\\singel_cell\\seurat")
pbmc.data <- Read10X(data.dir = "D:\\Desktop\\s_note\\data\\singel_cell\\seurat")
pbmc <- CreateSeuratObject(counts = pbmc.data,project = "pbmc.seurt",min.cells = 3,min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc,pattern = "^MT-")
pbmc <- subset(pbmc,subset = nFeature_RNA > 200 & nFeature_RNA <2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc,selection.method = "vst",nfeatures = 2000)
all.gene <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.gene)
pbmc <- RunPCA(pbmc)
pbmc <- RunTSNE(pbmc,dims = 1:10)
#计算Per umi，并写入meta.data
pbmc$GenesPerUMI <- pbmc$nFeature_RNA / pbmc$nCount_RNA
#提取tsne作图数据并作图
tsne_plot <- pbmc@reductions$tsne@cell.embeddings %>% 
  as.data.frame() %>% 
  mutate(GenesPerUMI=pbmc$GenesPerUMI) %>% 
  ggplot(tsne_plot,mapping = aes( tSNE_1 ,tSNE_2,color = GenesPerUMI)) + 
  geom_point()
tsne_plot
```

<img src="D:\Desktop\s_note\data\picture\image-20210716185319157.png" alt="image-20210716185319157" style="zoom:80%;" />

```r
#2、Generate an expression heatmap for the top 10 markers(10 each, total 30 genes) of Naive CD4+ T, CD8+ T and NK cells. (10)
pbmc <- FindNeighbors(pbmc,dims = 1:10)
pbmc <- FindClusters(pbmc,resolution = 0.5)
#注释cell
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono","NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.cluster.ids)
#注释信息写入meta.data
pbmc$cell.type <- Idents(pbmc)
#寻找marker基因
all.markers <- FindAllMarkers(pbmc,logfc.threshold = 0.25,min.pct = 0.25,test.use = "wilcox",only.pos = T)
#提取指定细胞的top10的基因
sub.markers <- 
  all.markers %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC) %>% filter(cluster %in% c("Naive CD4 T","CD8 T","NK"))
#得到指定细胞的top10的基因列表
top.gene.sun <- sub.markers$gene
#细胞信息的提取
top.sub <- pbmc@meta.data %>% filter(cell.type %in% c("Naive CD4 T","CD8 T","NK")) %>% select(cell.type)
#提取指定细胞和基因的表达矩阵
scale.data <- pbmc@assays$RNA@scale.data[top.gene.sun,rownames(top.sub)]
#绘制热图
colors <- colorRampPalette(c("purple","yellow"))(100)
pheatmap::pheatmap(scale.data,annotation_col = top.sub,cluster_cols = T,color = colors,cutree_cols = 3,show_colnames = F)
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210716222918331.png" alt="image-20210716222918331" style="zoom:80%;" />

### Q7 10X V2 chemistry has 10-bp long UMI but 12-bp in V3. What are the advantages and disadvantages of longer UMI? (10)

```r
更长的umi扩增时候可能会造成一定的噪音污染，同时会检测到更多的基因数量
```

### Q8 Write a script to subset 100 cells with the highest UMI count from pbmc3k data and save the matrix of these 100 cells in 10X format. (20)

```r
library(Matrix)
library(tidyverse)

#barcode为10x数据的barcode.tsv
#gene为10x数据的genes.tsv
#expr_mat为10x的matrix.mtx
#n_top: 提取前几歌
#sort_by：根据哪个指标提取nCount_RNA或nFeature_RNA
extract_top_cell <- function(barcode,gene,expr_mat,n_top=100,sort_by="nCount_RNA"){
  #对表达矩阵进行命名，列名为细胞名，行名为基因名称
  print("rename the expression matrix......")
  colnames(expr_mat) <- barcode[,1]
  rownames(expr_mat) <- gene[,1]
  
  #计算指标nCount_RNA
  print("caculate nCount_RNA......")
  nCount_RNA <- apply(expr_mat,2,function(cell){sum(cell)})
  
  #计算指标nFeature_RNA
  print("caculate nFeature_RNA......")
  nFeature_RNA <- apply(expr_mat,2,function(cell){sum(cell> 0) })
  
  #构建每个cell的count和feature信息
  print("create matrix......")
  meta.data <- data.frame(nCount_RNA = nCount_RNA,nFeature_RNA=nFeature_RNA)
  
  #排序meta.data并提取相应的cell
  print(str_c("extract cell top of ",n_top," by ",sort_by))
  meta.data.top <- meta.data %>% 
    arrange(desc(.[,sort_by])) %>% 
    head(n_top)
  
  #提取表达矩阵数据
  expr_mat_subset <- expr_mat[,rownames(meta.data.top)]
  barcode.file <- colnames(expr_mat_subset) %>% as.data.frame()
  gene.file <- rownames(expr_mat_subset) %>% as.data.frame()
  #保存数据
  print("Finish!")
  write_delim(barcode.file,".//bc/barcodes.tsv",col_names = F)
  write_delim(gene.file,".//bc/genes.tsv",col_names = F)
  writeMM(expr_mat_subset,".//bc/matrix.mtx")
}

############# 测 试 #############
setwd("D:\\Desktop\\s_note\\data\\singel_cell\\seurat")

bar <- read.delim("barcodes.tsv",header = F)
gen <- read.delim("genes.tsv",header = F)
MTX <- readMM("matrix.mtx")
#提取数据并保存为10x格式
extract_top_cell(barcode = bar,expr_mat = MTX,gene = gen,n_top = 100,sort_by ="nCount_RNA")
[1] "rename the expression matrix......"
[1] "caculate nCount_RNA......"
[1] "caculate nFeature_RNA......"
[1] "create matrix......"
[1] "extract cell top of 100 by nCount_RNA"
[1] "Finish!"
#测试是否可以创建10x对象
test <- Read10X("bc/",gene.column = 1)
test.1 <- CreateSeuratObject(counts = test)
#成功创建Seurat对象
> test.1
An object of class Seurat 
32738 features across 100 samples within 1 assay 
Active assay: RNA (32738 features, 0 variable features)
```

### Q9 Write a script to calculate the sequencing saturation of the CeleScope rna test data(https://github.com/singleron-RD/celescope_test_script/tree/master/rna). (20)

```
```

