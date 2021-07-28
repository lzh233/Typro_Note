# 单细胞-10X官网Q&A

## Q1 Why do I see high levels of Malat1 in my gene expression data?

[Why do I see high levels of Malat1 in my gene expression data? – 10X Genomics](https://kb.10xgenomics.com/hc/en-us/articles/360004729092-Why-do-I-see-high-levels-of-Malat1-in-my-gene-expression-data-)

```
Answer: Malat1 is frequently detected in poly-A captured RNA-seq data, independent of protocol. The level of detection appears to be cell type-specific and also shows some inverse correlation with cell health. We have observed that dead/dying cells have higher expression of Malat1.
```

## Q2 What is sequencing saturation?

[What is sequencing saturation? – 10X Genomics](https://kb.10xgenomics.com/hc/en-us/articles/115005062366-What-is-sequencing-saturation-)

```shell
#sequencing saturation-测序饱和度，用来表征如果需要发现新的转录本，还需要在cell中测得reads数，一般来说稀释曲线是区域平缓的，及随着reads数，即测序深度的增加发现新的转录本的数量也趋近于0，与文库的复杂程度(不同的细胞类型会有不同数量的RNA，因此在最终的文库中不同转录本的总数也会不同)和测序深度有关,sequencing saturation达到50%说明，测序测到每两条reads，就会有一个新的umi被发现，90%则为测到10条reads中，只有1个新的umi被发现
```

原回答: 

```she
Sequencing saturation is a measure of the fraction of library complexity that was sequenced in a given experiment. The inverse of the sequencing saturation can be interpreted as the number of additional reads it would take to detect a new transcript.

Sequencing saturation is dependent on the library complexity and sequencing depth. Different cell types will have different amounts of RNA and thus will differ in the total number of different transcripts in the final library (also known as library complexity). The figure below illustrates the median number of genes recovered from different cell types. As sequencing depth increases, more genes are detected, but this reaches saturation at different sequencing depths depending on cell type. 

Sequencing depth also affects sequencing saturation; generally, the more sequencing reads, the more additional unique transcripts you can detect. However, this is limited by the library complexity. 
```



![image-20210714171325645](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210714171325645.png)

## Q3 How much sequencing saturation should I aim for?

[How much sequencing saturation should I aim for? – 10X Genomics](https://kb.10xgenomics.com/hc/en-us/articles/115002474263-How-much-sequencing-saturation-should-I-aim-for-)

```R
#测序饱和度是衡量文库复杂度和cell捕获效率的重要指标，一般来说如果实验目标只是为了区分不同的cell就并不需要很高的测序饱和度，但是如果是为了研究一些很稀有的转录本，则需要更高的测序饱和度。有些cell本身的rna含量很低，需要很少的测序就能达到90%以上的测序深度
```

原回答

```R
Answer: Sequencing saturation is a measure of the fraction of library complexity that was captured during sequencing. Depending on the goals of the experiment, you may or may not need high sequencing saturation. If you are aiming to cluster cells into populations for downstream analysis, it is not necessary to detect every unique transcript (UMI count) in each cell and a lower sequencing saturation may be sufficient. However, if you are trying to recover very lowly expressed transcripts, higher sequencing saturation may be required to detect these transcripts.

Primary cells (e.g. PBMCs) generally have lower RNA content and may require less sequencing to achieve sequencing saturation rates of >90%.

The inverse of the sequencing saturation can be interpreted as roughly the number of new transcripts you expect to find with one new read. If sequencing saturation is at 50%, it means that every 2 new reads will result in 1 new UMI count (unique transcript) detected. In contrast, 90% sequencing saturation means that 10 new reads are necessary to obtain one new UMI count. If the sequencing saturation is high, additional sequencing would not recover much new information for the library.

Note: In earlier versions of the Cell Ranger pipeline, the sequencing saturation metric was referred to as cDNA PCR duplication rates. The previous term may be more intuitive to some people. We see a wide range of cDNA PCR duplication rates. These represent reads that map to the same reference gene, have the same cell barcode, AND the same transcript UMI. 
```



## Q4 What parameters are used for STAR index?

[What parameters are used for STAR index? – 10X Genomics](https://kb.10xgenomics.com/hc/en-us/articles/115004415786-What-parameters-are-used-for-STAR-alignment-)



## Q5 How to interpret the "Fraction Reads in Cells" metric?

[How to interpret the "Fraction Reads in Cells" metric? – 10X Genomics](https://kb.10xgenomics.com/hc/en-us/articles/360003919491-How-to-interpret-the-Fraction-Reads-in-Cells-metric-)

```r
#Fraction Reads in Cells:即来自细胞内的reads
#高的非内源rna可能是因为细胞死亡或裂解，但是有些情况比较特殊，2) The cell-calling heuristic did not apply. For example, there may be higher variation in RNA content than expected (more cells with lower RNA content). The current cell-calling heuristic assumes a ten-fold variation in RNA content.cellranger主要是通过细胞barcodes一定多于背景barcode来判断Fraction Reads in Cells
```

原回答

```
A low "Fraction Reads in Cells" value is typically explained by the following:

1) High ambient RNA (background) in your sample. This ambient RNA comes from lysed/dead cells in your sample. Cell Ranger is able to confidently align the reads from ambient RNA to the transcriptome but the reads are not associated with a valid cell-containing GEM.

2) The cell-calling heuristic did not apply. For example, there may be higher variation in RNA content than expected (more cells with lower RNA content). The current cell-calling heuristic assumes a ten-fold variation in RNA content.

Cell Ranger's algorithm for partitioning barcodes as cells versus background is based on the idea that barcodes for cells should have distinctly more transcript counts associated with them than the background barcodes. This can be visualized by the ranked barcode plot in the web_summary.html file. More details on the cell filtering algorithm can be found here.

If you suspect that Cell Ranger's cell calling algorithm did not work well for your sample, please re-run cellranger count again or cellranger reanalyze with --force-cells option to call the expected number of cells.
```



## Q6 How much of a single cell's genome is amplified? 

[How much of a single cell's genome is amplified? – 10X Genomics](https://kb.10xgenomics.com/hc/en-us/articles/360005108931-How-much-of-a-single-cell-s-genome-is-amplified-)

```
对于人类细胞而言，每个细胞的基因组大概被扩增了0.2倍，但是这些片段并不会被全部测序到，人的cell中有750000对reads大概有0.05可以被测到
```

原回答

```
Answer: For human samples, approximately 0.20x of each cell's genome is amplified. However, not all of these library fragments are sequenced. At 750,000 read pairs per cell (human), approximately 0.05x coverage per cell is sequenced. This is enough to perform CNV calling at 2 Mb resolution at the single cell level, but not to reliably detect smaller CNV events or SNPs.

To sequence all of the amplified fragments, a much higher sequencing depth (>>750,000 read pairs per cell) is required. There are diminishing returns to sequencing to higher depths because many more duplicate reads are sequenced as a fraction of total reads.
```



## Q7 Why do I see a high level of mitochondrial gene expression?

[Why do I see a high level of mitochondrial gene expression? – 10X Genomics](https://kb.10xgenomics.com/hc/en-us/articles/360001086611-Why-do-I-see-a-high-level-of-mitochondrial-gene-expression-)

```r
#细胞的质量较低，有较多的细胞发生了凋亡或溶解，也有可能是细胞本身原因所造成的，如肿瘤细胞的代谢和其坏死都会增加线粒体的比例，也有可能是细胞样品被死掉或溶解的细胞污染，导致线粒体的含量增加
```

原回答

```
Mitochondrial genes are expressed in most cells and their expression level is cell type-specific.

1 High expression levels of mitochondrial genes could be an indicator of:

2 Poor sample quality, leading to a high fraction of apoptotic or lysing cells.
Biology of the particular sample, for example tumor biopsies, which may have increased mitochondrial gene expression due to metabolic activity and/or necrosis.
```

## Q8 Does Cell Ranger automatically exclude doublets?

[Does Cell Ranger automatically exclude doublets? – 10X Genomics](https://kb.10xgenomics.com/hc/en-us/articles/360001074271-Does-Cell-Ranger-automatically-exclude-doublets-)

```
目前cellranger还没有方法可以自动识别一个barcode对应多的个cell的barcode，后续使用seurat，doublefinder等r-package可以进行过滤，
```

```
Answer: We currently don't have a method for computationally classifying whether a barcode contains more than one cell in single cell gene expression data from a single species.

At present, Cell Ranger software only detects doublets in the context of a barnyard or mixed-species experiment used for estimating multiplet rates.  

In a barnyard experiment, a 10x library is made from a mixture of cells from different species and mapped to a custom multi-species reference which combines the corresponding reference genomes. Multiplet rates are then estimated by counting the number of cells where transcripts are detected from multiple genomes. For more information, please see the "Estimating Multiplets" section in the Algorithms overview.
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview


While we currently do not have an official recommendation on this subject, a couple of methods might be explored. Please note that some of these methods are manual/ad-hoc or experimental and may filter out valid single cells.

Infer doublets in a single species case if there are known cell-type specific markers. For example, the presence of T cell and B cell specific markers coming from a single barcode may indicate a GEM with both T and B cells.
Evaluate the suggested workflow in the R package, Seurat, which has some functionality for flagging cells with a clear outlier number of UMIs or genes detected http://satijalab.org/seurat/pbmc3k_tutorial.html ("QC and selecting cells for further analysis.)
Evaluate the doublet identification method in the R package DoubletFinder.
```



## Q9 Are there methods for identifying multiplets?

[Are there methods for identifying multiplets? – 10X Genomics](https://kb.10xgenomics.com/hc/en-us/articles/360005165411-Are-there-methods-for-identifying-multiplets-)

```
At present, there is no explicit support for multiplet identification in the Cell Ranger gene expression analysis pipeline other than for the analysis of "barnyard" or mixed-species datasets mapped to a multiple-species reference.

Below are some third-party tools for detecting doublets/multiplets in single-cell gene expression datasets:

Scrublet:
https://www.biorxiv.org/content/early/2018/07/09/357368
https://github.com/AllonKleinLab/scrublet 
DoubletDecon:
https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31286-0 
https://github.com/EDePasquale/DoubletDecon 
DoubletFinder:
https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30073-0  
https://github.com/chris-mcginnis-ucsf/DoubletFinder
Solo:
https://www.cell.com/cell-systems/fulltext/S2405-4712(20)30195-2
https://github.com/calico/solo/blob/master/README.md 
Please contact the respective developers of each tool for more information.

Disclaimer: This article and code-snippet are provided for instructional purposes only. 10x Genomics does not support or guarantee the code.
```

## Q10 Can I remove ambient RNA contamination from cells in my gene expression data ?

[Can I remove ambient RNA contamination from cells in my gene expression data ? – 10X Genomics](https://kb.10xgenomics.com/hc/en-us/articles/360060415811-Can-I-remove-ambient-RNA-contamination-from-cells-in-my-gene-expression-data-1)

```
cellranger没得办法，以下的r包可以解决
```

- [SoupX](https://academic.oup.com/gigascience/article/9/12/giaa151/6049831)
- [DecontX](https://www.biorxiv.org/content/10.1101/704015v1)
- [DIEM](https://www.nature.com/articles/s41598-020-67513-5)
- [CellBender](https://www.biorxiv.org/content/10.1101/791699v1.full)

## Q11 How can I annotate the cell types from scRNA-seq data?

[How can I annotate the cell types from scRNA-seq data? – 10X Genomics](https://kb.10xgenomics.com/hc/en-us/articles/360050278851-How-can-I-annotate-the-cell-types-from-scRNA-seq-data-)

```
cellranger不可以
```

- Seurat's [find.markers()](https://satijalab.org/seurat/seurat_clustering_tutorial_part2.html) function
- [Human Cell Atlas Data Portal](https://data.humancellatlas.org/)
- [CellMarker](http://biocc.hrbmu.edu.cn/CellMarker/): a manually curated resource of cell markers in human and mouse
- [Identifying cell populations with scRNASeq](https://www.sciencedirect.com/science/article/pii/S0098299717300493)
- [CHETAH](https://academic.oup.com/nar/article/47/16/e95/5521789): a selective, hierarchical cell type identification method for single-cell RNA sequencing
- [Tabula Muris](https://tabula-muris.ds.czbiohub.org/): a compendium of single-cell transcriptome data from mice
- [Tabula Sapiens](https://tabula-sapiens-portal.ds.czbiohub.org/home): a growing repository of detailed gene expression profile in human tissues