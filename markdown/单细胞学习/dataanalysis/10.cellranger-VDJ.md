# Cellranger-VDJ

## 10X vdj分析

`id`: 结果输出目录名称与位置

`reference`: 指定参考序列位置

`fastqs`: 指定序列目录

`localcores`: cpu数

`localmem`: 内存

`lanes`： 指定lans

`chain`: `TR` for T cell receptors;  `IG` for B cell receptors;  `auto` for autodetection based on TR vs IG representation (default)

`denovo`: 是否根据参考序列拼接，如果使用`denovo`未指定`reference`，则需要指定`inner-enrichment-primers`

`inner-enrichment-primers`: 包含引物序列的文件的名称，每行一个，如果有两套引物，那应该是最内侧的反向PCR引物，与恒定区互补。例如，对于人类TCR，具有两行的文件：

AGTCTCTCAGCTGGTACACG

 TCTGATGGCTCAAACAGC

将等同于默认引物(如果未指定此选项)。

```shell
cellranger vdj --id=sample_test --reference=/Personal/liuzihao/cellranger_test/vdj_ref/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 --fastqs=/Personal/liuzihao/cellranger_test/vdj_demo/fastq/vdj_v1_hs_pbmc3_t_26x91_fastqs/data_002 --localcores=8 --localmem=64

#otuput
$ ll -h
total 1.3G
-rw-r-----. 1 liuzihao ssh.randd  18M Aug 16 18:03 airr_rearrangement.tsv
-rw-r-----. 1 liuzihao ssh.randd 2.3M Aug 16 18:03 all_contig_annotations.bed
-rw-r-----. 1 liuzihao ssh.randd 5.2M Aug 16 18:03 all_contig_annotations.csv
-rw-r-----. 1 liuzihao ssh.randd  60M Aug 16 18:03 all_contig_annotations.json
-rw-r-----. 1 liuzihao ssh.randd 1.2G Aug 16 18:02 all_contig.bam
-rw-r-----. 1 liuzihao ssh.randd 945K Aug 16 18:02 all_contig.bam.bai
-rw-r-----. 1 liuzihao ssh.randd 6.2M Aug 16 18:03 all_contig.fasta
-rw-r-----. 1 liuzihao ssh.randd 565K Aug 16 18:03 all_contig.fasta.fai
-rw-r-----. 1 liuzihao ssh.randd  13M Aug 16 18:03 all_contig.fastq
-rw-r-----. 1 liuzihao ssh.randd  99K Aug 16 18:03 cell_barcodes.json
-rw-r-----. 1 liuzihao ssh.randd 631K Aug 16 18:02 clonotypes.csv
-rw-r-----. 1 liuzihao ssh.randd 1.4M Aug 16 18:03 concat_ref.bam
-rw-r-----. 1 liuzihao ssh.randd 570K Aug 16 18:03 concat_ref.bam.bai
-rw-r-----. 1 liuzihao ssh.randd 6.9M Aug 16 18:03 concat_ref.fasta
-rw-r-----. 1 liuzihao ssh.randd 337K Aug 16 18:03 concat_ref.fasta.fai
-rw-r-----. 1 liuzihao ssh.randd 4.4M Aug 16 18:02 consensus_annotations.csv
-rw-r-----. 1 liuzihao ssh.randd 977K Aug 16 18:03 consensus.bam
-rw-r-----. 1 liuzihao ssh.randd 570K Aug 16 18:03 consensus.bam.bai
-rw-r-----. 1 liuzihao ssh.randd 3.9M Aug 16 18:02 consensus.fasta
-rw-r-----. 1 liuzihao ssh.randd 324K Aug 16 18:02 consensus.fasta.fai
-rw-r-----. 1 liuzihao ssh.randd 4.6M Aug 16 18:03 filtered_contig_annotations.csv
-rw-r-----. 1 liuzihao ssh.randd 4.1M Aug 16 18:03 filtered_contig.fasta
-rw-r-----. 1 liuzihao ssh.randd 8.0M Aug 16 18:03 filtered_contig.fastq
-rw-r-----. 1 liuzihao ssh.randd  919 Aug 16 18:04 metrics_summary.csv
-rw-r-----. 1 liuzihao ssh.randd  22M Aug 16 18:04 vdj_contig_info.pb
drwxr-x--x. 3 liuzihao ssh.randd   53 Aug 17 09:29 vdj_reference
-rw-r-----. 1 liuzihao ssh.randd 2.4M Aug 16 18:03 vloupe.vloupe
-rw-r-----. 1 liuzihao ssh.randd 2.3M Aug 16 18:04 web_summary.html
```

## cellreanger实现VDJ分析的算法

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210817133940573.png" alt="image-20210817133940573" style="zoom: 67%;" />

### Assembly Algorithm

reads被被输入后通过overlap进行组装，生成contig，这些contig代表对当前转录序列的最佳估计，如下图就是，一个barcode里的序列被assembly成两个contig，同时每个contig中的碱基都有一个质量值，同时也会得到每个contig对应的reads数量和umi数量

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210817134549152.png" alt="image-20210817134549152" style="zoom: 67%;" />

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210817135037811.png" alt="image-20210817135037811" style="zoom:67%;" />

生成contig是非常复杂的，因为数据中含有非常多的噪音，包括来自与胞外mRNA，cell doublets，errors in transcription in the cell, errors in reverse transcription to make cDNA, 测序错误.....

组装过程在某些地方使用参考序列, 在denovo模式下运行则不需要

https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/algorithms/glossary#full（名词解释）

#### Step of Assembly

https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/algorithms/assembly

**step1. subsample：**

抽平，将每个cell中的reads数量抽平到至多8000的水平，一个特定细胞的高测序深度可能是由于总体较高的测序深度，或者更常见的是由于BCR在浆细胞中的高表达。cell的覆盖率只增加了很少的信息，但降低了计算性能。因此，我们限制了每个条形码的reads。每个cell最多有80,000条reads。如果有超过这个数字，需要进行subsample。

**step2. Read trimming**

切除富集引物的序列

**step3. Graph formation**

构建`De Bruijn graph`(一种通过构建k-mer进行Assembly的算法)https://www.homolog.us/blogs/genome/2011/07/28/de-bruijn-graphs-i/，默认`k-mer = 20`,并将其转换为一个边为DNA序列的有向图

**step4. Reference-free Graph Simplification**

Simplify the graph by removing 'noise' edges.（通过取出”噪声“的边简化图）

```
A collection of heuristic steps are applied to simplify the graph. During this process we track and edit the read support on each edge. We describe here several examples of simplification steps.

Branch cleaning.

(I) For each branch in the graph, and for each UMI, if one branch has ten times more reads from the UMI supporting it than the other branch has from the UMI, we remove read support for the UMI from the other branch.

(II) Given two branches emanating from a vertex, this deletes the weak branch if all of the following hold: (a) there are at least twice as many reads on the strong branch; (b) the weak branch does not have more than 8 reads for any UMI; (c) for every UMI, the strong branch has at least as many reads, with at most one possible exception, where it can have just one read.

Path cleaning. For each UMI, define a single strongest path for the UMI. Then delete any graph edges that are not on such a strong path.

Component cleaning. For each UMI, if one graph component has ten times more reads supporting it from that UMI than another component, delete the read support for that UMI from the other components.
```

**step5. Reference-assisted graph simplification**

根据参考数据库对图进一步简化

```
If the pipeline is not run in denovo mode, we pop bubbles in the graph with the aid of the reference sequence. There are several heuristic tests, all of which require that both bubble branches have the same length. For example, if in addition, one branch is supported by at least three UMIs, and the other branch is supported by only one UMI, and the first branch has a kmer matching the reference, we delete the weak branch.
```

**step6. UMI filtering**

对UMI进行过滤，过滤的标准

```
As above, we define a single strongest path for each UMI. We consider only strong paths that either contain a reference kmer, or in the denovo case, match a primer.

Find good graph edges: edges that appear on one or more strong paths.

We find all the reads assigned to one or more good edge.

Find the UMIs for these reads.

Remove any UMI for which less than 50% of the kmers in the reads for the UMI are contained in good edges.

In the non-denovo case, if no strong path had a V segment annotation, remove all the UMIs.
```

**step7. 构建contig**

通过在graph中为每个UMI寻找最佳路径来创建contigs

**step8. Competitive deletion of contigs**

删除弱的barcode和可能是噪音的contig

**step9. Contig confidence**

定义高可信度的contigs，它可能代表与cell相关联的单个细胞的真实转录本。

**step10. Contig quality scores**

Assign a quality score to each base on each contig.

```
Each base in the assembled contigs is assigned a Phred-scaled Quality Value (QV), representing an estimate of the probability of an error at that base. The QV is computed with a hierarchical model that accounts for the errors in reverse transcription (RT), that will affect all reads with the same UMI, and sequencing errors which affect individual reads. The sequencing error model uses the reported sequencer QVs. At recommended sequencing depths, many reads per UMI are observed, so sequencing errors in individual reads are corrected rapidly. We estimate that the V(D)J RT reaction has an error rate of 1e-4 per base, so assembled bases that are covered by a single UMI will be assigned Q40, and bases covered by at least two UMIs will be assigned Q60.
```

#### About De Bruijn Graph

De Bruijn graph是一个展示符号序列之间重叠关系的有方向的示意图

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210817152556166.png" alt="image-20210817152556166" style="zoom:67%;" />

**k-mer**

把reads打断成长度为k的核酸片段（kmer）后，利用De Bruijn graph根据kmer间重复的部分构建graph，得到最优化路径从而拼接contig。

```
我们可以将一条reads切割成很多小的Kmer片段，从第一个碱基开始，每隔固定距离的碱基开始提取碱基。例如一条100bp长的reads，每隔一个位置取一个17mer的片段。也就是1-17取一个kmer出来，2-18取一个Kmer，3-19取一个Kmer。以此类推，最终84-100为1个kmer。那么最终将会生成将会生成100-17+1，也就是84个Kmer片段。原来一条长片段，变成了很多短的片段，碱基的数目也增加了很多倍。而且，每次取Kmer是同一条reads正反取两次，也就是对这条reads的反向互补序列再取一次Kmer
```

优点主要有两个：

> 我们将测序的reads切割成更短的kmer，形成了一个巨大的kmer的集合。所有的kmer片段长度相同，可以计算kmer出现的频数,绘制一张kmer频数的分布图。其中有一个峰值，这个峰值就是kmer的深度。可以利用这个值估计基因组大小。

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210817152848148.png" alt="image-20210817152848148" style="zoom:67%;" />

> 利用kmer可以用于去除测序错误碱基，得到更好的拼接效果。

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210817152916851.png" alt="image-20210817152916851" style="zoom:67%;" />

https://mp.weixin.qq.com/s/12F3MQaRoKm4qvVOMrpKcA

```python
#一个计算kmer的小脚本
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 09:06:06 2021

@author: liuzihao

create k-mer
"""

class k_mer:
    def __init__(self,seq,k_mer_len):
        self.seq = seq
        self.seq_len=len(seq)
        self.k_mer_len = k_mer_len
        self.k_mer_num = self.seq_len - self.k_mer_len + 1
        if self.k_mer_len > self.seq_len:
            print("Please input correct k-mer")
    def get_kmer(self):
        kmer_lst={}
        for base in range(self.k_mer_num):
            k = self.seq[base:base+self.k_mer_len]
            if k in kmer_lst:
               kmer_lst[k] =  kmer_lst[k] + 1
            else:
               kmer_lst[k] = 1
        #kmer_lst = pd.DataFrame(kmer_lst)
        return kmer_lst

k_lst = k_mer(seq = "ATTGGGCCCC",k_mer_len=3)
print(k_lst.get_kmer())
#输出
{'ATT': 1, 'TTG': 1, 'TGG': 1, 'GGG': 1, 'GGC': 1, 'GCC': 1, 'CCC': 2}
```



### VDJ Cell Calling Algorithm

在10x系统中，在液滴(GEMs)，一些含有细胞，而另一些含有T细胞或B细胞。T细胞和B细胞的发现依赖于VDJ区域的转录本的鉴定与计数，一些T和B细胞的VDJ区域表达量特别低，所以可能难以被检测到，相反，高水平的细胞外mRNA可能会导致一些cell被错误地识别为T或B细胞，因此，VDJ Cell Calling Algorithm的目标是近似包含T或B ，此算法是assembly的一部分

要识别为T或B细胞，barcode必须满足以下三个要求：

1. 其数量要高，高度可信，如果只用一个contig那么他的umi要大于1( (In the denovo case, we require only that there is a contig.), 虽然其他类型的细胞可以在TCR和BCR位点内表现出转录，但只有T和B细胞产生完全重排的转录本，其中既包含V片段，也包含C片段，因此数量较高的contig是认为该cell为T细胞和B细胞的证据之一，然而，转录本有可能是背景值--存在于cell外，而不是存在于完整的细胞中。因此需要这个contig有一个以上的UMI。

2. 必须至少有3个过滤的UMI，每个UMI至少有两个reads。这降低了仅仅根据背景转录本将细胞误认为T细胞或B细胞的可能性。

3. 计算每个barcode里的N50值，If for a given barcode, the maximum read pair count across filtered UMIs is less then 3% of this N50, do not call the barcode a cell**.(UMI过滤后数量最大的reads的长度小于N50的3%，则不认为该barcode为细胞？)** This provides some protection against transcripts arising from index hopping on an Illumina flowcell, and from other forms of cross-library contamination.

```
N50：Reads拼接后会获得一些不同长度的Contigs。将所有的Contigs长度相加，获得Contig总长度。然后将所有的Contigs按照从长到短进行排序，如获得Contig1，Contig2，Contig3…Contig25。将Contigs按照这个顺序依次相加，当相加的长度达到Contig总长度的一半时，最后一个加上的Contig长度即为Contig N50，可以作为基因组拼接的结果好坏的一个判断标准。
```

### Annotation Algorithm

注释的目的主要是将contig上的V、D、J区域进行区分，确定CDR3的序列，并且通过这些数据确定该contig是否有效(productive), 

#### Alignment to the V(D)J Reference

比对到VDJ的参考数据库

首先确定数据是TCR还是BCR(示例数据VDJ TCR - FASTQs为TCR数据)，后将这些contigs比对到参考序列（TCE/BCR），**In rare (mixed) cases contigs are aligned to both. Alignment is seeded on 12-mer perfect matches, followed by heuristic extension; we also search backward from C segment alignments for J segment alignments that do not have 12-mer perfect matches, as these will arise occasionally from somatic hypermutation.?**

V(D)J参考序列的选择可以是任意的，这取决于参考序列彼此有多相似。 For D segments, which are both short and more mutated, it is often not possible to find a confident alignment, and an alignment may not be shown.

#### Productive Contigs

有效的Contigs标准

**Full length requirement**：contig与V区域首先匹配，直到 J 区域结束

**Start requirement**：V的起始部分与contig上的起始密码子匹配。

**Nonstop requirement：** V区到J区之间没有起始密码子

**In-frame requirement：**The J stop minus the V start equals one mod three. This just says that the codons on the V and J segments are in frame.？？

**CDR3 requirement.** There is an annotated CDR3 sequence（必须有注释的CDR3序列）

**Structure requirement.** Let VJ denote the sum of the lengths of the V and J segments. Let len denote the J stop minus the V start, measured on the contig. Then VJ - len lies between -25 and +25, except for IGH, which must be between -55 and +25. This condition is imposed to preclude anomalous structure changes that are unlikely to correspond to functional proteins.

#### CDR3

对于每个重叠群，我们利用CDR3的侧翼保守序列来搜索CDR3序列。

我们比较了来自V和J参考片段的人和小鼠的基序，完全如下所示。在这里，字母代表一种特定的氨基酸，圆点代表任何氨基酸。

left flank: 左翼序列 

right flank右翼序列

CDR3序列

```
left flank   CDR3   right flank
LQPEDSAVYY   C...   LTFG.GTRVTV
VEASQTGTYF          LIWG.GSKLSI
ATSGQASLYL
```

我们要求CDR3序列的长度在5到27个氨基酸之间，以C开头，并且不包含终止密码子。将候选CDR3的侧翼序列与上述基序进行匹配，并为与一列中的一个条目匹配的每个位置打+1分。

```
LTY.... 
```

前三种氨基酸得分为2。(L与第一列中的条目匹配，因此分数为1。T与第二列中的条目匹配，因此分数为1。Y与第三列不匹配，因此不会影响分数。)。要将候选CDR3宣布为CDR3序列，其得分必须至少为10。此外，左侧必须至少贡献3分，右侧必须至少贡献4分。接下来，我们将找到重叠群上V片段末端的隐含停止位置。这是重叠群上V片段的起始位置，加上V片段的长度。然后我们要求CDR3序列在停止之前最多开始10个碱基，在V停止之后最多开始20个碱基(本段的条件不适用于重新开始的情况)。如果有多个CDR3序列，我们选择得分最高的序列。如果平局，我们选择重叠群上起始位置较晚的那个。如果仍然是平局，我们会选择较长的CDR3。

### Clonotype Grouping

在克隆型分组阶段，cellbarcode放置在称为克隆型的组中。根据计算，每个克隆由一个共同的祖先和来自于这个共同祖先的重排组成。在这个过程中，一些cellbarcode被标记为可能的噪音，并被过滤掉。

```
During the clonotype grouping stage, cell barcodes are placed in groups called clonotypes. Each clonotype consists of all descendants of a single, fully rearranged common ancestor, as approximated computationally. During this process, some cell barcodes are flagged as likely artifacts and filtered out, meaning that they are no longer called as cells.
```

**T细胞：**T细胞受体(TCR)缺乏体细胞高频突变，产生具有相同V(D)J转录本的生物克隆类型。技术误差(例如，在逆转录中出现)可能导致计算出的克隆类型具有差异。。

```
T cells. The lack of somatic hypermutation in T cell receptors (TCRs) yields biological clonotypes which have identical V(D)J transcripts. Technical artifacts (e.g. arising in reverse transcription) can result in the computed clonotypes having isolated differences. These are generally rare.
```

**B细胞：**完全重排的B细胞受体(BCR)可以发生体细胞高频突变(SHM)，从而增加抗原亲和力。因此，对于BCR，克隆类型中的V(D)J转录本在任何位置都可能不同，如下所示：由于SHM会产生大量的突变因此分辨B细胞克隆型不好搞 Cell Ranger version 5.0 and onwards accomplishes this by invoking a module for clonal analysis called `enclone`, which simultaneously groups cells into clonotypes and filters some out.

<img src="C:\Users\liuzihao\AppData\Roaming\Typora\typora-user-images\image-20210818104227965.png" alt="image-20210818104227965" style="zoom: 67%;" />

**具体分类算法**
https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/algorithms/clonotyping



