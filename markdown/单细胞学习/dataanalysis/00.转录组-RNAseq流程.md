# RNA-seq分析

## 软件下载

`fastqc 、muticqc 、Trimmomatic、 STAR 、rsem`，均使用**conda**安装

```shell
conda install fastqc
conda install muticqc 
conda install Trimmomatic
conda install STAR
conda install rsem
```

<u>***用于实验室服务器无法联网所以进行了离线创建***</u>

**1**、**首先在可以联网的服务器或虚拟机中出创建好相应的环境**

**2**、**打包所创建的环境：**

tar -tar -zcf [压缩文件名] [环境所在目录] &

**3**、**打包创建该环境所需的依赖：**

tar -zcf [压缩文件名] ~/miniconda2/pkgs &

**4**、**将两压缩包上传服务器并解压**

**5**、**（安装依赖包）**

**5-1** 将 ~/miniconda2/pkgs备份：mv ~/miniconda2/pkgs ~/miniconda2/pkgsbak

**5-2** 将解压后的新pkgs复制到~/miniconda2/

6**（创建环境）**conda create -n lzh --clone [指定步骤2中解压后目录的位置] --offline

7、测试环境没问题后，删除原来的包，rm -rf ~/miniconda2/pkgsbak/

## 参考基因组下载

**参考基因组的下载：**

1、ensembl数据库

2、JGI

.....

```shell
#参考基因组为
Araport11_GFF3_genes_transposons.201606.gtf
TAIR10_Chr.all.fasta
```

### gtf文件的解释

**gtf文件共有9列**，以`\t`分割

```shell
Chr1    Araport11    start_codon     7448    7450    .       -       .       transcript_id "AT1G01020.2"; gene_id "AT1G01020";
Chr1    Araport11       exon    7564    7649    .       -       .       transcript_id "AT1G01020.2"; gene_id "AT1G01020";
Chr1    Araport11       CDS     7564    7649    .       -       0       transcript_id "AT1G01020.2"; gene_id "AT1G01020";
Chr1    Araport11       exon    7762    7835    .       -       .       transcript_id "AT1G01020.2"; gene_id "AT1G01020";
Chr1    Araport11       CDS     7762    7835    .       -       2       transcript_id "AT1G01020.2"; gene_id "AT1G01020";
Chr1    Araport11       exon    7942    7987    .       -       .       transcript_id "AT1G01020.2"; gene_id "AT1G01020";
Chr1    Araport11       CDS     7942    7987    .       -       0       transcript_id "AT1G01020.2"; gene_id "AT1G01020";
Chr1    Araport11       exon    8236    8325    .       -       .       transcript_id "AT1G01020.2"; gene_id "AT1G01020";
Chr1    Araport11       CDS     8236    8325    .       -       0       transcript_id "AT1G01020.2"; gene_id "AT1G01020";
Chr1    Araport11       exon    8417    8464    .       -       .       transcript_id "AT1G01020.2"; gene_id "AT1G01020";
Chr1    Araport11       CDS     8417    8464    .       -       0       transcript_id "AT1G01020.2"; gene_id "AT1G01020";
Chr1    Araport11    stop_codon      8571    8573    .       -       .       transcript_id "AT1G01020.2"; gene_id "AT1G01020";

```

**seq_id：**序列的编号，一般为chr或者scanfold编号；

**source: **注释的来源，一般为数据库或者注释的机构，如果未知，则用点“.”代替；

**type: **注释信息的类型，比如Gene、cDNA、mRNA、CDS等

**start:**该基因或转录本在参考序列上的起始位置；

**end: **该基因或转录本在参考序列上的终止位置；

**score: **得分，数字，是注释信息可能性的说明，可以是序列相似性比对时的E-values值或者基因预测是的P-values值，“.”表示为空；

**strand: **该基因或转录本位于参考序列的正链(+)或负链(-)上;

**phase: **仅对注释类型为“CDS”有效，表示起始编码的位置，有效值为0、1、2(对于编码蛋白质的CDS来说，本列指定下一个密码子开始的位置。每3个核苷酸翻译一个氨基酸，从0开始，CDS的起始位置，除以3，余数就是这个值，表示到达下一个密码子需要跳过的碱基个数。该编码区第一个密码子的位置，取值0,1,2。0表示该编码框的第一个密码子第一个碱基位于其5'末端；1表示该编码框的第一个密码子的第一个碱基位于该编码区外；2表示该编码框的第一个密码子的第一、二个碱基位于该编码区外；如果Feature为CDS时，必须指明具体值。)；

**attributes:**一个包含众多属性的列表，格式为“标签＝值”（tag=value），标签与值之间以空格分开，且每个特征之后都要有分号；（包括最后一个特征），其内容必须包括gene_id和transcript_id。以多个键值对组成的注释信息描述，键与值之间用“=”，不同的键值用“；

![img](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/14977949-7e0452e2298ba318.png)

## 数据质控

### fastqc

使用**`fastqc`**软件对数据进行质量查看，基本用法如下,

| -o    | 指定结果文件的输出目录，目录应提前建好                       |
| ----- | ------------------------------------------------------------ |
| -f    | 指定序列格式，如`bam,sam`,`bam_mapped`,`sam_mapped`  以及`fastq`(默认) |
| -t    | 指定多线程                                                   |
| -help | 查看帮助，其他的选项                                         |

```shell
#对序列进行质控
mkdir ./fastqc
fastqc ./*.fastq.gz -o ./fastqc/ -t 32
#对应每个序列生成了一个html质控报告，zip里存放了所有html中的结果，详细看fastqc结果解读篇
(lzh2) [liuzihao@bogon fastqc]$ la | head
total 7.1M
-rw-r--r--. 1 liuzihao student 556K Apr 14 15:40 T1_1_R1_fastqc.html
-rw-r--r--. 1 liuzihao student 338K Apr 14 15:40 T1_1_R1_fastqc.zip
```

### mulitqc

由于序列较多，逐一查看麻烦所以使用**`multiqc`**对fastqc结果进行整合查看

```shell
#整合fastqc结果
multiqc ./fastqc/
#结果存放位置，html为报告，multiqc_data存放了所有html中的结果
(lzh2) [liuzihao@bogon rawdata]$ ll -d  m*
drwxr-xr-x. 2 liuzihao student     136 Apr 14 15:48 multiqc_data
-rw-r--r--. 1 liuzihao student 1195068 Apr 14 15:48 multiqc_report.html
```

## 数据过滤

### Trimmomatic

使用**`Trimmomatic`**对序列进行过滤与接头删除

### **基本用法**

![(Izh) [root@VM-€)-4-centos —]# trimmomatic  - -help  Usage :  PE [ -version] [ -threads <threads>] [ -phred331 -phred64] [ -trimlog [ -summar  y <statsSummaryFi1e>] [ -quiet] [ -validatePairs] [ -basein <inputBase> I <inputFi1e2  >] [ -baseout <outputBase> I <trimm  SE [ -version] [ -threads <threads>] [ -phred331 -phred64] [ -trimlog [ -summar  y [ -quiet] <trimmerl>. . .  -versxon ](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image001-1618391970026.png)

```shell
#对序列进行质控
trimmomatic PE -threads 16 \
./T1_1_R1.fastq.gz T1_1_R2.fastq.gz \
../cleandata/T1_1_R1_pair.fastq.gz \
../cleandata/T1_1_R1_unpair.fastq.gz \
../cleandata/T1_1_R2_pair.fastq.gz \
../cleandata/T1_1_R2_unpair.fastq.gz ILLUMINACLIP:/home/liuzihao/miniconda2/envs/lzh/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:2:30:10:1:true \
LEADIND:3 TRAILIND:3 \
SLIDINGWINDOW:4:20 MINLEN:50 TOPHRED33
```

### 选项参数

**`PE -threads 16`**: PE/SE指定双端还是单端序列，-thresds线程数

**`./T1_1_R1.fastq.gz T1_1_R2.fastq.gz`**: 指定输入文件

**`../cleandata/T1_1_R1_pair/unpair.fastq.gz`** : 指定输出文件，一组输出文件会有四个输出文件，分别是**T1_1_R1/2_pair.fastq.gz（表示质控过后双端数据全部保留的序列）**、**T1_1_R1/2_unpair.fastq.gz（表示质控过后双端数据未全部保留的序列，没有用了）**

**`ILLUMINACLIP:/home/liuzihao/miniconda2/envs/lzh/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:2:30:10:1:true `**: 

指定测序的接头，应该根据质控步骤中所显示的接头类型进行选择，**2:30:10:1意义见最后**

**`LEADIND:3 TRAILIND:3`**: `LEADING`从reads**开头**切除低于阈值的碱基；`TRAILING`从reads**结尾**切除低于阈值的碱基

**`SLIDINGWINDOW:4:20`** : 滑窗过滤，删除4个碱基平均质量低于20的滑窗

**`MINLEN:50`**: 序列剪切后，删除低于长度小于阈值的碱基

**`TOPHRED33`**: 指定质量体系，一般为33体系

### **结果查看**

T1-1为例，**有78.40%的序列被保留(Both Surviving)，即为paired的fastq文件**，**只有正向被保留的有14.26%（Forward Only Surviving）**，**只有反向被保留（Reverse Only Surviving）的有2.39%**，**扔掉的序列占(Dropped)4.94%**，最终得到想要过滤后以及删除接头的序列，是否需要删除接头应该从fastqc结果来看是否存在测序接头。`summary.txt 为下面脚本生成的`

```shell
(lzh2) [liuzihao@bogon rawdata]$ cat summary.txt 
T1_1
Input Read Pairs: 250000 Both SurvivinD: 196007 (78.40%) Forward Only SurvivinD: 35647 (14.26%) Reverse Only SurvivinD: 5987 (2.39%) Dropped: 12359 (4.94%)
T1_2
Input Read Pairs: 250000 Both SurvivinD: 213148 (85.26%) Forward Only SurvivinD: 23101 (9.24%) Reverse Only SurvivinD: 6381 (2.55%) Dropped: 7370 (2.95%)
······
```

### **批量质控**

```shell
#!/bin/bash
#文件命名规则：[处理名]_[重复编号]_R[1,2].fastq.gz
#Author:liuzihao
#Date:2021-4-14
#CentOS Linux release 7.8.2003 (Core)	
#Architecture:          x86_64
#CPU op-mode(s):        32-bit, 64-bit
#Byte Order:            Little Endian
#CPU(s):                144
#Model name:            Intel(R) Xeon(R) Gold 5220 CPU @ 2.20GHz
#获取处理名称
ls *.fastq.gz | cut -d "_" -f 1-2 | uniq > trname.txt
mkdir -p ./cleandata/pair
mkdir ./cleandata/unpair
function file_num(){
	count=0
	total=$(wc -l trname.txt | cut -d " " -f 1)
}
file_num
for i in $(cat trname.txt)
	do
		count=$(( $count + 1 ))
		echo -ne "共${total}个文件，正在质控第${count}个\r"
		trimmomatic PE -threads 16 ./${i}_R1.fastq.gz ${i}_R2.fastq.gz ./cleandata/pair/${i}_R1_pair.fastq.gz ./cleandata/unpair/${i}_R1_upair.fastq.gz ./cleandata/pair/${i}_R2_pair.fastq.gz ./cleandata/unpair/${i}_R2_upair.fastq.gz ILLUMINACLIP:/home/liuzihao/miniconda2/envs/lzh2/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:2:30:10:1:true LEADIND:3 TRAILIND:3 SLIDINGWINDOW:4:20 MINLEN:50 TOPHRED33 &> log.txt
		echo ${i} >> summary.txt
		grep -w "Input Read Pairs" log.txt >> summary.txt
	done
echo "质控完成！"
rm -rf log.txt
cat summary.txt
```

```shell
#输出目录查看
cleandata/
├── pair
│   ├── T1_1_R1_pair.fastq.gz
│   ├── T1_1_R2_pair.fastq.gz
└── unpair
    ├── T1_1_R1_upair.fastq.gz
    ├── T1_1_R2_upair.fastq.gz
······
```

## 序列比对

### 构建参考基因组索引（STAR）

使用`STAR`软件构建参考序列索引

```shell
STAR --runThreadN 32 --runMode genomeGenerate \
--genomeDir ../index_STAR/ \
--genomeFastaFiles TAIR10_Chr.all.fasta \
--sjdbGTFfile Araport11_GFF3_genes_transposons.201606.gtf \
--sjdbOverhang 149
#构建成功~WARNING似乎并无影响
Apr 14 22:34:21 ..... started STAR run
Apr 14 22:34:21 ... starting to generate Genome files
Apr 14 22:34:22 ..... processing annotations GTF
!!!!! WARNIND: --genomeSAindexNbases 14 is too large for the genome size=119667750, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 12
Apr 14 22:34:24 ... starting to sort Suffix Array. This may take a long time...
Apr 14 22:34:25 ... sorting Suffix Array chunks and saving them to disk...
Apr 14 22:34:33 ... loading chunks from disk, packing SA...
Apr 14 22:34:35 ... finished generating suffix array
Apr 14 22:34:35 ... generating Suffix Array index
Apr 14 22:35:06 ... completed Suffix Array index
Apr 14 22:35:06 ..... inserting junctions into the genome indices
Apr 14 22:36:05 ... writing Genome to disk ...
Apr 14 22:36:05 ... writing Suffix Array to disk ...
Apr 14 22:36:06 ... writing SAindex to disk
Apr 14 22:36:07 ..... finished successfully
```

**参数含义**

**`--runThreadN 32 --runMode  genomeGenerate `**: 指定线程数以及运行模式

**`--genomeDir`**: 指定输出目录     

**`--genomeFastaFiles / --sjdbGTFfile `**: 基因组文件以及注释文件

**`--sjdbOverhang`** : 指定可变剪切预测，一般设定为reads长度减1

### mapping（STAR）

使用`STAR`软件进行比对

```shell
STAR --runThreadN 8 \
--genomeDir ../index_g/ \
--readFilesCommand zcat \
--readFilesIn ./T1_1_R1_pair.fastq.gz T1_1_R2_pair.fastq.gz \
--outFileNamePrefix ../aligen/T1_1_
--outSAMtype BAM SortedByCoordinate \
--outBAMsortingThreadN 8 \
--quantMode TranscriptomeSAM GeneCounts
#比对成功
#Apr 15 20:49:37 ..... started STAR run
Apr 15 20:49:37 ..... loading genome
Apr 15 20:49:39 ..... started mapping
Apr 15 20:49:52 ..... finished mapping
Apr 15 20:49:52 ..... started sorting BAM
Apr 15 20:49:53 ..... finished successfully
Uniquely mapped reads
```

**`--genomeDir `**：指定索引目录

**`--readFilesCommand`**： 指定查看序列文件用的命令

**`--readFilesIn`**：指定输入文件，正向+反向

**`--outFileNamePrefix`**：指定输出文件位置和命名规则

**`--outSAMtype BAM SortedByCoordinate`**：将输出的sam文件转换为bam文件，并且基于位置进行排序

**`--outBAMsortingThreadN`** ：排序的线程数

**`--quantMode TranscriptomeSAM GeneCounts`** ：基于转录本进行比对，输出每一个基因上有几个reads匹配`T1_2_ReadsPerGene.out.tab`

**比对输出的文件：**

![image-20210415215135110](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210415215135110.png)

### 批量mapping

```shell
#!/bin/bash
#文件命名规则：[处理名]_[重复编号]_R[1,2]_pair.fastq.gz(上一个脚本的输出文件)
#Author:liuzihao
#Date:2021-4-14
#CentOS Linux release 7.8.2003 (Core)	
#Architecture:          x86_64
#CPU op-mode(s):        32-bit, 64-bit
#Byte Order:            Little Endian
#CPU(s):                144
#Model name:            Intel(R) Xeon(R) Gold 5220 CPU @ 2.20GHz
#脚本放在fastq文件同一目录，文件名称为：T_序列名_R1_pair.fastq.gz T_序列名_R2_pair.fastq.gz
#用法： ./aligan.sh [基因组位置]
#获取文件列表并mapping
ls *.fastq.gz | cut -d "_" -f 1-2 | uniq > ls.txt
for i in $(cat ls.txt)
	do
		STAR --runThreadN 8 --genomeDir $1 --readFilesCommand zcat --readFilesIn ./${i}_R1_pair.fastq.gz ${i}_R2_pair.fastq.gz --outFileNamePrefix ./aligan/${i}_ --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 8 --quantMode TranscriptomeSAM GeneCounts
	done
rm -rf ls.txt
#对输出文件分类
cd ./aligan
#比对日志文件
mkdir log
#比对结果的bam文件
mkdir ali_result
#其他文件
mkdir others
#进行分类
mv ./*final* ./log
mv ./*bam ./ali_result
mv ./T* ./others
#汇总基本比对结果，并展示
cd ./log
ls T* > ls.txt
for a in $(cat ls.txt)
	do
		echo ${a} >> summary.txt
		grep "Uniquely mapped reads" ${a}>> summary.txt
	done
rm -rf ls.txt
cat summary.txt
```

**基本比对情况查看**

```shell
T1_1_Log.final.out
                   Uniquely mapped reads number |	192956
                        Uniquely mapped reads % |	98.44%
T1_2_Log.final.out
                   Uniquely mapped reads number |	209905
                        Uniquely mapped reads % |	98.48%
T2_1_Log.final.out
                   Uniquely mapped reads number |	190153
                        Uniquely mapped reads % |	98.67%
T2_2_Log.final.out
                   Uniquely mapped reads number |	207634
                        Uniquely mapped reads % |	98.37%
```

### sam文件的解释

由于生成了bam文件，所以使用`samtools`进行一下转换，

```shell
samtools view -O sam -o ./T1_1_Aligned.sortedByCoord.out.sam ./T1_1_Aligned.sortedByCoord.out.bam
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210417195715040.png" alt="image-20210417195715040" style="zoom:150%;" />

**第1列：**reads的名称

**第2列：**表示比对的结果，以数字标识，类似于linux权限中的421**以83为例：（64+16+2+1）表示paired-end reads中的第一个reads比对到参考序列上了**,

[Explain SAM Flags (broadinstitute.github.io)](https://broadinstitute.github.io/picard/explain-flags.html), 通过这个网站可以轻松将比对的数字标识转换为具体比对情况

![](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image001-1618660783240.png)

**第3列：**表示参考序列的名称

**第4列：**表示比对的序列在基因组上的起始位置，没有匹配上则为0

**第5列：**比对的质量值，数字越大，特异性越好

**第6列：**比对的详细情况，如，`M`: 缺失；`I`: 代表插入；`D`: 代表缺失； `37M1D2M1I`: 37个匹配，1个参考序列上的删除，2个匹配，1个插入；具体如下

![ ](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image001-1618661142675.png)

**第7列：**表示双端测序的下一条reads比对到的参考序列，**=表示比对到了同一参考序列，*则代表没有比对到**

**第8列：**表示双端测序的下一条reads比对到的参考序列的位置，没有则是0

**第9列：**序列的长度

**第10列：**序列具体信息

**第11列：**序列质量值信息

**第12列：**如下，

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image001-1618661482930.png" style="zoom: 50%;" />

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image001-1618661510305.png" style="zoom: 80%;" />

## 基因定量（RSEM)

使用3个软件`RSEM`、 `kallisto`、`featurecounts`分别对基因进行定量

#### 构建RSEM索引

首先通过`rsem-prepare-reference`对参考基因组和GTF文件构建索引，

| --gtf ~/study/ref_seq/Araport11_GFF3_genes_transposons.201606.gtf | 指定基因组注释文件 |
| ------------------------------------------------------------ | ------------------ |
| ~/study/ref_seq/TAIR10_Chr.all.fasta                         | 指定基因组文件     |
| ./rsemindex/rsem_ref                                         | 指定输出目录       |

```shell
#创建索引存放位置
mkdir rsemindex
#构建索引
rsem-prepare-reference --gtf ~/study/ref_seq/Araport11_GFF3_genes_transposons.201606.gtf ~/study/ref_seq/TAIR10_Chr.all.fasta ./rsemindex/rsem_ref
#构建成功~，rsem_ref为索引名称
rsem-preref ./rsem_ref.transcripts.fa 1 ./rsem_ref
Refs.makeRefs finished!
Refs.saveRefs finished!
./rsem_ref.idx.fa is generated!
./rsem_ref.n2g.idx.fa is generated!
#查看索引文件结构rsem_ref为索引名称
tree rsemindex
rsemindex/
├── rsem_ref.chrlist
├── rsem_ref.grp
├── rsem_ref.idx.fa
├── rsem_ref.n2g.idx.fa
├── rsem_ref.seq
├── rsem_ref.ti
└── rsem_ref.transcripts.fa
```

**rsem_ref**：为索引的名称，并非目录，后续分析需要指定索引的名称，**而不是单单指定位置！！否则会报错，索引文件不存在**

#### 基因定量

通过`rsem-calculate-expression`对mapping结果定量

` --paired-end`  :指定数据类型 

` -q ~/study/rawdata/cleandata/pair/aligan/ali_result/T1_1_Aligned.toTranscriptome.out.bam` :指定序列比对结果，需要基于转录本转换过的bam文件 
` --alignments -p 32`   ：设置比对线程数
` ./rsemindex/rsem_ref`  ：索引**位置 **和**名称 **，位置为./rsemindex，名称：rsem_ref
` /home/liuzihao/study/aligen/bam/rsem/T1_1`  ： 指定输出文件与输出文件名前缀                   

```shell
#定量基因表达
rsem-calculate-expression --paired-end --no-bam-output --alignments -p 5 -q ~/study/rawdata/cleandata/pair/aligan/ali_result/T1_1_Aligned.toTranscriptome.out.bam ./rsemindex/rsem_ref ./count/T1_1
#运行成功
rsem-run-em ./rsemindex/rsem_ref 3 ./count/T1_1 ./count/T1_1.temp/T1_1 ./count/T1_1.stat/T1_1 -p 5 -q
Time Used for EM.cpp : 0 h 00 m 36 s
```

一共会生成**2个文件，1个目录**，

```shell
(lzh2) [liuzihao@bogon count]$ ll
total 5348
-rw-r--r--. 1 liuzihao student 2227781 Apr 17 21:12 T1_1.genes.results
-rw-r--r--. 1 liuzihao student 3245298 Apr 17 21:12 T1_1.isoforms.results
drwxr-xr-x. 2 liuzihao student      74 Apr 17 21:12 T1_1.stat
```

**T1_1.genes.results**：基于基因组的定量结果

```shell
gene_id transcript_id(s)        length  effective_length        expected_count  TPM     FPKM
AT1G01010       AT1G01010.1     1688.00 1533.47 3.00    11.68   10.53
AT1G01020       AT1G01020.1,AT1G01020.2,AT1G01020.3,AT1G01020.4,AT1G01020.5,AT1G01020.6 1087.00 932.47  2.00    12.81   11.54
AT1G01030       AT1G01030.1,AT1G01030.2 1870.50 1715.97 0.00    0.00    0.00
AT1G01040       AT1G01040.1,AT1G01040.2 5877.00 5722.47 2.00    2.09    1.88
AT1G01050       AT1G01050.1,AT1G01050.2 905.04  750.52  16.00   127.29  114.75
······
```

**T1_1.isoforms.results**：基于转录本的定向结果（后续以此为准）

```shell
transcript_id   gene_id length  effective_length        expected_count  TPM     FPKM    IsoPct
AT1G01010.1     AT1G01010       1688    1533.47 3.00    11.68   10.53   100.00
AT1G01020.1     AT1G01020       1329    1174.47 0.00    0.00    0.00    0.00
AT1G01020.2     AT1G01020       1087    932.47  2.00    12.81   11.54   100.00
AT1G01020.3     AT1G01020       1420    1265.47 0.00    0.00    0.00    0.00
······
```

**T1_1.stat**： 比对细节存放目录

`gene_id`: 基因ID

`transcript_id(s)`：每个基因对应的转录本id

`length`：基因长度

`effective_length`：有效长度

`expected_count`： 基因的原始counts

`RPKM`：RPKM值(标准化后的counts)

`FPKM`：FPKM值(标准化后的counts)

`IsoPct`: 该转录本表达量占基因总表达量的百分比

`TPM`: TPM值(标准化后的counts)

#### 关于各个丰度的表示方法



#### 批量定量脚本

```shell
#!/bin/bash
ls | grep "toTranscriptome" > ls.txt
mkdir count
for i in $(cat ls.txt)
	do
		name=$(echo "${i}" | cut -d "_" -f 1-2)
		rsem-calculate-expression --paired-end --no-bam-output --alignments -p 16  -q ${i} $1 ./count/${name} 
	done
#创建结果存放目录
cd count
mkdir stat
mkdir count_result
mv ./*stat ./stat
mv ./T* ./count_result
rm -rf ls.txt
```

#### 合并定量结果（以基于转录本定量结果为例子）

- **得到基于COUNTS数的矩阵**（基于基因组的定量结果）

  ```shell
  rsem-generate-data-matrix *.genes.results > expression.matrix
  ```

  **结果查看**

  ```shell
        "T1_1.genes.results"    "T1_2.genes.results"    "T2_1.genes.results"    "T2_2.genes.results"
  "AT1G01010"     3.00    0.00    3.00    4.00
  "AT1G01020"     2.00    1.00    1.00    2.00
  "AT1G01030"     0.00    0.00    0.00    0.00
  "AT1G01040"     2.00    6.00    2.00    10.00
  ······
  ```

- **得到基于FPKM/TPM的矩阵**（基于转录本定量结果）

  ```shell
  #修改rsem-generate-data-matrix
  vim ~/miniconda2/envs/lzh2/bin/rsem-generate-data-matrix
  #修改第11行，5-TPM; 6-FPKM; 7-IsoPct
  my $offsite = 5; 
  rsem-generate-data-matrix *.isoforms.results > expression.matrix2
  #结果查看
    "T1_1.isoforms.results" "T1_2.isoforms.results" "T2_1.isoforms.results" "T2_2.isoforms.results"
  "AT1G01010.1"   11.68   0.00    12.67   15.10
  "AT1G01020.1"   0.00    0.00    0.00    0.00
  "AT1G01020.2"   12.81   0.00    .00    0.00
  "AT1G01020.3"   0.00    0.00    0.00    0.00
  "AT1G01020.4"   0.00    0.00    0.00    0.00
  ······
  ```

## 基因定量（kallisto） 

https://pachterlab.github.io/kallisto/manual

kallisto功能如下，

```
kallisto 0.46.2

Usage: kallisto <CMD> [arguments] ..

Where <CMD> can be one of:

    index         Builds a kallisto index 
    quant         Runs the quantification algorithm 
    bus           Generate BUS files for single-cell data 
    pseudo        Runs the pseudoalignment step 
    merge         Merges several batch runs 
    h5dump        Converts HDF5-formatted results to plaintext
    inspect       Inspects and gives information about an index
    version       Prints version information
    cite          Prints citation information
```

使用**`kallisto`**可以直接对质控后的fastq文件进行定量，粗略

```shell
kallisto index -i kaillisto_index ./rsemindex/rsem_ref.transcripts.fa
```

`-i`: 指定索引名称

`./rsemindex/rsem_ref.transcripts.fa`:  指定基于转录本的序列

`-k`: 指定k-mer大小

```shell
#定量基因表达
kallisto quant -i ali_result/kaillisto_index -o ./K_T1_1 ../T1_1_R2_pair.fastq.gz ../T1_1_R1_pair.fastq.gz
```

`-i`: 指定索引名称

`-o`:输出目录

**结果查看**

```shell
#有3个输出文件
(lzh2) [liuzihao@bogon K_T1_1]$ tree ../K_T1_1/
../K_T1_1/
├── abundance.h5
├── abundance.tsv
└── run_info.json
```

`abundances.h5`是一个HDF5二进制文件，其中包含运行信息，bootstrap estimates, and transcript length information length可以通过sleuth读取此文件
`abundances.tsv`是丰度估计值的纯文本文件。第一行包含每列的标题，包括counts数，TPM，有效长度。
`run_info.json`是一个json文件，其中包含有关运行的信息

```
#kallisto
target_id       length  eff_length      est_counts      tpm
AT1G01010.1     1688    1533.74 3       11.0253
AT1G01020.1     1329    1174.74 0       0
AT1G01020.2     1087    932.744 3       18.1293
AT1G01020.3     1420    1265.74 0       0
AT1G01020.4     1397    1242.74 0       0
······
#RSEM
transcript_id   gene_id length  effective_length        expected_count  TPM     FPKM    IsoPct
AT1G01010.1     AT1G01010       1688    1533.47 3.00    11.68   10.53   100.00
AT1G01020.1     AT1G01020       1329    1174.47 0.00    0.00    0.00    0.00
AT1G01020.2     AT1G01020       1087    932.47  2.00    12.81   11.54   100.00
AT1G01020.3     AT1G01020       1420    1265.47 0.00    0.00    0.00    0.00
AT1G01020.4     AT1G01020       1397    1242.47 0.00    0.00    0.00    0.00
·······
```

## 基因定量（featurecounts）

https://www.jianshu.com/p/9cc4e8657d62

```shell
#基于基因组的定量
featureCounts -p -a ~/study/ref_seq/Araport11_GFF3_genes_transposons.201606.gtf -o counts.txt -B -C -g gene_id T*.sortedByCoord.out.bam
#与rsem结果基本一致
#输出文件
#比对细节
(lzh2) [liuzihao@bogon ali_result]$ head counts.txt.summary
Status	T1_1_Aligned.sortedByCoord.out.bam	T1_2_Aligned.sortedByCoord.out.bam	T2_1_Aligned.sortedByCoord.out.bam	T2_2_Aligned.sortedByCoord.out.bam
Assigned	185546	201269	181832	197073
Unassigned_Unmapped	0	0	0	0
Unassigned_Read_Type	0	0	0	0
Unassigned_Singleton	0	0	0	0
Unassigned_MappingQuality	0	0	0	0
#比对结果，如下图右下角
#与rsem结果基本一致
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210420215443345.png" alt="image-20210420215443345" style="zoom:67%;" />

`-p`: 双端测序

`-a < string >`	参考gtf文件名，支持Gzipped文件格式

`-o < string >`	输出文件的名字，输出文件的内容为read 的统计数目

`-t < string >`	设置feature-type，-t指定的必须是gtf中有的feature，同时read只有落到这些feature上才会被统计到，默认是“exon”

`-g < string >`	当参考的gtf提供的时候，我们需要提供一个id identifier 来将feature水平的统计汇总为meta-feature水平的统计，默认为gene_id，注意！选择gtf中提供的id identifier

`-B`	在-p选择的条件下，只有两端read都比对上的fragment才会被统计

`-C`	如果-C被设置，那融合的fragment（比对到不同染色体上的fragment）就不会被计数，这个只有在-p被设置的条件下使用

`-d < int >`	最短的fragment，默认是50
`-D < int >`	最长的fragmen，默认是600



