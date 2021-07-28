#  CeleScope

## 运行环境

```
Author:wangwe&liuzihao
Date:2021-4-9
CentOS Linux release 7.8.2003 (Core)	
Architecture:          x86_64
CPU op-mode(s):        32-bit, 64-bit
Byte Order:            Little Endian
CPU(s):                144
Model name:            Intel(R) Xeon(R) Gold 5220 CPU @ 2.20GHz
```

## 软件安装与数据下载

### 软件安装

```shell
git clone https://gitee.com/singleron-rd/celescope.git
cd CeleScope
#创建环境
conda create -n lzh
conda activate lzh
#安装依赖
conda install --file conda_pkgs.txt --channel conda-forge --channel bioconda --channel r --channel imperial-college-research-computing
#安装celescope
pip install -i https://pypi.mirrors.ustc.edu.cn/simple/ celescope
```

### 参考基因组与注释文件下载

```shell
#人参考基因组下载
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
#小鼠参考基因组下载
wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
```

### 参考基因组索引构建

使用`celescope rna mkref`进行参考基因组构建

`--genome_name` 基因组名称

`--fasta --gtf`: 指定fasta和gtf文件

`--genomeDir` ：指定输出目录，默认当前目录

`--thread` ：线程数，默认6

`--dry_run`：仅仅输出config文件

`--mt_gene_list`: 线粒体基因列表文件。它是每行一个基因的纯文本文件。如果没有提供，将使用`mt-`和`mt -`确定线粒体基因。

```shell
#构建小鼠的参考基因组索引
celescope rna mkref --genome_name Mus_musculus_ensembl_99 --fasta ./Mus_musculus.GRCm38.dna.primary_assembly.fa --gtf Mus_musculus.GRCm38.99.gtf &>> log.txt &
#构建人类参考基因组序列
celescope rna mkref --genome_name Homo_sapiens_ensembl_99 --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa --gtf Homo_sapiens.GRCh38.99.gtf &>> log.txt &
```

## 单细胞转录组流程

使用`multi_rna`构建分析脚本

`--mapfile` : 指定mappfile的位置

`--genomeDir`: 指定参考基因组位置

`--thread 16 --mod shell`: 线程数与输出脚本格式

```shell
#构建分析脚本
multi_rna --mapfile ../mapp.file --genomeDir ../ref/human/hg38_index --thread 16 --mod shell
```

构建mappfile

```shell
$ ls ../../test_data_rna/
test1_1.fq.gz  test1_2.fq.gz  test2_1.fq.gz  test2_2.fq.gz
#第一列为文件名称前缀，第二列文件目录，第三列分组信息用_1和_2区分两个方向
$ cat ../../mapp.file 
test1	/Personal/liuzihao/data/test_data_rna	sample1
test2	/Personal/liuzihao/data/test_data_rna	sample1
```

分析

```shell
$ chmod 755 sample1.sh 
$ ./sample1.sh 
```





