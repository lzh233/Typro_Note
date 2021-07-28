#  Cellranger

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

需要**`SRAtool`**、**`Aspera`**、**`cellranger`**、**`bcl2fastq`**、**`tree`**

### 软件安装

```shell
#SRAtool安装v-2.9.2
wget  https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.2/sratoolkit.2.9.2-centos_linux64.tar.gz
#下载后解压直接修改环境变量即可
```

```shell
#Aspera安装v-3.7.4
wget http://download.asperasoft.com/download/sw/connect/3.7.4/aspera-connect-3.7.4.147727-linux-64.tar.gz
tar -zxvf aspera-connect-3.7.4.147727-linux-64.tar.gz
#安装，安装结束后家目录如果有.aspera证明安装成功
./aspera-connect-3.7.4.147727-linux-64.sh
```

```shell
#cellranger安装直接安装最新版的6.0.0，软件官网：https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

wget -O cellranger-6.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.0.0.tar.gz?Expires=1618175827&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2MTgxNzU4Mjd9fX1dfQ__&Signature=BZ1xzJU1NGM-eikz5V59t0JIM2o6m2bXxM9N5Ce3GfuBkxLI2YJZskHeSFIBE9B9xiJpG22yFUBxSQOIAE2ABuyf4BvbKOBm4iRnPvCiePumiOQG9D-ocvCmhIQkFvFVkqR0700OazZSkdw-IyF6ywXNDdKuWBMq~mKYJnVq~nSC-cftxb5jco1wJvYD2PtC6vEYefY7G18XsSYrRif6HTa13cgGvK5MBsT6qachS6FOGafAalxlY5EGj428FruAVk2v-lnJXIDz97iIegyNPR6qfo4Qr9YwYrWw-zYRwmhWpWpW~HwyI-rOqioe6IzTQIjJ1mOAwAT7HT5bNFMmgw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
#解压
tar -zxvf cellranger-6.0.0.tar.gz
#修改环境变量后即可使用
#测试，可以输出系统配置信息
cellranger sitecheck
```

```shell
#bcl2fastq安装
#将下载好的rpm包进行解压
unzip bcl2fastq2-v2-20-0-linux-x86-64.zip
#因为没有root权限所以要直接提取rpm包中的文件
rpm2cpio bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm |cpio -idvm
#提取出来后当前目录会出现usr目录，修改环境变量（将bin目录加入）
/home/liuzihao/study/scrna/usr/local/bin:（加入环境变量）
#测试
(base) [liuzihao@bogon scrna]$ bcl2fastq -v
BCL to FASTQ file converter
bcl2fastq v2.20.0.422
Copyright (c) 2007-2017 Illumina, Inc.
```

```shell
#tree安装
wget ftp://mama.indstate.edu/linux/tree/tree-1.6.0.tgz
tar -zxvf tree-1.6.0.tgz
#修该Makefile，指定安装路径(第19行)
19 prefix = /home/liuzihao/tree
#安装
make -j 4 && make install
```

## 数据下载

### 参考基因组索引下载

```shell
#10x官方已经构建好的人基因组索引
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
```

### 测序数据下载（GSE117988）

```shell
#prefech无法调用ASPERA未解决，所以用了迅雷替代，只下载了一半，另一半实在是太大，网费emmmmm
prefetch -t ascp -a "~/aspera/connet/bin/ascp|/~/aspera/connet/etc/asperaweb_id_dsa.openssh" --option-file ./data/list.txt
#最终得到以下数据6个SRA文件，后续进行拆分与定量
-rw-r--r--. 1 liuzihao student  2.6G Apr 11 10:37 SRR7722937
-rw-r--r--. 1 liuzihao student  3.7G Apr 11 11:12 SRR7722938
-rw-r--r--. 1 liuzihao student  1.8G Apr 11 10:54 SRR7722939
-rw-r--r--. 1 liuzihao student  2.0G Apr 11 11:02 SRR7722940
-rw-r--r--. 1 liuzihao student 1016M Apr 11 10:38 SRR7722941
-rw-r--r--. 1 liuzihao student  1.9G Apr 11 11:04 SRR7722942
```

## SRA文件拆分

### fastq-dump工具

使用**`fastq-dump`**命令拆分所有SRA文件，

- **`--gzip / --bzip2 `**:输出文件压缩，输出fastq.gz

- **`--split-file`**：将双端测序分为两份,放在不同的文件,但是对于一方有而一方没有的reads直接丢弃https://www.biostars.org/p/186741/

- **`--split-spot`**：将双端测序分为两份,放在不同的文件,但是对于一方有而一方没有的reads会单独放在一个文件夹里

- **`--split-3`**：将双端测序分为两份,放在不同的文件,但是对于一方有而一方没有的reads会单独放在一个文件夹里（优先）

- **`--split-file`**：

- **`-A`**：指定输出文件的名字

- **`-I | --readids`**:双端reads加.1和.2区分 不需要

- **`--defline-seq <fmt> `**:指定序列名称的显示格式，如, `--defline-qual '+'`第三行则只显示加号

- **`--defline-qual <fmt>`**：指定质量信息名称的显示格式

  **`<fmt>`**具体格式：?

  - `$ac:序列编号`

  - `$si: spot id`

  - `$sG: spot group`

  - `$sl:碱基长度`

  - `$ri:reads编号，$rn:reads名称，$rl:reads长度`

    ```shell
    #具体格式
    #<fmt> is string of characters and/or 
    #                                   variables. The variables can be one of: $ac 
    #                                   - accession, $si spot id, $sn spot 
    #                                   name, $sg spot group (barcode), $sl spot 
    #                                   length in bases, $ri read number, $rn 
    #                                  read name, $rl read length in bases. '[]' 
    #                                  could be used for an optional output: if 
    #                                   all vars in [] yield empty values whole 
    #                                   group is not printed. Empty value is empty 
    #                                   string or for numeric variables. Ex: 
    #                                   @$sn[_$rn]/$ri '_$rn' is omitted if name 
    #                                  is empty
    #举例子
    fastq-dump --gzip --split-files --defline-qual '+' --defline-seq '>\$ac-\$si/\$ri' ./SRR7722941 
    #查看
    zless -SN SRR7722941_1.fastq.gz | head
    
    >\SRR7722941-\1/\1
    CTAGGTGA
    +
    GGGGGIII
    >\SRR7722941-\2/\1
    CTAGGTGA
    +
    GGGGGIII
    >\SRR7722941-\3/\1
    CTAGGTGA
    ```

### 批量拆分脚本

**1小时左右**，脚本名：fastq.sh

```shell
#!/bin/bash
######################################################
#Author:liuzihao
#Date:2021-4-11
#CentOS Linux release 7.8.2003 (Core)	
#Architecture:          x86_64
#CPU op-mode(s):        32-bit, 64-bit
#Byte Order:            Little Endian
#CPU(s):                144
#Model name:            Intel(R) Xeon(R) Gold 5220 CPU @ 2.20GHz
######################################################
function file_num(){
	count=0
	total=$(wc -l ls.txt | cut -d " " -f 1)
}

ls S* > ls.txt
file_num
for i in $(cat ls.txt)
	do
		count=$(( $count + 1 ))
		echo -ne "共${total}个文件，正在拆分第${count}个\r"
		fastq-dump --gzip --split-files -A ${i} ${i} &>> log2.txt && echo "${i}_Success" >> ./log.txt || echo "${i}_Fail" >> ./log.txt
	done
rm -rf ls.txt
#创建fastq文件存放位置
mkdir fastqdata
mv ./*.gz ./fastqdata

```

### 拆分结果查看

```shell
(base) [liuzihao@bogon fastqdata]$ la | head -n 6
total 19G
-rw-r--r--. 1 liuzihao student 417M Apr 11 17:07 SRR7722937_1.fastq.gz
-rw-r--r--. 1 liuzihao student 889M Apr 11 17:07 SRR7722937_2.fastq.gz
-rw-r--r--. 1 liuzihao student 2.4G Apr 11 17:07 SRR7722937_3.fastq.gz
```

每个`SRA`文件被拆分成3个fastq文件，通过**`zless -NS <FILE>`**查看文件内容,`zless`与`less`用法一致

```shell
#<file>_1文件中每个序列均只有8bp可能为barcode或UMI或index（详细见后）
(base) [liuzihao@bogon fastqdata]$ zless -SN SRR7722937_1.fastq.gz | head -n 6
@SRR7722937.1 SN367:911:HKMNCBCXY:1:1103:1119:1866 length=8
TTTCATGA
+SRR7722937.1 SN367:911:HKMNCBCXY:1:1103:1119:1866 length=8
GGGGGIII
@SRR7722937.2 SN367:911:HKMNCBCXY:1:1103:1182:1935 length=8
TTTCATGA

#<file>_2文件中每个序列均只有26bp可能为barcode或UMI或index（详细见后）
(base) [liuzihao@bogon fastqdata]$ zless -SN SRR7722937_2.fastq.gz | head -n 6
@SRR7722937.1 SN367:911:HKMNCBCXY:1:1103:1119:1866 length=26
AGCAGCCGTGACTACTGTATTGCGAC
+SRR7722937.1 SN367:911:HKMNCBCXY:1:1103:1119:1866 length=26
AGGGGIIIIIGGIIIIIIIIIIIIGG
@SRR7722937.2 SN367:911:HKMNCBCXY:1:1103:1182:1935 length=26
GCTCCTAAGACACTAAGGCCTGTACC

#<file>_3文件中每个序列均只有98bp，应该为测序reads（详细见后）
(base) [liuzihao@bogon fastqdata]$ zless -SN SRR7722937_3.fastq.gz | head -n 6
@SRR7722937.1 SN367:911:HKMNCBCXY:1:1103:1119:1866 length=98
NNNNNCTGTAATCCCAGCCAGGAGGACTGCTTGAACCCGGGAGGCAGAGGTTTCAGTGAGCTGAGTGCCACTGCACTCCAGCCTGGGTGACAGAGTGA
+SRR7722937.1 SN367:911:HKMNCBCXY:1:1103:1119:1866 length=98
#####<.<G<GGAAGGA.AGGGGGGIIIGAAGGGAAGG<AAAAAGGGAAA<GGGGGGGAAGGGAGGGA<GGGGIGGGAAAGGGGGGGGGGAGG..AAA
@SRR7722937.2 SN367:911:HKMNCBCXY:1:1103:1182:1935 length=98
GCGAAAAGTCCTAATAGTAGAAGAACCCTCCATAAACCTGGAGTGACTATATGGATGCCCGCCACCCTACCACACATTCGAAGAGCCCGTATACATAA
```

### 10x文库特征

**（文库特征与白洋的土壤细菌宏培养组设计相似，index区分96孔板，barcode区分单菌16s）**

<img src="G:\Desktop\s_note\data\picture\image-20210412092554946.png" alt="单细胞学习-cellranger流程" style="zoom:110%;" />

**`P5/P7接头`**:与Flowcell上的p5p7接头互补，使reads可以固定于Flowcell上进行桥式PCR扩增

**`index`**：用于区分不同样本（文库）的序列，96孔板的每个孔中都加入了**4种不同的index oligos混合**，后续拆分序列可根据**index**进行区分，后续将具有相同4种oligo的fq文件组合在一起，表示来自同一个样本（文库）。

**`barcode`**：用于区分同样本（文库）中的cell

**`UMI`**：UMI就是**分子标签（Unique Molecular Identifier）**，由4-10个随机核苷酸组成，在mRNA反转录后，进入到文库中，每一个mRNA随机连上一个UMI，结果可以计数不同的UMI，最终统计mRNA的数量，

```shell
#UMI原理: 就是给每一条原始DNA片段**加上一段**特有的**标签序列，经文库构建及PCR扩增后一起进行测序。这样，根据**不同的标签序列**我们就可以区分**不同来源
#的DNA模板，分辨哪些是PCR扩增及测序过程中的随机错误造成的假阳性突变。
```

**`Poly(dT)VN`**: 随机引物，mRNA可以与其结合序列：TTTTTTTTT

![image-20210412103453346](G:\Desktop\s_note\data\picture\image-20210412103453346.png)

- 首先，**1-26个cycle就是测序得到了26个碱基，先是16个Barcode碱基，然后是10个UMI碱基；**
- 然后，**27-34这8个cycle得到了8个碱基，就是i7的sample index；**
- 最后，**35-132个cycle得到了98个碱基，就是转录本reads**

```shell
#根据reads长度，上述所拆分出的3个fastq文件分别为，index、barcode、reads
#<file>_1: 8bp>>>index
#<file>_2: 26bp>>>barcode(16bp)+UMI(10bp)
#<file>_3: 98bp >>>reads
```

### fastq重命名

**根据10x官方说明，需要对文件名进行修改**

![image-20210412104846914](G:\Desktop\s_note\data\picture\image-20210412104846914.png)

- `I1`：**index序列**
- `R1`：**barcode序列**
- `R2`：**reads序列**

如，

`SRR7722937_S1_L001_R1_001.fastq.gz `

`SRR7722937_S1_L001_R2_001.fastq.gz `

`SRR7722937_S1_L001_I1_001.fastq.gz`

**使用以下脚本批量重命名：**(rename.sh)

```shell
#!/bin/bash
######################################################
#Author:liuzihao
#Date:2021-4-12
#CentOS Linux release 7.8.2003 (Core)	
#Architecture:          x86_64
#CPU op-mode(s):        32-bit, 64-bit
#Byte Order:            Little Endian
#CPU(s):                144
#Model name:            Intel(R) Xeon(R) Gold 5220 CPU @ 2.20GHz
######################################################
#获取各个序列类型的文件列表
echo "或取文件列表"
ls -l | awk  '{print $9}' | fgrep "_2" | cut -d "_" -f 1 > ls2.txt
ls -l | awk  '{print $9}' | fgrep "_1" | cut -d "_" -f 1 > ls1.txt
ls -l | awk  '{print $9}' | fgrep "_3" | cut -d "_" -f 1 > ls3.txt
echo "开始重命名...."
#重名名
for a in $(cat ls1.txt)
	do
		mv ${a}_1.fastq.gz ${a}_S1_L001_I1_001.fastq.gz
	done

for b in $(cat ls2.txt)
	do
		mv ${b}_2.fastq.gz ${b}_S1_L001_R1_001.fastq.gz
	done

for c in $(cat ls3.txt)
	do		
		mv ${c}_3.fastq.gz ${c}_S1_L001_R2_001.fastq.gz
	done
echo "命名结束，删除临时文件...."
rm -rf ls*.txt

#重命名结果如下
(base) [liuzihao@bogon fastqdata]$ la S* | head -n 3
-rw-r--r--. 1 liuzihao student 417M Apr 11 17:07 SRR7722937_S1_L001_I1_001.fastq.gz
-rw-r--r--. 1 liuzihao student 889M Apr 11 17:07 SRR7722937_S1_L001_R1_001.fastq.gz
-rw-r--r--. 1 liuzihao student 2.4G Apr 11 17:07 SRR7722937_S1_L001_R2_001.fastq.gz
```

## fastqc查看数据质量情况

```shell
#使用fastqc对R1和R2碱基质量查看，主要为R2以为其为真实reads
fastqc ./*R2* -o ./fastqc/ -t 64
#multicqc整合fastqc结果
multicqc ./fastqc/
```

### 结果查看

查看**`muticqc_report.html`**

## cellranger使用

### 软件自检

```shell
cellranger testrun --id=tiny
#以下信息自检成功
Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!
```

### mkfastq

```shell
#官方资料
#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq#example_data
```

#### mkfastq基本原理与特点

`mkfastq`可以将混合测序的样品通过`index`信息将Illumina 测序仪的下机数据转换为fastq文件，本功能类似于Illumina 的`bcl2fastq`，特性如下，

- 将10X 样本index名称与四种寡核苷酸对应起来，比如A1孔是样本`SI-GA-A1`，然后对应的寡核苷酸是`GGTTTACT, CTAAACGG, TCGGCGTC, and AACCGTAA` ，那么程序就会去index文件中将存在这四种寡核苷酸的fastq组合到A1这个样本
- 提供质控结果，包括barcode 质量、总体测序质量如Q30、R1和R2的Q30碱基占比、测序reads数等
- 可以使用10X简化版的样本信息表

**不同文库在同一条Lane测序**

![image-20210412153942463](G:\Desktop\s_note\data\picture\image-20210412153942463.png)

**同一文库在不同Lane测序**

![image-20210412154140496](G:\Desktop\s_note\data\picture\image-20210412154140496.png)

#### mkfastq命令选项

**`--run`**：BCL文件的位置（路径）

**`--id`**：mkfastq创建的输出目录的名称

**`--samplesheet`**：Illumina测序数据文件的相关信息文件的位置，包括样品名、index信息、lane信息，**`--sample-sheet`**与上述一致

**`--csv`**：仅包括样品名、index信息、lane等信息的简单版samplesheet，**`--simple-csv`**与上述一致

**`--filter-dual-index`**: 仅拆分i5和i7确定的reads

**`--filter-single-index`**：仅拆分i7确定的reads

**`--lanes`**：指定拆分哪些lanse的reads

**`--output-dir`**：设定输出目录，代替默认的`flowcell_id/outs/fastq_path`

**`--localcores`**：限定使用cpu数量

**`--localmem`**：限定使用内存数量，GB

......

![image-20210413092645774](G:\Desktop\s_note\data\picture\image-20210413092645774.png)

#### mkfastq用法

**下载官网提供的3个示例文件，**

```shell
#BCL文件
wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
#samplesheet，10x简化版
wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv
#samplesheet，Illumina版
wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-samplesheet-1.2.0.csv
```

**samplesheet，Illumina版**

![image-20210412161415412](G:\Desktop\s_note\data\picture\image-20210412161415412.png)

**10x简化版**

![image-20210412161531775](G:\Desktop\s_note\data\picture\image-20210412161531775.png)

```shell
#使用简化版samplesheet对bcl进行拆分
cellranger mkfastq --id=tin-bcl --run=./cellranger-tiny-bcl-1.2.0 --csv=./cellranger-tiny-bcl-simple-1.2.0.csv
#运行结束后会输出各个文件的存放位置
Outputs:
- FASTQ output folder:   /home/liuzihao/study/scrna/exbcl/tin-bcl/outs/fastq_path
- Interop output folder: /home/liuzihao/study/scrna/exbcl/tin-bcl/outs/interop_path
- Input samplesheet:     /home/liuzihao/study/scrna/exbcl/tin-bcl/outs/input_samplesheet.csv

(base) [liuzihao@bogon fastq_path]$ la *.gz
-rw-r--r--. 1 liuzihao student  20M Apr 13 09:17 Undetermined_S0_L001_I1_001.fastq.gz
-rw-r--r--. 1 liuzihao student  50M Apr 13 09:17 Undetermined_S0_L001_R1_001.fastq.gz
-rw-r--r--. 1 liuzihao student 146M Apr 13 09:17 Undetermined_S0_L001_R2_001.fastq.gz
```

**Illumina版**

```shell
cellranger mkfastq --id=illumia --run=./cellranger-tiny-bcl-1.2.0 --samplesheet=cellranger-tiny-bcl-samplesheet-1.2.0.csv
#运行结束后会输出各个文件的存放位置
Outputs:
- FASTQ output folder:   /home/liuzihao/study/scrna/exbcl/illumia/outs/fastq_path
- Interop output folder: /home/liuzihao/study/scrna/exbcl/illumia/outs/interop_path
- Input samplesheet:     /home/liuzihao/study/scrna/exbcl/illumia/outs/input_samplesheet.csv

(base) [liuzihao@bogon fastq_path]$ la *.gz
-rw-r--r--. 1 liuzihao student  20M Apr 13 09:17 Undetermined_S0_L001_I1_001.fastq.gz
-rw-r--r--. 1 liuzihao student  50M Apr 13 09:17 Undetermined_S0_L001_R1_001.fastq.gz
-rw-r--r--. 1 liuzihao student 146M Apr 13 09:17 Undetermined_S0_L001_R2_001.fastq.gz
```

### counts

使用**`cellranger count`**进行定量（整个定量过程包含了比对（STAR）、质控、定量）

#### counts选项参数

https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count

`--id`：指定输出文件目录的名称

`--transcriptome`：指定参考基因组索引位置

`--fastqs`：指定fastq文件位置

`--expect-cells`：预期得到cell数目，默认3000

`--sample`：指定样品名称，fastq文件的前缀中的sample保持一致，作为软件识别的标志，不指定则分析所有fastq文件

`--nosecondary`：不进行下游降维聚类等分析

`--r1-length / --r2-length`：将输入的R2/R1序列截短到指定长度

`--localcores / --localmem`：指定处理器数目和内存使用量

#### counts用法

```shell
cellranger count --id=count --transcriptome=/home/liuzihao/study/scrna/ref/refdata-gex-GRCh38-2020-A --fastqs=/home/liuzihao/study/scrna/exbcl/test --expect-cells=3000 --nosecondary

##输出文件说明与位置
Outputs:
- Run summary HTML:                         /home/liuzihao/study/scrna/exbcl/test/count/outs/web_summary.html
- Run summary CSV:                          /home/liuzihao/study/scrna/exbcl/test/count/outs/metrics_summary.csv
- BAM:                                      /home/liuzihao/study/scrna/exbcl/test/count/outs/possorted_genome_bam.bam
- BAM index:                                /home/liuzihao/study/scrna/exbcl/test/count/outs/possorted_genome_bam.bam.bai
- Filtered feature-barcode matrices MEX:    /home/liuzihao/study/scrna/exbcl/test/count/outs/filtered_feature_bc_matrix
- Filtered feature-barcode matrices HDF5:   /home/liuzihao/study/scrna/exbcl/test/count/outs/filtered_feature_bc_matrix.h5
- Unfiltered feature-barcode matrices MEX:  /home/liuzihao/study/scrna/exbcl/test/count/outs/raw_feature_bc_matrix
- Unfiltered feature-barcode matrices HDF5: /home/liuzihao/study/scrna/exbcl/test/count/outs/raw_feature_bc_matrix.h5
- Secondary analysis output CSV:            /home/liuzihao/study/scrna/exbcl/test/count/outs/analysis
- Per-molecule read information:            /home/liuzihao/study/scrna/exbcl/test/count/outs/molecule_info.h5
- CRISPR-specific analysis:                 null
- Loupe Browser file:                       /home/liuzihao/study/scrna/exbcl/test/count/outs/cloupe.cloupe
- Feature Reference:                        null
- Target Panel File:                        null
Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!

#最重要的三个文件
(base) [liuzihao@bogon filtered_feature_bc_matrix]$ la
total 2.8M
-rw-r--r--. 1 liuzihao student  13K Apr 13 14:24 barcodes.tsv.gz
-rw-r--r--. 1 liuzihao student 326K Apr 13 14:24 features.tsv.gz
-rw-r--r--. 1 liuzihao student 2.5M Apr 13 14:24 matrix.mtx.gz
```

- `web_summary.html`：官方说明 summary HTML file 
- `metrics_summary.csv`：CSV格式数据摘要
- `possorted_genome_bam.bam`：比对文件
- `possorted_genome_bam.bam.bai`：索引文件
- **`filtered_gene_bc_matrices`：*是重要的一个目录，下面又包含了 barcodes.tsv.gz、features.tsv.gz、matrix.mtx.gz，是下游Seurat、Scater、Monocle等分析的输入文件**
- `filtered_feature_bc_matrix.h5`：过滤掉的barcode信息HDF5 format
- `raw_feature_bc_matrix`：原始barcode信息
- `raw_feature_bc_matrix.h5`：原始barcode信息HDF5 format
- `analysis`：数据分析目录，下面又包含聚类clustering（有graph-based & k-means）、差异分析diffexp、主成分线性降维分析pca、非线性降维tsne
- `molecule_info.h5`：**下面进行aggregate使用的文件**
- `cloupe.cloupe`：官方可视化工具Loupe Cell Browser 输入文件

### mkgtf和mkref

使用**`cellranger mkgtf`**和**`cellranger mkref`**筛选GTF文件（构建新的gtf文件）以及构建参考基因组索引，也可以使用`awk`等命令自行提取

- 参考序列只能有很少的 overlapping gene annotations，因为reads比对到多个基因会导致流程检测的分子数更少(它只要uniquely mapped的结果)

- FASTA与GTF比对和STAR兼容，GTF文件的第三列（feature type）必须有exon

- 可以自己增加基因注释信息，按照gtf的格式添加，gtf具体信息参考`RNA-seq`分析

  eg：mylocus   annotation   exon    100    200   .    +   .   gene_id "mygene"; transcript_id "mygene";

#### mkgtf

```shell
#参考基因组下载
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#参考基因组注释下载
wget ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz
```

使用**`cellranger mkgtf`**对GTF进行筛选：`cellranger mkgtf <input_gtf> <output_gtf> [--attribute=KEY:VALUE...]`

```shell
#!/bin/bash
cellrangerr mkgtf Homo_sapiens.GRCh38.ensembl.gtf Homo_sapiens.GRCh38.ensembl.filtered.gtf \
                   --attribute=gene_biotype:protein_coding \
                   --attribute=gene_biotype:lincRNA \
                   --attribute=gene_biotype:antisense \
                   --attribute=gene_biotype:IG_LV_gene \
                   --attribute=gene_biotype:IG_V_gene \
                   --attribute=gene_biotype:IG_V_pseudogene \
                   --attribute=gene_biotype:IG_D_gene \
                   --attribute=gene_biotype:IG_J_gene \
                   --attribute=gene_biotype:IG_J_pseudogene \
                   --attribute=gene_biotype:IG_C_gene \
                   --attribute=gene_biotype:IG_C_pseudogene \
                   --attribute=gene_biotype:TR_V_gene \
                   --attribute=gene_biotype:TR_V_pseudogene \
                   --attribute=gene_biotype:TR_D_gene \
                   --attribute=gene_biotype:TR_J_gene \
                   --attribute=gene_biotype:TR_J_pseudogene \
                   --attribute=gene_biotype:TR_C_gene
```

筛选结果查看

```shell
#假定只筛选TR_C_gene
cellranger mkgtf Homo_sapiens.GRCh38.84.gtf filtered.gtf --attribute=gene_biotype:TR_C_gene
#查看筛选后的GTF
cat filtered.gtf |grep -v "#" |awk -v FS='gene_biotype ' 'NF>1{print $2}'|awk -F ";" '{print $1}'|sort | uniq -c
#只有1种
125 "TR_C_gene"

#查看原gtf文件
cat Homo_sapiens.GRCh38.84.gtf |grep -v "#" |awk -v FS='gene_biotype ' 'NF>1{print $2}'|awk -F ";" '{print $1}'|sort | uniq -c
    128 "3prime_overlapping_ncrna"
  45662 "antisense"
     24 "bidirectional_promoter_lncrna"
    213 "IG_C_gene"
     33 "IG_C_pseudogene"
    152 "IG_D_gene"
     76 "IG_J_gene"
......
```

#### mkref

使用**`cellranger mkref`**构建参考基因组索引

`--genome`: 输出的索引目录名称

`--fasta`：参考基因组序列

`--genes`：参考基因组GTF文件

`--nthreads`：线程数

`--memgb`：内存限制数量

`--ref-version`：Reference version string to include with reference.？

```shell
# 构建参考基因组index
cellranger mkref --genome=testout --fasta=Homo_sapiens.GRCh38.dna.primary_assembly.fa --genes=Homo_sapiens.GRCh38.ensembl.filtered.gtf --nthreads=32
#以下提示构建成功
>>> Reference successfully created! <<<

#构建的索引目录结构
testout/
├── fasta
│   ├── genome.fa
│   └── genome.fa.fai
├── genes
│   └── genes.gtf.gz
├── reference.json
└── star
    ├── chrLength.txt
    ├── chrNameLength.txt
    ├── chrName.txt
    ├── chrStart.txt
    ├── exonGeTrInfo.tab
    ├── exonInfo.tab
    ├── geneInfo.tab
    ├── Genome
    ├── genomeParameters.txt
    ├── SA
    ├── SAindex
    ├── sjdbInfo.txt
    ├── sjdbList.fromGTF.out.tab
    ├── sjdbList.out.tab
    └── transcriptInfo.tab
#自建索引测试
cellranger count --id=count37 --fastqs=./srr37/ --transcriptome=/home/liuzihao/study/scrna/refseq/testout --expect-cells=3000 --nosecondary &> log2.txt &
#成功
Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!

2021-04-13 23:20:55 Shutting down.
```

**多物种整合**

```shell
#没有gtf可以自己根据NCBI上的序列信息自己造一个https://mp.weixin.qq.com/s/_VGFGmYBJmYm_4KLc9zamg
cellranger mkref --genome=hg19 --fasta=hg19.fa --genes=hg19-filtered-ensembl.gtf \
                   --genome=mm10 --fasta=mm10.fa --genes=mm10-filtered-ensembl.gtf
```

### **aggr**

#### aggr选项参数

`--id`：输出结果目录名 

`--csv`：样品信息文件，详细见下

当处理多个生物学样本或者一个样本存在**多个重复/文库时**，最好的操作就是先分别对每个文库进行单独的count定量，然后将定量结果利用`cellranger aggr`组合起来

```shell
#分别进行3个定量的流程（下载的SRA数据，37 38 39 分别跑）
cellranger count --id=count37 --fastqs=./srr37/ --transcriptome=/home/liuzihao/study/scrna/ref/refdata-gex-GRCh38-2020-A --expect-cells=3000 --nosecondary &>log.txt &
cellranger count --id=count38 --fastqs=./srr38/ --transcriptome=/home/liuzihao/study/scrna/ref/refdata-gex-GRCh38-2020-A --expect-cells=3000 --nosecondary &>log2.txt &
cellranger count --id=count39 --fastqs=./srr39/ --transcriptome=/home/liuzihao/study/scrna/ref/refdata-gex-GRCh38-2020-A --expect-cells=3000 --nosecondary &>log3.txt &
```

构建分组信息，**csv文件**，包括`sample_id`和`molecule_h5`文件的位置

```shell
(base) [liuzihao@bogon fastqdata]$ cat agg.csv 
sample_id,molecule_h5
count37,/home/liuzihao/study/scrna/rawdata/fastqdata/count37/outs/molecule_info.h5
count38,/home/liuzihao/study/scrna/rawdata/fastqdata/count38/outs/molecule_info.h5
count39,/home/liuzihao/study/scrna/rawdata/fastqdata/count39/outs/molecule_info.h5
```

合并文件

```shell
cellranger aggr --id=counts --csv=agg.csv

#提示以下信息则代表合并成功
Outputs:
- Copy of the input aggregation CSV: /home/liuzihao/study/scrna/rawdata/fastqdata/count/outs/aggregation.csv
- Aggregation metrics summary HTML:  /home/liuzihao/study/scrna/rawdata/fastqdata/count/outs/web_summary.html
- count:
    Aggregation metrics summary JSON:         /home/liuzihao/study/scrna/rawdata/fastqdata/count/outs/count/summary.json
    Secondary analysis output CSV:            /home/liuzihao/study/scrna/rawdata/fastqdata/count/outs/count/analysis
    Filtered feature-barcode matrices MEX:    /home/liuzihao/study/scrna/rawdata/fastqdata/count/outs/count/filtered_feature_bc_matrix
    Filtered feature-barcode matrices HDF5:   /home/liuzihao/study/scrna/rawdata/fastqdata/count/outs/count/filtered_feature_bc_matrix.h5
    Unfiltered feature-barcode matrices MEX:  /home/liuzihao/study/scrna/rawdata/fastqdata/count/outs/count/raw_feature_bc_matrix
    Unfiltered feature-barcode matrices HDF5: /home/liuzihao/study/scrna/rawdata/fastqdata/count/outs/count/raw_feature_bc_matrix.h5
    Loupe Browser file:                       /home/liuzihao/study/scrna/rawdata/fastqdata/count/outs/count/cloupe.cloupe
- vdj_t:                             null
- vdj_b:                             null
- V(D)J reference:                   null

Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!
```

