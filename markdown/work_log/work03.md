# work03-9.3

## 描述

增加snp的测试数据到github仓库

```apl
仓库：https://github.com/singleron-RD/celescope_test_script
数据：/SGRNJ03/randd/RD20051303_Panel/20210901_2/D0817_PHA6567_SR_11LGR_TS

https://github.com/singleron-RD/CeleScope/blob/master/scripts/extract_read.py  
提取出5个有较多突变（尽量包含T790M和L858R突变）的细胞的fq1和fq2，
用这5个细胞的fq1和fq2创建测试数据；match_dir可以参考tag和vdj的match_dir
```

## 解决流程

### step1. 首先进行具有较多突变和包含T790M和L858R突变的barcode的提取

使用以下脚本进行提取

输入文件

- snp分析结果中的`<sample_name>__count_tsne.tsv`，从这个文件中得到CID和barcode的对应关系以及每个barcode的突变数量

```python
	CID	alt_count	VID	barcode	cluster	tSNE_1	tSNE_2	Gene_Counts
0	1	2.0	[189, 194]	AAACATCGAAACATCGCCTCTATC	1	-13.115079104195699	-16.3938413442029	2227
3	4	2.0	[194, 225]	AAACATCGACGTATCAAGAGTCAA	1	-17.2405509512317	-9.58486363185642	2454
4	5	3.0	[194, 195, 214]	AAACATCGATGCCTAAAGTCACTA	4	-6.312338659492969	15.6323557309553	1384
5	6	3.0	[194, 195, 214]	AAACATCGCAAGGAGCCATCAAGT	5	5.23383980122655	26.6677560190516	311
```

- snp分析结果中的`<sample_name>_variant_table.tsv`，此文件得到突变的类型

```shell
Chrom	Pos	Alleles	VID	CID	Gene	mRNA	Protein	COSMIC	nCell
1	114704958	C-T	1	731	NRAS	NM_002524:*3136G>A			1
1	114707701	G-C	2	152	NRAS	NM_002524:*393C>G			1
1	114713969	T-G	3	437	NRAS	exon3:A121C	R41R		1
1	114713974	G-A	4	384	NRAS	exon3:C116T	S39F		1
1	114716052	C-G	5	('385', '470', '497', '681', '768', '859')	NRAS	exon2:G109C	E37Q		6
1	114716056	G-A	6	679	NRAS	exon2:C105T	T35T		1
.......
```

脚本，

**位置：/Personal/liuzihao/work/script/extract_snp_barcode/**

```python
import pandas as pd
import argparse

class Extract_Data:
    def __init__(self,count_file,variant_table):
        #<sample_name>__count_tsne.tsv
        count = pd.read_table(count_file)
        
        #<sample_name>_variant_table.tsv
        variant = pd.read_table(variant_table).fillna("no_pro")
        #barcode和cid的对应关系
        self.cid2barcode = count.set_index(['CID'])['barcode'].to_dict()
        #每个barcode的突变数量
        self.barcode_vcount = count.set_index(['barcode'])['alt_count'].to_dict()
        
        #筛选包含T790M和L858R的突变的细胞(CID)
        gene_T790M = variant[variant['Protein'].str.contains("T790M")].loc[:,"CID"]
        gene_L858R = variant[variant['Protein'].str.contains("L858R")].loc[:,"CID"]
        
        #得到cid列表，删除cid的括号 并转换为cid列表
        gene_T790M_cid = set(list(gene_T790M)[0].strip("(").strip(")").split(","))
        gene_L858R_cid = set(list(gene_L858R)[0].strip("(").strip(")").split(","))
        
        #共有T790M和L858R的barcode(cid)
        self.cid = gene_T790M_cid.intersection(gene_L858R_cid)
       
    #计算共有T790M和L858R的barcode的各个barcode的vid数
        self.barcode = {}
        for i in self.cid:
            #cid转换为barcode，cid有两层引号
            key = self.cid2barcode[eval(eval(i))]
            #得到barcao的突变数量
            value = self.barcode_vcount[key]
            self.barcode[key] = value

    def get_brcode_test(self):
        #将count数排序
        count_lst = sorted(list(self.barcode.values()),reverse = True)
        #挑选计数位点最多的barcode前五名
        barcode_test = [bar for bar,value in self.barcode.items() if value in count_lst[0:5]]
        if len(barcode_test) > 5:
            barcode_test = barcode_test[0:5]
        
        with open("./test_data/barcode_file.txt","w") as fd:
            for line in barcode_test:
                fd.write(f"{line}\n")
        #return barcode_test

test = Extract_Data(count_file="D0817_PHA6567_SR_11LGR_TS_count_tsne.tsv",
                    variant_table="D0817_PHA6567_SR_11LGR_TS_variant_table.tsv")
test.get_brcode_test()

"""
def main():
    parser = argparse.ArgumentParser(description='extract barcode to make test data')
    parser.add_argument("--count_file", help="count_data", required=True)
    parser.add_argument("--variant_table", help="variant_table", required=True)
    #parser.add_argument("--cell_number",help = "number of cells wants to get")

    args = parser.parse_args()

    test = Extract_Data(count_file=str(args.count_file),variant_table=str(args.variant_table))
    
    test.get_brcode_test()

if __name__ == "__mian__":
    main()

"""

```

输出文件，

```shell
$ cat barcode_file.txt 
AAGGTACACCTAATCCACCACTGT
AGAGTCAACCAGTTCACATCAAGT
ACCTCCAAAACAACCAACGCTCGA
TCTTCACAACAAGCTAACACAGAA
ACGTATCACCTCTATCACAGCAGA
```

### step2.  使用extract_read.py提取数据

```shell
  --barcode_file BARCODE_FILE
                        file with barcode to extract (default: None)
  --match_dir MATCH_DIR
                        Match celescope scRNA-Seq directory. (default: None)
  --R1_read R1_READ     R1 read path. (default: None)
  --outdir OUTDIR       Output directory. (default: ./)
```

输入文件

`--barcode_file`：提供的需要提取的`barcode file`

`--match_dir`: 原始数据结果的的`01.barcode`结果**(snp的结果),** 相当于提供R2

```shell
$ tree match_dir/
match_dir/
└── 01.barcode
    ├── D0817_PHA6567_SR_11LGR_TS_2.fq
    └── stat.txt
```

`--R1_read`: R1的数据，通过mappfile可以 找到(mappfile位置如下，第二列是位置 然后根据样品名找到对应的R1)

![image-20210913171640672](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210913171640672.png)

提取出的数据为`test1_1.fq` `test1_2.fq`

### step3. 构建测试数据

- 构建fastq, 上一步已经构建
- 构建match_dir, match_dir的位置在原始数据的mapfile中

![image-20210913171901751](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210913171901751.png)

```shell
$ ll 
total 7272
drwxr-xr-x. 2 zhuwenqi ssh.randd      26 Aug 31 16:13 00.sample
drwxr-xr-x. 2 zhuwenqi ssh.randd      74 Aug 31 17:03 01.barcode
drwxr-xr-x. 2 zhuwenqi ssh.randd     110 Aug 31 17:35 02.cutadapt
drwxr-xr-x. 2 zhuwenqi ssh.randd     511 Aug 31 19:12 03.star
drwxr-xr-x. 2 zhuwenqi ssh.randd     270 Aug 31 19:44 04.featureCounts
drwxr-xr-x. 4 zhuwenqi ssh.randd     306 Aug 31 20:00 05.count
drwxr-xr-x. 3 zhuwenqi ssh.randd     241 Aug 31 20:36 06.analysis
-rw-r--r--. 1 zhuwenqi ssh.randd 6142757 Aug 31 20:38 D0817_PHA6567_SR_11LGR_ZL_report.html
```

构建snp测试数据的match_dir如下

```shell
$ tree snp_match_dir/
snp_match_dir/
├── 05.count
│   └── snp_matrix_10X
│       └── barcodes.tsv
└── 06.analysis
    ├── snp_markers.tsv
    └── snp_tsne_coord.tsv
```

- 构建mappfile

```shell
$ cat snp.mapfile 
snp	../celescope_test_data/fastqs	snp	../celescope_test_data/snp_match_dir
```

- 构建`ANNOVA`软件配置文件

```shell
$ cat annovar.config 
[ANNOVAR]
dir = /Public/Software/annovar/
db = /SGRNJ/Database/script/database/annovar/humandb
buildver = hg38
protocol = refGene,cosmic70
operation = g,f
```

- ### 构建genelist

```shell
$ cat gene_list.tsv 
EGFR
PIK3CA
BRAF
......
```

### step4. 运行测试

```shell
$ cat run_shell.sh 
multi_snp --mapfile ./snp.mapfile \
          --gene_list ../celescope_test_data/gene_list.tsv \
          --genomeDir /SGRNJ03/randd/public/genome/rna/hs_ensembl_99 \
          --annovar_config ./annovar.config \
          --mod shell \
```

