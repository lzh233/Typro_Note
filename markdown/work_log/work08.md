# work-9.27

##  10X tag数据拆分

### 数据

三个文件分别为

```
10X-HCC19-Index_S1_L001_R1_001.fastq.gz
10X-HCC19-Index_S1_L001_R2_001.fastq.gz
HCC19-10x-valid-barcode.txt
```

### 数据内容查看

10X-HCC19-Index_S1_L001_R1_001.fastq.gz

**barcode + UMI序列**	C16U12

```shell
$ zcat 10X-HCC19-Index_S1_L001_R1_001.fastq.gz | head -n 8
@A00358:556:H3WYYDSX2:3:1101:1787:1110 1:N:0:TAAGGCGA+AAGGCTAT
TCAGCCTAGGCGAAGGACTCTCCTCGAG
+
FF:FFFFFFFFFFFFFFFFFFFFFFFFF
@A00358:556:H3WYYDSX2:3:1101:2003:1110 1:N:0:TAAGGCGA+AAGGCTAT
CAGCGTGTCTACCCACGCCTGTCCTGCA
+
FFFFFFF:FFFFFFFFFFFFFFFFFFFF
```

10X-HCC19-Index_S1_L001_R2_001.fastq.gz

**linker + tag序列**	L25C15

```shell
$ zcat 10X-HCC19-Index_S1_L001_R2_001.fastq.gz | head -n 8
@A00358:556:H3WYYDSX2:3:1101:1787:1110 2:N:0:TAAGGCGA+AAGGCTAT
GTTGTCAAGATGCTACCGTTCAGAGTAAGAGCCCGGCAAGAAAAAAAAAAAAAAAAAAAAAAAAAACTCGAGGAGATTACTTCGCCTAGGC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFF:FFFF,,::,::,:,,,,:FF:,F,,:FF,
@A00358:556:H3WYYDSX2:3:1101:2003:1110 2:N:0:TAAGGCGA+AAGGCTAT
GTTGTCAAGATGCTACCGTTCAGAGTACGAGCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGCAGGAAAGGCGGGGGTAGACCCGGGG
+
FFFFFFFFFFFFFFFFFFFFFFF:FFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,F:,,F,:,,FF,F,FF::,,,,,,,,,,
```

HCC19-10x-valid-barcode.txt 

**barcode序列**

```shell
$ head HCC19-10x-valid-barcode.txt 
barcode
AAACCCAAGAGGGTCT
AAACCCAAGCATCAGG
AAACCCAAGCTACTAC
AAACCCAAGTACTGTC
AAACCCAAGTGCGACA
AAACCCACACCTATCC
AAACCCACAGAGATTA
AAACCCAGTAGGTAGC
AAACCCAGTATGGGAC
```

### tag相关数据

 tag_barcode.fasta 

```shell
$ head tag_barcode.fasta 
>CLindex_TAG_1
CGTGTTAGGGCCGAT
>CLindex_TAG_2
GAGTGGTTGCGCCAT
>CLindex_TAG_3
AAGTTGCCAAGGGCC
>CLindex_TAG_4
TAAGAGCCCGGCAAG
>CLindex_TAG_5
TGACCTGCTTCACGC
```

tag_linker.fasta 

```shell
$ cat tag_linker.fasta 
>CLindex_LINKER
GTTGTCAAGATGCTACCGTTCAGAG
```
### 需求

使用celescope tag流程拆分。
read1 是10X的文库结构, 16bp barcode + 12bp UMI
read2 是Clindex的文库结构。L25C15

### 拆分命令和map文件

```shell
#!/bin/bash
multi_tag \
    --mapfile ./mapfile.tag.map\
    --barcode_fasta ./tag_seq/tag_barcode.fasta\
    --linker_fasta ./tag_seq/tag_linker.fasta\
    --chemistry customized\
    --pattern C16U12\
    --fq_pattern L25C15\
    --mod shell\
    --thread 8
```

```shell
$ cat mapfile.tag.map 
10X-HCC19-Index_S1_L001	./fasta	10X-HCC19	./match_dir_tag
```

match_dir 结构 先跑一下barcode步骤得到`10X-HCC19_2.fq`,  `barcodes.tsv`**是提供的**

```shell
$ tree ./match_dir_tag/
./match_dir_tag/
├── 01.barcode
│   └── 10X-HCC19_2.fq
├── 05.count
│   └── 10X-HCC19_match_dir_matrix_10X
│       └── barcodes.tsv
└── 06.analysis
    ├── markers.tsv
    └── tsne_coord.tsv
```

### celescope相关文档

`--mod` mod, sjm or shell

`--mapfile` tsv file, 4 columns:
                1st col: LibName;
                2nd col: DataDir;
                3rd col: SampleName;
                4th col: optional;

`--rm_files` remove redundant fq.gz and bam after running

`--steps_run` Steps to run. Multiple Steps are separated by comma.

`--outdir` Output directory.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--chemistry` Predefined (pattern, barcode whitelist, linker whitelist) combinations. Can be one of:  

- `auto` Default value. Used for Singleron GEXSCOPE libraries >= scopeV2 and automatically detects the combinations.  
- `scopeV1` Used for legacy Singleron GEXSCOPE scopeV1 libraries.  
- `customized` Used for user defined combinations. You need to provide `pattern`, `whitelist` and `linker` at the 
  same time.

`--pattern` The pattern of R1 reads, e.g. `C8L16C8L16C8L1U12T18`. The number after the letter represents the number 
        of bases.  

- `C`: cell barcode  
- `L`: linker(common sequences)  
- `U`: UMI    
- `T`: poly T

`--whitelist` Cell barcode whitelist file path, one cell barcode per line.

`--linker` Linker whitelist file path, one linker per line.

`--lowQual` Default 0. Bases in cell barcode and UMI whose phred value are lower than lowQual will be regarded as low-quality bases.

`--lowNum` The maximum allowed lowQual bases in cell barcode and UMI.

`--nopolyT` Outputs R1 reads without polyT.

`--noLinker` Outputs R1 reads without correct linker.

`--allowNoPolyT` Allow valid reads without polyT.

`--allowNoLinker` Allow valid reads without correct linker.

`--gzip` Output gzipped fastq files.

`--adapter_fasta` Addtional adapter fasta file.

`--minimum_length` Default `20`. Discard processed reads that are shorter than LENGTH.

`--nextseq_trim` Default `20`. Quality trimming of reads using two-color chemistry (NextSeq). 
Some Illumina instruments use a two-color chemistry to encode the four bases. 
This includes the NextSeq and the NovaSeq. 
In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. 
However, dark cycles also occur when sequencing “falls off” the end of the fragment.
The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.

`--overlap` Default `10`. Since Cutadapt allows partial matches between the read and the adapter sequence,
short matches can occur by chance, leading to erroneously trimmed bases. 
For example, roughly 0.25 of all reads end with a base that is identical to the first base of the adapter. 
To reduce the number of falsely trimmed bases, the alignment algorithm requires that 
at least {overlap} bases match between adapter and read.

`--insert` Default `150`. Read2 insert length.

`--fq_pattern` Required. R2 read pattern. The number after the letter represents the number of bases.         
`L` linker(common sequences)  
`C` tag barcode

`--barcode_fasta` Required. Tag barcode fasta file. It will check the mismatches between tag barcode 
sequence in R2 reads with all tag barcode sequence in barcode_fasta. 
It will assign read to the tag with mismatch < len(tag barcode) / 10 + 1. 
If no such tag exists, the read is classified as invalid.

```
>tag_0
GGGCGTCTGTGACCGCGTGATACTGCATTGTAGACCGCCCAACTC
>tag_1
TTCCTCCAGAGGAGACCGAGCCGGTCAATTCAGGAGAACGTCCGG
>tag_2
AGGGCTAGGCGTGTCATTTGGCGAGGTCCTGAGGTCATGGAGCCA
>tag_3
CACTGGTCATCGACACTGGGAACCTGAGGTGAGTTCGCGCGCAAG
```

`--linker_fasta` Optional. If provided, it will check the mismatches between linker sequence in R2 reads 
with all linker sequence in linker_fasta. If no mismatch < len(linker) / 10 + 1, the read is classified as invalid.

`--UMI_min` Default='auto'. Minimum UMI threshold. Cell barcodes with valid UMI < UMI_min are classified as *undeterminded*.

`--dim` Default=1. Tag dimentions. Usually we use 1-dimentional tag.

`--SNR_min` Default='auto'. Minimum signal-to-noise ratio. 
Cell barcodes with UMI >=UMI_min and SNR < SNR_min are classified as *multiplet*.

`--combine_cluster` Conbine cluster tsv file.

`--coefficient` Default=0.1. If `SNR_min` is 'auto', minimum signal-to-noise ratio is calulated as 
`SNR_min = max(median(SNRs) * coefficient, 2)`. 
Smaller `coefficient` will cause less *multiplet* in the tag assignment.

`--split_fastq` If used, will split scRNA-Seq fastq file according to tag assignment.

`--split_matrix` If used, will split scRNA-Seq matrix file according to tag assignment.

