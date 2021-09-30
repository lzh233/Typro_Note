# work-9.15

## 描述

输入数据：
--bam bam：/SGRNJ03/randd/RD20102301_DZH/Blood_panel/20210907/B0831_R1_V2FJ/06.target_metrics/B0831_R1_V2FJ_filtered_sorted.bam
--bed bed: /SGRNJ03/randd/test_rd/dxh/data/Blood/blood_Aligned.sortedByCoord.out.bed
--bp_min 50 (可选，默认为50bp)
--zl zl：/SGRNJ03/randd/RD20102301_DZH/Blood_panel/20210906/B0831_R1_V2ZL/06.analysis/B0831_R1_V2ZL_tsne_coord.tsv

需求：
  1，统计bam文件中在指定位置（根据bed文件）的reads和UMI的cover情况，输出统计表格
  position：bed文件中的位置对应的名称
  cell:在该位置cover的细胞数
  reads_number:在该位置的reads数（要求read 在指定位置bp至少大于bp_min）
  umi_number：在该位置的umi种类数
  mean_reads:对每个位置的每个cell的read数进行排序，reads的平均数
  medium_reads：对每个位置的每个cell的read数进行排序，reads的中位数
  mean_umi：对每个位置的每个cell的umi数进行排序，umi的平均数
  medium_umi:对每个位置的每个cell的umi数进行排序，umi的中位数

  \##example
  position  cell  reads_number  umi_number  mean_reads  medium_reads  mean_umi  medium_umi
  位置1  
  位置2
  位置3
  ...

2，根据position的cell barcode和zl文件，整合输出position ，cell barcode list，cluster 表格

\##example
position  cell barcode list  cluster
位置1    AAACATCGAAACATCGCCGAAGTA  1
位置2
位置3...

## 解决流程

脚本位置

/SGRNJ03/randd/user/liuzihao/statistics_bam

主要问题

`-bp_min`的计算

使用了pysam.fetch(contig, start, end) 抓取了在指定位置匹配的reads，在通过get_overlap() 方法获得overlap的长度，`-bp_min`即其阈值

实现：

[pysam文档](https://pysam.readthedocs.io/_/downloads/en/latest/pdf/)

```python
....
for chr,start,end in zip(position_star_end[0],position_star_end[1],position_star_end[2]):
                for readr in samfile_sample.fetch(str(chr),int(start),int(end)):
                    length_overlap = readr.get_overlap(int(start),int(end))
                    #print(readr.get_refernece_name())
                    #length_overlap = readr.reference_length
                    #length_overlap  = readr.query_alignment_length
                    #length_overlap = readr.infer_query_length()
                    #length_overlap = readr.infer_read_length()
                    """
                    https://pysam.readthedocs.io/_/downloads/en/latest/pdf/ 
                    page 14
                    get_overlap(self, uint32_t start, uint32_t end)
                    return number of aligned bases of read overlapping the interval start and end on the        
                    reference sequence.
                    Return None if cigar alignment is not available.
                    """
...
```

bed文件

```shell
$ cat data/sample1/blood_Aligned.sortedByCoord.out.bed 
#第一列染色体编号
#第二列起始position
#结束position
11	32396278	32399968	WT1_PCR2	255	-
11	32396281	32399970	WT1_PCR3	255	-
11	108319960	108321352	Blood-ATM6491_2轮扩增-1	255	+
12	25227401	25250767	Blood-KRAS-2轮扩增-2	255	-
13	28018481	28023391	Blood-FLT3-2轮扩增-2	255	-
13	28028176	28033924	Blood-FLT3-2轮扩增-1	255	-
13	28028245	28033992	Blood-FLT3-2轮扩增-3	255	-

```



## 脚本使用

```python
python 3.6 or later
require: pysam samtools
```

## options

`--data_dir`: The directory which saved the data, make sure this directory only save the result which you need, if not, please use  `--manual_sample_list` and `--sample_list_dir` to provide a sample name list.

 `--results_dir` : Output directory.

`--manual_sample_list`: If provide this option, the script will not create samlple name list by auto.

`--sample_list_dir`: When use `--manual_sample_list` you must use  `--sample_list_dir` to provide a sample name list file.

`--bp_min`: Provide a min length to filter data (default: 50)

## data_dir  sample list

### data_dir

**The bam file may not have index file,  the script will run `samtools index <bam file>` by auto to make a bai file**

```shell
$ ls data/
test1  test2  test3
$ ls ./data/test1/
B0831_R1_V2FJ_filtered_sorted.bam  B0831_R1_V2FJ_filtered_sorted.bam.bai  B0831_R1_V2ZL_tsne_coord.tsv  blood_Aligned.sortedByCoord.out.bed
```

### sample  list

**if you provide a manual sample list ,it should not have duplicate name, if have, the script will drop it by auto**

```shell
test1
test2
```

## usage

```shell
#Quick start
python get_cover_condition.py --data_dir ./data/ --results_dir ./results
```

```shell
# If you want set a min length of reads use option of --bp_min n
python  get_cover_condition.py --data_dir ./data/ --results_dir ./results --bp_min 100
```

```shell
# If you want to manual make sample name list
python get_cover_condition.py --data_dir ./data/ --results_dir ./results --manual_sample_list --sample_list_dir ./sample_list.txt
```
