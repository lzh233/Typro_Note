# work-9.9

## 描述

把https://github.com/singleron-RD/CeleScope/pull/59 这个功能抽出来单独变成一个脚本，这样已经跑完的样本直接就可以出结果----将`work04`的功能单独写成一个脚本, 统计已经完成的样本，**并且增加每个细胞的cluster来源信息 (9.14)**

**脚本修改(9.29):  **

**问题：  正常情况不会出现no_cluster ,结果却出现了**

**将cluster信息的来源改为转录组的数据，解决`no_cluster`问题**

**出问题原因是snp的tsne的cluster信息是仅仅包含了有突变信息的细胞，并非全部细胞，所以出现了`no_cluster`**

**已经修复**，更换tsne文件

![image-20210930170052718](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210930170052718.png)

## 脚本位置

/SGRNJ03/randd/user/liuzihao/vid_sumarrize 或 

/Personal/liuzihao/work/vid_sumarrize

## 使用方法

## environment

```python
python 3.6 or later
```

```apl
Script: /SGRNJ03/randd/user/liuzihao/vid_sumarrize/
Test data: /SGRNJ03/randd/user/liuzihao/vid_sumarrize/snp_data
```

## options

`--data_dir`: The directory which saved the celescope snp results, make sure this directory only save the data of snp results, if not please use  `--manual_sample_list` and `--sample_list_dir` to provide a sample name list.

 `--results_dir` : Output directory.

`--detect_share`: If provide this option, the script will find the vid which meanwhile in alt_count and ref_count. This result will save in output file's `both_ref_and_variant` column.

`--manual_sample_list`: If provide this option, the script will not create samlple name list by auto.

`--sample_list_dir`: When use `--manual_sample_list` you must use  `--sample_list_dir` to provide a sample name list file.

## data_dir  and sample list

### data_dir

```shell
$ ls ./snp_data
snp1  snp2  snp3
$ ls ./snp_data/snp1
00.sample  01.barcode  02.cutadapt  03.consensus  04.star  05.featureCounts  06.target_metrics  07.variant_calling  08.analysis_snp  snp_report.html
```

### sample  list

```shell
snp1
snp2
snp3
```

## usage

```shell
#Quick start
python get_vid_summarize.py --data_dir ./snp_data/ --results_dir ./results
```

```shell
# If you want to find the vid which meanwhile in alt_count and ref_count
python get_vid_summarize.py --data_dir ./snp_data/ --results_dir ./results --detect_share
```

```shell
# If you want to manual make sample name list
python get_vid_summarize.py --data_dir ./snp_data/ --results_dir ./results --manual_sample_list --sample_list_dir ./sample_list.txt
```

## results

```python
"""
Format
{sample name}_summarize_capture_vid.tsv 
"""
$ cat snp1_summarize_capture_vid.tsv 
VID	nCell_with_read_count	with_ref_read	with_variant_read protein cluster_1 cluster_2 ...
1	3	1	2	protein1 1	2 ...
2	3	2	1	protein2 1	2 ...
3	2	1	1	protein3 1	1 ...
4	2	1	1	protein4 1	1 ...
......
```

## log.txt

```python
#You can check log.txt to konw is it successful
$ cat ./log.txt
snp1	Success
snp2	Success
snp3	Success
```