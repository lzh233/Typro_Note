# work-9.6

## 描述

统计每个突变（VID）检测到该位点的细胞数

一个细胞ref_count > 0或者alt_count>0 即认为检测到了这个位点。

celescope.snp.variant_calling 增加输出文件：对于每个VID，检测到该位点的细胞数，用于统计捕获效率。

## 解决流程

PR过程

https://github.com/singleron-RD/CeleScope/pull/59/files

### step1. 查看此步骤的输出文件`07.variant_calling/snp_variant_count.tsv`

```shell
VID	CID	ref_count	alt_count
1	185	0	2
1	231	6	0
1	916	0	2
2	185	2	0
2	231	7	0
2	916	0	2
3	213	1	0
3	231	0	4
4	213	0	1
......
```

### step2. 计算`nCell_with_read_count`、`with_ref_read `,`with_variant_read`

```python
.......
#df_filter为上述文件
#判断: 如果ref_count和alt_count和为零 表示为检测到 丢弃
df_filter.loc[:,"vid_judge"] = df_filter.loc[:,"ref_count"] + df_filter.loc[:,"alt_count"]
df_filter = df_filter[df_filter.loc[:,"vid_judge"] > 0]
#summarize
vid_summarize = {}
#得到vid
vid_summarize["VID"] = list(set(df_filter.loc[:,"VID"]))
#计算每个突变（VID）检测到该位点的细胞数
vid_summarize["nCell_with_read_count"] = list(df_filter.groupby("VID")["vid_judge"].count())

#计算该突变位点来自哪里 ref还是alt，以及多少是来自ref多少是alt，将那两列转换为01即可
variant_count =  (df_filter.loc[:,"alt_count"] != 0).astype(int)
ref_count =  (df_filter.loc[:,"ref_count"] != 0).astype(int)
variant_count["VID"] = df_filter.loc[:,"VID"]
ref_count["VID"] = df_filter.loc[:,"VID"]

vid_summarize["with_ref_read"] = list(ref_count.groupby("VID").sum())
vid_summarize["with_variant_read"] = list(variant_count.groupby("VID").sum())

vid_summarize = pd.DataFrame(vid_summarize)

vid_summarize.to_csv(self.summarize_capture_vid, sep='\t', index=False)
......
```

### step3. 输出文件

```shell
$ cat snp_summarize_capture_vid.tsv 
VID	nCell_with_read_count	with_ref_read	with_variant_read
1	3	1	2
2	3	2	1
3	2	1	1
4	2	1	1
5	2	1	1
```

