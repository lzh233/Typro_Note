# work-9 9.28

## 描述

1. 在07.variant_calling中添加一个输出文件

{sample}_variant_ncell.tsv with 4 columns

col1 : VID

col2: ncell_cover: number of cells with read count

col3: ncell_alt: number of cells with variant read count only

col4: ncell_ref: number of cells with reference read count only

同时将08.analysis_snp/{sample}_variant_table.tsv中的nCell 替换成 ncell_cover和ncell_alt两列。

2. 取ncell_alt 从高到低排序前5的variant， 画一个venn图， 输出到08.analysis_snp中

可以用https://github.com/tctianchi/pyvenn 这个包

## 解决流程

需求1 与work04 05相似 

PR流程https://github.com/singleron-RD/CeleScope/pull/72

venn图使用的**venn**包绘制, 新增的包需要写入`requirement.txt`

具体使用：

https://github.com/LankyCyril/pyvenn

绘制方法

```python
....
"""
cid_lst = [(1,2,3),(3,5,6),(1,3)......]
vid_lst = [1,2,3....]
cid_lst and vid_lst have same length
"""
 def get_venn_plot(self):
        df_top_5 = self.get_df_table().sort_values(by = "ncell_alt",ascending=False).iloc[:5,:]
        plot = {}
        cid_lst = df_top_5.loc[:,"CID"].to_list()
        vid_lst = df_top_5.loc[:,"VID"].to_list()
        for cid,vid in zip(cid_lst,vid_lst):
            plot[f"VID_{vid}"] = set(cid)
        #venn plot
        set_cid = list(plot.values())
        set_name = list(plot.keys())
        labels = generate_petal_labels(set_cid)
        plot = draw_venn(
                         petal_labels=labels, 
                         dataset_labels=set_name,
                         hint_hidden=False,
                         colors=generate_colors(n_colors=5), 
                         figsize=(8, 8),
                         fontsize=14, 
                         legend_loc="best",
                         ax=None
                         )
        fig = plot.get_figure()
        fig.savefig(f'{self.outdir}/{self.sample}_variant_top5.jpg',dpi = 600)
....
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210930160821595.png" alt="image-20210930160821595" style="zoom: 50%;" />

celescope增加选项参数ncell_file

修改需要添加参数的脚本`analysis_snp.py`

```python
....
def get_opts_analysis_snp(parser, sub_program):
    parser.add_argument('--annovar_config', help='ANNOVAR config file.', required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help=HELP_DICT['match_dir'], required=True)
        parser.add_argument('--filter_vcf', help='filter vcf file.', required=True)
        parser.add_argument('--CID_file', help='CID_file.', required=True)
        parser.add_argument('--filter_variant_count_file', help='filter variant count file.', required=True)
        #增加需要添加的参数 --ncell_file
        parser.add_argument('--ncell_file', help='filter cell count file.', required=True)
```

修改`multi_snp.py`

```python
def analysis_snp(self, sample):
    step = 'analysis_snp'
    filter_vcf = f'{self.outdir_dic[sample]["variant_calling"]}/{sample}_filter.vcf'
    CID_file = f'{self.outdir_dic[sample]["variant_calling"]}/{sample}_CID.tsv'
    filter_variant_count_file = f'{self.outdir_dic[sample]["variant_calling"]}/{sample}_filter_variant_count.tsv'
    cmd_line = self.get_cmd_line(step, sample)
    #指定参数的值是多少
    ncell_file = f'{self.outdir_dic[sample]["variant_calling"]}/{sample}_variant_ncell.tsv'
    cmd = (
            f'{cmd_line} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--filter_vcf {filter_vcf} '
            f'--CID_file {CID_file} '
            f'--filter_variant_count_file {filter_variant_count_file} '
            f'--ncell_file {ncell_file}'
    )
    self.process_cmd(cmd, step, sample, m=8, x=1)
```

**结果**

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210930165406731.png" alt="image-20210930165406731" style="zoom: 80%;" />

