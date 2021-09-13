# work02-9.1

## 描述

target metrics 统计指标 UMI改成read

```python
#位置
celescope.tools.target_metrics

全部改成统计read
        self.add_metric(
            name="Total UMIs",
            value=total_UMIs,
        )

        self.add_metric(
            name="Enriched UMIs",
            value=enriched_UMIs,
            total=total_UMIs,
        )
        self.add_metric(
            name="Enriched UMIs in Cells",
            value=enriched_UMIs_in_cells,
            total=total_UMIs,
        )
        self.add_metric(
            name="Median Enriched UMIs per Cell",
            value=np.median(enriched_UMIs_per_cell_list),
        )
```



## 解决流程

PR过程

https://github.com/singleron-RD/CeleScope/pull/57/files

### step1.  定位到文件位置, 将计算的UMI数量改为reads数量

celescope.tools.target_metrics

需要修改的代码以及修改后如下，

```python
......
def parse_count_dict_add_metrics(self):
    total_reads = 0
    enriched_reads = 0
    enriched_reads_in_cells = 0
    enriched_reads_per_cell_list = []
    #count_dict : 为上一步通过pysam对bam文件迭代得来，基本结构如下，
    """
    count_dict结构字典套字典 ({barcode1:{gene_name1:{umi:1,umi2:2,umi3:3....}}},.....)
    """
#读取barcode
    for barcode in count_dict:
    #每个cell富集的umi数量
        cell_enriched_reads = 0
    #取出每个barcode中包含的基因名
        for gene_name in count_dict[barcode]:
        #每个基因中对应的umi个数 所以用len(),方法为len()
               #gene_UMI = len(count_dict[barcode][gene_name])
        #计算每个reads被富集的数量
               gene_reads = sum(count_dict[barcode][gene_name].values())
        #每次总的umi 累加
               total_reads += gene_reads
        #判断基因是否在需要的基因列表中，如果该基因在，则认为该基因为要的基因 对应的umi个数(len())为富集的umi数量
               if gene_name in gene_list:
                   enriched_reads += gene_reads
        #进一步判断这个barcode是不是在barcode列表中，如果在 那么这个为富集在细胞内的umi
                if barcode in match_barcode:
                   enriched_reads_in_cells += gene_reads
                   cell_enriched_reads += gene_reads

      if barcode in match_barcode:
           enriched_reads_per_cell_list.append(cell_enriched_reads)
        
        
 self.add_metric(
        name="Enriched Reads",
        value=enriched_reads,
        total=total_reads,
        )
self.add_metric(
       name="Enriched Reads in Cells",
       value=enriched_reads_in_cells,
       total=total_reads,
        )
self.add_metric(
       name="Median Enriched Reads per Cell",
       value=np.median(enriched_reads_per_cell_list),
        )
......
```

### step2. 修改相应的html文件，更改输出report

```html
......
    <div class="box">
      <div class="description" style="display: none;">
          <!--修改描述和名称-->
        <p><b>Number of Target Genes</b> : number of target genes in gene_list.</p>
        <p><b>Enriched UMIs</b> : number of UMIs mapped to target genes.</p>
        <p><b>Enriched UMIs in Cells</b> : number of enriched cell UMIs.</p>
        <p><b>Median Enriched UMIs per Cell</b> : median of enriched UMIs per cell(cell with Enriched UMIs > 0).</p>
        <p><b>Enriched Reads</b> : number of reads mapped to target genes.</p>
        <p><b>Enriched Reads in Cells</b> : number of enriched cell reads.</p>
        <p><b>Median Enriched Reads per Cell</b> : median of enriched reads per cell(cell with Enriched reads > 0).</p>
      </div>
      <table style="float: left; margin-left: 0%; margin-right:3%; width: 47%">
        {% for item in target_metrics_summary %}
 ......
```

