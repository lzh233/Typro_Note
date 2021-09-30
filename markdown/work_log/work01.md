# work1- 8.27

## 描述

report添加检出病毒的细胞数，和每个cluster含有的病毒的细胞数

**9.29 multi_capture_virus  bug 百分比小100倍 原因是百分比计算后忘记*100**

已经修复 但是没有PR

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210930170541785.png" alt="image-20210930170541785" style="zoom:80%;" />

## 解决流程

PR过程

https://github.com/singleron-RD/CeleScope/pull/56/files

### step 1. 首先根据分析结果去计算检出病毒的细胞数量

**tag列带有EBV标签为检出病毒的cell**

使用分组求和计算即可

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210913141220805.png" alt="image-20210913141220805" style="zoom:80%;" />

```python
        #get the sum of virus
        total_virus = virus_tsne_df.loc[:,"tag"].count()
        #get the sum of virus in each cluster
        cluster_virus = virus_tsne_df.groupby("cluster")["tag"].count()
        #Total cell each cluster
        total_cell_cluster = virus_tsne_df.groupby("cluster")["barcode"].count()
        #add sum of virus to metric
        self.add_metric(
                        name = 'Number of Cells with Virus',
                        value = f'{total_virus}({round(total_virus/sum(total_cell_cluster),2)}%)'
                        )
        #add sum of virus in each cluster to metric
        cluster = 1    
        for num_virus in cluster_virus:
            total_cell = total_cell_cluster[cluster]
            self.add_metric(
                            name = f'Cluster_{cluster}',
                            value = f'{num_virus}({round(num_virus/total_cell,2)}%)'
                            )
            cluster += 1
```

### step 2. 使用step类中的 add_metric方法将结果写入

```python
  self.add_metric(
                  name = 'Number of Cells with Virus',
                  value = f'{total_virus}({round(total_virus/sum(total_cell_cluster),2)}%)'
                  )
  self.add_metric(
                  name = f'Cluster_{cluster}',
                  value = f'{num_virus}({round(num_virus/total_cell,2)}%)'
                  )
```

实际写入的是`.metrics.json`文件

查看一下`.metrics.json`, 每一项的名称都是`<step_name>_summarry`, 如`analysis_capture_virus_summary`

```apl
    ......
    "analysis_capture_virus_summary": {
        "otsu UMI threshold": 1,
        "Total Virus Cell": 53,
        "Total Virus Cell Fraction": 0.0005,
        "Cluster_1": 19,
        "Cluster_1 Fraction": 0.0004,
        "Cluster_2": 15,
        "Cluster_2 Fraction": 0.0004,
        "Cluster_3": 19,
        "Cluster_3 Fraction": 0.0006
    }
}
```

### step 3. 修改templates目录

在template目录中，每个方法都会对应了一个html模板，每个报告都是一个独立的`html`文件

<img src="C:\Users\liuzihao\AppData\Roaming\Typora\typora-user-images\image-20210913142929519.png" alt="image-20210913142929519" style="zoom:80%;" />

- 首先增加一个`virus_summary.html`，可以将别的html文件复制一下，然后**重新修改3个地方**即可

```html
<div class="abc" style="float: left; margin-left: 15%; margin-right:15%; width: 70%" >
    <!-- Cells 为输出报告的名称-->
    <h2>Cells    <i class="fa fa-question-circle" onClick="toggle1(this)" style="cursor:pointer;"></i></h2>
    <div class="box">
      <div class="description" style="display: none;">
          <!--描述信息-->
        <p><b>Number of Cells with Virus</b> : The sum of virus cells.</p>
        <p><b>Cluster_n(n = 1,2,3...)</b> : The number virus cells in each cluster.</p>
      </div>
      <table style="float: left; margin-left: 0%; margin-right:3%; width: 47%">
          <!-- 指定描述哪一个项目，需要和.metrics.json中一致，如本例为"analysis_capture_virus_summary-->
        {% for item in analysis_capture_virus_summary %}
          {% if loop.index <= (loop.length+1)/2 %}
          <tr>
            {% for i in item %} 
            <td>{{ i|e }}</td>
            {% endfor %}
          </tr>
          {% endif %}
        {% endfor %}
      </table>

      <table style="float: left; margin-left: 3%; margin-right:0%; width: 47%">
           <!-- 指定描述哪一个项目，需要和.metrics.json中一致，如本例为"analysis_capture_virus_summary-->
        {% for item in analysis_capture_virus_summary %}
          {% if loop.index > (loop.length+1)/2 %}
          <tr>
            {% for i in item %} 
            <td>{{ i|e }}</td>
            {% endfor %}
          </tr>
          {% endif %}
        {% endfor %}
      </table>

      <div class="clear" ></div>

    </div>
  </div>
```

- 然后下一步修改`base.html`, 增加以下内容

```html
<!--可以将其他的复制过来修改相应地方即可，第一行要和.metrics.json中一致，第二行调用的文件为上一步创建的html文件-->      
{% if analysis_capture_virus_summary is defined %}
{% include "html/capture_virus/virus_summary.html"%}
{% endif %}
```

### step4. 查看report

![image-20210913144101429](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210913144101429.png)

