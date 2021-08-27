# Veen图和upset图

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-7-22
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)
```

### 2021-7-23Veen图和upset图

### veen图

#### 使用ggVennDiagram绘制veen图

```r
#install.package("ggVennDiagram")
library(tidyverse)
library(ggVennDiagram)
library(ggsci)
#构建数据
genes <- paste0("gene",1:1000)
set.seed(20210502)
#将维恩图的每个圈里的数据放到一个list中
gene_list <- list(A = sample(genes,100),
                  B = sample(genes,200),
                  C = sample(genes,300),
                  D = sample(genes,200))
#绘制韦恩图
ggVennDiagram(gene_list, 
              #指定每个list的名字
              category.names = c("A","B","C","D"),
              #指定每个名字的颜色
              set_color = c("red1","red2","red3","red4"),
              #指定每个名字的大小
              set_size = 10,
              #指定线的类型，具体线条类型见下
              edge_lty = "dashed",
              #指定线的粗细
              edge_size = 1,
              #指定标签的内容: counts(仅仅计数结果),percent(百分比),both(全部显示),none
              label = "both",
              #指定label的颜色
              label_color = "black",
              #标签的透明度
              label_alpha = 0.5) +
  #具体配色代码见下方
  scale_fill_distiller(palette = "RdBu") +
  #自行指定外框的颜色
  scale_color_manual(values = c(rep("black",4)))
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/Typro_imgimage-20210723214549290.png" alt="image-20210723214549290" style="zoom: 67%;" />



#### R的线条类型与点类型以及调色板

```r
#install.packages("ggpubr")
a <- ggpubr::show_point_shapes()
b <- ggpubr::show_line_types()

cowplot::plot_grid(a,b)
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/Typro_imgimage-20210723214926407.png" alt="image-20210723214926407" style="zoom: 50%;" />

```r
#调色板
RColorBrewer::display.brewer.all()
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/Typro_imgimage-20210723215145679.png" alt="image-20210723215145679" style="zoom:50%;" />

### Upset图

韦恩图不适合展示维度过高的数据，个人认为最多四个维度，当维度更高时，用***Upset***图会更好

#### 基本使用

- 基本用法, 使用数据集为package里自带的`movies`数据集

```r
#install.package("UpSetR")
library(UpSetR)
#载入数据集
movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=TRUE, sep=";" )
#简单查看数据集，数据集第一列为电影的名称，第二列为发行年份，第三列往后为电影类型，1代表为对应类型，
> head(movies)
                                Name ReleaseDate Action Adventure Children Comedy Crime Documentary Drama Fantasy Noir Horror
1                   Toy Story (1995)        1995      0         0        1      1     0           0     0       0    0      0
2                     Jumanji (1995)        1995      0         1        1      0     0           0     0       1    0      0
3            Grumpier Old Men (1995)        1995      0         0        0      1     0           0     0       0    0      0
4           Waiting to Exhale (1995)        1995      0         0        0      1     0           0     1       0    0      0
5 Father of the Bride Part II (1995)        1995      0         0        0      1     0           0     0       0    0      0
6                        Heat (1995)        1995      1         0        0      0     1           0     0       0    0      0
  Musical Mystery Romance SciFi Thriller War Western AvgRating Watches
1       0       0       0     0        0   0       0      4.15    2077
2       0       0       0     0        0   0       0      3.20     701
3       0       0       1     0        0   0       0      3.02     478
4       0       0       0     0        0   0       0      2.73     170
5       0       0       0     0        0   0       0      3.01     296
6       0       0       0     0        1   0       0      3.88     940
#绘制Upset图
upset(data = movies, 
      sets = c("Drama", "Comedy", "Action", "Thriller", "Western", "Documentary"))
```

- `data = `: 指定数据集
- `sets = `: 指定作图的数据

```r
#图表的左下角为指定的六个数据集合的大小，即所有上映的电影中属于各个类型电影的数量
> movies %>% select(c("Drama", "Comedy", "Action", "Thriller", "Western", "Documentary")) %>% colSums()
      Drama      Comedy      Action    Thriller     Western Documentary 
       1603        1200         503         492          68         127
#右上角的图是对应右下方的集合中的元素的个数，如第一根柱子为1184，代表仅属于Daram类型的电影有1184部，第七根柱子则代表同时属于Drama和Comedy的电影为211部
```





<img src="https://aironi.oss-cn-beijing.aliyuncs.com/Typro_imgimage-20210724163050394.png" alt="image-20210724163050394" style="zoom:70%;" />

#### 展示指定部分

- 展示指定部分**（*intersections*）**

  使用`queries`选项来对指定的数据进行展示，默认提供`intersections`和` elements`两种展示方式，首先是`intersections`，对指定的交集进行展示

  命令格式, 

  `query`为指定展示数据的方式, 数据类型为列表, 

  `params`, 指定展示哪个交集，数据类型为列表，列表内的元素为需要展示部分的数据名称，

  `color`为指定展示的颜色，

  `active`为是否对柱子着色，`=T`代表上色，`=F`则在柱子上用一个与`color`指定的一致颜色的小三角表示

  **命令格式**

  `queries = list(query = intersects, params = list(a,b,c,d), color = "", active = T)`

```r
#展示"Action","Comedy", "Drama"三种类型的交集
#展示"Drama", "Comedy"两种类型的交集
#展示"Western","Drama"两种类型的交集，但是柱子不填充颜色
#展示"Action"特有的部分
upset(data = movies,sets = c("Drama", "Comedy", "Action", "Thriller", "Western", "Documentary"),
      queries = list(list(query = intersects,params = list("Action","Comedy", "Drama"),color = "red1",active = T),
                     list(query = intersects,params = list("Drama", "Comedy"),color = "blue2",active = T),
                     list(query = intersects,params = list("Western","Drama"),color = "green2",active = F)))
```

![image-20210724165132800](https://aironi.oss-cn-beijing.aliyuncs.com/Typro_imgimage-20210724165132800.png)

#### 展示指定部分在某因子上的水平

- 展示某部分在某因子上的水平**（*elements*）**

  **命令格式**

`queries = list(query = elements, params = list("factors",b,c,d), color = "", active = T)`

```r
#展示各个交集在分布在1998和1995年的数量
upset(data = movies,sets = c("Drama", "Comedy", "Action", "Thriller", "Western", "Documentary"),
      queries = list(list(query = intersects,params = list("Action","Comedy", "Drama"),color = "red1",active = T),
                     list(query = intersects,params = list("Drama", "Comedy"),color = "blue2",active = T),
                     list(query = intersects,params = list("Western","Drama"),color = "green2",active = F),
                     list(query = intersects,params = list("Action"),color = "pink",active = T),
                     list(query = elements, params = list("ReleaseDate", 1998,1995),color = "purple",active = T)))
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/Typro_imgimage-20210724184452922.png" alt="image-20210724184452922" style="zoom:80%;" />

#### 增加箱线图

- 可以通过指定因子来展示在不同因子上，每个交集的分布情况，如，通过箱线图展示在1920-2000年间，不同的intersects在每所有年份的分布情况

  `boxplot.summary = "factors"`: 可以指定多个因子

```r
upset(data = movies,sets = c("Drama", "Comedy", "Action", "Thriller", "Western", "Documentary"),
      queries = list(list(query = intersects,params = list("Action","Comedy", "Drama"),color = "red1",active = T),
                     list(query = intersects,params = list("Drama", "Comedy"),color = "blue2",active = T),
                     list(query = intersects,params = list("Western","Drama"),color = "green2",active = F),
                     list(query = intersects,params = list("Action"),color = "pink",active = T),
                     list(query = elements, params = list("ReleaseDate", 1998,1995),color = "purple",active = T)),
      boxplot.summary = "ReleaseDate")
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/Typro_imgimage-20210724184353583.png" alt="image-20210724184353583" style="zoom:80%;" />

#### 增加其他图表

- 通过`attribute.plots`可以增加其他的图表，如直方图，点图，箱线图等等......

  选项用法：`attribute.plots = list()`，可以通过构建函数的方法进行绘图，每张图片放在一个`list`中，`gridrows`为图的高度，`ncols`用于控制列数

使用方法

```r
#添加一个直方图，表示每个年份的电影数量情况
#首先构建一个绘制直方图的function
plot1 <- function(data,x){
  plot.his <- (ggplot(data,aes_string(x=x,fill = "color"))) + 
    geom_histogram(bins = 30)+  
    scale_fill_identity() + 
    labs(y = "Number of movies",x = "Years") +theme_bw()
}
#将直方图写入attribute.plots，其中x= "ReleaseDate"为传给x的参数，queries = TRUE表示将上面展示的数据同样展示在直方图中，
upset(data = movies,sets = c("Drama", "Comedy", "Action", "Thriller", "Western", "Documentary"),
      queries = list(list(query = intersects,params = list("Action","Comedy", "Drama"),color = "red1",active = T),
                     list(query = intersects,params = list("Drama", "Comedy"),color = "blue2",active = T),
                     list(query = intersects,params = list("Western","Drama"),color = "green2",active = F),
                     list(query = intersects,params = list("Action"),color = "pink",active = T),
                     list(query = elements, params = list("ReleaseDate", 1998,1995),color = "purple",active = T)),
      attribute.plots = list(gridrows = 55,
                             plots = list(list(plot = plot1,x= "ReleaseDate", queries = TRUE)),ncols = 1))
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/Typro_imgimage-20210725213819483.png" alt="image-20210725213819483" style="zoom: 80%;" />

```r
#添加一个散点图（如，Avgrating和观看人数的关系），同样是先构建绘图函数
plot2 <- function(data,x,y){
  polt.point <- (ggplot(data,aes_string(x=x,y=y,color = "color"))) +
    geom_point() +
    scale_color_identity() +
    theme_bw()
}
##将散点图图写入attribute.plots
upset(data = movies,sets = c("Drama", "Comedy", "Action", "Thriller", "Western", "Documentary"),
      queries = list(list(query = intersects,params = list("Action","Comedy", "Drama"),color = "red1",active = T),
                     list(query = intersects,params = list("Drama", "Comedy"),color = "blue2",active = T),
                     list(query = intersects,params = list("Western","Drama"),color = "green2",active = F),
                     list(query = intersects,params = list("Action"),color = "pink",active = T),
                     list(query = elements, params = list("ReleaseDate", 1998,1995),color = "purple",active = T)),
      attribute.plots = list(gridrows = 55,
                             plots = list(list(plot = plot1,x= "ReleaseDate", queries = TRUE),
                                          list(plot = plot2,x = "AvgRating",y = "Watches",queries = TRUE)),ncols = 2))
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/Typro_imgimage-20210725215101963.png" alt="image-20210725215101963" style="zoom:70%;" />

```r
#增加一个抖动图
  plot3 <- function(data,x,y){
    polt.jit <- (ggplot(data,aes_string(x=x,y=y,color = "color"))) +
      geom_jitter() +
      scale_color_identity() +
      theme_bw()
  }
##将抖动图写入attribute.plots
upset(data = movies,sets = c("Drama", "Comedy", "Action", "Thriller", "Western", "Documentary"),
      queries = list(list(query = intersects,params = list("Action","Comedy", "Drama"),color = "red1",active = T),
                     list(query = intersects,params = list("Drama", "Comedy"),color = "blue2",active = T),
                     list(query = intersects,params = list("Western","Drama"),color = "green2",active = F),
                     list(query = intersects,params = list("Action"),color = "pink",active = T),
                     list(query = elements, params = list("ReleaseDate", 1998,1995),color = "purple",active = T)),
      attribute.plots = list(gridrows = 55,
                             plots = list(list(plot = plot1,x= "ReleaseDate", queries = TRUE),
                                          list(plot = plot2,x = "AvgRating",y = "Watches",queries = TRUE),
                                          list(plot = plot3,x = "Watches",y = "AvgRating")),ncols = 3))

```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/Typro_imgimage-20210725220526150.png" alt="image-20210725220526150" style="zoom:80%;" />



