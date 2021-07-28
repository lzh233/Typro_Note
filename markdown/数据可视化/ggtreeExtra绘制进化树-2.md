# ggtreeExtra绘制发育树

## R版本与运行环境信息

```R
https://www.bilibili.com/read/cv8074214?share_medium=iphone&share_plat=ios&share_source=WEIXIN&share_tag=s_i&timestamp=1618403699&unique_k=Aga8BO
Author:liuzihao
Date:2021-4-14
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)
```

## 载入相关包

```R
library("ggtree")
library("tidyverse")
library("treeio")
library("ggtreeExtra")
library("ggsci")
#ggtreeExtraA安装
#install.packages("remotes")
#remotes::install_github("YuLab-SMU/ggtreeExtra")
```

## 设置工作目录

```R
setwd("D:\\Desktop\\s_note\\data\\ggtreeextract")

#node.label = "support"`**: 表示节点可以增加其他信息
#**`split()`**: 按照分组因子，把向量，矩阵和数据框进行适当的分组。它的返回值是一个列表，代表分组变量每个水平的观测。
```

## 绘制进化树

### 构建进化树

```R
#读入数据
tree <- read.tree("tree.nwk")
#使用ggtree构建系统发育树
p1 <- ggtree(tree,layout = "fan",open.angle =10,size = 0.5)
###常用选项
#layout：树形状，常见的又fan、circual（圆形）、rectangular（正常的竖）
#open.angle：开口角度
#size: 线条粗细
#branch.length：none则忽略枝条长度
#root.position：指定根节点，默认0为根节点
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210606210356674.png" alt="image-20210606210356674" style="zoom:50%;" />

### 读入第一圈数据

```R
dat1 <- read.csv("tree_tippoint_bar.csv",header = T)
#查看第一圈数据情况
> dim(dat1)
[1] 230   5
> head(dat1)
                ID Location    Length Group Abundance
1 DE0655_HCMC_2001       HK 0.1786629   Yes 12.331055
2 MS0111_HCMC_1996       HK 0.2105236   Yes  9.652052
3 MS0063_HCMC_1995       HK 1.4337983   Yes 11.584822
4 DE0490_HCMC_2000       HK 0.3823731   Yes  7.893231
5 DE0885_HCMC_2001       HK 0.8478901   Yes 12.117232
6 DE0891_HCMC_2001       HK 1.5038646   Yes 10.819734
#加入第一圈数据，构建新的发育树，加入星星图，使用geom_star函数
##geom_fruit选项
#geom: 指定加入图的类型，如，geom_bar,geom_col,geom_boxplot,geom_violin,geom_tile、geom_star....
#mapping ：构建数据映射，注意x和y！！
#starstroke: 外框的宽度
#size：星星大小
#alpha: 透明度
#starshape ：指定形状，本例为1：星星；15：点
#show.legend = T：是否显示图例
#offset：图距树的距离

##guide对图例进行调整
#keywidth/keyheight ：图例的宽和长（整体包括字）
#direction：图例的排序，horizontal(水平)、vertical(垂直)
#title.positio：图例标题位置，left right top .....
#title： 图例标题
#order：图例处于第几个

p2 <- p1 + 
  geom_fruit(data = dat1,
            geom = geom_star,
            mapping = aes(y = ID,fill = Location, starshape = Group),
            position = "identity",
            starstroke = 0.2,
            size = 2) + 
  scale_starshape_manual(values = c(1,15),
                         guide  = guide_legend(keywidth = 0.5,keyheight = 0.5,order = 3,direction = "horizontal",title.position = "top")) + 
  scale_fill_npg(guide = guide_legend(keywidth = 0.5,keyheight = 0.5,order = 2,override.aes=list(starshape=15)))
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210606214145224.png" alt="image-20210606214145224" style="zoom:50%;" />

### 读入第二圈数据

```r
#读入第二圈数据
dat2 <- read.csv("first_ring_discrete.csv",header = T)
> dim(dat2)
[1] 1024    3
> head(dat2)
                ID Pos  Type
1 DE0846_HCMC_2001   8 type2
2 MS0034_HCMC_1995   8 type2
3 EG1017_HCMC_2009   6 type2
4   KH18_HCMC_2009   5 type2
5  10365_HCMC_2010   7 type2
6 EG1021_HCMC_2009   1 type1
#####构建第二圈，热图
p3 <-p2 +  new_scale_fill() + 
  geom_fruit(data = dat2,
             geom = geom_tile,
             mapping = aes(y = ID ,x = Pos,fill = Type),pwidth = 0.2,offset = 0.08) + 
  scale_fill_d3(guide = guide_legend(keywidth = 0.5,keyheight = 0.5,order = 1,direction = "horizontal",title.position = "top"))
p3
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210606214617470.png" alt="image-20210606214617470" style="zoom:50%;" />

### 读入第三圈数据

```r
dat3 <- read.csv("second_ring_continuous.csv",header = T)
> dim(dat3)
[1] 200   3
> head(dat3)
                ID Type2     Alpha
1 MS0004_HCMC_1995    p3 0.2256195
2 DE1150_HCMC_2002    p2 0.2222086
3 MS0048_HCMC_1995    p2 0.1881510
4  HUE57_HCMC_2010    p3 0.1939088
5 DE1486_HCMC_2002    p2 0.2018493
6 DE1165_HCMC_2002    p3 0.1812997

#构建第三圈数据
p4 <- p3 +
  new_scale_fill()+
  geom_fruit(data = dat3,
             geom = geom_tile,
             mapping = aes(y = ID,x = Type2,fill = Type2),pwidth = 0.1,offset = 0.05,axis.params = list(axis = "x",color = "black",text.angle = 45,hjust = 0.5,vjust = 1))+
  scale_fill_brewer(direction = -1,palette = 5,guide = guide_legend(keywidth = 0.5,keyheight = 0.5,order = 5,direction = "vertical",title.position = "top"))
p4

#使用axis.params设定坐标轴
#axis：:"x","y","xy"
#title: 坐标轴标题
#title.angle/size/height/color/angel: 字体角度/大小/距离/颜色/角度
#hjust/vjust: 调节字的位置
#nbreak：坐标轴分为几份
#line.size/color/alpha: 坐标轴的大小/颜色/透明度
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210607111633220.png" alt="image-20210607111633220" style="zoom:67%;" />

### 读入第四圈数据

```r
dat5 <- read.csv("dat5.csv",header = T)

> dim(dat5)
[1] 690   5
> head(dat5)
                ID Location    Length Group Abundance
1 DE0655_HCMC_2001       HK 0.1786629   Yes 12.331055
2 MS0111_HCMC_1996       HK 0.2105236   Yes  9.652052
3 MS0063_HCMC_1995       HK 1.4337983   Yes 11.584822
......

#构建第四圈数据，加入一个箱线图
p6 <- p5 + new_scale_fill() + geom_fruit(data = dat5,
                                        geom = geom_boxplot,
                                        mapping = aes(y= ID, x= Abundance,fill = Location),
                                        size = 0,
                                        grid.params = list()) + scale_fill_npg() 
```

![image-20210607155059022](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210607155059022.png)
