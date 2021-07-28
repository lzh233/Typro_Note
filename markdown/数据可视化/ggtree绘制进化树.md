# ggtree绘制进化树

## R版本与运行环境信息

```R
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

## 设置工作目录以及输入文件查看

```R
setwd("G:\\Desktop\\s_note\\data\\ggtree")
dir()
#tree_group.txt为每个节点对应的分组与值
> dir()
[1] "tree.nwk"       "tree_group.txt"
```

## 绘制进化树

### 读入数据

**`node.label = "support"`**: 表示节点可以增加其他信息

**`split()`**: 按照分组因子，把向量，矩阵和数据框进行适当的分组。它的返回值是一个列表，代表分组变量每个水平的观测。

```R
#读入数据
tree <- read.newick("tree.nwk",node.label = "support")
group <- read.delim("tree_group.txt",header = T)
#构建节点分组信息
groupinf <- split(group$Sample,group$Group)
#将节点分组信息加入进化树中
tree <- groupOTU(tree,.node = groupinf)
```

### 绘制初步进化树

```R
p1 <- ggtree(tree,branch.length = "none",layout = "circular",linetype = 1,size=1,ladderize = F,aes(color=group))
p1
```

![image-20210414224604650](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210414224604650.png)

0