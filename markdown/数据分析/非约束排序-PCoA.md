# 非约束排序-PCoA

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-3-28
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)
```

```R
#工作路径设定
setwd("G:\\Desktop\\s_note\\data\\16s\\beta_diversity\\PCOA")
rm(list = ls())
```

## 数据读入

**载入相关包**，**读入数据**

```R
library("vegan")
library("tidyverse")
gc()
#读入数据
qu_otu <- read.csv("quhua_otu.csv",header = T,row.names = 1,encoding = "UTF-8")
group <- read.csv("quhua_group.csv",header = T)
```

## PCoA分析简介

是一种基于距离的非度量多维度分析方法，与PCA、CA（DCA）不同，后者主要是基于原始数据所进行的分析但是，本本质上还是基于欧氏距离和卡方距离的进行降维分析。若PCoA输入的距离矩阵是欧式距离矩阵，则会与PCA分析结果一致,卡方距离CA也类似(但是CA会考虑到物种丰度，会略有不同)
**算法简述**
（1）基于变量组成，计算x个对象间的距离矩阵，计算时根据实际情况选择合适的距离测度。
（2）按行和列中心化矩阵。
（3）对中心化后的距离矩阵执行特征分解。
（4）新的特征向量即表征了主坐标。
**轴：**
与PCA CA类似但是PCoA分析的过程可能会出现负特征根，但是如果只在后几个轴则可以忽略，如果负值绝对值较大，则应该进行数据转换(可能会丢失差异，慎重)：
将非欧式距离测度转化为欧式几何性质，
    Bray-curtis距离是在群落分析中典型的非欧式距离测度，可对其平方根转化，即可使它获得了欧式距离属性；或者直接在距离的平方基础上加个常数（Lingoes校正）；或是距离本身直接加个常数（Cailliez校正）等。
**缺点是PCoA会使物种信息丢失，但是物种变量可以通过和PCoA轴进行相关分析来投影到PCoA排序图中，进而表明某一物种对该样方的贡献程度。**

## 距离计算

```R
#计算样品间距离：bray-crutis
qu_dis <- as.matrix(vegdist(t(qu_otu),method = "bray"))
```

## PCoA分析

```R
#pcoa分析(默认只有前两轴，参数2指定轴的数目,一般为样品数减1（'k'必需是{1, 2, ..  n - 1}中的一个），eig为特征值是否输出，F则只输出排序结果)
pcoa <- cmdscale(qu_dis,k = nrow(qu_dis)-1,eig = T)
```

## 特征值查看与解释变量计算

```R
#特征值图查看
pcoa$eig
plot(pcoa$eig,type = "b")
abline(h = mean(pcoa$eig),col = "blue")
abline(h = 0,col = "gray")
```

<img src="C:\Users\lzh233\AppData\Roaming\Typora\typora-user-images\image-20210328190830820.png" alt="image-20210328190830820" style="zoom:67%;" />

```R
#解释变量计算
pcoa_eig <- as.data.frame(pcoa$eig)
pcoa1 <- round((pcoa_eig[1,]/sum(pcoa_eig)*100),2)
pcoa2 <- round((pcoa_eig[2,]/sum(pcoa_eig)*100),2)
pcoa3 <- round((pcoa_eig[3,]/sum(pcoa_eig)*100),2)
pcoa4 <- round((pcoa_eig[4,]/sum(pcoa_eig)*100),2)

> pcoa1
[1] 65.79
> pcoa2
[1] 20.63
```

## 特征轴提取

```R
pcoa_plot <- as.data.frame(pcoa$points[,c(1:4)])
names(pcoa_plot) <- c("PCoA1","PCoA2","PCoA3","PCoA4")
```

## 合并分组信息并可视化

```R
pcoa_plot <- cbind(pcoa_plot,group)
names(pcoa_plot) <- c("PCoA1","PCoA2","PCoA3","PCoA4","TR","GR")
#可视化
ggplot(pcoa_plot,mapping = aes(PCoA1,PCoA2,color = TR)) + geom_point(size = 3) + 
  geom_vline(xintercept = 0,color = "gray") +
  geom_hline(yintercept = 0,color = "gray") +
  labs(x = paste("PCoA1",pcoa1,"%"),y = paste("PCoA1",pcoa2,"%")) + 
  theme_bw()
```

<img src="C:\Users\lzh233\AppData\Roaming\Typora\typora-user-images\image-20210328200447168.png" alt="image-20210328200447168" style="zoom:67%;" />