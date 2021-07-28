# 非约束排序-PCA

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-3-25
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)
```

```R
#工作路径设定
setwd("G:\\Desktop\\s_note\\data\\16s\\beta_diversity\\PCA")
rm(list = ls())
```

## 数据读入

**载入相关包**，**读入数据**

```R
library("vegan")
library("tidyverse")
gc()
#读入数据
lh_data <- as_tibble(read.csv("LH.csv",header = T,row.names = 1))
group <- as_tibble(read.csv("LH_group.csv",header = T))
names(group) <- "group"
```

## PCA分析简介

**PCA分析**
主成分分析（Principal Component Analysis，PCA）
是基于**特征向量**的线性无约束排序方法。能够将大量相关变量转化为一组很少的不相关变量，这些无关变量称为主成分（Principal Component，PC）可用于替代原始的大量相关变量，进而简化分析过程。

**数据特征**
对于非物种组成数据而言，如**环境数据，应对数据进行标准化后进行分析**
对于物种组成的数据来说，首先要确定**数据的分布要是均匀的和少零的**，因为，PCA本质还是表征了对象之间的**欧几里得距离**。
一般对于物种数据而言要首先进行**Hellinger转化标准化数据**作为pca的标准数据，但是对于物种数据而言，PcoA更加合适（Bray-curtis距离、Unifrac距离等）

**主成分轴的选择**
首选绘制排序图，观测对象是否在主要的排序轴（如前2-3轴）中具有明显的分布特征
**判断模型**：计算各个主成分贡献度的均值，选择大于这一均值的主成分即可，一般为前两轴
**-例外-**：如果数据在其他的轴上分布具有很强的规律性，并且这两个轴的解释变量没有特别低，此时是可以选择这两个轴的
**I型标尺与II型标尺**
	**I型标尺**：解释对象之间的关系，特征向量被标准化为单位长度，关注的是对象之间的关系。双序图中对象之间的距离近似于多维空间内的欧式距离，代表变量箭头之间的角度没有意义。
	**II型标尺**：解释变量和变量之间的关系，II型标尺可以根据一定的比例进行缩放，同样平衡贡献圈的半径也要进行相应的缩放

**平衡贡献圈**
可以用于评估变量的相对重要性。
$$
r=\sqrt(\frac{d}{p})
$$
d：所展示的PCA轴数（通常d = 2)
p：变量个数
平衡贡献圈的半径代表变量的向量长度对排序的平均贡献率，如果某个变量的箭头长度大于圆的半径，代表它对这个排序空间的贡献大于所有变量的平均贡献。

## PCA分析

 使用vegan包中的**`rda()`**函数进行PCA分析

```R
#scale=T表示对数据标准化
pca_lh <- rda(lh_data,scale = T)

#summary()函数查看结果
#Eigenvalue:每个轴的贡献度
#Proportion Explained：解释变量
#Cumulative Proportion：累积变量
#Species scores变量的排序坐标
#Site scores 对象排序坐标
> summary(pca_lh)
Call:
rda(X = lh_data, scale = T) 

Partitioning of correlations:
              Inertia Proportion
Total              13          1
Unconstrained      13          1

Eigenvalues, and their contribution to the correlations 

Importance of components:
                         PC1    PC2    PC3     PC4     PC5     PC6     PC7     PC8     PC9     PC10      PC11     PC12      PC13
Eigenvalue            4.8305 3.6780 1.5638 0.93934 0.61612 0.50126 0.36821 0.22123 0.14514 0.120945 0.0103937 0.005071 2.096e-06
Proportion Explained  0.3716 0.2829 0.1203 0.07226 0.04739 0.03856 0.02832 0.01702 0.01116 0.009303 0.0007995 0.000390 1.613e-07
Cumulative Proportion 0.3716 0.6545 0.7748 0.84705 0.89444 0.93300 0.96132 0.97834 0.98951 0.998810 0.9996098 1.000000 1.000e+00
```

## 主成分数量的判断

通过碎石图进行判断，使用**`pca_lh$CA$eig`**提取各个主成分的贡献度，利用**`plot()`**函数绘制碎石图判断主成分数量

```R
#type为绘图类型
plot(pca_lh$CA$eig,type = "b")
#添加均值线，可以得出前3个主成分对变量的解释贡献度较大
abline(h=mean(pca_lh$CA$eig),lwd=1,col="blue")
```

<img src="C:\Users\lzh233\AppData\Roaming\Typora\typora-user-images\image-20210325094352901.png" alt="image-20210325094352901" style="zoom:67%;" />

## 提取各个轴的解释变量

```R
PC1 <- paste((round((as.data.frame(pca_lh$CA$eig)[1,]/sum(as.data.frame(pca_lh$CA$eig))),4)*100),"%")
PC2 <- paste((round((as.data.frame(pca_lh$CA$eig)[2,]/sum(as.data.frame(pca_lh$CA$eig))),4)*100),"%")
PC3 <- paste((round((as.data.frame(pca_lh$CA$eig)[3,]/sum(as.data.frame(pca_lh$CA$eig))),4)*100),"%")
#查看各个轴的解释变量
> PC1
[1] "37.16 %"
> PC2
[1] "28.29 %"
> PC3
[1] "12.03 %"
```

## 排序坐标提取

提取前3轴的排序坐标

```R
#处理坐标提取，并合并分组信息
scaling_I <- as.data.frame(pca_I$sites)[,c(1:3)]
scaling_I <- data.frame(scaling_I,Label=group)
#变量坐标提取,并且合并变量信息
scaling_II <- as.data.frame(pca_II$species)[,c(1:3)]
scaling_II <- data.frame(scaling_II,vars=rownames(scaling_II))
```

## 计算平衡贡献圈的长度

```R
#平衡贡献圈半径，可利用这个筛选变量
r_pca <- round(((2/13)^0.5),4)
> r_pca
[1] 0.3922

#计算变量的半径,即箭头的长度（前两轴）
r <- as_tibble(((scaling_II$PC1)^2 + (scaling_II$PC2)^2)^0.5) %>% mutate(vars=rownames(scaling_II))

> head(r)
# A tibble: 6 x 2
  value vars       
  <dbl> <chr>      
1 0.836 Temperature
2 1.02  pH         
3 1.19  TOC        
4 1.09  IC         
5 1.15  TC         
6 1.17  TN   
```

## 结果可视化

```R
ggplot(scaling_I,mapping = aes(PC1,PC2)) + 
  geom_point(scaling_I,mapping = aes(shape = group,color = group)) +
  geom_segment(scaling_II,mapping = aes(x = 0,y = 0,xend = PC1,yend = PC2),
               arrow = arrow(length = unit(0.3, 'cm')), size = 0.3, color = 'blue',alpha = 0.5) +
  geom_vline(xintercept = 0,color = "gray") +
  geom_hline(yintercept = 0,color = "gray") +
  geom_text(scaling_II,mapping = aes(label = Label)) + 
  labs(x = paste("PC1:",PC1), y = paste("PC2:",PC2)) + 
  theme_bw() 
```

![image-20210325100823058](C:\Users\lzh233\AppData\Roaming\Typora\typora-user-images\image-20210325100823058.png)

