# 非约束排序-CA和DCA

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
setwd("G:\\Desktop\\s_note\\data\\16s\\beta_diversity\\CA_DCA")
rm(list = ls())
```

## 数据读入

**载入相关包**，**读入数据**

```R
library("vegan")
library("tidyverse")
gc()
#读入数据
#读取门水平丰度表top10
phy <- read.csv("phy2.csv",header = T,row.names = 1)
#分组读取
group <- read.csv("group.csv",header = T)
```

## CA和DCA分析简介

**可用于判断采用PCA或CA（DCA）分析**
CA分析即对应分析（Correspondence Analysis，CA），在群落分析中常用于物种组成数据
排序结果展示的是样方间χ2距离。χ2距离不受双零问题的影响，适用于原始物种多度数据（无需像PCA那样要进行数据转换），CA分析对象**必须是频度或类频度、同量纲的非负数据，如物种个体计数或二元数据。**

**CA分析的过程**

原始计算方法：迭代计算
1、给每个样方给予随机得分
2、物种的得分：每一个样方中物种丰度为加权计算平均得分
3、样方得分：样方中每个物种丰度为加权计算平均得分
CA轴的选择与判断与PCA相似，CA排序轴承载的总变差不是用总方差来表示，而是通过一个叫总惯量（total inertia）的指标表征

I型标尺与II型标尺排序图
**I型标尺**
PCA类似，样方与样方之间的关系，样方得分的计算方法是样方中出现的物种得分的平均值，并通过物种丰度加权
物种通常分布在样方范围之外的原因
排序图内样方之间的距离近似于它们的**χ2距离**
1、排序图中两个样方点越近，代表这些样方内的物种组成越相似
2、一个样方点靠近一个物种点，表示该物种对于该样方的贡献比较大

II型标尺
PCA类似，物种与物种之间的关系，物种得分的计算方式是出现该物种的所有样方得分的平均值，并按出现该物种的所有样方中该物种丰度加权
样方通常分布在物种范围之外
1、两个物种点越近，代表它们的相对多度沿样方分布越相似
2、一个物种点靠近一个样方点，表示该物种在该样方内存在的可能性很大
靠近（0,0）的物种则在整个生态系统中分布均匀

**DCA分析**
去趋势对应分析，解决了CA分析的弓型效应，DCA还可用于估算梯度长度（gradient length），上述提到DCA将CA轴重标定为物种更替标准差单位（以**SD**表示），可用于表征环境异质性或物种β多样性特征。
**SD≤3可表明物种沿排序轴更有可能是线性分布适合PCA分析，SD≥4表明物种沿排序轴更有可能是单峰分布，适合DCA/CA分析**

## CA分析

使用**`cca()`**函数进行CA分析

```R
phy_ca <- cca(phy)
#Eigenvalues:各个CA轴的贡献度，CA排序轴承载的总变差不是用总方差来表示，而是通过一个叫总惯量（total inertia）的指标表征，即Inertia
>phy_ca
Call: cca(X = phy)

              Inertia Rank
Total          0.6311     
Unconstrained  0.6311   10
Inertia is scaled Chi-square 

Eigenvalues for unconstrained axes:
    CA1     CA2     CA3     CA4     CA5     CA6     CA7     CA8     CA9    CA10 
0.28884 0.23923 0.04914 0.03957 0.00662 0.00468 0.00202 0.00051 0.00037 0.00011
```

## 轴数目的确定

```R
#查看轴数目,与PCA类似，本例前两轴
plot(phy_ca$CA$eig,type = "b")
abline(h = mean(phy_ca$CA$eig))
```

<img src="C:\Users\lzh233\AppData\Roaming\Typora\typora-user-images\image-20210325162422270.png" alt="image-20210325162422270" style="zoom:67%;" />

## I和II型标尺的提取

```R
#提取I型标尺
cca_I <- summary(phy_ca,scaling = 1)
#提取II型标尺
cca_II <- summary(phy_ca,scaling = 2)

#提取I型标尺的前两轴
plot_I <- as.data.frame(cca_I$sites[,c(1:2)])
```

## 贡献度的计算

```R
CA1 <- paste((round((as.data.frame(phy_ca$CA$eig)[1,]),4)*100),"%")
CA2 <- paste((round((as.data.frame(phy_ca$CA$eig)[2,]),4)*100),"%")

> CA1
[1] "28.88 %"
> CA2
[1] "23.92 %"
```

## I型标尺的可视化

```R
#合并分组信息
plot_I <- cbind(plot_I,group$gr)
names(plot_I) <- c("CA1","CA2","Group")
#绘图（I型标尺）
ggplot(plot_I,mapping = aes(CA1,CA2,color = Group,shape = Group)) + geom_point() + theme_classic()
```

<img src="C:\Users\lzh233\AppData\Roaming\Typora\typora-user-images\image-20210325162954936.png" alt="image-20210325162954936" style="zoom:67%;" />

## II型标尺的提取与可视化

```R
#提取II型标尺前两轴并合并门水平信息
plot_II <- as.data.frame(cca_II$species[,c(1,2)])
plot_II <- cbind(plot_II,rownames(plot_II))
names(plot_II) <- c("CA1","CA2","Phylum")
#作图（II型标尺）
ggplot(plot_II,mapping = aes(CA1,CA2)) + 
  geom_point(color = "red") + 
  geom_text(aes(label = Phylum)) + 
  theme_classic()
```

<img src="C:\Users\lzh233\AppData\Roaming\Typora\typora-user-images\image-20210325163352411.png" alt="image-20210325163352411" style="zoom:67%;" />

## 双序图

```R
#双序图,注意数据不要放在主函数中！
ggplot() + 
  geom_point(plot_I,mapping = aes(CA1,CA2,color = Group,shape = Group)) + 
  geom_point(plot_II,size = 2,mapping = aes(CA1,CA2)) + 
  geom_text(plot_II,mapping = aes(CA1,CA2,label = Phylum)) + 
  theme_classic()
```

<img src="C:\Users\lzh233\AppData\Roaming\Typora\typora-user-images\image-20210325163625620.png" alt="image-20210325163625620" style="zoom:67%;" />

## DCA分析(去趋势分析)

DCA分析可以用于判断物种沿样方分布的模型，单峰型还是线性

```R
#读入数据与分组信息
dca <- as.data.frame(t(read.csv("otutable.csv",header = T,row.names = 1)))
group_dca <- read.csv("group_dca.csv",header = T)

> dca[1:4,5:10]
     1000757 1000876 100100 1001013 100104 1001564
CK_1       0       1      0       0      0       0
CK_2       0       1      0       0      0       0
CK_3       1       1      0       0      0       0
CK_4       0       0      0       0      0       0
```

## 进行CA分析并可视化

```R
#CA分析，弓型效应,样品点过于聚集
a <- cca(dca)
#I型标尺提取
b <- summary(a,scaling = 1)
#提取前两轴
c <- as.data.frame(b$sites[,c(1,2)])
#合并分组信息
c <- cbind(c,group_dca)
names(c) <- c("CA1","CA2","G")
ggplot(c,mapping = aes(CA1,CA2)) + geom_point() + geom_text(aes(label = G))
#产生了弓型效应,样品点过于聚集
```

<img src="C:\Users\lzh233\AppData\Roaming\Typora\typora-user-images\image-20210325164224879.png" alt="image-20210325164224879" style="zoom:67%;" />

## 消除弓形效应

使用**`decorana()`**函数消除弓形效应，进行 DCA分析

```R
dca_otu <- decorana(dca)
```

## 轴的提取与分析方法判断

```R
#判断物种分布模型，根据Axis lengths（物种更替标准差单位（standard deviation units of species turnover，简称SD））判断
#SD<=3 线性分布，更适合PCA
#介于3，4二者皆可
#SD>=4 单峰分布，CA
#提取DCA1和DCA2（DCA分析一般提取Axis lengths最大的前两轴）
#提取前两轴信息
dca_I <- summary(dca_otu,scaling = 1)
dca_II <- summary(dca_otu,scaling = 2)

> dca_I <- summary(dca_otu,scaling = 1)

Call:
decorana(veg = dca) 

Detrended correspondence analysis with 26 segments.
Rescaling of axes with 4 iterations.

                  DCA1    DCA2    DCA3    DCA4
Eigenvalues     0.1118 0.03742 0.03224 0.01794
Decorana values 0.1121 0.05665 0.02614 0.01653
Axis lengths    1.2726 0.94871 0.78649 0.75094
```

## 可视化

```R
#提取前两轴并合并分组信息
plot_dca1 <- as.data.frame(dca_I$site.scores[,c(1,2)])
plot_dca2 <- as.data.frame(dca_II$spec.scores[,c(1,2)])
names(plot_dca1) <- c("DCA1","DCA2","G")
#可视化
ggplot(plot_dca1,mapping = aes(DCA1,DCA2)) + 
  geom_point(aes(color = G)) + 
  geom_text(aes(label = G,color = G)) + 
  theme_classic()
```

<img src="C:\Users\lzh233\AppData\Roaming\Typora\typora-user-images\image-20210325165141966.png" alt="image-20210325165141966" style="zoom:67%;" />