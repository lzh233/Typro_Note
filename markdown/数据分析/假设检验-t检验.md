# 假设检验：t-test

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-3-24
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)
```

```R
#工作路径设定
setwd("G:\\Desktop\\s_note\\data\\hypothesis_testing")
rm(list = ls())
```

## 数据读入

**载入相关包**，**读入数据**

```R
library("car")
gc()
t_data <- read.csv("ttest.csv",header = T)
```

## 正态检验

```R
qqPlot(lm(Value~Group,data = t_data), simulate = T, main = 'QQ Plot', labels = T)
#没有明显的离群点，说明数据符合正态分布
```

<img src="C:\Users\lzh233\AppData\Roaming\Typora\typora-user-images\image-20210324214703579.png" alt="image-20210324214703579" style="zoom:67%;" />

## 方差齐性分析

使用**`shapiro.test`**进行方差齐性分析

```R
shapiro <- tapply(t_data$Value, t_data$Group,shapiro.test)
shapiro
#p大于0.05说明方差相等
> shapiro
$F

	Shapiro-Wilk normality test

data:  X[[i]]
W = 0.96766, p-value = 0.7049


$Y

	Shapiro-Wilk normality test

data:  X[[i]]
W = 0.93901, p-value = 0.2296

```

## t检验

使用**`t.test()`**进行t检验

```R
#paired = F：指定是独立样本t检验 paired = T非独立样本
#altemative = "two.sided" ：指定双尾，less；greater分别为单尾
#var.equal= T ：假定方差相等，默认方差不一致，根据齐性检验确定
#conf.level = 0.95(默认)
t_test <- t.test(Value~Group,t_data,paired = F,altemative = "two.sided",var.equal= T,conf.level = 0.95)

#结果查看
> t_test

	Two Sample t-test

data:  Value by Group
t = 27.637, df = 38, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 19.38466 22.44890
sample estimates:
mean in group F mean in group Y 
       35.79433        14.87755 

#可以直接提取p-value
t_test[["p.value"]]

> t_test[["p.value"]]
[1] 9.069093e-27
```

