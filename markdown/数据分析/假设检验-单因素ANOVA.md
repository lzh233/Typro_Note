# 假设检验-单因素ANOVA

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
ano <- read.csv("anova.csv",header = T)
```

## 正态检验

使用**`qqPlot()`**函数进行正态检验，同t检验

```R
qqPlot(lm(Value~Group,ano))
```

<img src="C:\Users\lzh233\AppData\Roaming\Typora\typora-user-images\image-20210324221055306.png" alt="image-20210324221055306" style="zoom:67%;" />

## 方差齐性检验

使用**`bartlett.test()`**进行方差齐性检验

```R
bartlett.test(Value~Group,ano)

#p>0.05方差齐性
> bartlett.test(Value~Group,ano)

	Bartlett test of homogeneity of variances

data:  Value by Group
Bartlett's K-squared = 1.7958, df = 2, p-value = 0.4074
```

## 方差分析

使用**`aov()`**函数进行方差分析,**`summary()`**函数查看结果

```R
anova <- aov(Value~Group,ano)
summary(anova)
#可以看出整体数据差异显著
> summary(anova)
            Df Sum Sq Mean Sq F value  Pr(>F)   
Group        2 1.2121  0.6061   8.452 0.00859 **
Residuals    9 0.6453  0.0717                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

## 多重比较

以`duncan`和`HSD`为例子,**`agricolae::duncan.test()`**函数进行duncan法事后比较，**`TukeyHSD()`**函数使用HSD法事后比较

```R
#duncan法
duncan <- agricolae::duncan.test(anova,"Group",alpha=0.05)

#比较结果查看
> duncan[["groups"]]
      Value groups
F  4.605126      a
NS 4.127301      b
Y  3.833950      b
#均值等结果
> duncan[["means"]]
      Value       std r      Min      Max      Q25      Q50      Q75
F  4.605126 0.3566586 4 4.269298 4.918298 4.310681 4.616454 4.910898
NS 4.127301 0.1489123 4 3.963860 4.324476 4.063662 4.110434 4.174073
Y  3.833950 0.2563868 4 3.551964 4.076358 3.650419 3.853739 4.037270

#HSD法
turk <- TukeyHSD(anova,conf.level = 0.95)
#查看比较结果
> turk$Group
           diff        lwr         upr       p adj
NS-F -0.4778248 -1.0064835  0.05083394 0.075912736
Y-F  -0.7711758 -1.2998346 -0.24251710 0.007062042
Y-NS -0.2933510 -0.8220098  0.23530770 0.314952790
#可视化比较结果,超过虚线则无差异
plot(turk)
```

<img src="C:\Users\lzh233\AppData\Roaming\Typora\typora-user-images\image-20210324222047075.png" alt="image-20210324222047075" style="zoom:67%;" />