# 可视化-boxplot

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
library("tidyverse")
```

## 绘制热图

```R
setwd("")
library("ggplot2")
tem <- read.csv("Temperature.csv",header = T)
names(tem) <- c("Sample_id","group","Tempreture")
ggplot(tem,mapping = aes(group,Tempreture,color = group)) +
  geom_boxplot() + scale_fill_manual(values = c("red","gray","blue","pink"))
################
#增加均值等信息fun：mean/median等均可
#+ stat_summary(fun=mean, geom="point", shape=23, size=4)
#增加抖动图，0.2为抖动程度，越大越稀疏
#+ geom_jitter(shape=16, position=position_jitter(0.2))
#+ geom_dotplot(binaxis = "y",stackdir = "center",position = "jitter",dotsize = 0.4)
#横向图
#+ coord_flip()
#横轴的顺序&只有几组
#+ scale_x_discrete(limits=c("T3","T2","CK","T1"))
#+ scale_x_discrete(limits=c("T3","T2"))
###################
#箱线图参数
# geom_boxplot()
#开口箱线图：geom_boxplot(notch=T)
#单色填充
#geom_boxplot(fill='#A4A4A4', color="black")
#异常值的设置：
#outlier.colour="red", outlier.shape=8, outlier.size=4
#手工颜色设置(填充和边框色)：
#+ scale_fill_manual(values = c("red","gray","blue","pink"))
#+ scale_color_manual(values = c("red","gray","blue","pink"))
#添加误差线
#+ stat_boxplot(geom = "errorbar" , width = 0.1)
#修改箱子宽度
# + geom_boxplot(width = 0.1)
```

