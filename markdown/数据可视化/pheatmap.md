# pheatmap绘制热图

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
library("pheatmap")
library("Hmisc")
```

## 绘制热图

```R
#填上路径
setwd("")
#填文件名
data <- as.matrix(read.csv("data_normalized.csv",header = T,row.names = 1))
pheatmap(data,cluster_rows = F,
         cluster_cols = F)
#相关系数热图
cor_num <- cor(data)
#保存相关系数文件
write.csv(cor_num,"cor_data.csv")
pheatmap(cor_num,cluster_rows = F,
         cluster_cols = F)


##############################################
#cor()函数method选项：格式method=""
#pearson/kendall/spearman三个相关系数，默认pearson
#pheatmap选项
#fontsize = 字体大小
#border = #边界大小或者存在情况
#cellwidth = #方框高度
#cellheight = #方框宽度
#display = #是否显示数值，可以显示显著性
#number_color = "black" #数值颜色
#treeheight_row = 列树高
#treeheight_col = 行树高
#分组
#annotation_col =annotation_row 跟上分组文件（列名必须指定）
#annotation_row =annotation_row 跟上分组文件（列名必须指定）
##############################################


#显著性检验
#Hmisc包，rcorr包分析，得出p-vaule
mysor <- rcorr(as.matrix(data),type = "spearman")
#相关系数
R <- mysor[["r"]]
P <- mysor[["P"]]
P <- as.vector(P)
#校正p-value,校正后的p值应该转换为原来的矩阵格式
q_value <- p.adjust(P,"fdr")
q_value <- as.matrix(q_value)
```

