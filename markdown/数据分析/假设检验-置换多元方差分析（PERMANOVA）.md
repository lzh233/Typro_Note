# 置换多元方差分析（PERMANOVA）

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-3-31
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)
version: Seurat_4.0.1
```

```R
#设定工作目录
setwd("G:\\Desktop\\s_note\\data\\hypothesis_testing")
```

## 数据读入载入包

```R
set.seed(123)
library("vegan")
library("ggplot2")
otu <- read.delim("PERMANOVA_otu_table.txt",sep = "\t",header = T,row.names = 1)
group <- read.delim("PERMANOVA_group.txt",sep = "\t",header = T)
otu <- as.data.frame(t(otu))
```

## 计算距离矩阵

```R
#计算距离矩阵
bray <- vegdist(otu,method = "bray")
```

## PERMANOVA分析

### 全局比较

使用**`adonis()`**函数进行**PERMANOVA分析**

```R
#PERMANOVA分析(全局分析)
adonis <- adonis(bray~site,group,permutations = 999)

#在site水平上整体差异显著
#Df，自由度，其值=所比较的分组数量-1
#SumsOfSqs，即Sums of squares，总方差，又称离差平方和
#MeanSqs，即Mean squares，均方（差）
#F.Model，F检验值
#R2，即Variation (R2)，
#方差贡献，表示不同分组对样品差异的解释度，即分组方差与总方差的比值
#p值，默认p<0.05即存在显著差异
> adonis
Call:
adonis(formula = bray ~ site, data = group, permutations = 999) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
site       4   0.40572 0.101430  4.9804 0.36273  0.001 ***
Residuals 35   0.71280 0.020366         0.63727           
Total     39   1.11852                  1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### 组间两两比较

```R
#流程：根据分组名称提取otu表信息→→→计算距离矩阵→→→PERMANOVA分析并保存结果
#提取处理信息
group_name <- unique(group$site)
#循环数确定
a <- length(group_name) - 1
b <- length(group_name) + 0
#设定结果储存变量
result <- NULL
#两两比较
for (i in 1:a){
  for (j in (i+1):b) {
    print(paste("开始比较:",group_name[i],"/",group_name[j],"√"))
    #得到分组信息（两个方法）
    group_ij <- subset(group,site %in% c(group_name[i],group_name[j]))
    #group_ij <- subset(group,site == group_name[i] | site==group_name[j])
    #根据分组信息提取otu表格信息
    otu_ij <- otu[group_ij$names,]
    #计算距离矩阵
    bray_ij <- vegdist(otu_ij,method = "bray")
    #保存矩阵
    bray_tmp <- as.matrix(bray_ij)
    write.csv(bray_tmp,paste(group_name[i],"&",group_name[j],".csv"))
    #adonis分析
    result_ij <- adonis(bray_ij~site,group_ij,permutations = 9999)
    result <- rbind(result,c(paste(group_name[i],group_name[j],sep = "/"),
                             "bc-distance",
                             unlist(data.frame(result_ij$aov.tab, check.names = FALSE)[1, ])))
  }
}
#输出结果
    V1          V2 Df          SumsOfSqs            MeanSqs          F.Model                 R2 Pr(>F)       p.adj
1  s1/s2 bc-distance  1 0.0191970029355519 0.0191970029355519 1.51974196714995 0.0979231465553183 0.1552 0.172444444
2  s1/s3 bc-distance  1 0.0811658139120037 0.0811658139120037 4.93221327369211  0.260519634043307 0.0013 0.002166667
3  s1/s4 bc-distance  1  0.174591586154343  0.174591586154343 8.44104613075124  0.376143165589075  4e-04 0.000800000
```

## p值校正与结果保存

```R
result <- as.data.frame(result)
result <- cbind(p.adj=result,p.adjust(result$"Pr(>F)",method = "BH"))
write.csv(result,"PERMANOVA_result.csv")
```

