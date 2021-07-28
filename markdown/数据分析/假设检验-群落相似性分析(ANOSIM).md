# ANOSIM群落相似性分析

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
otu <- read.delim("anosim_otu_table.txt",sep = "\t",header = T,row.names = 1)
group <- read.delim("anosim_group.txt",sep = "\t",header = T)
otu <- as.data.frame(t(otu))
```

## 计算距离矩阵

```R
#计算距离矩阵
bray <- vegdist(otu,method = "bray")
```

## ANOSIM群落相似性分析

### 全局比较

使用**`anosim()`**进行分析

```R
anosim_res <- anosim(bray,group$site,permutations = 999)
summary(anosim_res)
#statistic R:r>0,表示组间差异大于组内差异，反之组内差异大于组间差异，其绝对值越大差异越大
#Significance:p值
> summary(anosim_res)

Call:
anosim(x = bray, grouping = group$site, permutations = 999) 
Dissimilarity: bray 

ANOSIM statistic R: 0.3483 
      Significance: 0.001 

Permutation: free
Number of permutations: 999

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99% 
0.0548 0.0724 0.0933 0.1123 

Dissimilarity ranks between and within classes:
        0%    25%   50%    75% 100%   N
Between  8 231.75 415.5 607.25  780 640
s1       1  17.75 166.5 227.00  360  28
s2       2  60.50 145.5 257.25  578  28
s3      17 100.75 238.5 438.25  761  28
s4      26 290.75 462.0 595.50  725  28
s5       4 123.00 437.5 544.25  680  28
```

### 两两比较

```R
#得到分组名称
group_name <- unique(group$site)
#创建一个空数据集用于存储结果
anosim_result_two <- NULL
#循环次数的确定
a <- length(group_name)-1
b <- length(group_name) + 0
#组间两两比较
for (i in 1:a){
  for (j in (i+1):b){
    print(paste("开始比较:",group_name[i]," & ",group_name[j],"√"))
    #获得分组
    group_ij <- subset(group,site %in% c(group_name[i],group_name[j]))
    #group_ij <- subset(group,site == group_name[i] | site == group_name[j])
    #根据分组提取otu数据,按行提取
    otu_ij <- otu[group_ij$names,]
    #计算bray距离
    bray_ij <- vegdist(otu_ij,method = "bray")
    #保存矩阵
    bray <- as.matrix(bray_ij)
    write.csv(bray,paste(group_name[i],"&",group_name[j],".csv"))
    #比较组件差异,方法1可以省略计算距离的代码
    #anosim_ij <- anosim(otu_ij,group_ij$site,permutations = 999,distance = "bray")
    anosim_ij <- anosim(bray_ij,group_ij$site,permutations = 999)
    #保存每一次的结果，表格的形式
    anosim_result_two <- rbind(anosim_result_two,
                               c(paste(group_name[i],group_name[j],sep = "/"),
                                 "bray",
                                 anosim_ij$signif,
                                 anosim_ij$statistic))
  }
}
```

## p值校正与结果保存

```R
#保存比较结果
anosim_result_two <- as.data.frame(anosim_result_two)
names(anosim_result_two) <- c("sample","distance","p_value","R")
anosim_result_two <- cbind(anosim_result_two,p.adj=p.adjust(anosim_result_two$p,method = "BH"))
write.csv(anosim_result_two,"result_between.csv")

> anosim_result_two
   sample distance p_value                  R      p.adj
1   s1/s2     bray   0.134          0.0859375 0.14888889
2   s1/s3     bray   0.001          0.3828125 0.00200000
3   s1/s4     bray   0.003  0.598214285714286 0.00500000
4   s1/s5     bray   0.001  0.534040178571429 0.00200000
5   s2/s3     bray   0.001  0.473772321428571 0.00200000
6   s2/s4     bray   0.001  0.586495535714286 0.00200000
7   s2/s5     bray   0.001  0.565848214285714 0.00200000
8   s3/s4     bray    0.07         0.13671875 0.08750000
9   s3/s5     bray   0.041  0.181919642857143 0.05857143
10  s4/s5     bray   0.201 0.0491071428571428 0.20100000
```

