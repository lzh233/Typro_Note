# 约束排序-RDA分析

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-1-30
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)
```

```R
#设定工作目录
setwd("G:\\Desktop\\s_note\\data\\16s\\beta_diversity\\RDA")
```

## 数据读入载入包

```R
set.seed(123)
library("vegan")
library("ggplot2")
```

## 读入数据

```R
otu <- read.csv("phylum_table.csv",header = T,row.names = 1)
group <- read.csv("group.csv",header = T,row.names = 1)
env_table <- read.delim("env_table.txt",sep = "\t",row.names = 1)
```

## RDA分析简介

**冗余分析**
RDA是响应变量矩阵与解释变量矩阵之间多元多重线性回归的拟合值矩阵的PCA分析
**RDA算法简述**
Y-中心化后的响应变量矩阵，X-中心化后的解释变量矩阵，如环境变量矩阵
1、通过多重线性回归对矩阵Y和X进行拟合，得到两个矩阵，分别为拟合值矩阵A和残差值矩阵B
**残差：拟合值与真实值之间的差值**
2、对矩阵A进行PCA分析，得到典范特征向量（eigenvectors）矩阵U使用矩阵U计算两套样方的排序得分

## 数据转换

```R
#由于counts数会有很多0值，因此进行数据转换,hellinger默认按行
otu_tb <- decostand(otu,method = "hellinger",MARGIN = 1)
```

## RDA分析

使用**`rda()`**函数进行RDA分析，需要进一步校正与变量筛选

```R
#得到初步的rda模型.~代表使用全部环境因子
rda_ori <- rda(otu_tb~.,env_table,scale = F)
#查看初步分析的结果
summary(rda_ori)
#可以看出模型解释了0.7338的方差(r2)，02662的没有解释，RDA1解释了0.53445的方差....
> summary(rda_ori)

Call:
rda(formula = otu_tb ~ pH + TC + DOC + SOM + TN + NO3 + NH4 +      AP + AK, data = env_table, scale = F) 

Partitioning of variance:
              Inertia Proportion
Total         0.02686     1.0000
Constrained   0.01971     0.7338
Unconstrained 0.00715     0.2662

Eigenvalues, and their contribution to the variance 

Importance of components:
                         RDA1     RDA2     RDA3      RDA4      RDA5      RDA6      RDA7      RDA8
Eigenvalue            0.01435 0.003181 0.001172 0.0004967 0.0002044 0.0001612 7.278e-05 4.314e-05
Proportion Explained  0.53445 0.118465 0.043633 0.0184940 0.0076120 0.0060011 2.710e-03 1.606e-03
Cumulative Proportion 0.53445 0.652913 0.696546 0.7150401 0.7226520 0.7286531 7.314e-01 7.330e-01
```

## r2的校正

```R
#查看并校正r2
r <- RsquareAdj(rda_ori)
r_adj <- r$adj.r.squared
r_adj
#校正后r2降低
> r_adj
[1] 0.5928396
```

## 变量的选择

**选择变量，剔除对模型贡献度较小以及共线性明显的变量，重新构建模型**
***方差膨胀因子***分析查看变量vif值,如果存在vif较为大的变量则需要后续对某些变量剔除**（删除离群的vif值）**
并非直接删除vif较大的的变量，https://zhuanlan.zhihu.com/p/56468729

``vif``：**对因子之间的线性相关关系进行检验**,得出的值，为了检验因子之间的线性相关关系，我们可以通过OLS对单一因子和解释因子进行回归，然后如果其R^2较小，说明此因子被其他因子解释程度较低，线性相关程度较低。**VIF越高解释变量和因变量之间线性相关性就越强**
$$
VIF=\frac{1}{1-R^2}
$$

使用**`vif.cca()`**函数计算各个变量的VIF值，可以结合后续分析去剔除共线性严重以及对模型贡献不大的变量

```R
sort(vif.cca(rda_ori))
> sort(vif.cca(rda_ori))
      NO3        pH       NH4       DOC        TN        AP        TC        AK       SOM 
 1.777526  2.995094  3.109630 18.924057 20.626384 26.854763 45.040644 46.166302 68.852428
#可以通过boxplot()查找离群值对利群变量进行删除
boxplot(vif)$out
```

使用**`ordiR2step()`**对变量执行选择，分为前向选择、后向选择以及逐步分析**（一般步骤）**

（1）首先，运行包含所有解释变量的全模型置换检验，当且仅当置换检验显示显著性后再执行变量的前向选择。

（2）根据解释变量的高低对所有的变量进行降序排序

（3）运行置换检验检验各个解释变量的显著性，如果显著则将该变量加入到模型中，否则终止选择

（4）将此前已经选择的解释变量作为协变量，执行偏RDA分析，计算剩余未被选择的各解释变量（每个单独）解释的变差

（5）根据解释变差的高低（该变差代表了每个变量的局部效应，partial effect），对这些解释变量再次降序排序后，选择其中的最佳解释变量（降序排序后排在第一位的）并检验其显著性。若显著则将该解释变量添加在模型中，否则终止选择。

（6）若步骤（6）中的置换检验结果显著，则继续执行步骤（5）和（6），直到终止选择为止（即剩余解释变量中，最佳变量所解释的变差不显著为止）。

由上述步骤可知，变量的显著性是选择停止的规则之一，即如果加入新变量的偏RDA置换检验显著性p值不显著，选择过程即被终止。然而，这个标准过于宽松，有时会选择显著但不包含任何变量的模型（夸大I类错误），或选择包括过多变量的模型（夸大被解释方差的量）。Blanchet等（2008）关注这两个问题并提出解决方案，使选择停止规则达到全模型的校正后R2为准：

为了防止夸大I类错误 ，首先运行包含所有解释变量的全模型置换检验，当且仅当置换检验显示显著性后，再执行变量的前向选择。

为了减少纳入太多变量的风险，首先计算包含所有解释变量的全模型的R2adj，将R2adj作为第二个终止原则。如果备选变量的偏RDA置换检验不显著或当前模型的R2adj超过全模型的R2adj，前向选择即被终止。

```R
#指定rda分析的对象？
#scope：变量范围
#R2scope：指定r2阈值，若改变某变量超过该阈值则停止筛选
#direction：指定变量选择方法c("both", "backward", "forward")
#permutations置换检验的次数
#Pin：解释变量的p值小于该值则加入模型
#Pout：解释变量p值大于某值则删除该变量
rda_adj <- ordiR2step(object=rda(otu_tb~1, env_table, scale = F), 
                      scope = formula(rda_ori), 
                      R2scope = r_adj, 
                      direction = 'both', 
                      permutations = 9999,
                      Pin = 0.05)

Step: R2.adj= 0 
Call: otu_tb ~ 1 
 
                 R2.adjusted
<All variables>  0.592839620
+ pH             0.132370285
+ DOC            0.072061204
+ NO3            0.060902955
+ AP             0.060410543
+ NH4            0.046565232
+ AK             0.040344442
+ SOM            0.004455678
+ TC             0.002502013
<none>           0.000000000
+ TN            -0.006905803

     Df     AIC      F Pr(>F)   
+ pH  1 -99.578 4.9667 0.0098 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Step: R2.adj= 0.1323703 
Call: otu_tb ~ pH 
 
                R2.adjusted
<All variables>   0.5928396
+ AK              0.2807891
+ AP              0.2512868
+ NH4             0.2237985
+ DOC             0.2091295
+ TN              0.1542857
+ SOM             0.1534109
+ TC              0.1443998
<none>            0.1323703
+ NO3             0.1267164

     Df     AIC      F Pr(>F)   
+ AK  1 -103.75 6.1591 0.0036 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Step: R2.adj= 0.2807891 
Call: otu_tb ~ pH + AK 
 
                R2.adjusted
<All variables>   0.5928396
+ DOC             0.3503766
+ AP              0.3309766
+ SOM             0.3282698
+ TC              0.3202571
+ TN              0.3131210
+ NH4             0.2924521
+ NO3             0.2829532
<none>            0.2807891

      Df     AIC      F Pr(>F)  
+ DOC  1 -105.64 3.5709 0.0259 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Step: R2.adj= 0.3503766 
Call: otu_tb ~ pH + AK + DOC 
 
                R2.adjusted
<All variables>   0.5928396
+ SOM             0.4663563
+ TN              0.4615550
+ TC              0.4423237
+ NH4             0.4175384
+ AP              0.3554377
<none>            0.3503766
+ NO3             0.3309934

      Df     AIC      F Pr(>F)   
+ SOM  1 -110.15 5.9987 0.0028 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Step: R2.adj= 0.4663563 
Call: otu_tb ~ pH + AK + DOC + SOM 
 
                R2.adjusted
<All variables>   0.5928396
+ NH4             0.5228360
+ AP              0.4987266
+ TN              0.4685682
<none>            0.4663563
+ NO3             0.4506937
+ TC              0.4483916

      Df     AIC     F Pr(>F)  
+ NH4  1 -112.43 3.604 0.0207 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Step: R2.adj= 0.522836 
Call: otu_tb ~ pH + AK + DOC + SOM + NH4 
 
                R2.adjusted
<All variables>   0.5928396
+ AP              0.5814848
+ NO3             0.5360746
<none>            0.5228360
+ TN              0.5111369
+ TC              0.5091855

     Df     AIC      F Pr(>F)  
+ AP  1 -115.29 3.9428 0.0131 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Step: R2.adj= 0.5814848 
Call: otu_tb ~ pH + AK + DOC + SOM + NH4 + AP 
 
                R2.adjusted
+ NO3             0.6153910
<All variables>   0.5928396
<none>            0.5814848
+ TC              0.5704847
+ TN              0.5684003
```

```R
#查看剔除某些变量后校正后的rda模型,只有6个变量最终被选择
summary(rda_adj)
sum.rda <- summary(rda_adj)
Biplot scores for constraining variables


        RDA1    RDA2    RDA3     RDA4      RDA5     RDA6
pH   0.50787  0.4921  0.4381 -0.02456  0.448065 -0.32657
AK  -0.05717  0.7382  0.5521  0.31174  0.215610  0.05716
DOC  0.15568  0.8823  0.4036  0.18293 -0.009949 -0.03054
SOM  0.14416  0.1426  0.8530  0.32282  0.334745  0.12254
NH4  0.13381 -0.7516 -0.4022  0.47888  0.152667 -0.05228
AP  -0.01394  0.8968  0.3127  0.18401  0.184616  0.17260
#查看新的vif值，比之前要好一些，可以考虑删除AK，本例不删除
vif.cca(rda_adj)
> vif.cca(rda_adj)
       pH        AK       DOC       SOM       NH4        AP 
 2.815855 43.866758 15.810348  7.418528  2.297437 24.547635
```

## 校正优化后的RDA模型

```R
#对优化后的模型r2进行校正,可以看出进行变量选择后r2与之前随略减小（肯定减少）且相差不大
r2 <- RsquareAdj(rda_adj)

>r2$adj.r.squared
[1] 0.5814848
#校正前的r2
> r_adj
[1] 0.5928396
```

## 约束轴检验与显著性检验

置换检验检验约束轴和变量的显著性，使用**`anova()`**函数（跟自带的方差分析函数冲突可以用**`anova.cca()`**代替）

```R
#对约束轴进行检验（整体检验）
rda_adj_test <- anova(rda_adj,permutations = 999)
#模型整体显著
> rda_adj_test
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = otu_tb ~ pH + AK + DOC + SOM + NH4 + AP, data = env_table, scale = F)
         Df  Variance      F Pr(>F)    
Model     6 0.0182101 7.0207  0.001 ***
Residual 20 0.0086458                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#对约束轴逐一检验
rda_adj_test_axis <- anova(rda_adj,by = "axis",permutations = 999)
#p值矫正
rda_adj_test_axis$`Pr(>F)` <- p.adjust(rda_adj_test_axis$`Pr(>F)`,method = "fdr")
#RDA1和RDA2显著，因此后续分析可以选择这两轴，即这两轴可以解释的变量显著高于其他轴
> rda_adj_test_axis
Permutation test for rda under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

Model: rda(formula = otu_tb ~ pH + AK + DOC + SOM + NH4 + AP, data = env_table, scale = F)
         Df  Variance       F Pr(>F)   
pH        1 0.0044511 10.2965 0.0060 **
AK        1 0.0045755 10.5843 0.0060 **
DOC       1 0.0023961  5.5427 0.0105 * 
SOM       1 0.0033066  7.6489 0.0060 **
NH4       1 0.0017763  4.1091 0.0240 * 
AP        1 0.0017045  3.9428 0.0240 * 
Residual 20 0.0086458                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#对变量逐一检验
rda_adj_test_var <- anova(rda_adj,by = "term",permutations = 999)
#p值矫正
rda_adj_test_var$`Pr(>F)` <-  p.adjust(rda_adj_test_var$`Pr(>F)`,method = "fdr")
#检验对模型显著相关的变量，p越小则代表该变量对模型的贡献越大
> rda_adj_test_var
Permutation test for rda under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

Model: rda(formula = otu_tb ~ pH + AK + DOC + SOM + NH4 + AP, data = env_table, scale = F)
         Df  Variance       F Pr(>F)   
pH        1 0.0044511 10.2965 0.0060 **
AK        1 0.0045755 10.5843 0.0060 **
DOC       1 0.0023961  5.5427 0.0120 * 
SOM       1 0.0033066  7.6489 0.0060 **
NH4       1 0.0017763  4.1091 0.0240 * 
AP        1 0.0017045  3.9428 0.0204 * 
Residual 20 0.0086458                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

## 根据R2对解释变量进行矫正

```R
#解释变量校正，本例取前两轴
#r2_adj :校正后的r2
#r2：校正前的r2
exp <- ((sum.rda$cont$importance[,c(1:2)][2,])*r2_adj)/r2
rda1 <- exp[1]*100
rda2 <- exp[2]*100

> rda1
    RDA1 
42.90547 
> rda2
    RDA2 
9.796678
```

## 提取排序结果

使用**`scores()`**函数提取结果

### 提取样方排序结果

```R
weight_site <- scores(rda_adj, choices = 1:2, scaling = 1, display = 'wa')
write.csv(weight_site,"weight_site.csv")
```

### 提取物种排序图

```R
#提取物种排序
species <- sum.rda$species[,c(1:2)]
write.csv(species,"species_scores.csv")
```

### 解释变量得分提取

```R
#解释变量得分提取
#explain <- sum.rda$biplot[,c(1:2)]
explain <- scores(rda_adj, choices = 1:2, scaling = 1, display = 'bp')
write.csv(explain,"explain_score.csv")
```

### 回归系数提取

```R
#回归系数的提取(斜率):前两轴
k <- coef(rda_adj)[,c(1:2)]
write.csv(k,"rda_coef.csv")
```

### r2提取

```R
write.csv(RsquareAdj(rda_adj),"r2.csv")
```

### 约束轴显著性检验结果保存

```R
#约束轴校验结果保存
write.csv(rda_adj_test,"rda_adj_test.csv")
write.csv(rda_adj_test_axis,"rda_adj_test_axis.csv")
```

### 残差特征值提取

未解释的方差的信息

```R
###提取残差特征值###
pca_eig <- rda_adj$CA$eig
pca_eig_axis <- pca_eig[pca_eig>mean(pca_eig)]
barplot(pca_eig)
abline(h=mean(pca_eig))
```

## 可视化

```R
#ggplot2可视化
#作图文件，分组信息
plot_file <- data.frame(weight_site,group)
plot_explain_file <- data.frame(explain)

ggplot(plot_file,aes(RDA1,RDA2)) + 
  geom_point(aes(color = group),size = 2) +
  geom_segment(plot_explain_file,mapping = aes(x = 0,y = 0,xend = RDA1,yend = RDA2),
               arrow = arrow(length = unit(0.2,"cm")), size = 0.3, color = 'blue') + 
  geom_text(plot_explain_file,mapping = aes(label = rownames(plot_explain_file),vjust = "inward",hjust = "inward")) +
  stat_ellipse(aes(fill = group),type = "norm", geom = "polygon",alpha= 0.3) + 
  geom_vline(xintercept = 0,color = "gray")+
  geom_hline(yintercept = 0,color = "gray")+
  labs(x=paste("RDA1:",round(exp[1]*100,2),"%"),y=paste("RDA2:",round(exp[2]*100,2),"%")) + 
  theme_bw()
```

<img src="C:\Users\lzh233\AppData\Roaming\Typora\typora-user-images\image-20210331095707185.png" alt="image-20210331095707185" style="zoom:67%;" />