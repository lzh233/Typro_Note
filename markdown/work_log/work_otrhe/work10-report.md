## Q1 

```
```



数据基本情况

| group   | total_reads | total_UMI | max_reads | max_UMI | min_reads | min_UMI | mean_reads  | mean_UMI    | median_reads | median_UMI |
| ------- | ----------- | --------- | --------- | ------- | --------- | ------- | ----------- | ----------- | ------------ | ---------- |
| drug_ts | 2431134     | 99106     | 3765      | 75      | 1         | 1       | 259.7365385 | 10.58824786 | 151          | 8          |
| drug_zl | 3818519     | 3269340   | 1150      | 972     | 78        | 74      | 388.1792213 | 332.3513266 | 385          | 330        |
| pt_ts   | 614785      | 29401     | 4550      | 46      | 1         | 1       | 94.40801597 | 4.514895577 | 34           | 3          |
| pt_zl   | 3177427     | 2755937   | 1109      | 962     | 71        | 67      | 347.2597814 | 301.1953005 | 344          | 298        |

**UMI**

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210927140643868.png" alt="image-20210927140643868" style="zoom: 67%;" />

**reads**

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210927135047216.png" alt="image-20210927135047216" style="zoom:50%;" />

相关分析

**drug_S_H1975_203**

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210927143259500.png" alt="image-20210927143259500" style="zoom:50%;" />

**polyT_H1975**

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210927144102534.png" alt="image-20210927144102534" style="zoom: 50%;" />

```R

library(tidyverse)
#---Q1 -- ----------------------------------------------------------------------

#读如四组数据的counts

setwd("D:\\Desktop\\data\\Q1\\counts")
pt_zl <- read_tsv("polyT_H1975_ZL_counts.txt") %>% 
  filter(mark == "CB") %>% 
  select(c("Barcode","readcount","UMI")) %>% 
  mutate(group = rep("pt_zl",nrow(pt_zl)))

pt_ts <- read_tsv("polyT_H1975_TS_count.tsv") 

drug_zl <- read_tsv("drug_S_H1975_203_zl_counts.txt") %>% 
  filter(mark == "CB") %>% 
  select(c("Barcode","readcount","UMI"))%>% 
  mutate(group = rep("drug_zl",nrow(drug_zl)))

drug_ts <- read_tsv("drug_S_H1975_203_TS_count.tsv") 

#计算reads
pt_ts_df1 <- pt_ts %>% group_by(barcode) %>% summarise(reads = sum(read_count))
#计算UMI
pt_ts_df2 <- pt_ts %>% 
  group_by(barcode,UMI) %>% 
  count() %>% 
  group_by(barcode) %>% 
  summarise(umi = sum(n)) 

pt_ts <- pt_ts_df1 %>% left_join(pt_ts_df2,by = "barcode")
pt_ts <- pt_ts %>% mutate(group = rep("pt_ts",nrow(.)))


drug_ts_df1 <- drug_ts %>% group_by(barcode) %>% summarise(reads = sum(read_count))
drug_ts_df2 <- drug_ts %>% 
  group_by(barcode,UMI) %>% 
  count() %>% 
  group_by(barcode) %>% 
  summarise(umi = sum(n)) 

drug_ts <- drug_ts_df1 %>% left_join(drug_ts_df2,by = "barcode")
drug_ts <- drug_ts %>% mutate(group = rep("drug_ts",nrow(.)))

col.name <- c("barcode","reads","umi" ,"group")

rm(drug_ts_df1,drug_ts_df2,pt_ts_df1,pt_ts_df2)

names(drug_ts) <- col.name
names(pt_ts) <- col.name
names(drug_zl) <- col.name
names(pt_zl) <- col.name

write_tsv(drug_ts,"drug_ts.tsv")
write_tsv(drug_zl,"drug_zl.tsv")
write_tsv(pt_ts,"pt_ts.tsv")
write_tsv(pt_zl,"pt_zl.tsv")
# -- ----------------------------------------------------------------------

df_tmp <- list()
df_tmp[[1]] <- drug_ts
df_tmp[[2]] <- drug_zl
df_tmp[[3]] <- pt_ts
df_tmp[[4]] <- pt_zl

df <-  do.call('rbind',df_tmp)
write_tsv(df,"reads_umi.tsv")

#Violin
plot_theme <- 
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 13,color = "black"))

ggplot(df,aes(group,reads,fill = group,color = group)) + 
  geom_violin() + labs(x = "")  + plot_theme


ggplot(df,aes(group,umi,fill = group,color = group)) + 
  geom_violin() + labs(x = "") + plot_theme

df_summary <- df %>% 
  group_by(group) %>% summarise(
    total_reads = sum(reads),
    total_UMI = sum(umi),
    max_reads = max(reads),
    max_UMI = max(umi),
    min_reads = min(reads),
    min_UMI = min(umi),
    mean_reads = mean(reads),
    mean_UMI = mean(umi),
    median_reads = median(reads),
    median_UMI = median(umi)
  )

write_tsv(df_summary,"df_summary.tsv")

df1 <- df %>%filter(group == "drug_ts")
df2 <- df %>%filter(group == "drug_zl")

(df1$barcode %in% df2$barcode) %>% table()

df1 <- df1 %>% full_join(df2,by = "barcode")

lm.fix <- lm(umi.x~umi.y,data = df1)
r <- summary(lm.fix)
r_drug <- round(r$r.squared,4)

p1 <- ggplot(df1,aes(umi.x,umi.y)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(x = "UMI:TS", y = "UMI:ZL")+
  ggtitle("drug_S_H1975_203",subtitle = str_glue("r2 = {r_drug}"))+
  plot_theme

lm.fix <- lm(reads.x~reads.y,data = df1)
r2 <- summary(lm.fix)
r_drug_read <- round(r2$r.squared,4)

p2 <- ggplot(df1,aes(reads.x,reads.y)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(x = "reads:TS", y = "reads:ZL")+
  ggtitle("drug_S_H1975_203",subtitle = str_glue("r2 = {r_drug_read}"))+
  plot_theme
cowplot::plot_grid(p1,p2)


# -- ----------------------------------------------------------------------



df1 <- df %>%filter(group == "pt_ts")
df2 <- df %>%filter(group == "pt_zl")

df1 <- df1 %>% full_join(df2,by = "barcode")

lm.fix <- lm(umi=.x~umi.y,data = df1)
r <- summary(lm.fix)
r_drug <- round(r$r.squared,4)

p1 <- ggplot(df1,aes(umi.x,umi.y)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(x = "UMI:TS", y = "UMI:ZL")+
  ggtitle("polyT_H1975",subtitle = str_glue("r2 = {r_drug}"))+
  plot_theme

lm.fix <- lm(reads.x~reads.y,data = df1)
r2 <- summary(lm.fix)
r_drug_read <- round(r2$r.squared,4)

p2 <- ggplot(df1,aes(reads.x,reads.y)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(x = "reads:TS", y = "reads:ZL")+
  ggtitle("polyT_H1975",subtitle = str_glue("r2 = {r_drug_read}"))+
  plot_theme
cowplot::plot_grid(p1,p2)
```

