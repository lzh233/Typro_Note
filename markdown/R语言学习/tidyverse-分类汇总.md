# 数据表分类汇总

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-4-20
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)
```

## 载入相关包

```R
library(tidyverse)
library(nycflights13)
```

## 数据查看

```r
#载入测试数据
> flights
# dep_time arr_time ：实际出发与到达时间
#sched_dep_time, sched_arr_time：预计到达与出发时间
#dep_delay, arr_delay：出发与到达的延迟/提前时间
......
# A tibble: 336,776 x 19
    year month   day dep_time sched_dep_time dep_delay arr_time sched_arr_time arr_delay carrier flight tailnum origin
   <int> <int> <int>    <int>          <int>     <dbl>    <int>          <int>     <dbl> <chr>    <int> <chr>   <chr> 
 1  2013     1     1      517            515         2      830            819        11 UA        1545 N14228  EWR   
 2  2013     1     1      533            529         4      850            830        20 UA        1714 N24211  LGA   
......
```

## 使用summarise对数据进行汇总计算

```r
###使用summarize汇总数据
flights %>% 
  group_by(year,month) %>% 
  summarise(dep_delay=mean(dep_delay,na.rm = T),
            dep_sum=sum(dep_delay,na.rm = T),
            time_spend=mean((arr_time - dep_time) %% 60,na.rm = T)) %>% head()
# A tibble: 6 x 5
# Groups:   year [1]
   year month dep_delay dep_sum time_spend
  <int> <int>     <dbl>   <dbl>      <dbl>
1  2013     1      10.0    10.0       29.7
2  2013     2      10.8    10.8       29.7
3  2013     3      13.2    13.2       30.2
4  2013     4      13.9    13.9       30.3
5  2013     5      13.0    13.0       31.2
6  2013     6      20.8    20.8       30.7
```





