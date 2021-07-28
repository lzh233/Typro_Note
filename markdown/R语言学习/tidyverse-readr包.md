# 数据读入(readr包)

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-6-14
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)
```

## 载入相关包

```R
library(tidyverse)
setwd("G:\\Desktop\\s_note\\data\\data2clean")
```

## 使用readr读入数据

### readr的优点

`reader`包的优点，速度快要比R基础包的函数读取速度快10倍，在读取大文件时有很大的优势，读取的数据会转换为`tibble`不会使用行名称也不会更改列名称，也不会将字符向量自动转换为因子

### 基本数据的读入

使用以下函数可以实现对不同分隔符的数据进行读入

`read_csv()`：读入以逗号分隔的`csv`文件

`read_csv2()`：读入以分号`;`分隔的文件

`read_tsv()`：读入以制表符`\t`为分隔符号的文件

`read_delim()`：读入**任意分隔符**的文件，与`read.delim()`不同，后者默认读入文件分隔符为制表符

```R
#基本的文件读入
df <- read_csv("otutable.csv")
#读入后readr_csv函数会返回数据表的总结信息，如下，包含了数据中每一列的数据类型
-- Column specification ---------
cols(
  OTU_ID = col_character(),
  CK_1 = col_double(),
  CK_2 = col_double(),
  CK_3 = col_double(),
  CK_4 = col_double(),
  CK_6 = col_double(),
  CK_7 = col_double(),
  CK_8 = col_double(),
......

```

**常见的选项参数**(适用于一些系列的readr函数)

- `skip = n`: 可以通过该选项来指定跳过前几行
- `comment = "#"`: 可以用于指定跳过以某些字符开头的行, 如，读入数据时候跳过以`#`开头的行
- `skip_empty_rows=TURE`: 跳过数据框中的空格
- `col_name = TURE`: 指定第一行为列名，相当于`header = TURE`, 也可以指定一组与列数相等的向量来作为列名
- `na = `: 表示文件的**某个符号**应作为`NA`处理

```R
#跳过前两行读取数据，数据中有#开头的行，因此使用两种方法跳过
df2 <- read_csv("otutable2.csv",skip = 2)
df2 <- read_csv("otutable2.csv",comment = "#")
#均读入成功
#读入后readr_csv函数会返回数据表的总结信息，如下，包含了数据中每一列的数据类型
-- Column specification ---------
cols(
  OTU_ID = col_character(),
  CK_1 = col_double(),
  CK_2 = col_double(),
 .....
    
#使用readr_csv()构建数据，并指定列名
df3 <- read_csv("1,2,3\n4,5,6",col_names = c("A","B","C"))
df3
> df3
# A tibble: 2 x 3
      A     B     C
  <dbl> <dbl> <dbl>
1     1     2     3
2     4     5     6

#如将点替换为NA
df3 <- read_csv("1,2,3\n4,5,.",col_names = c("A","B","C"),na = ".")
> df3
# A tibble: 2 x 3
      A     B     C
  <dbl> <dbl> <dbl>
1     1     2     3
2     4     5    NA
```

### readr读取数据的原理与过程

### parse_*族函数



