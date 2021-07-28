# stringr

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-5-18
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)
```

## 载入相关包

```R
library("tidyverse")
```

## stringr函数简要列表

简要归纳为6个大功能：**拆、连、取、删、替、筛**

## R语言正则表达式

### 正则表达式

| 表达式 | 含义                                                         |
| :----: | ------------------------------------------------------------ |
|   \n   | 换行符（unix）                                               |
|   \r   | 回车符                                                       |
|   \t   | 制表符                                                       |
|   \f   | 换页符                                                       |
|  [\b]  | 退格符                                                       |
|  \r\n  | 换行符（windows）                                            |
|   .    | 可以匹配任何单个的字符字母数字甚至.字符本身。同一个正则表达式允许使用多个.字符。但不能匹配换行符 |
|   *    | 表示*前的字符出现至少0次，*.*组合代表任意字符                |
|   +    | 表示+前的字符至少出现一次                                    |
|   \    | 转义符                                                       |
|   ^    | 表示开头，如 ^a                                              |
|   $    | 表示结尾, 如a$                                               |
|   ()   | 表示匹配括号内的字符串/表达式                                |
|   []   | 表示匹配[]中的任意字符如，[a-z] [0-9] [AaCcSs]               |
|   {}   | {}前的字符出现次数，{5}表示5次，{1,3}表示出现1到3次，{2,} 出现两次以上 |
|   ?    | 匹配零个或一个字符                                           |
|   \|   | 逻辑运算符“或”                                               |
|  [^]   | 表示括号里的匹配内容取反                                     |
|   \d   | 匹配任意一个数字，使用时应该在加一个"\", 下同                |
|   \D   | 匹配任意一个非数字                                           |
|   \w   | 匹配任意一个小写字母（等价于[a-zA-Z0-9]）                    |
|   \W   | 匹配任意一个大写字母（等价于`[^a-zA-Z0-9]`）                 |
|   \s   | 任何一个空白字符（等价于[\f\n\r\t\v]）                       |
|   \S   | 任何一个非空白字符（等价于`[^\f\n\r\t\v]`）                  |

### POSIX字符类

| [:alnum:]  | 任何一个字母或数字（等价于[a-ZA-Z0-9]）                      |
| ---------- | ------------------------------------------------------------ |
| [:alpha:]  | 任何一个字母（等价于[a-ZA-Z]）                               |
| [:blank:]  | 空格或制表符（等价于[\t ]）  注:t后面有一个空格              |
| [:cntrl:]  | ASCII控制字符（ASCII 0到31，再加上ASCII 127）                |
| [:digit:]  | 任何一个数字（等价于[0-9])                                   |
| [:graph:]  | 和[:print:]一样，但不包括空格                                |
| [:lower:]  | 任何一个小写字母（等价于[a-z])                               |
| [:print:]  | 任何一个可打印字符                                           |
| [:punct:]  | 既不属于[:alnum:]，也不属于[:cntrl:]的任何一个字符           |
| [:space:]  | 任何一个空格字符，包括空格（等价于[f\n\r\t\v ] 注:v后面有一个空格 |
| [:upper:]  | 任何一个大写字母（等价于[A-Z])                               |
| [:xdigit:] | 任何一个十六进制数字(等价于[a-fA-F0-9])                      |

## string函数使用方法

- **`str_c()`**: 用于连接字符串成一个向量（paste/paste0），要连接的字符串在开头则在在被连接的前面，反之则在后方

  选项：

  `seq=""`: 指定连接的分隔符  

  `collapse =""`: 将多个字符串连接成一个向量时的连接符

```r
#例如生成一系列处理
str_c("T",c(1:8),sep = "_")
#输出结果
[1] "T_1" "T_2" "T_3" "T_4" "T_5" "T_6" "T_7" "T_8"
#另一种连接方式
> str_c("T_",c(1:8))
[1] "T_1" "T_2" "T_3" "T_4" "T_5" "T_6" "T_7" "T_8"
#同样可以指定多个字符串同时连接，例如生成一系列带有重复的处理
> str_c("T",c(1:8),"_1")
[1] "T1_1" "T2_1" "T3_1" "T4_1" "T5_1" "T6_1" "T7_1" "T8_1"
#例如生成一组处理名
> a <- str_c("T",c(1:3),"_1")
> b <- str_c("T",c(1:3),"_2")
> c <- str_c("T",c(1:3),"_3")
> data.frame(Treatment = c(a,b,c))
  Treatment
1      T1_1
2      T2_1
3      T3_1
4      T1_2
5      T2_2
6      T3_2
7      T1_3
8      T2_3
9      T3_3
```

使用`collapse =""`将多个字符串连接成一个向量

```r
> str_c(letters, collapse = "::")
[1] "a::b::c::d::e::f::g::h::i::j::k::l::m::n::o::p::q::r::s::t::u::v::w::x::y::z"
```

使用负索引可以对数据循环连接

```r
> str_c(letters[-26], " comes before ", letters[-1])
 [1] "a comes before b" "b comes before c" "c comes before d" "d comes before e" "e comes before f"
 [6] "f comes before g" "g comes before h" "h comes before i" "i comes before j" "j comes before k"
[11] "k comes before l" "l comes before m" "m comes before n" "n comes before o" "o comes before p"
......
```

- **`str_conv()`**: 修改字符串编码格式

```r
str_conv(string, encoding)
```

- **`str_count()`**: 计算某一字符出现的频次，**支持正则表达式**

```r
fruit <- c("apple", "banana", "pear", "pineapple")
fruit2 <- c("aaa","abb")
#计算"a"在fruit中出现的频次
> str_count(fruit, "a")
[1] 1 3 1 1
#计算"a"在fruit和fruit2中出现的频率
> str_count(c(fruit,fruit2), "a")
[1] 1 3 1 1 3 1
#统计指定字符在指定元素中的出现频次
> str_count(fruit, c("a", "b", "p", "p"))
[1] 1 1 1 3
#统计向量中每个元素的长度
> str_count(c("abc", "bac", "cba1a"), ".")
[1] 3 3 5

###正则表达式简单举例
#统计"b"在开头的频次
> str_count(fruit, pattern = "^b")
[1] 0 1 0 0
#统计"a"在结尾的频次
> str_count(fruit, pattern = "a$")
[1] 0 1 0 0
#统计每个元素中，"a"或"b"或"e"出现的频次
> str_count(fruit, pattern = "[a,b,e]")
[1] 2 4 2 3
```

- **`str_detect()`**: 判断字符串中是否含有某元素，**支持正则**，与`grep()`类似，但是后者默认显示index，可以用于筛选，**`negate = TRUE`代表取反**

```r
#查看"p"是在每个元素中是否存在
> str_detect(fruit,"p")
[1]  TRUE FALSE  TRUE  TRUE
#与grep做个比较
> grep("^p",fruit )
[1] 3 4
#对结果取反，即判断p是不是不在各个元素中
> str_detect(fruit,"p",negate = T)
[1] FALSE  TRUE FALSE FALSE
```

- **`str_dup()`**: 将指定字符串重复，与`rep()`函数类似，`str_dup(string, times)`

```r
fruit <- c("apple", "pear", "banana")
#将每个元素重读3次
> str_dup(fruit,3)
[1] "appleappleapple"    "pearpearpear"       "bananabananabanana"
#每个元素分别重复，2次 4次 6次
> str_dup(fruit,c(2,4,6))
[1] "appleapple"                           "pearpearpearpear"                    
[3] "bananabananabananabananabananabanana"
#某字符串分别重复0，1，2，3，4，5次
> str_dup("na", 0:5)
[1] ""           "na"         "nana"       "nanana"     "nananana"   "nanananana"
```

- **`str_start()/str_ends()`**: 判断某字符是否在字符串的开头/结尾，类似于`str_detect(str,pattern="^str")/str_detect(str,pattern="str$")` ,`negate = TRUE`代表取反

```r
fruit <- c("apple", "pear", "banana")
#判断是每个元素是否是以“a”开头
> str_starts(fruit,"a")
[1]  TRUE FALSE FALSE
#判断是否每个元素是以“a“结尾
> str_ends(fruit,"a")
[1] FALSE FALSE  TRUE
#取反
> str_starts(fruit,"a",negate = T)
[1] FALSE  TRUE  TRUE
#简单正则判断
> str_starts(fruit,"[a,b]",negate = T)
[1] FALSE  TRUE FALSE
```

- **`str_flatten()`**: 将字符串连接，`collapse =""`可以指定连接符

```
> str_flatten(letters,collapse = "_")
[1] "a_b_c_d_e_f_g_h_i_j_k_l_m_n_o_p_q_r_s_t_u_v_w_x_y_z"
```

- **`str_glue()`** /**` str_glue_data()`**:  字符串的参数传递，**类似于函数的形参和实参的传递**

```r
#字符串的参数传递
> name <- "bill"
> age <- "50"
> anniversary <- as.Date("1991-10-12")
> str_glue(
+   "my name is {name} " ,
+   "my age next {age}",
+   "my data is {format(anniversary, '%A, %B %d, %Y')}"
+ )
my name is bill my age next 50my data is 星期六, 十月 12, 1991

#从数据框中进行参数传递
> mtcars %>% str_glue_data("{rownames(.)} has {hp} hp") %>% head()
Mazda RX4 has 110 hp
Mazda RX4 Wag has 110 hp
Datsun 710 has 93 hp
Hornet 4 Drive has 110 hp
Hornet Sportabout has 175 hp
Valiant has 105 hp
```

- **`str_length()`**: 统计字符串的长度，或向量中每个元素的长度

```r
#统计单一字符串的长度
> leter <- "abcdaaa"
> str_length(leter)
[1] 7
#统计向量中每个元素的字符串长度，如 统计mtcars中每种车的名字长度
> str_length(rownames(mtcars))
 [1]  9 13 10 14 17  7 10  9  8  8  9 10 10 11 18 19 17  8 11 14 13 16 11 10 16  9 13 12 14 12 13 10
```

- **`str_locate()`** / **`str_locate_all()`**: 匹配指定字符串**第一次/全部**出现的位置，支持正则

```R
fruit <- c("apple", "banana", "pear", "pineapple")
#匹配字母“ap”第一次出现的位置
>  str_locate(fruit,"ap")
     start end
[1,]     1   2
[2,]    NA  NA
[3,]    NA  NA
[4,]     5   6
#匹配特殊字符，如结尾$, end列为每个元素结尾的位置 If the match is of length 0, (e.g. from a special match like $) end will be one character less than start.
>  str_locate(fruit,"$")
     start end
[1,]     6   5
[2,]     7   6
[3,]     5   4
[4,]    10   9
#统计“a”出现的全部位置，返回对象为list
>  str_locate_all(fruit,"a")
[[1]]
     start end
[1,]     1   1

[[2]]
     start end
[1,]     2   2
[2,]     4   4
[3,]     6   6

[[3]]
     start end
[1,]     3   3

[[4]]
     start end
[1,]     5   5
```

- **`str_order()`** / **`str_sort()`**: 对字符串排序，前者返回值是排序后的**index值**，后者返回**实际值**, 二者参数一致，

  `decreasing=`: 默认FALSE从大到小，反之从小到大

  `na_last=`: 默认`TURE`, 即NA排在最后，`FALSE`, NA排在开头，` NA` 丢弃NA

  `locale=`: 指定语言习惯， 默认en（英语）即可

  `numeric=`: 是否按照字符串中的数字排序字符串，默认`FALSE`

```r
#以排列处理为例子
> a <- str_c("T",c(1:3),"_1")
> b <- str_c("T",c(1:3),"_2")
> c <- str_c("T",c(1:3),"_3")
> d <- c(a,b,c)
> d
[1] "T1_1" "T2_1" "T3_1" "T1_2" "T2_2" "T3_2" "T1_3" "T2_3" "T3_3"
> str_sort(d)
[1] "T1_1" "T1_2" "T1_3" "T2_1" "T2_2" "T2_3" "T3_1" "T3_2" "T3_3"
#按照字符串中的数字排序字符串
> x <- c("100a10", "100a5", "2b", "2a")
> str_sort(x, numeric = TRUE)
[1] "2a"     "2b"     "100a5"  "100a10"
```

- **`str_pad()`**: 字符串补齐功能，选项

  `width=`: 指定补齐的长度

  `side=c("left", "right", "both")`: 指定补齐的字符位于哪里，左中右

  `pad=""`: 指定补齐的字符串时使用的符号

```r
> str_pad(c("a", "abc", "abcdef"), width = 11,side = "both",pad = "+")
[1] "+++++a+++++" "++++abc++++" "++abcdef+++"
```

- **`str_trunc()`**: 字符串截齐

  `width=`: 指定截齐的长度

  `side=c("left", "right", "both")`: 指定截齐方式，左中右

  `ellipsis=` : 指定被截去字符的替代符，默认`···`

```r
> x <- "This string is moderately long"
> rbind(
+     str_trunc(x, 20, "right"),
+     str_trunc(x, 20, "left"),
+     str_trunc(x, 20, "center",ellipsis = "++++")
+ )
     [,1]                  
[1,] "This string is mo..."
[2,] "...s moderately long"
[3,] "This str++++ely long"
```
- **`str_remove()`**/**`str_remove_all()`**： 删除字符串中的指定字符（**第一次出现/全部**），**支持正则表达式**

```r
> fruits <- c("one apple", "two pears", "three bananas")
#删除元音字母（仅仅第一次出现）
> str_remove(fruits, "[aeiou]")
[1] "ne apple"     "tw pears"     "thre bananas"
#删除所有处理结尾表示重复的数字与下划线
> d
[1] "T1_1" "T2_1" "T3_1" "T1_2" "T2_2" "T3_2" "T1_3" "T2_3" "T3_3"
> str_remove_all(d,"_[1-3]$")
[1] "T1" "T2" "T3" "T1" "T2" "T3" "T1" "T2" "T3"
```

- **`str_replace()/str_replace_all()/str_replace_na()`**： **第一次/全部出现**的指定字符串的替换，支持正则，最后一个指将NA由数值转换为**普通字符串**

```r
#将第一次出现的元音字母替换为__
> str_replace(fruits,pattern = "[aeiou]",replacement = "___")
[1] "___ne apple"     "tw___ pears"     "thr___e bananas
#将重复编号全部转换为a，b，c
> str_replace_all(d,c("1$" = "a","2$" = "b" ,"3$" = "c"))
[1] "T1_a" "T2_a" "T3_a" "T1_b" "T2_b" "T3_b" "T1_c" "T2_c" "T3_c"
#指定字符转换为大/小写
> str_replace_all(fruits, "[aeiou]", toupper)
[1] "OnE ApplE"     "twO pEArs"     "thrEE bAnAnAs"
> str_replace_all(fruits, "[AEIOU]", tolower)
[1] "one apple"     "two pears"     "three bananas"
#NA转换为普通的字符串"NA"
> str_replace_na(c(NA, "abc", "def"))
[1] "NA"  "abc" "de f"
```

- **`str_lower()/str_upper()/str_title()/str_to_sentence ()`**: 字符串大小写转换

```r
#大小写转换
> str_to_upper("i")
[1] "I"
> str_to_lower("I")
[1] "i"
#所有首字母大小写转换
> str_to_title(dog)
[1] "The Quick Brown Dog"
#句首字母大写
> str_to_sentence("the quick brown dog")
[1] "The quick brown dog"
```

- **`str_split()/str_split_fixed()`**: 分割字符串，默认返回值类型为**列表/矩阵**

  `simplify=FALSE`: 默认FALSE，返回值为列表（适用于`str_split()`）,TURE返回值为矩阵

  `n`: 将字符串分割成几部分，当n大于最大可以分割的数量时，会自动以空字符串补齐

```r
> fruits
[1] "apples and oranges and pears and bananas"

#将处理以“_”分割后可直接接索引
> str_split(d, "_",simplify = T,n = 2)
      [,1] [,2]
 [1,] "T1" "1" 
 [2,] "T2" "1" 
 [3,] "T3" "1" 
 [4,] "T1" "2" 
 [5,] "T2" "2" 
......
#simplify = F 返回列表
> str_split(fruits,"and",simplify = F)
[[1]]
[1] "apples "   " oranges " " pears "   " bananas" 
#str_split_fixed()可以直接返回矩阵
> str_split_fixed(fruits,"and",n = 3)
     [,1]      [,2]        [,3]                
[1,] "apples " " oranges " " pears and bananas"
```

- **`str_trim()/str_squish()`**: 删除字符串指定位置（开头、结尾、开头和结尾的多余的空格）/删除所有多余的空格

  `side =c("left","right","both")`: 分别指定开头、结尾、开头和结尾的多余的空格

```r
#指定位置删除空格
> str_trim("  String with trailing and leading white space        asda",side = "left")
[1] "String with trailing and leading white space        asda"
> str_trim("  String with trailing and leading white space        asda    ",side = "both")
[1] "String with trailing and leading white space        asda"
> str_trim("  String with trailing and leading white space        asda    ",side = "right")
[1] "  String with trailing and leading white space        asda"
#删除全部多余空格
> str_squish(c("  String   with trailing and leading white space        asda    ","   sa   adad  asda     adad  ada "))
[1] "String with trailing and leading white space asda"
[2] "sa adad asda adad ada"  
```

- **`str_sub()/str_sub() <- `**: 根据**索引提取/替换**字符串

`start= end=`: 指定要提取的字符串开头和结尾的index

`omit_na=`: 默认FALSE, 是否忽略NA

```r
hw <- "Hadley Wickham"
#提取前6个字符
> str_sub(hw,1,6)
[1] "Hadley"
#仅指定start/end提取
##从第二个字符开始取
> str_sub(hw,start = 2)
[1] "adley Wickham"
##取到第二个字符
> str_sub(hw,end  = 2)
[1] "Ha"
#使用符索引进行提取
##提取倒数第四到最后一个字符
> str_sub(hw,-4,-1)
[1] "kham"
#如果输入字符串为向量，则对每个元素均提取
> a
[1] "a"    "ab"   "abc"  "abcd"
> str_sub(a,1,2)
[1] "a"  "ab" "ab" "ab"
#可以直接进行赋值替换
> str_sub(hw,1,6) <- "REPLACE"
[1] "REPLACE Wicreplace"
```

- **`str_subset()/str_which()`**: 分别包含指定字符的**字符串/index**， 支持正则表达式以及`negate`选项

```r
> fruit
[1] "apple"    "banana"   "pear"     "pinapple"
> str_subset(fruit,"p")
[1] "apple"    "pear"     "pinapple"
> str_which(fruit,"p")
[1] 1 3 4
#显示除了NA外所有的字符/index
> str_subset(c("a", NA, "b"), ".")
[1] "a" "b"
> str_which(c("a", NA, "b"), ".")
[1] 1 3
```

- **`str_extract()/str_extract_all()`**: 用于提取字**符串中符合条件**的部分（第一次/全部）

```r
 shopping_list <- c("apples x4", "bag of flour", "bag of sugar", "milk x2")
#返回含有a或p的部分，仅仅第一次出现
> str_extract(shopping_list, "[ap]")
[1] "a" "a" "a" NA 
#返回全部符合条件的部分
> str_extract_all(shopping_list, "[ap]",simplify = T)
     [,1] [,2] [,3]
[1,] "a"  "p"  "p" 
[2,] "a"  ""   ""  
[3,] "a"  "a"  ""  
[4,] ""   ""   "" 
```


- **`str_wrap() `**: 用于控制字符串的输出格式

  `width`: 输出每行中字符串的个数

  `indent`: 控制每段首行缩进字符数

  `exdent`: 悬挂缩进字符数

- **`str_view()/str_view_all()`**：用于正则表达式的匹配情况查看