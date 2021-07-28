# awk学习

## awk介绍

AWK 是一种处理文本文件的语言，是一个强大的文本分析工具。

之所以叫 AWK 是因为其取了三位创始人 Alfred Aho，Peter Weinberger, 和 Brian Kernighan 的 Family Name 的首字符。

## awk基本命令格式与用法

```shell
#基本用法
awk '{[pattern] action}' {filenames} 
#举例
awk '{print $0}' ls.txt

#$0 表示输出全部内容；$1...n:表示输出第1....n列
root:x:0:0:root:/root:/bin/bash
daemon:x:1:1:daemon:/usr/sbin:/usr/sbin/nologin
bin:x:2:2:bin:/bin:/usr/sbin/nologin
sys:x:3:3:sys:/dev:/usr/sbin/nologin
sync:x:4:65534:sync:/bin:/bin/sync
games:x:5:60:games:/usr/games:/usr/sbin/nologin
man:x:6:12:man:/var/cache/man:/usr/sbin/nologin
lp:x:7:7:lp:/var/spool/lpd:/usr/sbin/nologin
mail:x:8:8:mail:/var/mail:/usr/sbin/nologin
news:x:9:9:news:/var/spool/news:/usr/sbin/nologin
```

### 指定分隔符

使用**`awk -F`**指定分隔符

```shell
awk -F ":" '{print $1}' ls.txt | head -n 3
root
daemon
bin
#可以同时指定多个分隔符
awk -F [":","/"] '{print $1}' ls.txt | head -n 3
```

### 内置变量

```shell
#输出行号
awk -F ":" '{print NR ")" $0}' ls.txt
1)root:x:0:0:root:/root:/bin/bash
2)daemon:x:1:1:daemon:/usr/sbin:/usr/sbin/nologin
3)bin:x:2:2:bin:/bin:/usr/sbin/nologin
awk -F ":"  'NR == 2 {print $NF}' ls.txt
/usr/sbin/nologin

#FILENAME：当前文件名
#FS：字段分隔符，默认是空格和制表符。
#RS：行分隔符，用于分割每一行，默认是换行符。
#OFS：输出字段的分隔符，用于打印时分隔字段，默认为空格。
#ORS：输出记录的分隔符，用于打印时分隔记录，默认为换行符。
#OFMT：数字输出的格式，默认为％.6g
#NF：表示当前行有多少字段
```

### 内置函数

```shell
#tolower()：字符转为小写。
#length()：返回字符串长度。
#substr()：返回子字符串。
#sin()：正弦。
#cos()：余弦。
#sqrt()：平方根。
#rand()：随机数。

awk -F ":" '{print toupper($1) ":" $2 ":" $3}' ls.txt
ROOT:x:0
DAEMON:x:1
BIN:x:2
SYS:x:3
SYNC:x:4
```

### 条件判断

```shell
#所有的逻辑判断符都可以用
#awk '条件 动作' 文件名
#只输出包含某字符的行
awk -F ":" '/root/ {print $0}' ls.txt
root:x:0:0:root:/root:/bin/bash
#输出偶数/奇数行的第一个字段
awk -F ":" 'NR%2 == 0 {print $1}' ls.txt
awk -F ":" 'NR%2 != 0 {print $1}' ls.txt
#根据指定字段条件输出
awk -F ":" '$1 == "root" {print $0}' ls.txt
```

### if函数判断

```shell
awk -F ":" '{if($1 == "root") print $0;else print $2}' ls.txt
awk -F ":" '{if(NR == 1) print $0;else print $2}' ls.txt
```

### awk用法举例

```shell
#筛选所有序列长度大于等于1000的rRNA
awk -F "\t" 'BEGIN{OFS="\t"; print "accessnum","source","type","Start","End","length"}{OFS="\t";if($5-$4 >=1000 && $3 == "rRNA") print $1,$2,$3,$4,$5,$5-$4} ' Bac_meg.gff | head
```

