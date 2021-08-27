# awk学习-day1

## awk介绍

AWK 是一种处理文本文件的语言，是一个强大的文本分析工具。

之所以叫 AWK 是因为其取了三位创始人 Alfred Aho，Peter Weinberger, 和 Brian Kernighan 的 Family Name 的首字符。

### 安装gawk

gawk是awk的最新版，功能更强大

```shell
wget --no-check-certificate https://mirrors.tuna.tsinghua.edu.cn/gnu/gawk/gawk-4.2.0.tar.gz
tar -zxvf gawk-4.2.0.tar.gz
cd gawk-4.2.0/
./configure --prefix=/usr/local/gawk
make -j 4 && make install
ln -fs /usr/local/gawk4.2/bin/gawk /usr/bin/awk
```

### 测试文件

```shell
wget https://aironi.oss-cn-beijing.aliyuncs.com/test.txt
```

## awk基本用法

```shell
awk 'awk_program' a.txt
```

`a.txt`: 是awk要读取的文件

`''`: 单引号内是代码块（尽量使用单引号，双引号内的$会被shell解析），**代码块位于大括号**中用内部代码`;`隔开，不同的代码块放于不同的大括号`{}`中

```shell
#awk读取所有行
awk '{print $0}' test.txt | head -n 5
#输出
ID name gender age email phone
1 Bob male 28 abc@qq.com 18023394012
2 Alice female 24 def@gmail.com 18084925203
3 Tony male 21 aaa@163.com 17048792503
4 Kevin male 21 bbb@189.com 17023929033
......
#多个代码块以及同一个代码块中多个语句
awk '{print $0}{print "hello";print "world"}' test.txt  | head -n 7
#printf 为数据字符时不带换行符
#---output----
ID name gender age email phone
hello world
1 Bob male 28 abc@qq.com 18023394012
hello world
2 Alice female 24 def@gmail.com 18084925203
hello world
3 Tony male 21 aaa@163.com 17048792503
.......
```

### awk读取文件过程

```shell
awk '{print $0}' test.txt
```

- `awk`：读取文件是按行进行读取的，`awk`每次读取文件将每一行（记录）的储存在`$0`，最后读取完成后，`print`输出，退出`awk`，
- 类似于`while / for`, 只不过`awk`隐藏了这个过程

## BEGIN与END

awk的所有代码大多是写在语句块中的

```shell
awk '{print $0}' test.txt
awk '{print $0}{print "hello";print "world"}' test.txt
```

`awk`的每个语句块前还可以跟上`pattern  `格式如下，一般语句块分为三类分别为`BEGIN` `END`  `main`，其中前两者的格式都是固定的，`BEGIN{...}`和`END{...}`，**而main语句块是一种统称， 它的pattern部分没有固定格式， 也可以省略， main代码块是在读取文件的每一行的时候都执行的代码块**

```shell
pattern1{statement1}pattern2{statement3;statement4;...}
```

### BEGIN代码块

在读取文件前所运行的代码块, BEGIN中无法使用`$0`这样的特殊变量

```shell
awk 'BEGIN{print "BEGIN运行!";print "test!"}{print $0}' test.txt | head -n 4
#----output----
BEGIN运行!
test!
ID name gender age email phone
1 Bob male 28 abc@qq.com 18023394012
```

### main代码块

读取文件时，对每一行所运行的代码，每读取一行， 就执行一次main代码块 ，同时main代码块可以有多个

### END代码块

在读取文件结束的时候运行，可以使用`$0`这样的特殊变量，但是有些变量只会储存最后一轮循环的值

```shell
awk '{print $0}END{print "END test!!"}' test.txt | tail -n 4
#---output----
8 Peter male 20 bax@qq.com 17729348758
9 Steven female 23 bc@sohu.com 15947893212
10 Bruce female 27 bcbd@139.com 13942943905
END test!!
```

## awk语法结构

`awk`的语法结构，即`main`代码块内部的代码结构，一般来说`awk`的语法结构如下，`pattern{action}`成为`awk rules`

```shell
awk pattern{action} file
```

举例子`/Alice/`即为`pattern`部分, `{print $0}`即为`action`部分

```shell
awk '/Alice/{print $0}' test.txt
#----output----
2 Alice female 24 def@gmail.com 18084925203
```

- 多个 `pattern{action}` 可以直接连接  

```shell
awk '/Alice/{print $0} /Peter/{print $0}' test.txt
#----output----
2 Alice female 24 def@gmail.com 18084925203
8 Peter male 20 bax@qq.com 17729348758
```

- `action`中多个语句如果写在同一行， 则需使用分号分隔  

- **`pattern`部分用于筛选行， `action`表示在筛选通过后执行的操作**  

- `pattern`和`action`都可以省略

  **省略`pattern`，等价于对每一行数据都执行action**

```shell
awk '{print $0}' test.txt | head -n 3
#----output----
ID name gender age email phone
1 Bob male 28 abc@qq.com 18023394012
2 Alice female 24 def@gmail.com 18084925203
```

​	 **省略代码块 `{action} `， 等价于 `{print}` 即输出所有行  **

```shell
awk '/Alice/' test.txt
#----output----
2 Alice female 24 def@gmail.com 18084925203
```

​	 **省略代码块中的` action `， 表示对筛选的行什么都不做**  

```shell
awk '/Alice/{}' test.txt
#----output----
None
```

​	 **全部省略**，什么也不干

```shell
awk '' test.txt
```

### pattern和action

- `pattern`可以使用以下模式

```shell
# 特殊pattern
BEGIN
END

#bool代码块
/regular expression/	# 正则匹配成功与否 /a.*ef/{action}
relational expression	# 即等值比较、 大小比较 3>2{action}
pattern && pattern		# 逻辑与 3>2 && 3>1 {action}
pattern || pattern		# 逻辑或 3>2 || 3<1 {action}
! pattern				# 逻辑取反 !/a.*ef/{action}
(pattern)				# 改变优先级
```

`action`部分， 可以是任何语句， 例如`print`

## awk读取文件详细过程

### 换行符RS

`awk`读取输入文件时， 每次读取一条记录`(record)`(默认情况下按行读取， 所以此时记录就是行)。 每
读取一条记录， 将其保存到`$0` 中， 然后执行一次`main`代码段。  

```shell
awk '{print $0}' test.txt
```

如果是空文件， 则因为无法读取到任何一条记录， 将导致直接关闭文件， 而不会进入`main`代码段  

```shell
touch test2.txt
awk '{print $0}' test2.txt
#---output---
None
```

可设置表示输入**记录分隔符的预定义变量RS(Record Separator)**来改变每次读取的记录的模式 ,**RS通常设置在BEGIN代码块中， 因为要先于读取文件就确定好RS分隔符  **

`RS`指定输入记录分隔符时， 所读取的记录中是不包含分隔符字符的。 例如 `RS="a"` ， 则 `$0` 中一定不可能出现字符`a`  ，`RS`一般为单个字符，当为多个字符时，会将其当作正则对待

```shell
awk 'BEGIN{RS="\n"}{print $0}' test.txt | head -n 4
##----output----##
ID name gender age email phone
1 Bob male 28 abc@qq.com 18023394012
2 Alice female 24 def@gmail.com 18084925203
3 Tony male 21 aaa@163.com 17048792503

#如果指定其他字符为分隔符，例如指定a为分隔符，awk在读取文件时，一旦遇到a就会把a当作换行符
#原来的表头：ID name gender age email phone
awk 'BEGIN{RS="a"}{print $0}' test.txt | head -n 4
##----output----##
ID n
me gender 
ge em
il phone

#指定多个字符当作分隔符(a或b)
awk 'BEGIN{RS="[ab]"}{print $0}' test.txt | head -n 6
##----output----##
ID n
me gender 
ge em
il phone
1 Bo
 m
```

**特殊的`RS`值用来解决特殊读取需求**  

`RS="" `： 按段落读取

`RS="\0"` ： 一次性读取所有数据， 但有些特殊文件中包含了空字符 \0

`RS="^$" `： 真正的一次性读取所有数据， 因为非空文件不可能匹配成功

`RS="\n+"` ： 按行读取， 但忽略所有空行  

### 行号 (每条记录的编号)NR和FNR

在读取每条记录之后， 将其赋值给$0， 同时还会设置`NR`、 `FNR`、 `RT`

`NR`： 所有文件的行号计数器

`FNR`： 是各个文件的行号计数器  

如，输出时候显示行号

```shell
awk '{print NR "\t"  $0}' test.txt | head -n 3
#最左侧为行号
1	ID name gender age email phone
2	1 Bob male 28 abc@qq.com 18023394012
3	2 Alice female 24 def@gmail.com 18084925203
```

## awk字段分割 

### 记录与字段的关系

```shell
1 Bob male 28 abc@qq.com 18023394012
2 Alice female 24 def@gmail.com 18084925203
3 Tony male 21 aaa@163.com 17048792503
4 Kevin male 21 bbb@189.com 17023929033
5 Alex male 18 ccc@xyz.com 18185904230
6 Andy female 22 ddd@139.com 18923902352
7 Jerry female 25 exdsa@189.com 18785234906
8 Peter male 20 bax@qq.com 17729348758
9 Steven female 23 bc@sohu.com 15947893212
10 Bruce female 27 bcbd@139.com 13942943905
#awk在读这个表时，会将每一行作为一个记录，每一行的Bob male 28....则为每一个条记录的字段，类似于数据库
```

### 引用字段的方式  

`awk`读取每一条记录之后， 会将其赋值给` $0 `， 同时还会对这条记录按照预定义变量`FS`划分字段， 将划分好的各个字段分别赋值给 `$1 $2` `$3` `$4`...`$N` ， 同时将划分的字段数量赋值给预定义变量`NF`  

`$N` 引用字段：

`N=0` ： 即 `$0` ， 引用记录本身

`0<N<=NF` ： 引用对应字段

`N>NF` ： 表示引用不存在的字段， 返回空字符串

`N<0` ： 报错 

```shell
#显示每个记录的第五个字段
awk '{n = 5;print $n}' test.txt | head -n 3
##----output----
email
abc@qq.com
def@gmail.com

#显示每条记录的第四个字段
awk '{print $(2+2)}' test.txt | head -n 3
##----output----
age
28
24
awk '{print $(NF-3)}' test.txt | head -n 3
##----output----
gender
male
female

#显示最后一个字段
awk '{print $NF}' test.txt | head -n 4
##----output----
phone
18023394012
18084925203
17048792503
```

###  分割字段的方式（FS）

读取`record`之后， 将使用预定义变量`FS`、 `FIELDWIDTHS`或`FPAT`中的一种来分割字段。 分割完成之后， 再进入`main`代码段(所以， 在`main`中设置`FS`对本次已经读取的`record`是没有影响的， 但会影响下次读取)  , 即`FS`为用于分割记录的分隔符

使用`FS` 或者` -F` ： 指定字段的分隔符号，常见的如`:` `\t`  ` ` 

- **`FS`为单个字符时， 该字符即为字段分隔符**  
- `FS`为多个字符时， 则采用正则表达式模式作为字段分隔符  
- **默认情况**：  `FS`为单个空格时， 将以连续的空白（ 空格、 制表符、 换行符） 作为字段分隔符
- **特殊的，` FS`为空字符串`""`时， 将对每个字符都进行分隔， 即每个字符都作为一个字段**    
- 设置预定义变量`GNORECASE`为非零值， 正则匹配时表示忽略大小写(只影响正则， 所以`FS`为单字时无影响)  
- **如果`record`中无法找到FS指定的分隔符(例如将`FS`设置为`\n`)， 则整个记录作为一个字段，即 `$1`和 `$0` 相等**  

```shell
#以passwd文件为例
awk 'BEGIN{FS=":"}{print $1}' passwd | head -n 4
#第二种方式
awk -F ":" '{print $1}' passwd | head -n 4
##----output----
root
bin
daemon
adm
#以正则的形式匹配
awk 'BEGIN{FS="@|+"}{print $1 $2 $3 $4 $5 $6}' test.txt | head -n 3
##----output----
ID name gender age email phone
1 Bob male 28 abcqq.com 18023394012
2 Alice female 24 defgmail.com 18084925203
```

### 按照FIELDWIDTHS  对字段进行划分

指定预定义变量`FIELDWIDTHS`按字符宽度分割字段， 这是`gawk`提供的高级功能。 在处理某字段缺失时非常好用(仅限`gawk`)  

简单举例

- `FIELDWIDTHS = "3 5 6 9"`: 表示第一个字段3个字符，第二个字段5字符，第三个字段6字符，第四个字段9字符

- `FIELDWIDTHS = "8 1:5 6 2:33`: 第一个字段8个字符，**然后跳过1个字段在读5个字符**作为第2字段，读取6个字符作为第3字段，跳过

  两个字符在读取33个字符作为最后一个字段（不足33个则默认读到最后一个字符）

- `FIELDWIDTHS = "3 5 *"`: 读取3个字符作为第一个字段，读取5个字符作为第二个字段，***代表剩下的所有字符作为最后一个字段, 只可以用于结尾**

## awk数据筛选示例  

### 根据行号进行筛选

```shell
#筛选得到第4行数据
awk 'NR==4' test.txt
#筛选第8行后的所有数据
awk 'NR>=8' test.txt 
#输出行号
awk '{print NR $0}' test.txt
```

### 根据正则表达式筛选行

```shell
#输出包含qq.com的行
awk '/qq.com/{print $0}' test.txt
awk '$0 ~ /qq.com/' test.txt
#输出不包含@的行
awk '/^[^@]+$/{print $0}' test.txt
awk '!/@/{print $0}' test.txt
```

### 根据字段筛选行

```shell
#筛选所有男性的数据
awk '$3=="male"{print $0}' test.txt
#筛选大于21岁的数据
awk '$4>21 {print $0}' test.txt
#输出第五字段包含qq.com的记录
awk '$5 ~ /qq.com/ {print $0}' test.txt
```

### 多条件筛选

```shell
#筛选第3到7行(第2到第7条记录)
awk 'NR>=3 && NR<=7' test.txt
#筛选男性并且手机号以170开头的
awk '$3=="male && $6 ~ /^170/"{print $0}' test.txt
#筛选男性或手机号为170开头的记录
 awk '$3=="male" || $6 ~ /^170/' test.txt
```

### 按照范围筛选

```shell
#筛选第2到7个记录
awk 'NR==2,NR==7' test.txt
#筛选第二条记录直到手机号以170开头的记录为止
awk 'NR==2, $6 ~ /^170/' test.txt
##output
1 Bob male 28 abc@qq.com 18023394012
2 Alice female 24 def@gmail.com 18084925203
3 Tony male 21 aaa@163.com 17048792503
#小例子 筛选ip地址
ifconfig | awk '/inet/ && $2 ~ /^172/  {print $2}'
ifconfig | awk '/inet/ && $2 ~ /^[0-9]/ && $2 != "127.0.0.1" {print $2}'
```

