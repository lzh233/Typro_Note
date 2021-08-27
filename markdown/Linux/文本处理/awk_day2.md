# awk学习-day2

## awk的输出-print

awk可以通过print、 printf将数据输出到标准输出或重定向到文件。  逗号分隔要打印的字段列表， 各字段都会自动转换成字符串格式， 然后通过预定义变量OFS(outputfield separator)的值(其默认值为空格)连接各字段进行输出。  

```shell
#举例子
awk 'BEGIN{print "hello","world"}'
#output,默认以空格为连接符号
hello world
#指定连接符OFS
awk 'BEGIN{OFS = "++";print "hello","world"}
#output 
hello++world
```

`print`：使用`print`输出时，默认的换行符(`ORS`)为`\n`，可以自行指定

```shell
awk 'BEGIN{ORS = "++\n";print "hello","word"}'
#output
hello world++
```

### print

`print`在输出数据时， 总是会先转换成字符串再输出 

对于数值而言， 可以自定义转换成字符串的格式， 例如使用`sprintf()`进行格式化。print在自动转换数值（ 专指小数） 为字符串的时候， 采用预定义变量`OFMT（ Output format）` 定义的格式按照`sprintf()`相同的方式进行格式化。 `OFMT`默认值为 `%.6g` ， 表示有效位(整数部分加小数部分)最多为`6`   

```shell
awk 'BEGIN{print 1.12345678}'
#output
1.12346
```

**可以修改`OFMT`， 来自定义数值转换为字符串时的格式**  

```shell
#转换为整数类型
awk 'BEGIN{OFMT = "d%";print 1.234}'
#output
1
#保留三位小数
awk 'BEGIN{OFMT = "%.3f";print 1.12345678}'
#output
1.123
```

### 常见的格式化字符的方法

`c%`: 将ASCII码转换为字符

`d%`：转换为整数

`%f`：转换为小数，%.3f 保留三位小数

`%g `：保留多少有效数字`%.6g` ， 表示有效位(整数部分加小数部分)最多为`6`

`%s`：转换为字符串

## awk变量**

### awk变量类型

`awk`的变量是动态变量， 在使用时声明。所以`awk`变量有3种状态：

-  未声明状态： 称为`untyped`类型

- 引用过但未赋值状态：`unassigned`类型

- 已赋值状态

- 引用未赋值的变量， 其默认初始值为空字符串或数值0。

  在awk中未声明的变量称为`untyped`， 声明了但未赋值(**只要引用了就声明了**)的变量其类型为`unassigned`。  

### awk变量赋值

awk中的变量赋值语句也可以看作是一个有返回值的表达式。例如， a=3 赋值完成后返回3， 同时变量a也被设置为3。
基于这个特点， 有两点用法：

- 可以 x=y=z=5 ， 等价于 z=5 y=5 x=5
- 可以将赋值语句放在任意允许使用表达式的地方  

```shell
 awk 'BEGIN{a=5;print a}'
 #output 
 5
```

### awk中声明变量的位置  

- 在`BEGIN`或`main`或`END`代码段中直接引用或赋值  

```shell
 awk 'BEGIN{a=5;print a}'
 awk '{a=5;print a}'
```

- **使用`-v`选项定义变量与接收外部参数**

  使用 `-v var=val` 选项， 可定义多个， 必须放在`awk`代码的前面

  它的变量声明早于`BEGIN`， 所以`-v`选项声明的变量在`BEGIN{}`、 `END{}`和`main`代码
  段中都能直接使用  

  普通变量： `awk -v age=123 'BEGIN{print age}'`

  使用shell变量赋值： `awk -v age=$age 'BEGIN{print age}`  

```shell
awk -v age=123 'BEGIN{print age}
#当变量来自外部
a=1
b=3
c=23
awk -v aa=$a -v bb=$b -v cc=$c 'BEGIN{print aa,bb,cc}'
#output
1 3 23
```

### awk中的数据类型

gawk有两种基本的数据类型： 数值和字符串。 在gawk 4.2.0版本中， 还支持第三种基本的数据类型：正则表达式类型。  

数据是什么类型在使用它的上下文中决定： 在字符串操作环境下将转换为字符串， 在数值操作环境下将转换为数值。  

通过`typeof()`函数可以查看数据的类型

```shell
 awk 'BEGIN{a="123a";print a,typeof(a)}'
 #output
 123a string
 awk 'BEGIN{a="123a"+1;print a,typeof(a)}'
 #output
 124 number
```

### awk字面量  

awk中有3种字面量： 字符串字面量、 数值字面量和正则表达式字面量  

- 数值字面量  

算术运算  

```shell
++ -- 自增、 自减， 支持i++和++i或--i或i--
^ 幂运算(**也用于幂运算)
+ - 一元运算符(正负数符号)
* / % 乘除取模运算
+ - 加减法运算
```

- 字符串字面量  

字符串连接（ 串联） ： awk没有为字符串的串联操作提供运算符， 可以直接连接或使用空格连接。  

- **正则表达式字面量**

`/[0-9]+/`

匹配方式： `"str" ~ /pattern/` 或 `"str" !~ /pattern/`

匹配结果返回值为0(匹配失败)或1(匹配成功)

任何单独出现的 `/pattern/` 都等价于 `$0 ~ /pattern/`

`if(/pattern/)` 等价于 `if($0 ~ /pattern/)  `

## awk内置变量与函数***

### awk字符串类内置函数

**基本函数**

`sprintf(format, expression1, ...) `： 返回格式化后的字符串  

`length() `： 返回字符串字符数量、 数组元素数量、 或数值转换为字符串后的字符数量  

`tolower(str) `： 转换为小写

`toupper(str)` ： 转换为大写  

`index(str,substr)` ： 从str中搜索substr(子串)， 返回搜索到的索引位置(索引从1开始)， 搜索不到则返回0  

**awk substr()  **

`substr(string,start[,length]) `： 从string中截取子串 ，如果start值小于1， 则将其看作为1对待， 如果start大于字符串的长度， 则返回空字符串

```shell
awk 'BEGIN{s="ATTGCATTGCGGGAGG";print substr(s,2,5)}'
#output
TTGCA
#可以应用于截取fasta文件中的序列或者基因名等等
#以剪切序列为例子
```

### awk split()和patsplit()  

`split(string, array [, fieldsep [, seps ] ])` ： 将字符串分割后保存到数组`array`中，数组索引从1开始存储。 并返回分割得到的元素个数  ,其中`fieldsep`指定分隔符， 可以是正则表达式方式的。 如果不指定该参数， 则默认使用`FS`作为分隔符， 而`FS`的默认值又是空格  

```shell
```



