# Perl

## 运行环境

```
Architecture:        x86_64
CPU op-mode(s):      32-bit, 64-bit
CPU(s):              4
Thread(s) per core:  1
Core(s) per socket:  4
Socket(s):           1
Model name:          Intel(R) Core(TM) i5-7300HQ CPU @ 2.50GHz
```

## perl版本

```perl
lzh@DESKTOP-S5MC5J5:~$ perl -v
This is perl 5, version 26, subversion 1 (v5.26.1) built for x86_64-linux-gnu-thread-multi
(with 67 registered patches, see perl -V for more detail)
```

## 标量变量

**标量的计算**：`+ - * / % **....` ps：取模仅对整数起作用

**字符串操作符**：

`.`: 字符串连接

`x`: 字符串重复出现多少次

**数字与字符串的比较操作符**：

| 数值 | 字符串 |
| ---- | ------ |
| ==   | eq     |
| !=   | ne     |
| <    | lt     |
| >    | gt     |
| <=   | le     |
| >=   | ge     |

**标量变量的命名与其他语言基本一致**

变量的定义：

`$variate =value;`

`undef`表示未定义变量，如`$_=undef`

换行符：Linux:`\n`  Mac: `\r` windows: `\r\n`  `cat -A`可以查看换行符，使用`dos2unix`系列工具可以实现转换。

**`chomp函数`: 去掉结尾的**换行符** 

**`chop`**: 去掉结尾的**一个字符**

```perl6
#用法举例
#!/usr/bin/perl -w
$dna="ATGC\n";
print $dna;

#去掉结尾换行符
chomp($dna);
print $dna;

#删除结尾最后一个字符
chop($dna);
print "\n";
print $dna;

###输出结果
lzh@DESKTOP-S5MC5J5:~$ perl ./per2.pl
ATGC
ATGC
ATGlzh@DESKTOP-S5MC5J5:~$
```

## 列表和数组

**`列表(list)`**：标量的有序集合

**`数组(array)`**：储存列表的变量

```perl6
#定义一个数组
@array=(1,4,2,5,6,4);
@array=(1,2,"abc",undef,$a,23);
```

列表的索引，类似于python，使用`[]`进行索引，索引下标为`0`开始，也可以负索引`-1:最后一个标量`，

**`$#`**：可以获得数组的最后一个元素的**`index`**值

`..`: 范围操作符，生成一串连续的数值标量。

**`qw`**: 省略逗号和引号

```perl6
#定义一个数组
@array=(1..10);
#数组索引
$array[0]=1
$array[1]=1
$array[-1]=10
#获得最后一个index
$#array=9
#数组元素个数计算
$#array + 1;

###举例

#!/usr/bin/perl -w
@abc=(1..100);
print "$abc[0]\n";
print $#abc+1;
print "\n";
print $#abc;
print "\n";

##输出结果
lzh@DESKTOP-S5MC5J5:~$ perl ./per2.pl
1
100
99
```

## 构建列表

使用`()`、`!!`、`<>`、`//`、均可以构建列表

`@list=(array)`

```perl6
@num1=(1,2,3,4,5,6);
@num=(1..6);
#使用qw操作符构建列表,省略逗号和引号
@num3=qw (1 2 3 4 6);
@str=qw (a b c d r f);

#也可使用()外的符号构建列表
#!/usr/bin/perl -w
@str2=qw !a b c d!;
@str3=qw (a b c c d);
@str4=qw /a b c d e/;
@str5=qw {a b c d};
@str=qw <a s d d 3>;
print "@str2\n";
print "@str3\n";
print "@str4\n";
print "@str5\n";
print "@str\n";
##输出结果
lzh@DESKTOP-S5MC5J5:~$ ./per2.pl
a b c d
a b c c d
a b c d e
a b c d
a s d d 3
```

**可以利用数组对多个标量同时赋值**

```perl6
#!/usr/bin/perl -w
($a,$b,$c,$d)= ("aaa",233,"111a","ddd");
print "$a\n";
print "$b\n";
print "$c\n";
print "$d\n";
#输出结果
lzh@DESKTOP-S5MC5J5:~$ ./per2.pl
aaa
233
111a
ddd
```

## 数组的操作

使用以下操作符可以对数组进行操作

**`split`**:将**标量**按照指定的分割符号进行分割，返回一个**数组**

```perl
#!/usr/bin/perl -w
$d="a:b:c:d";
@d=split /:/,$d;
print "@d\n";
##输出结果
lzh@DESKTOP-S5MC5J5:~$ ./per2.pl
a b c d
```

**`join`**: 将**数组**按照指定分隔符连接成一个**标量**

```perl
#!/usr/bin/perl -w
$d="a:b:c:d";
@d=split /:/,$d;
print "@d\n";

$e=join "_",@d;
print "$e\n"
##输出结果
lzh@DESKTOP-S5MC5J5:~$ ./per2.pl
a b c d
a_b_c_d
```

**`pop`**: 删除数组的最后的一个标量

**`push`**:将指定标量添加进数组

```perl6
#!/usr/bin/perl -w
@num=qw /a b c 1 2 3/;
print "@num\n";
#删除最后一个标量
pop @num;
print "@num\n";
#末尾增加一个标量
push @num,"233";
print "@num\n";
###输出结果
lzh@DESKTOP-S5MC5J5:~$ ./per2.pl
a b c 1 2 3
a b c 1 2
a b c 1 2 233
```

**`shift`**：数组开头**减少**一个标量

**`unshift`**：数组开头**增加**一个标量

```perl
#!/usr/bin/perl -w
@num=qw /a b c 1 2 3/;
print "@num\n";
unshift @num,"233";
print "@num\n";
shift @num;
print "@num\n";
####输出结果
lzh@DESKTOP-S5MC5J5:~$ ./per2.pl
a b c 1 2 3
233 a b c 1 2 3
a b c 1 2 3
```

**`sort`**：对数组中的标量进行排序，默认按照ASCII码排序，从小到大

**`reverse`**：从大到小排序

```perl6
#!/usr/bin/perl -w
@num=qw (10 9 8 a 1 12);
@num2=sort @num;
@num3=sort reverse @num;
print "@num2\n";
print "@num3\n"
###输出结果
lzh@DESKTOP-S5MC5J5:~$ ./per2.pl
1 10 12 8 9 a
1 10 12 8 9 a
```

**`foreach`**：遍历数组

```perl6
#!/usr/bin/perl -w
@num=qw (10 9 8 a 1 12);
foreach $numb (@num) {
        print "$numb\n";
}
#输出结果
lzh@DESKTOP-S5MC5J5:~$ ./per2.pl
10
9
8
a
1
12
```

如果未指定变量，则变量默认保存在`$_`

```perl
#!/usr/bin/perl -w
@num=qw (10 9 8 a 1 12);
foreach (@num) {
        print"$_\n";
}
###输出结果
lzh@DESKTOP-S5MC5J5:~$ ./per2.pl
10
9
8
a
1
12
```

