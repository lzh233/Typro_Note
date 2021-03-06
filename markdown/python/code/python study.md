# Python

## argparse包的使用

### 用途

通过命令行控制输入变量

### 使用方法

```python
import python
"""
Creating index for reference from fasta
Aligling fastq to index

lunzi 1.0_only for pairs-end
"""
import os
import argparse

from pandas.io.parsers import count_empty_vals 
#定义一个Bowtie2类
class Bowtie2:
    def __init__(self,reference_path,index_path,bam_path):
        self.reference_path = reference_path
        self.index_path = index_path
        self.bam_path = bam_path
    
    def build_index(self):
        cmd1 = f"bowtie2-build {self.reference_path} {self.index_path}"
        return cmd1

    def align_pairs(self,fasta_path,fasta_path_R):
        cmd2 = f"bowtie2  -x {self.index_path} -1 {fasta_path} -2 {fasta_path_R} -S {self.bam_path}"
        return cmd2


def main():
    #创建argparse对象，通过ArgumentParser函数
    parser = argparse.ArgumentParser(description='allign by Bowtie2')
    #通过add_argument方法可以增加选项，帮助信息为help中的内容，ruquire为是否必须选项，每个argparse对象都有一个固定的选项-h，可以输出明命令的帮助信息
    parser.add_argument('--fq1', help='fastq', required=True)
    parser.add_argument('--fq2', help='fastq_r.', required=True)
    parser.add_argument('--reference', help='reference path', required=True)
    parser.add_argument('--index_out', help='output path of indeies', required=True)
    parser.add_argument('--bam_out', help='output path of bam_file', required=True)
    #parser.add_argument('--to_sam', help='to sam?', required=True)
    #将选项信息储存在args中
    args = parser.parse_args()
   #选项的使用
#使用args.<选项的全名(--指定的名称)> = 选项的值
#如args.index_out
#  args.reference
b = Bowtie2(index_path = str(args.index_out),
                reference_path=str(args.reference),
                bam_path = str(args.bam_out))
    
    cmd1 = b.build_index()
    cmd2 =b.align_pairs(fasta_path=args.fq1,
                        fasta_path_R=args.fq2)

    file_handle=open('run.sh',mode='w')
    file_handle.writelines([f'{cmd1}\n',f'{cmd2}'])
    file_handle.close()

if __name__ == '__main__':
    main()

```

```shell
 ##使用方法, 相应的参数就会传到对应的位置
python allign.py --fq1 ./seq/fq1 --fq2 ./seq/fq2 --index_out ./result/out_index --bam_out ./result/bam_out --reference /vdj_ref
```

```python
#输出信息
bowtie2-build /vdj_ref ./result/out_index
bowtie2  -x ./result/out_index -1 ./seq/fq1 -2 ./seq/fq2 -S ./result/bam_out
```

## python文件读取与写入

### 打开和关闭文件

```python
#使用open函数打开文件，创建一个file对象, 可以对文件进行操作
#file object = open(file_name [, access_mode][, buffering])
#常用的mode如下，
```

`t`: 文本模式，默认

`x`: 写模式，新建一个文件，如果该文件已存在则会报错

`b`: 二进制模式

`+`: 打开一个文件进行更新(可读可写)

`r`: 以只读方式打开文件。文件的指针将会放在文件的开头

`w`:  打开一个文件只用于写入。如果该文件已存在则打开文件，并从开头开始编辑，**原有内容会被删除**。如果该文件不存在，创建新文件

`a`: 打开一个文件用于追加。如果该文件已存在，**文件指针将会放在文件的结尾, 新的内容将会被写入到已有内容之后。**如果该文件不存在，创建新文件进行写入。

`a+`:  打开一个文件**用于读写**。如果该文件已存在，文件指针将会放在文件的结尾。文件打开时会是追加模式。如果该文件不存在，创建新文件用于读写

```python
#打开一个文件,追加模式，会在文件末尾进行添加
fd = open(fo = open("01_save_kmerww.txt",mode="a"))
print(f"文件名称: {fd.name}")
print(f"文件访问模式: {fd.mode}")
print(f"是否关闭: {fd.close}")
##output
#文件名称: 01_save_kmerww.txt
#文件访问模式: a
#是否关闭: <built-in method close of _io.TextIOWrapper object at 0x01D60E38>
```

### close方法

```python
# 关闭打开的文件
fo.close()
```

### write方法

```python
fd.write("test")
fd.write("\ntest2")
#--output
#testtest
#test2
```

## python进阶

### 生成器和迭代器

迭代器是一个可以让程序员**遍历一个容器（尤其是列表）的对象**，然而一个**迭代器在遍历啊一个容器的数据元素时候，并不会执行一个迭代**

- 可迭代对象（Iterable）：Python中任意的对象，只要它定义了可以返回一个迭代器的 `__iter__` 方法，或者定义了可以支持下标索引的 `__getitem__` 方法，那么它就是一个可迭代对象

  ```python
  #以列表为例
  a = [1,2,3]
  print(dir(a))
   '__getitem__', '__iter__'
      ......
  ```

- 迭代器（Iterator）：任意对象，只要定义了 `next`（Python2） 或者 `__next__` 方法，它就是一个迭代器

  ```python
  #以字符串为例
  a = "aa"
  #将a转换为迭代器
  a = iter(a)
  print(dir(a))
  '__class__', '__delattr__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__'
  ......
  ```

- 迭代（Iteration）:从某个地方（比如一个列表）取出一个元素的过程。当我们使用一个循环来遍历某个东西时，这个过程本身就叫迭代

#### 生成器(Generators)

**生成器也是一种迭代器，但是你只能对其迭代一次。这是因为它们并没有把所有的值存在内存中，而是在运行时生成值。**你通过遍历来使用它们，要么用一个 “for” 循环，要么将它们传递给任意可以进行迭代的函数和结构。大多数时候生成器是以函数来实现的。然而，它们并不返回一个值，而是 `yield` (暂且译作“生出”)一个值

```python
#一个简单的生成器
def genera_function():
    for i in range(3):
        yield i
test = genera_function()
for a in test:
    print(a)
print(test)
##output
0
1
2
<generator object genera_function at 0x0000020CAAFAF580>
```

生成器最佳应用场景是：你不想同一时间将所有计算出来的大量结果集分配到内存当中，特别是结果集里还包含循环

**`next()`**函数：它允许我们获取一个序列的下一个元素。

```python
def generator_function():
    for i in range(3):
        yield i

gen = generator_function()
print(next(gen))
# Output: 0
print(next(gen))
# Output: 1
print(next(gen))
# Output: 2
print(next(gen))
# Output: Traceback (most recent call last):
#            File "<stdin>", line 1, in <module>
#         StopIteration
```

在 `yield` 掉所有的值后，`next()` 触发了一个 `StopIteration` 的异常。这个异常告诉我们，所有的值都已经被 `yield` 完了

关于`for`循环：

- 先判断对象是否为可迭代对象，不是的话直接报错，抛出TypeError异常，是的话，调用 iter方法，返回一个迭代器

- 不断地调用迭代器的next方法，每次按序返回迭代器中的一个值

- 迭代到最后，没有更多元素了，就抛出异常 StopIteration，这个异常 python 自己会处理，不会暴露给开发者

#### 迭代器

```python
my_string = "Yasoob"
next(my_string)
# Output: Traceback (most recent call last):
#      File "<stdin>", line 1, in <module>
#    TypeError: str object is not an iterator
```

这个异常说那个 `str` 对象不是一个迭代器。**它是一个可迭代对象，而不是一个迭代器**。这意味着它**支持迭代**，但我们不能直接对其进行

迭代操作。那我们怎样才能对它实施迭代呢？是时候学习下**另一个内置函数，`iter`。它将根据一个可迭代对象返回一个迭代器对象。**这

里是我们如何使用它：

```python
my_string = "Yasoob"
my_iter = iter(my_string)
next(my_iter)
# Output: 'Y'
```

### map  /  Filter / Reduce 

#### 匿名函数`lambda`

**形式**

**`argument_list`是参数列表，它的结构与Python中函数`(function)`的参数列表是一样的**

**`expression`是一个关于参数的表达式。表达式中出现的参数需要在`argument_list`中有定义，并且表达式只能是单行的。**

```python
 lambda argument_list: expression
```

**特性**

lambda函数有如下特性：

- `lambda`函数是匿名的：所谓匿名函数，通俗地说就是没有名字的函数。`lambda`函数没有名字。

- `lambda`函数有输入和输出：输入是传入到参数列表`argument_list`的值，输出是根据表达式`expression`计算得到的值。

- `lambda`函数一般**功能简单**：单行expression决定了`lambda`函数不可能完成复杂的逻辑，只能完成非常简单的功能。由于其实现的功能一目了然，甚至不需要专门的名字来说明

**用法**

- **将lambda函数赋值给一个变量，通过这个变量间接调用该lambda函数**

```python
test = lambda x,y:x+y
test(2,3)
#output
5
```

- 将lambda函数赋值给其他函数，从而将其他函数用该lambda函数替换

-  将lambda函数作为其他函数的返回值，返回给调用者

-  将lambda函数作为参数传递给其他函数

#### map

`Map` 会将一个函数映射到一个输入列表的所有元素上。这是它的规范：

大多数时候，我们使用匿名函数`lambdas`来配合 `map`

```python
map(function_to_apply, list_of_inputs)
```

例如实现列表中每个元素+1

```python
#一般方法
items = [1,2,3,4,5]
for i in range(len(items)):
    items[i] += 1
#map方法
items = list(map(lambda x: x+1,items))
print(items)
#output
[2, 3, 4, 5, 6]
```

#### Filter

`filter` 过滤列表中的元素，并且返回一个由所有符合要求的元素所构成的列表，**符合要求**即函数映射到该元素时返回值为True。

```python
num_list = [-5,-4,-3,-2,-1,0,1,2,3,4,5]
#将num_list中的每个参数传个x，然后判断x是否小于0
less_zero_lst = filter(lambda x : x < 0,num_list)
print(list(less_zero_lst))
#output
[-5,-4,-3,-2,-1]
```

#### Reduce

当需要对一个列表进行一些计算并返回结果时，`Reduce` 是个非常有用的函数。举个例子，当你需要计算一个整数列表的乘积时

```python
#常规方法:循环
num_list = [-5,-4,-3,-2,-1]
a = 1
for i in num_list:
    a = a*i
print(a)

#reduce方法
from functools import reduce
#过程:前两个数相乘然后构成一个新的n-1项的列表，重复上述过程，直到返回一个值
#[2,3,4]
#第一次[6,4]
#第二次
# 24
a = reduce( (lambda x,y:x*y),num_list)
print(a)
#output
-120
```

### Set 集合

`set`（集合）是一个非常有用的数据结构。它与列表（`list`）的行为类似，区别在于 `set` 不能包含重复的值

如**重复值**提取

```python
#普通办法
some_list = ['a', 'b', 'c', 'b', 'd', 'm', 'n', 'n']
dup =[]
for i in some_list:
    if some_list.count(i) > 1 and i not in dup:
        dup.append(i)
print(dup)
#output
['b','n']
#用set的方法(也可以用于提取出现了指定次数以上的值)
dup2 = set([x for x in some_list if some_list.count(x) > 1])
print(dup2)
#output
{'b','n'}
```

#### 交集

```python
set1 = {"a","b","c"}
#set1 = {["a","b","c"]}
set2 = {"a","b"}
#set2 = {["a","b"]}
print(set1.intersection(set2))
#output
{'b', 'a'}
```

#### 差集

```python
set1 = {"a","b","c"}
set2 = {"a","b"}
print(set1.difference(set2))
#output
{'c'}
```

### 装饰器

装饰器（Decorators）是 Python 的一个重要部分。简单地说：他们是修改其他函数的功能的函数。他们有助于让我们的代码更简短，

#### 内置装饰器

`property` 装饰器用于类中的函数，使得我们可以像访问属性一样来获取一个函数的返回值。

```python
class XiaoMing:
    first_name = '明'
    last_name = '小'

    @property
    def full_name(self):
        return self.last_name + self.first_name

xiaoming = XiaoMing()
print(xiaoming.full_name)
```

`staticmethod` 装饰器同样是用于类中的方法，这表示这个方法将会是一个静态方法，意味着该方法可以直接被调用无需实例化，但同样意味着它没有 `self` 参数，也无法访问实例化后的对象。

`classmethod` 依旧是用于类中的方法，这表示这个方法将会是一个类方法，意味着该方法可以直接被调用无需实例化，但同样意味着它没有 `self` 参数，也无法访问实例化后的对象。相对于 `staticmethod` 的区别在于它会接收一个指向类本身的 `cls` 参数。

```python
class XiaoMing:
    name = '小明'

    @classmethod
    def say_hello(cls):
        print('同学你好， 我是' + cls.name)
        print(cls)

XiaoMing.say_hello()
```

https://www.zhihu.com/question/26930016/answer/1904166977

### Gllobe和retuen

```python
#如果要在函数外部访问其内部变量，可以使用globe关键字将其设置为全局变量(实际情况能不这样就不这样)
def add(value1,value2):
    global result
    result = value1 + value2

add(3,5)
print(result)
# Output: 8
```

### 多个return值

那如果你想从一个函数里返回两个变量而不是一个呢？ 。最著名的方法，是使用 `global` 关键字。**（绝对不能这样）**

```python
#错误示范
def profile():
    global name
    global age
    name = "Danny"
    age = 30
```

返回一个包含多个值的 `tuple`(元组)，`list`(列表)或者 `dict`(字典)，来解决这个问题，或者直接return两个值

```python
def test():
    name = "aa"
    age ="12"
    return name,age
#这样默认返回的是一个元组对象
a = test()
print(a)
print(type(a))
#output
('aa', '12')
<class 'tuple'>
```

### 对象变动**

曾经踩过的坑......

可变（mutable）意味着"可以被改动"，而不可变（immutable）的意思是“常量（constant）

```python
#对象的可变
a = [1,2,3]
b = a
b[0] = 999
print(a)
print(b)
#output
[999, 2, 3]
[999, 2, 3]
```

对象可变性（**mutability**）在作怪。每当你将一个变量赋值为另一个可变类型的变量时，对这个数据的任意改动会同时反映到这两个变量上去。新变量只不过是老变量的一个别名而已。这个情况只是针对可变数据类型

#### 解决办法1

```python
def add_to(num, target=[]):
    target.append(num)
    return target

print(add_to(1))
print(add_to(2))
print(add_to(3))
#output
[1]
[1, 2]
[1, 2, 3]
##解决
def add_to(element, target=None):
    if target is None:
        target = []
    target.append(element)
    return target
#output
[1]
[2]
[3]
```

#### **解决办法2**

**通过`copy`包解决**

```python
"""
函数每次运行时，会通过copy.deepcopy创建一个target的副本，每次在这个副本中添加元素
"""
import copy

def add_to(element, target=[]):
    target = copy.deepcopy(target)
    target.append(element)
    return target

print(add_to(1))
print(add_to(2))
print(add_to(3))
```

### __slots__魔法

在 Python 中，每个类都有实例属性。默认情况下 Python 用一个字典来保存一个对象的实例属性它会消耗掉很多内存。 不过还是有一个方法来规避这个问题。这个方法需要使用 `__slots__` 来告诉 Python 不要使用字典，而且只给一个固定集合的属性分配空间

```python
class test:
    __slots__ = ["name","age"]
    def __init__(self,name,age):
        self.name = name
        self.age = age

a  = test(name="tom",age = 12)
print(a.name)
#output
tom
```

### enumerate

枚举（`enumerate`）是 Python 内置函数。可以将列表转成`[(index:[list[0]]),....]`

```python
a = ["a","b","c"]
for k,v in enumerate(a):
    print(k,v)
#output
0 a
1 b
2 c
#也可以指定索引开始的数字
1 a
2 b
3 c
```

### 对象自省

自省`（introspection）`，在计算机编程领域里，是指在运行时来判断一个对象的类型的能力。

#### dir

返回对象所有的方法和属性

```python
a = ["a","b","c"]
print(dir(a))
#output
['__add__', '__class__', '__contains__', '__delattr__', '__delitem__', '__dir__', ...... 'append', 'clear', 'copy', 'count', 'extend', 'index', 'insert', 'pop', 'remove', 'reverse', 'sort']
```

#### type

返回对象的类型

```python
print(type(a))
#output
<class 'list'>
```

### 推导式 Comprehension

#### 列表解析

结构是在一个中括号里包含一个表达式，然后是一个 `for` 语句，然后是0个或多个 `for` 或者 `if` 语句。那个表达式可以是任意的，意思是你可以在列表中放入任意类型的对象。返回结果将是一个新的列表

```python
lst = [a**2 for a in range(5) if a % 3 == 0]
lst2 = [a**2 for a in range(5)]
print(f"{lst}\n{lst2}")
#print
[0, 9]
[0, 1, 4, 9, 16]
```

#### 集合解析

与列表一致，`[]`改为`{}`即可

```python
lst = {a**2 for a in range(5) if a % 3 == 0}
lst2 = {a**2 for a in [1,1,1,1,1,1]}
print(f"{lst}\n{lst2}")
#output
{0, 9}
{1}
```

#### 字典解析

用的少，与列表 集合基本一致，只不过字典要注意输出的值为`k:v`形式

```python
#键值对互换
a = {"Tom":"name",12:"age","coder":"work"}
b = {v:k for k,v in a.items()}
print(b)
#output
{'name': 'Tom', 'age': 12, 'work': 'coder'}

#对每个值加10
c = {"A":1,"B":2,"C":3}
c = {k:v+10 for k,v in c.items()}
print(c)
#output
{'A': 11, 'B': 12, 'C': 13}
```

### __getitem__

在类中创建`__getitem__`方法，将允许使用索引的方式访问属性

```python
class test:
    def __init__(self):
        self.Tom = {
            "name" : "Tom",
            "age" : 12
        }
    def __getitem__(self,i):
        return self.Tom[i]

a = test()
print(a['name'])
#output
Tom
```

### with

```python
f = open('photo.jpg', 'r+')
jpgdata = f.read()
f.close()
```

`open` 的返回值是一个文件句柄，从操作系统托付给你的 Python 程序。一旦你处理完文件，你会想要归还这显式地调用 `close` 关闭了这个文件句柄，但前提是只有在 read 成功的情况下。如果有任意异常正好在 `f = open(...)` 之后产生，`f.close()` 将不会被调用个文件句柄，只有这样你的程序不会超出一次能打开的文件句柄的数量上限。

为了确保不管异常是否触发，文件都能关闭，我们将其包裹成一个 `with` 语句:

```python
with open('photo.jpg', 'r+') as f:
    jpgdata = f.read()
```

`mode` 打开模式见前

### for-else语句

`else`从句会在循环**正常结束**的时候执行, `else`可以判断循环结束是因为以下哪些原因导致的

- 当一个元素被找到，`break` 被触发。

- 是循环结束。  

```python
for n in range(2, 10):
    for x in range(2, n):
        if n % x == 0:
            print(n, 'equals', x, '*', n / x)
            break
    else:
        # loop fell through without finding a factor
        print(n, 'is a prime number')
#output
2 is a prime number
3 is a prime number
4 equals 2 * 2.0
5 is a prime number
6 equals 2 * 3.0
7 is a prime number
8 equals 2 * 4.0
9 equals 3 * 3.0
```

###  hypot()

返回欧几里德范数 sqrt(x*x + y*y)  (计算直角三角形边长)

```python
import math

math.hypot(x, y)
math.hypot(3, 4)
#output 
5
```

