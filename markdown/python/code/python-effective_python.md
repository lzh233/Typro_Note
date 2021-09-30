# Effective python

### 1. 使用生成器表达式改写数据量较大的列表推导式

对于较大的文件来说，如果直接采用列表推导式来做，会大量的消耗内存，导致程序崩溃，python提供了生成器来解决这个问题，他是对列表推导和生成器的一种泛化，生成器在运行的时候，并不会把整个输入序列都展示出来，而是会估值为迭代器，每次可以根据迭代器的表达式生成数据，**但是生成器只能使用一次**，**生成器的好处之一还有可以结合使用, 而且结合在一起的生成器的运行速度非常快**

```python
#将列表推导式转换为生成器，仅仅需要将[]转换为()即可
test1 = [i for i in open("./test.txt")] 
#转换成生成器
test2 = (i for i in open("./test.txt"))
#与其他的生成器结合使用
test3 = (x**4 for x in test1)
#可以使用循环使用生成器或者next()
#循环
for i in test1:
    pass
#next
next(test1)
```

### 2. 异常捕获 try、except、finally

- finally块

  既想将异常向上传播(传播给调用者)，又想在异常发生时执行异常清理工作，就可以使用`try / finally`结构

  ```python
  #eg:
  a = "error run!"
  try:
      print(b)
  finally:
      print(a)
  #既然运行了 finally语句块也将异常向上传播(代码正确运行的时候也会运行finally的语句块)
  error run!
  ......
  NameError: name 'b' is not defined
  #应用: 无论文件是否正确打开, 均能关闭文件句柄
  fd = open("test.txt")
  try:
      data = fd.read()
  except:
      fd.close()
  ```

- else块

  `try / except / else` 结构可以清晰的描述出哪些异常会由代码自己去处理，哪些异常需要传播到上一级，如果`try`代码块没有发生异常，那么就 执行`else`

  ```python
  a = "error run!"
  c = "correct!"
  
  try:
      print(d)
  except NameError as e:
      raise KeyError from e
  else:
      print(f"{a} var a run correct!")
  #output
  Traceback (most recent call last):
    File "d:/Desktop/test/a.py", line 5, in <module>
      print(d)
  NameError: name 'd' is not defined
  
  The above exception was the direct cause of the following exception:
  
  Traceback (most recent call last):
    File "d:/Desktop/test/a.py", line 7, in <module>
      raise KeyError from e
  KeyError
  ```

  <img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210914105150142.png" alt="image-20210914105150142" style="zoom:70%;" />

### 3. 使用zip函数同时对多个迭代器实现迭代

例如，在两个列表中的数据可能是**关联**的时候，可以利用zip函数对两个列表进行迭代，从而避免`for`循环与各种下标的嵌套

```python
name = ["Tom","Jerry","Agua"]
age = [5,6,7]
food = ["FISH","Fruit","beef"]

for name,age,food in zip(name,age,food):
    print(name,age,food)
#output
Tom 5 FISH
Jerry 6 Fruit
Agua 7 beef
```

`python`自带的zip函数有一定的缺陷，当列表的长度不同时，会以最短的列表为为准，

```python
name = ["Tom","Jerry","Agua","Achang"]
age = [5,6,7]
food = ["FISH","Fruit","beef"]

for name,age,food in zip(name,age,food):
    print(name,age,food)
#output, name的最后一个元素并没有生成
Tom 5 FISH
Jerry 6 Fruit
Agua 7 beef
```

为了解决上述缺点可以使用`itertools`库中的`zip_longest()`解决

```python
import itertools

name = ["Tom","Jerry","Agua","Achang"]
age = [5,6,7]
food = ["FISH","Fruit","beef"]
#fillvalue 指定None值为其他值，不指定则为None
for name,age,food in itertools.zip_longest(name,age,food,fillvalue="no_detect"):
    print(name,age,food)
#output
Tom 5 FISH
Jerry 6 Fruit
Agua 7 beef
Achang no_detect no_detect
```

### 4. 可以使用生成器来改写直接返回列表的函数

​	如果函数需要产生一系列的记过，一般的做法就是将函数的结果放到列表中，并将其返回给调用者，但是如果数据非常大的话以及让代码更加简洁，**可以使用生成器**进行改写，即，使用`yield`表达式的函数，调用生成器函数的时候并不会真正的运行，而是返回一个迭代器，每次会运行迭代器内置的`next()`函数，生成器产`yield`产生的值都会由迭代器传给调用者

```python
#eg:统计首字母位置
a = "One two three four"
def test(txt):
    if txt:
        yield 0
    for index,letter in enumerate(txt):
        if letter == " ":
            yield index + 1 
aa = list(test(txt = a))
#可以利用迭代获取和操作数据
#for loc in test(txt = a):
    ......
print(aa)
#output
[0, 4, 8, 14]
```

该函数如果接受的文本非常大，如果不使用生成器改写，那对内存的消耗是非常大的，如果改写为生成器函数，每次的内存消耗则为读入的每行数据占用的内存，而不是所有文本

```python
#对比
import sys
a = "One two three four" * 1000000

def test(txt):
    if txt:
        yield 0
    for index,letter in enumerate(txt):
        if letter == " ":
            yield index + 1 

aa = test(txt = a)
print(sys.getsizeof(aa))

def test2(txt):
    result = []
    if txt:
        result.append(0)
    for index,letter in enumerate(txt):
        if letter == " ":
            if letter == " ":
                result.append(index + 1)
    return result 

bb = test2(txt = a)
print(sys.getsizeof(bb))

#output
size of aa is 112
size of bb is 25105984
```

