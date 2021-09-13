# Effective python

### 1. 使用生成器表达式改写数据量较大的列表推导式

对于较大的文件来说，如果直接采用列表推导式来做，会大量的消耗内存，导致程序崩溃，python提供了生成器来解决这个问题，他是对列表推导和生成器的一种泛化，生成器在运行的时候，并不会把整个输入序列都展示出来，而是会估值为迭代器，每次可以根据迭代器的表达式生成数据，**但是生成器只能使用一次**

```python
```

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

