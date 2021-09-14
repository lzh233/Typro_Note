# itertools

Python的内建模块`itertools`提供了非常有用的用于操作迭代对象的函数。

## count

count在没有中止条件的时候会无限输出所有的数字，可以指定步长

```python
import itertools

a = itertools.count(1)
for i in a:
    print(i)
#output
1
2
3
....

#指定步长
a = itertools.count(1,10)
for i in a:
    print(i)
1
11
21
...
```

## cycle

`cycle()`会把传入的一个序列无限重复下去, 列表则会将每个元素无限循环

```python
cs = itertools.cycle('ABC') # 注意字符串也是序列的一种
for c in cs:
     print(c)
'A'
'B'
'C'
'A'
'B'
'C'
...
```

## repeat

`repeat()`负责把一个元素无限重复下去，不过如果提供第二个参数就可以限定重复次数

```python
ns = itertools.repeat('A', 3)
for n in ns:
     print(n)
A
A
A
```

## chain

`chain()`可以把一组迭代对象串联起来，形成一个更大的迭代器：

```python
for c in itertools.chain('ABC', 'XYZ'):
     print(c)
# 迭代效果：'A' 'B' 'C' 'X' 'Y' 'Z'
```

## groupby

`groupby()`把迭代器中相邻的重复元素挑出来放在一起

```
>>> for key, group in itertools.groupby('AAABBBCCAAA'):
...     print(key, list(group))
...
A ['A', 'A', 'A']
B ['B', 'B', 'B']
C ['C', 'C']
A ['A', 'A', 'A']
```

实际上挑选规则是通过函数完成的，只要作用于函数的两个元素返回的值相等，这两个元素就被认为是在一组的，而函数返回值作为组的key。如果我们要忽略大小写分组，就可以让元素`'A'`和`'a'`都返回相同的`key`：

```python
>>> for key, group in itertools.groupby('AaaBBbcCAAa', lambda c: c.upper()):
...     print(key, list(group))
...
A ['A', 'a', 'a']
B ['B', 'B', 'b']
C ['c', 'C']
A ['A', 'A', 'a']
```