# 容器 Collections

Python 附带一个模块，它包含许多容器数据类型，名字叫作 `collections`

```
defaultdict
counter
deque
namedtuple
enum.Enum
```

## defaultdict

作用：使用普通的字典时，用法一般是`dict={}`,添加元素的只需要`dict[element] =value`即，调用的时候也是如此，`dict[element] = xxx`,但前提是`element`字典里，如果不在字典里就会报错

```python
strings = ('puppy', 'kitten', 'puppy', 'puppy',
           'weasel', 'puppy', 'kitten', 'puppy')
counts = {}
#不设定默认值 会报错
for kw in strings:
    counts[kw] += 1
#KeyError
```

如果想要不报错，就要使用字典的`setdefault`方法设定默认值后再进行后续操作

```python
for kw in strings:
    counts.setdefault(kw,0)
    counts[kw] += 1
```

python标准库中提供了`collection`包中的`defaultdict`本身提供了默认值的功能

```python
dd = collections.defaultdict(list)
dd["foo"].append(1)
dd["foo"].append(2)
dd["bar"].append(3)
dd["bar"].append(4)
print(dd)
#----output----不会再报错
defaultdict(<class 'list'>, {'foo': [1, 2], 'bar': [3, 4]})

#其中,
#collections.defaultdict后面可以跟的参数，如，list int string dict  .... 设定了不同的参数后，
#defaultdic中values的数据类型则为设置的那些，同时可以调用相应的属性与方法
#list默认为[]   int默认0    string默认""".....
```

回到第一个例子，如果用`default`实现

```python
cc = collections.defaultdict(int)
for kw in strings:
    cc[kw] += 1
print(cc)
#----output----
defaultdict(<class 'int'>, {'puppy': 5, 'kitten': 2, 'weasel': 1})
```

可以用一个函数的返回值作为默认值，函数应该不带参数

```python
def test1():
    return [1,3,4]
dd = collections.defaultdict(test1)
print(dd)
print(dd["a"])
#----output----
defaultdict(<function test1 at 0x7f9bad6ed488>, {})
[1, 3, 4]
##使用lambda
dd = collections.defaultdict(lambda: [1,2,3])
print(dd)
print(dd["a"])
#----output----
defaultdict(<function <lambda> at 0x7efd3d592510>, {})
[1, 2, 3]
```

## counter

Counter 是一个计数器，它可以帮助我们针对某项数据进行计数，避免循环。

```python
aa = ["a","s","a",1,2,2,2,2]
#创建一个Counter
fav = collections.Counter(aa)
print(fav)
Counter({2: 4, 'a': 2, 's': 1, 1: 1})

#计数碱基分布
cc = "AATTGGC"
fav2 = collections.Counter(cc)
print(fav2)
Counter({'A': 2, 'T': 2, 'G': 2, 'C': 1})
```

## deque

deque 提供了一个双端队列，你可以从头/尾两端添加或删除元素	`appendleft`和`popleft`，列表其他的方法均可以用

```python
import collections
a = [1,2,3]
#将普通列表转换为deque
d = collections.deque(a)
print(d)
deque([1, 2, 3])

#左边添加元素
d.appendleft("left")
print(d)
deque(['left', 1, 2, 3])

#左边删除元素
d.popleft()
print(d)
deque([1, 2, 3])
```

我们也可以限制这个列表的大小，当超出你设定的限制时，数据会从对队列另一端被挤出去（pop）

```python
d = deque(maxlen=30)
```

## namedtuple

一个元组是一个不可变的列表，你可以存储一个数据的序列，**它和命名元组（`namedtuples`）非常像，但有几个关键的不同**

主要相似点是都不像列表，你不能修改元组中的数据。为了获取元组中的数据，你需要使用整数作为索引

**`namedtuple`**: 把元组变成一个针对**简单任务的容器**。你**不必使用整数索引来访问一个 `namedtuples` 的数据**。你可以像**字典**（`dict`）一样访问 `namedtuples`，但 `namedtuples` 是不可变的。

```python
import collections
#几种不同的创建方式 喜欢第二种
message = collections.namedtuple("person_mesage","name age hobby work code")
#message1 = collections.namedtuple("person_mesage",['name','age','hobby','work','code'])
#message3 = collections.namedtuple("person_mesage","name,age,hobby,work,code")

#举例子
#首先创建一个命名元组
Cat_message = collections.namedtuple('Cat_message',field_names = ["name","age","id"])#Cat_message与后面的'Cat_message'需要一致
cat1 = Cat_message(name=["AChangQiuJun"],age = 2,id = "003")
```

由上述例子可以知道**一个命名元组（`namedtuple`）有两个必需的参数。它们是元组名称和字段名称**，属性值在 `namedtuple` 中是不可变的，既使用整数索引，也可以使用名称来访问

```python
##访问属性
print(cat1.age)
print(cat1.id)
#----output----
2
003
#修改属性，如果是列表，可以直接利用append方法添加值
cat1.name.append("XueDing")
cat1.name.append("ChangZai")
print(cat1)
#----output----
['AChangQiuJun', 'XueDing', 'ChangZai']
#其他的属性没法直接修改，所以需要用_replace()方法
cat1 = cat1._replace(age = 5)
print(cat1.age)
#----output----
5
```

可以将一个命名元组转换为字典(`OrderedDict`), `_asdict()`

```python
print(cat1._asdict())
OrderedDict([('name', ['AChangQiuJun', 'XueDing', 'ChangZai']), ('age', 5), ('id', '003')])
```

### 