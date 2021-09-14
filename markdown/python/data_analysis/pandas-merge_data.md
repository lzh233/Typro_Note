# Pandas-Merge data

## 基本的数据增删查改



## 数据合并与连接

### 数据合并

使用`pd.concat`实现简单合并`Series`对象

```python
#构建两个Series
ser1 = pd.Series(['A', 'B', 'C'], index=[1, 2, 3])
ser2 = pd.Series(['D', 'E', 'F'], index=[7, 8, 9])
#合并
pd.concat([ser1,ser2])
#output
1    A
2    B
3    C
7    D
8    E
9    F
#通过修改axis参数实现按列合并，默认按行合并
pd.concat([ser1,ser2],axis = 1)
#output
     0    1
1    A  NaN
2    B  NaN
3    C  NaN
7  NaN    D
8  NaN    E
9  NaN    F
```

使用`pd.concat()`合并`DataFrame`

```python
#构建两个DataFrame
df1 = pd.DataFrame({"grammer":["python","C","Go"],"score":[1,2,3]})
df2 = pd.DataFrame({"grammer":["Java","R","C++"],"score":[4,6,21]})
print(df1)
print(df2)
#output
  grammer  score
0  python      1
1       C      2
2      Go      3

  grammer  score
0    Java      4
1       R      6
2     C++     21

#合并数据 默认按行合并, ignore_index = True 表示忽略原index 重新构建index
print(pd.concat([df1,df2],ignore_index=True))
#output
  grammer  score
0  python      1
1       C      2
2      Go      3
3    Java      4
4       R      6
5     C++     21
#按列合并
print(pd.concat([df1,df2],axis = 1,ignore_index=True))
#output
        0  1     2   3
0  python  1  Java   4
1       C  2     R   6
2      Go  3   C++  21

#增加多级索引
print(pd.concat([df1,df2],ignore_index=False,keys=["A","B"]))
#output
    grammer  score
A 0  python      1
  1       C      2
  2      Go      3
B 0    Java      4
  1       R      6
  2     C++     21
```

### 数据连接(关系型数据合并)

关系型数据的合并主要通过`pd.merge`实现
pd.merge() 函数实现了三种数据连接的类型：
一对一、 多对一和多对多。 这三种数据连接类型都通过 `pd.merge() `接口进行调用， 根据不同的数据连接需求进行不同的操作。

#### 一对一连接

 关系型数据中，某列的列名完全和另一个表的列名与内容完全对应，对这样的数据merge是一对一连接

```python
df1 = pd.DataFrame({'employee': ['Bob', 'Jake', 'Lisa', 'Sue'],'group': ['Accounting', 'Engineering', 'Engineering', 'HR']})
df2 = pd.DataFrame({'employee': ['Lisa', 'Bob', 'Jake', 'Sue'],'hire_date': [2004, 2008, 2012, 2014]})
print(f"{df1}\n\n{df2}")
#output
  employee        group
0      Bob   Accounting
1     Jake  Engineering
2     Lisa  Engineering
3      Sue           HR

  employee  hire_date
0     Lisa       2004
1      Bob       2008
2     Jake       2012
3      Sue       2014

#数据合并,会自动发现两组数据的共有列employee,
"""
pd.merge() 方法会发现两个 DataFrame 都有“employee”列， 并会自动以这列作为键进行连接。 两个输入的合并结果是一个新的DataFrame。 
需要注意的是， 共同列的位置可以是不一致的。
例如在这个例子中， 虽然 df1 与 df2 中“employee”列的位置是不一样的， 但是 pd.merge() 函数会正确处理这个问题。 
另外还需要注意的是， pd.merge() 会默认丢弃原来的行索引， 不过也可以自定义
"""
df3 = pd.merge(left = df1, right = df2)
#output
  employee        group  hire_date supervisor
0      Bob   Accounting       2008      Carly
1     Jake  Engineering       2012      Guido
2     Lisa  Engineering       2004      Guido
3      Sue           HR       2014      Steve
```

#### 多对一连接

在需要连接的两个列中， 有一列的值有重复。如本例子，将df4的信息（supervisor）加入到df3中,**类似物种注释**

```python
df4 = pd.DataFrame({'group': ['Accounting', 'Engineering', 'HR'],'supervisor': ['Carly', 'Guido', 'Steve']})
print(f"{df3}\n\n{df4}")
#output
  employee        group  hire_date
0      Bob   Accounting       2008
1     Jake  Engineering       2012
2     Lisa  Engineering       2004
3      Sue           HR       2014

         group supervisor
0   Accounting      Carly
1  Engineering      Guido
2           HR      Steve
#合并数据
pd.merge(left=df3,right=df4)
#output
  employee        group  hire_date supervisor
0      Bob   Accounting       2008      Carly
1     Jake  Engineering       2012      Guido
2     Lisa  Engineering       2004      Guido
3      Sue           HR       2014      Steve
```

#### 多对多连接

如果左右两个输入的共同列都包含重复值， 那么合并的结果就是一种多对多连接。要合并的列 同一个名称对应多个数据, 被并表格对应的列也是存在重复，下面的例子， 里面有一个 DataFrame 显示不同岗位人员的一种或多种能力。

```python
df5 = pd.DataFrame({'group': ['Accounting', 'Accounting','Engineering', 'Engineering', 'HR', 'HR'],
                    'skills': ['math', 'spreadsheets', 'coding', 'linux','spreadsheets', 'organization']})
print(f"{df1}\n\n{df5}")
#output
  employee        group
0      Bob   Accounting
1     Jake  Engineering
2     Lisa  Engineering
3      Sue           HR

         group        skills
0   Accounting          math
1   Accounting  spreadsheets
2  Engineering        coding
3  Engineering         linux
4           HR  spreadsheets
5           HR  organization
#合并数据
pd.merge(left=df1,right=df5)
#output
  employee        group        skills
0      Bob   Accounting          math
1      Bob   Accounting  spreadsheets
2     Jake  Engineering        coding
3     Jake  Engineering         linux
4     Lisa  Engineering        coding
5     Lisa  Engineering         linux
6      Sue           HR  spreadsheets
7      Sue           HR  organization
```

#### 不同列名的数据连接

合并两个列名不同的数据集， 例如前面的员工信息表中有一个字段不是`employee`而是`name`。 在这种情况下， 就可以用
`left_on` 和 `right_on` 参数来指定列名

```python
df1 = pd.DataFrame({'employee': ['Bob', 'Jake', 'Lisa', 'Sue'],'group': ['Accounting', 'Engineering', 'Engineering', 'HR']})
df6 = pd.DataFrame({'name': ['Bob', 'Jake', 'Lisa', 'Sue'],'salary': [70000, 80000, 120000, 90000]})
print(f"{df1}\n\n{df6}")
#output
  employee        group
0      Bob   Accounting
1     Jake  Engineering
2     Lisa  Engineering
3      Sue           HR

   name  salary
0   Bob   70000
1  Jake   80000
2  Lisa  120000
3   Sue   90000
#数据连接，使用left_on和right_on指定连接键,合并后会出现多余的一列，使用drop()方法去掉
pd.merge(left=df1,left_on="employee",right = df6,right_on="name").drop('name',axis=1)
#output
  employee        group  salary
0      Bob   Accounting   70000
1     Jake  Engineering   80000
2     Lisa  Engineering  120000
3      Sue           HR   90000
```

