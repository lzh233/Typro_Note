# interval

用于判断一个数是否在一个区间内的 操作

https://www.cnblogs.com/cotyb/p/5256303.html

```python
from interval import Interval
#创建区间默认左闭右闭
test1 = Interval(2,5)
print(1 in test1)
False

print(2 in test1)
#output
True

#左开右闭区间
test2 = Interval(2,5,lower_closed = False)
print(2 in test2)
#output
False

#区间合并，判断是否交叉重复及是否比邻
zoom_1_3 = Interval(1,3,lower_closed=False)
zoom_3_5 = Interval(3,5,lower_closed=False)
#合并
zoom_1_5 = zoom_1_3.join(zoom_3_5)
print(zoom_1_5)
#output
#(1..5]

#判断交集
test4 = Interval(1,10)
test5 = Interval(2,3)
print(test4.overlaps(test5))
#output
#True

#比邻是指：1、两个区间相邻的边界值相等 2、同时这两个边界一个为开区间，一个为闭区间，即两个区间没有交集
test6 = Interval(0, 1, upper_closed=False)
test7 = Interval(1, 3, lower_closed=False)
print(test6.adjacent_to(test7))
#output
#False
```

