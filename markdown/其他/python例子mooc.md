# python例子(mooc)

###### 题目：判断一个整数是否是回文数。回文数是指正序（从左向右）和倒序（从右向左）读都是一样的整数。示例1：输入: 121输出: true示例2：输入: -121输出: false解释: 从左向右读, 为 -121 。从右向左读, 为 121- 。因此它不是一个回文数示例3：输入: 10输出: false解释: 从右向左读, 为 01 。因此它不是一个回文数

```python
class solution:
    #def __init__(self,number):
        #self.number=number
    
    #判断是否是回文数
    def res_num(self,number):
        res_number=""
        for i in range(len(str(number))):
            res_number = res_number + str(number)[-1-i]
        if eval(res_number) == number:
            print("{}是回文数".format(number))
        else:
            print("{}不是是回文数".format(number))     
                
test = solution()
test.res_num(123)
test.res_num(121212121)
```

