# os模块

[os模块](https://www.jianshu.com/p/6b699e9cea1f)

```python
import os

#显示当前系统类型
print(os.name)

#显示当前工作目录
print(os.getcwd())

#改变工作目录
os.chdir("/Personal/liuzihao/vscode_use")

#删除文件
os.remove("aaaa.py")

#创建目录
os.mkdir("test")

#删除空目录
os.rmdir("test/")

#运行系统命令
os.system("ls -la")

#分割路径和文件，默认最后一个为文件(没有/)
print(os.path.split("/Personal/liuzihao/vscode_use/"))

#判断给定的是否是目录
print(os.path.isdir(os.getcwd()))

#判断给定的是否是文件
print(os.path.isfile('log.py'))

#判断路径是否存在
print(os.path.exists('C:\\Python25\\abc.txt'))
```

[os.path.dirname](https://www.cnblogs.com/annawong/p/10417257.html)

```python
import os

path1=os.path.abspath(__file__)
print(path1)#当前文件的绝对路径

path2=os.path.dirname(os.path.abspath(__file__))
print(path2)#当前文件的上一层目录的绝对路径

path3=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print(path3)#当前文件的上上层目录的绝对路径

curpath = os.path.dirname(os.path.realpath(__file__)) # 当前文件夹的路径

#实例
import os
path2=os.path.dirname(os.path.realpath(__file__))
print(path2)

"""
语法：os.path.dirname(path) 

功能：去掉文件名，返回目录 

 

__file__ 为内置属性，表示当前文件的path

os.path.dirname((__file__) ：指的是，得到当前文件的绝对路径，是去掉脚本的文件名，只返回目录。

os.path.dirname(os,path.realname(__file__))：指的是，获得你刚才所引用的模块所在的绝对路径
"""
```



