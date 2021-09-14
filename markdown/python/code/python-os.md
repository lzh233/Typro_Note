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

