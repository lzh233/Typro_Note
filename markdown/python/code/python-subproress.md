# subproress

[来源](https://www.jianshu.com/p/2eb33b491024)

从`Python 2.4`开始，Python引入`subprocess`模块来管理子进程，以取代一些旧模块的方法：如 `os.system`、`os.spawn`、`os.popen`、`popen2.`、`commands.`不但可以调用外部的命令作为子进程，而且可以连接到子进程的`input/output/error`管道，获取相关的返回信息

运行python的时候，我们都是在创建并运行一个进程。像Linux进程那样，一个进程可以fork一个子进程，并让这个子进程exec另外一个程序。在Python中，我们通过标准库中的`subprocess`包来`fork`一个子进程，并运行一个外部的程序。 `subprocess`包中定义有数个创建子进程的函数，这些函数分别以不同的方式创建子进程，所以我们可以根据需要来从中选取一个使用。

#### subprocess.call()

```python
父进程等待子进程完成
返回退出信息(returncode，相当于Linux exit code)
```

#### subprocess.check_call()

```shell
父进程等待子进程完成
返回0(运行正确)
检查退出信息，如果returncode不为0，则举出错误subprocess.CalledProcessError，该对象包含有returncode属性，可用try…except…来检查
```

举例

```python
import subprocess
subprocess.check_call('ls -l',shell=True)
##
subprocess.check_call('ls -1wd',shell=True)
#输出
subprocess.CalledProcessError: Command 'ls -1wd' returned non-zero exit status 2.
```
