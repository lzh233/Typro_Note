# logging模块

[logging详解](https://blog.csdn.net/pansaky/article/details/90710751)

`basicConfig`: 日志输出的基本设置，`level`为输出级别（输出级别从高到低为：FATAL、ERROR、WARN、INFO、DEBUG）级别设定好以后，只会输出设定级别之上的日志，`format`为输出格式

常见的`format`内容如下

```python
"""
%(levelno)s：打印日志级别的数值
%(levelname)s：打印日志级别的名称
%(pathname)s：打印当前执行程序的路径，其实就是sys.argv[0]
%(filename)s：打印当前执行程序名
%(funcName)s：打印日志的当前函数
%(lineno)d：打印日志的当前行号
%(asctime)s：打印日志的时间
%(thread)d：打印线程ID
%(threadName)s：打印线程名称
%(process)d：打印进程ID
%(message)s：打印日志信息
"""
```

具体日志级别

```python
"""
日志等级：使用范围
FATAL：致命错误
CRITICAL：特别糟糕的事情，如内存耗尽、磁盘空间为空，一般很少使用
ERROR：发生错误时，如IO操作失败或者连接问题
WARNING：发生很重要的事件，但是并不是错误时，如用户登录密码错误
INFO：处理请求或者状态变化等日常事务
DEBUG：调试过程中使用DEBUG等级，如算法中每个循环的中间状态
""""
```

基本的应用：

```python
import logging
logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

logger = logging.getLogger("TEST!!!")
logger.info("Start print log")
logger.debug("Do something")
logger.warning("Something maybe fail.")
logger.info("Finish")
#--------=output---------
2021-09-03 09:40:31,363 - TEST!!! - INFO - Start print log
2021-09-03 09:40:31,364 - TEST!!! - WARNING - Something maybe fail.
2021-09-03 09:40:31,364 - TEST!!! - INFO - Finish
```

将日志输出到文件: 

设置`logging`，创建一个`FileHandler`，并对输出消息的格式进行设置，将其添加到`logger`，然后将日志写入到指定的文件中

```python
#创建一个日志
logger = logging.getLogger(__name__)
#设定等级
logger.setLevel(level = logging.INFO)
#增加FileHandler
handler = logging.FileHandler("log.txt")
#设定等级与格式
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
#将设定好的FileHandler添加进logger主对象中
logger.addHandler(handler)
#运行，就会输出到log.txt中
logger.info("Start print log")
logger.debug("Do something")
logger.warning("Something maybe fail.")
```

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210903095457896.png" alt="image-20210903095457896" style="zoom:80%;" />

将日志输出到文件和屏幕:

logger中同时添加`StreamHandler`，可以将日志输出到屏幕上


```python
logger = logging.getLogger(__name__)
logger.setLevel(level = logging.INFO)

handler = logging.FileHandler("log.txt")
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

#屏幕输出设定
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logger.addHandler(console)
```

