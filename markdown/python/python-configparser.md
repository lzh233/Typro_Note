# configparser

ConfigParser模块在python中用来读取配置文件，配置文件的格式跟windows下的`ini`配置文件相似，可以包含一个或多个节(section), 每个节可以有多个参数（键=值），可以使程序更灵活。

```shell
#ini风格的配置文件
[MainInfo]
ProtoSum=5
ProtoId.0=8021
ProtoId.1=1100
......
[8021]
SID=X1
ProtoName=802.1X
;2052 -- 语言ID，代表“简体中文”
Desciption_2052=关于802.1X协议的描述信息在此
.....
```

## 使用方法

示例文件如下,

```shell
[ANNOVAR]
dir = /Public/Software/annovar/
db = /SGRNJ/Database/script/database/annovar/humandb
buildver = hg38
protocol = refGene,cosmic70
operation = g,f
[DATA]
path = ./test_file/
size = 120
```

使用

```python
import configparser
#创建ConfigParser对象
conf = configparser.ConfigParser()
#读入配置文件
conf.read("annovar.config")
#配置文件中的[]内的内容称为section，使用section()方法显示配置文件中所有的section
conf.section()		#['ANNOVAR', 'DATA']

#配置文件section下的内容称为option，使用option()显示某section下的option
conf.options("ANNOVAR")		#返回option列表: ['dir', 'db', 'buildver', 'protocol', 'operation']

#获取指定section的每个option的键值对，使用item()方法
conf.items("ANNOVAR")		#[('dir', '/Public/Software/annovar/'), ('db', '/SGRNJ/Database/script/database/annovar/humandb'), ('buildver', 'hg38'), ('protocol', 'refGene,cosmic70'), ('operation', 'g,f')]

#获取所有指定section的option信息
conf.get("ANNOVAR","db")		#/SGRNJ/Database/script/database/annovar/humandb

```

写一个配置文件

**所有`option`的值必须转换为字符串格式**

```python
import configparser
config = configparser.ConfigParser()
config['library'] = {}
config['library']['lib1'] = "path1"
config['library']['lib2'] = "path2"

config['DATA'] = {}
config['DATA']['size'] = str(1024)
config['DATA']['dir']  = "test/dir/"

with open("configure.conf","w") as config_data:
    config.write(config_data)
```

查看输出的配置文件

```python
[library]
lib1 = path1
lib2 = path2

[DATA]
size = 1024
dir = test/dir/
```

