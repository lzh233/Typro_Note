# conda环境变量自动修改

## 运行环境

```
Author:liuzihao
Date:2021-4-9
CentOS Linux release 7.8.2003 (Core)	
Architecture:          x86_64
CPU op-mode(s):        32-bit, 64-bit
Byte Order:            Little Endian
CPU(s):                144
Model name:            Intel(R) Xeon(R) Gold 5220 CPU @ 2.20GHz
```

## 说明

```shell
##############
#任意位置直接运行
##############
```

## 脚本

```shell
#!/bin/bash
#2020-01-11
#author：lzh
#获取用户所属组和家目录
grup=$(whoami | groups)

#检查
function mk(){
	if [ $? -eq 0 ] 
	then
		echo "环境变量修改成功！"
	else
		echo "修改失败！" && exit
	fi
}
#修改
function env_change(){
	cp -a ${HOME}/.bash_profile ${HOME}/.bash_profile.bak
	if [ "${grup}" = "student" ]
		then
			sed -i '/^PATH/{s/$/:\/home\/student\/miniconda2\/bin\//}' ${HOME}/.bash_profile
	elif [ "${grup}" = "teacher" ]
		then
			sed -i '/^PATH/{s/$/:\/home\/teacher\/miniconda2\/bin\//}' ${HOME}/.bash_profile
	else
		echo "您不属于student或teacher用户组，如有需要请手动修改.bash_profile文件" && exit
	fi
}

#修改和生效以及删除脚本
echo "您所属用户组为${grup},家目录为${HOME}"
env_change
mk
source ${HOME}/.bash_profile
mk
echo "配置文件已生效！(如果输入conda命令提示command not found，请输入：source $HOME/.bash_profile，回车运行命令)"	
rm -rf ${HOME}/env_change.sh
```