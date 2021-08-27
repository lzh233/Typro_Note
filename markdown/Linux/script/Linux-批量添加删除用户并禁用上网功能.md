# 批量添加删除用户并禁用上网功能

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
#添加用户
#用户名和所属组信息写在userlist.txt中，与脚本位于同一目录，运行即可，用户信息格式如下
AAAA1:student
AAAA2:teacher
#名称：add.sh
#批量删除用户
#用户名和所属组信息写在userlist_del.txt中，与脚本位于同一目录，运行即可，用户信息格式如下
AAAA1:student
AAAA2:teacher
#名称：del.sh
##############
```

## 脚本（添加用户并禁用上网功能）

```shell
#!/bin/bash
#20210131lzh
function mk_add(){
	if [ $? -eq 0 ] 
	then
		echo "添加用户${username}成功！"
	else
		echo "添加户${username}失败！" 
		continue
	fi
}
function mk_passwd(){
	if [ $? -eq 0 ] 
	then
		echo "用户${username}初始密码设置成功！"
	else
		echo "添加用户${username}失败！" 
		continue
	fi
}
function mk_net(){
	if [ $? -eq 0 ] 
	then
		echo "已禁用用户${username}上网功能！"
	else
		echo "禁用用户${username}上网功能失败！"	
		continue
	fi
}

function add(){
	for username in $(cat userlist.txt | cut -d : -f 1)
		do
			useradd ${username} -g $(grep "${username}" userlist.txt | cut -d : -f 2 )
			mk_add
			echo "zhiwuyingyang609" | passwd --stdin ${username} &> /dev/null
			mk_passwd
			iptables -I OUTPUT 1 -m owner --uid-owner ${username} -j DROP
			mk_net
		done
}
#防火墙信息备份
iptables-save > /root/\.iptable.bak
add
```

## 脚本（删除用户并初始化防火墙）

```shell
#!/bin/bash
#20210131lzh

function mk_del(){
	if [ $? -eq 0 ] 
	then
		echo "删除用户${username}成功！"
	else
		echo "删除${username}失败！"
		continue
	fi
}


function mk_net(){
	if [ $? -eq 0 ] 
	then
		echo "已初始化防火墙"
	else
		echo "防火请初始化失败" && exit
	fi
}

function del(){
	for username in $(cat userlist_del.txt)
		do
			userdel -r ${username}
			mk_del
		done
}
del
iptables -F
mk_net
```

