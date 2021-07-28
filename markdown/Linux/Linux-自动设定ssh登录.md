# linux自动设定ssh登录

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
#名称：sshd.sh
##############
```

## 脚本

```shell
#!/bin/bash
#DATE:20210115
#Author:lzh
function mk_add(){
	if [ $? -eq 0 ] 
	then
		echo "公钥信息写入成功！"
	else
		echo "写入失败！请联系管理员！" && exit
	fi
}
cd $HOME
echo "直接输入回车键即可！（3次）
"
#公钥生成
ssh-keygen -t rsa
#写入公钥信息
cd $HOME/.ssh
cat id_rsa.pub >> $HOME/.ssh/authorized_keys
mk_add
chmod 600 authorized_keys
#备份公钥和私钥以及公钥信息文件
cp -a authorized_keys ../\.authorized_keys.bak
cp -a id_rsa ../\.id_rsa.bak
cp -a id_rsa.pub ../\.id_rsa.pub.bak
cp -a id_rsa ../id_rsa_$(whoami)
```

