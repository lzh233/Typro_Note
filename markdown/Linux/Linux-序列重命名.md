# 序列重命名

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
#0、保证名称和序列是对应的，三个文件放在同一个文件夹下面，每条序列、名称保证为单独一行
#1、名称放在name.txt文件内
#2、序列放在seq.fasta中，notepad或者记事本都可以
#4、seq_name.fasta是结果文件
#5、序列不能太多（5000左右）
#名称：seq_rename
##############
```

## 脚本

```shell
#!/bin/bash
#Date:20201020
#Author:lzh
#定义行号变量，line_num 
#当前序列替换条数，seq_count
#当前序列替换条数，seq_count
#总条数，seq
#序列名，seq_name
#输入文件：seq.fasta、name.txt
#输出文件：seq_name.fasta

line_num=1
seq_count=0
cp -a ./seq.fasta ./seq_name.fasta
cp -a ./name.txt ./name.txt.tmp
seq=$(wc -l seq_name.fasta | awk '{print $1}')

#名称开头加大于号
function format(){
	sed -i 's/^/\>/g' name.txt.tmp
}

#批量命名，按照行号1 3 5 7 ....插入
function rename(){
	for seq_name in $(cat name.txt.tmp)
		do
			sed -i "${line_num} i ${seq_name}" seq_name.fasta
			line_num=$(($line_num+2))
			seq_count=$(($seq_count+1))
			echo -ne "共${seq}条序列，已经命名${seq_count}条\r"
		done
}
#若名称提前加入“>”,在语句format前加#---> #format
#format
rename
rm -rf name.txt.tmp
```

