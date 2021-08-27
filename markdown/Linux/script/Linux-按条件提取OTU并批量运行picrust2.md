# 按条件提取OTU并批量运行picrust2

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
#门水平为例子，任意水平即可（或组合），只需给出信息
#提取指定门水平otu信息，建立每个门水平的文件

#需要以下3/4个文件
#$1 指定分类水平信息文件
#$2 物种注释表
#$3 otu丰度表  
#$4 OTU代表序列（需要运行picrust2时指定）
#4/3个文件放置于同一目录
##############
```

## 脚本

```shell
#!/bin/bash
#定义命令检查功能
function mk(){
	if [ $? -eq 0 ] 
	then
		echo "${i}Success！" >> log2.txt
	else
		echo "${i}FAIL！" >> log2.txt
		exit
	fi
}
#定义自动化picrust2
function picrust2() {
    cd ./otuinfo
    #文件列表获取
    ls otu_* > ls_pic.txt
    #批量picrust2
    for i in $(cat ls_pic.txt)
        do
            picrust2_pipeline.py -s $4 -i ${i} -o pic_${i} -p 16 &>> log.txt
            mk
        done
    rm -rf ls.txt
}
picrust2

#从物种注释文件中提取相应的OTUid fgrep为了关闭正则
echo "提取OTU_ID信息......"
for i in $(cat $1)
	do
		fgrep "${i}" $2 | awk '{print $1}' > ./${i}.txt
		mk
	done

#可以对门对应OTU合并，得到总丰度表
#cat f_* > totalotu.txt

#从根据上一步中的OTU_ID提取丰度信息
echo "提取丰度信息......"
#根据id从丰度表中提取
ls f_* > ls1.txt
for i in $(cat ls1.txt)
	do
		for otu in $(cat ${i})
			do
				grep -w "${otu}" $3 >> ./otu_${i}	
			done	
	done

#对为匹配到的空白行进行删除，同时增加处理信息
echo "删除空白行与合并分组信息......"
ls otu_* > ls2.txt
for file in $(cat ls2.txt)
	do
		 sed -i '/^$/d' ${file}
	done
fir=$(head -n 1 oturare.txt)
for i in $(cat ls2.txt)
	do
		sed -i "1i ${fir}" ${i}
	done
#是否后续进行picrust2，取消注释即可运行自动化picrust2
#picrust2

#建立结果存放文件并删除临时文件
echo "结果文件生成与临时文件删除......"
mkdir finfo
mkdir otuinfo
mv ./f_* ./finfo
mv ./otu_* ./otuinfo
rm -rf ls*
```

