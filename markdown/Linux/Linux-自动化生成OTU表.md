# 自动化生成OTU表

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
#默认使用grengene数据库，可以自行指定，数据库放在当前目录database_greengene/97_otus.fasta 
#0.97聚类
#默认16线程
#支持自动化安装usearch10和vsearch ！默认仅centos适用！
#原始数据解压到ori_data目录下，返回上一级目录运行otu.sh
#名称：otu_tool
##############
```

## 脚本

```shell
#!/bin/bash
#当前目录：/root/data
#定义目录
dir=$(pwd)
date=$(date +%y%m%d)
function mk(){
	if [ $? -eq 0 ] 
	then
		echo "命令执行成功！"
		sleep 1
	else
		echo "安装失败！"
		exit
	fi
}
#usearch检查
function check_usearch(){
	usearch &> /dev/null
	if [ $? -eq 0 ] 
		then
			echo "usearch检查通过"
		else
			echo "开始部署！"
			cp -a ${dir}/soft/usearch* /usr/bin/usearch
			if [ $? -eq 0 ] 
				then
					chmod 755 /usr/bin/usearch
					usearch
					mk
					echo "usearch部署成功！"
				else
					echo "usearch部署失败！"
					exit
			fi
	fi
}

function format_fastq(){
	ls *.fastq > ls.txt
	echo "开始替换格式"
	file_num
	for seq_format in  $(cat ls.txt)
		do
			sed -i 's/^+*/+/g' ${seq_format}
			count=$(( $count + 1 ))
		done
}

#vsearch检查
function check_vsearch(){
	vsearch &> /dev/null
	if [ $? -eq 0 ] 
		then
			echo "usearch检查通过"
		else
			echo "开始部署！"
			echo "开始安装依赖！autoconf libtool automake gcc"
			yum -y install autoconf libtool automake gcc*
			mk
			cd ${dir}/soft/vsearch-2.7.1 
			./autogen.sh
			mk
			./configure
			mk
			make
			mk
			make install
			mk
			vsearch &> /dev/null
			mk

	fi
}

#定义功能函数——序列拼接/质控/聚类/代表序列/注释/结果文件存放/临时文件删除
function mergepairs(){
	usearch \
	-fastq_mergepairs ${seq} \
	-fastqout ${dir}/result_${date}/mergepairs_result/${seq} \
	-relabel @ &>> ${dir}/result_${date}/ana_log
}
function fastq_filter(){
	usearch \
	-fastq_filter ${seq_merg} \
	-fastq_trunclen 250 \
	-fastq_maxee 0.5 \
	-fastaout ${dir}/result_${date}/filter_result/${seq_merg} &>> ${dir}/result_${date}/ana_log
}
function otu97(){
	vsearch \
	-usearch_global ${dir}/result_${date}/filter_${date} \
	--db ${dir}/database_greengene/97_otus.fasta \
	--id 0.97 \
	--otutab ${dir}/result_${date}/otu_filter/otutable_0.97.txt \
	--threads 16
}
function rep_seq(){
	usearch \
	-fastx_getseqs ${dir}/database_greengene/97_otus.fasta \
	-labels ${dir}/result_${date}/otu_filter/otu_id.txt \
	-fastaout ${dir}/result_${date}/otu_filter/rep_seq.fa
}
function note(){
	awk \
	'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $1,a[$1]}' \
	${dir}/database_greengene/97_otu_taxonomy.txt \
	${dir}/result_${date}/otu_filter/otu_id.txt > ${dir}/result_${date}/otu_filter/otu.tax
}

function save_result(){
	if  [ ! -d ${dir}/result_${date}/ ]
		then
			mkdir -p ${dir}/result_${date}/mergepairs_result/
			mkdir -p ${dir}/result_${date}/filter_result/
			mkdir -p ${dir}/result_${date}/otu_filter
			touch ${dir}/result_${date}/ana_log
			echo "开始分析！"
			sleep 1
		else
			echo "结果存放目录已存在！开始分析！"
			sleep 1
	fi
}

function tmp_file_del(){
	rm -rf otu_final_tmp.txt
	rm -rf otu_tmp.tax8
	rm -rf otu_id.txt
	rm -rf otu.tax
	rm -rf ${dir}/result_${date}/filter_${date}
	rm -rf ${dir}/ori_data/ls.txt
	rm -rf ${dir}/result_${date}/mergepairs_result/ls.txt
	rm -rf ${dir}/result_${date}/filter_result/ls.txt
}
function file_num(){
	count=0
	total=$(wc -l ls.txt | cut -d " " -f 1)
}

#检查分析环境
check_usearch
check_vsearch
#创建结果存放目录
save_result

#获取文件列表
cd ${dir}/ori_data
ls *R1.fastq > ls.txt
#序列格式更改
#format_fastq
#序列拼接
cd ${dir}/ori_data
file_num
for seq in $(cat ls.txt)
	do
		mergepairs
		count=$(( $count + 1 ))
		echo -ne "序列拼接：共${total}个文件，已完成${count}个\r"
	done
mk
sleep 2

#质控
cd ${dir}/result_${date}/mergepairs_result/
ls *.fastq > ls.txt
file_num
for seq_merg in $(cat ls.txt)
	do
		fastq_filter
		count=$(( $count + 1 ))
		echo -ne "序列质控：共${total}个文件，已完成${count}个\r"
	done
mk

#合并
cd ${dir}/result_${date}/filter_result/
ls *.fastq > ls.txt 
echo "
开始合并序列..."
for fa in $(cat ls.txt)
	do
		cat ${fa} >> ../filter_${date}
	done
mk

#改处理格式
sed -i 's/-/z/g' ${dir}/result_${date}/filter_${date}
cd ${dir}/result_${date}

echo "
开始聚类！相似度97%
"
#聚类
otu97
echo "物种注释"
#物种注释
#得到otuid
cd ./otu_filter
cut -f 1 ./otutable_0.97.txt | tail -n +2 > ${dir}/result_${date}/otu_filter/otu_id.txt

echo "获得代表序列"
#代表序列
rep_seq

echo "匹配物种信息"
#匹配物种信息
note
cd ${dir}/result_${date}/otu_filter/

#调整结构
sed 's/; /\t/g' otu.tax | sed '1 i ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' > otu.tax8
cut -f 2-8 otu.tax8 > otu_tmp.tax8

echo "合并物种注释和otu表，并除去叶绿体，线粒体"
#合并物种注释和otu表，并除去叶绿体，线粒体
paste otutable_0.97.txt otu_tmp.tax8 > otu_final_tmp.txt
grep -vi 'Archea' otu_final_tmp.txt | grep -vi 'mitochondria' | grep -vi 'chloroplast' > otu_0.97final.txt
echo "删除临时文件"
tmp_file_del
echo "分析完毕！请查看结果目录result_${date}"
exit
```

