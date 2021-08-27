# 批量Picrust1

## 运行环境

```
Author:wangwe&liuzihao
Date:2021-4-9
CentOS Linux release 7.8.2003 (Core)	
Architecture:          x86_64
CPU op-mode(s):        32-bit, 64-bit
Byte Order:            Little Endian
CPU(s):                144
Model name:            Intel(R) Xeon(R) Gold 5220 CPU @ 2.20GHz
```

## 脚本

```shell
#Author: liuzihao
#Date: 20210604
#Discreption: Fast Picrust1

#!/bin/bash

#check
function mk(){
	if [ $? -eq 0 ] 
	then
		echo "${otu}.txt Success！" >> ./result/summary.txt
	else
		echo "${otu}.txt FAIL！" >> ./result/summary.txt
		exit
	fi
}

function getpic1(){
	#convert
	biom convert -i ${otu}.txt -o ./result/${otu}/otu_table2.biom --table-type="OTU table" --to-json &>> ./result/log.txt
	mk
	
	#normalize
	normalize_by_copy_number.py -i ./result/${otu}/otu_table2.biom -o ./result/${otu}/normalized_otus.biom &> ./result/log.txt
	mk
	
	#predict
	predict_metagenomes.py -i ./result/${otu}/normalized_otus.biom -o ./result/${otu}/metagenome_predictions.biom &> ./result/log.txt
	mk
	
	#biom2tsv
	biom convert -i ./result/${otu}/metagenome_predictions.biom -o ./result/${otu}/metagenome_predictions.txt --table-type="OTU table" --to-tsv &> ./result/log.txt
	mk
	
	#kegg note
	categorize_by_function.py -i ./result/${otu}/metagenome_predictions.biom -c KEGG_Pathways -l 3 -o ./result/${otu}/metagenome_predictions.L3.biom &> ./result/log.txt
	mk
	
	#get result
	biom convert -i ./result/${otu}/metagenome_predictions.L3.biom -o ./result/${otu}/metagenome_predictions.L3.txt --table-type="OTU table" --to-tsv &> ./result/log.txt
	mk
}

function file_num(){
	count=0
	total=$(wc -l ls.txt | cut -d " " -f 1)
}


####picrust1#####

#get file list

ls *.txt >  ls.txt
cat ls.txt | cut -d "." -f 1 > ls2.txt

mkdir result

#file number get
file_num

for otu in $(cat ls2.txt)
	do
		mkdir ./result/${otu}
		count=$(( $count + 1 ))
		echo -ne "：Total: ${total} files Finish ${count} \r"
		getpic1		
	done
```

