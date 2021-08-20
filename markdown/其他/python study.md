# Python

## argparse包的使用

### 用途

通过命令行控制输入变量

### 使用方法

```python
import python
"""
Creating index for reference from fasta
Aligling fastq to index

lunzi 1.0_only for pairs-end
"""
import os
import argparse

from pandas.io.parsers import count_empty_vals 
#定义一个Bowtie2类
class Bowtie2:
    def __init__(self,reference_path,index_path,bam_path):
        self.reference_path = reference_path
        self.index_path = index_path
        self.bam_path = bam_path
    
    def build_index(self):
        cmd1 = f"bowtie2-build {self.reference_path} {self.index_path}"
        return cmd1

    def align_pairs(self,fasta_path,fasta_path_R):
        cmd2 = f"bowtie2  -x {self.index_path} -1 {fasta_path} -2 {fasta_path_R} -S {self.bam_path}"
        return cmd2


def main():
    #创建argparse对象，通过ArgumentParser函数
    parser = argparse.ArgumentParser(description='allign by Bowtie2')
    #通过add_argument方法可以增加选项，帮助信息为help中的内容，ruquire为是否必须选项，每个argparse对象都有一个固定的选项-h，可以输出明命令的帮助信息
    parser.add_argument('--fq1', help='fastq', required=True)
    parser.add_argument('--fq2', help='fastq_r.', required=True)
    parser.add_argument('--reference', help='reference path', required=True)
    parser.add_argument('--index_out', help='output path of indeies', required=True)
    parser.add_argument('--bam_out', help='output path of bam_file', required=True)
    #parser.add_argument('--to_sam', help='to sam?', required=True)
    #将选项信息储存在args中
    args = parser.parse_args()
   #选项的使用
#使用args.<选项的全名(--指定的名称)> = 选项的值
#如args.index_out
#  args.reference
b = Bowtie2(index_path = str(args.index_out),
                reference_path=str(args.reference),
                bam_path = str(args.bam_out))
    
    cmd1 = b.build_index()
    cmd2 =b.align_pairs(fasta_path=args.fq1,
                        fasta_path_R=args.fq2)

    file_handle=open('run.sh',mode='w')
    file_handle.writelines([f'{cmd1}\n',f'{cmd2}'])
    file_handle.close()

if __name__ == '__main__':
    main()

```

```shell
 ##使用方法, 相应的参数就会传到对应的位置
python allign.py --fq1 ./seq/fq1 --fq2 ./seq/fq2 --index_out ./result/out_index --bam_out ./result/bam_out --reference /vdj_ref
```

```python
#输出信息
bowtie2-build /vdj_ref ./result/out_index
bowtie2  -x ./result/out_index -1 ./seq/fq1 -2 ./seq/fq2 -S ./result/bam_out
```

## python文件读取与写入

### 打开和关闭文件

```python
#使用open函数打开文件，创建一个file对象, 可以对文件进行操作
#file object = open(file_name [, access_mode][, buffering])
#常用的mode如下，
```

`t`: 文本模式，默认

`x`: 写模式，新建一个文件，如果该文件已存在则会报错

`b`: 二进制模式

`+`: 打开一个文件进行更新(可读可写)

`r`: 以只读方式打开文件。文件的指针将会放在文件的开头

`w`:  打开一个文件只用于写入。如果该文件已存在则打开文件，并从开头开始编辑，**原有内容会被删除**。如果该文件不存在，创建新文件

`a`: 打开一个文件用于追加。如果该文件已存在，**文件指针将会放在文件的结尾, 新的内容将会被写入到已有内容之后。**如果该文件不存在，创建新文件进行写入。

`a+`:  打开一个文件**用于读写**。如果该文件已存在，文件指针将会放在文件的结尾。文件打开时会是追加模式。如果该文件不存在，创建新文件用于读写

```python
#打开一个文件,追加模式，会在文件末尾进行添加
fd = open(fo = open("01_save_kmerww.txt",mode="a"))
print(f"文件名称: {fd.name}")
print(f"文件访问模式: {fd.mode}")
print(f"是否关闭: {fd.close}")
##output
#文件名称: 01_save_kmerww.txt
#文件访问模式: a
#是否关闭: <built-in method close of _io.TextIOWrapper object at 0x01D60E38>
```

### close方法

```python
# 关闭打开的文件
fo.close()
```

### write方法

```python
fd.write("test")
fd.write("\ntest2")
#--output
#testtest
#test2
```

