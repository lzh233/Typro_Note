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

