# pysam

Pysam是htslib C-API的包装器，提供读写功能`SAM/BAM/VCF/BCF/BED/GFF/GTF/FASTA/FASTQ`文件，以及访问`samtools`和`bcftools`包的命令行功能。该模块支持通过索引进行压缩和随机访问。

## 安装

```python
pip install pysam
```

## 基本使用

### 读入数据

`samfile`是一个迭代器

```python
import pysam
#创建AlignmentFile对象
samfile = pysam.AlignmentFile("ex1.bam", "rb")
#通过index_filename可以指定索引
samfile = pysam.AlignmentFile("ex1.bam", index_filename = "ex1.bam.bai","rb")
```

## 内置方法

### fetch()

使用`fetch()`方法可以抓取到匹配到指定区域的所有`reads`, 返回一个迭代器，迭代器内包含了匹配到该区域的`reads`的所有信息

```python
for read in samfile.fetch('chr1', 100, 120):
    print(read)
```

### write()

使用`write()`方法可以将数据写入`bam`文件

```python
import pysam
samfile = pysam.AlignmentFile("ex1.bam", "rb")
pairedreads = pysam.AlignmentFile("allpaired.bam", "wb", template=samfile)

for read in samfile.fetch():
    if read.is_paired:
        pairedreads.write(read)

pairedreads.close()
samfile.close()
```

### pileup()

```python
import pysam
samfile = pysam.AlignmentFile("ex1.bam", "rb" )
for pileupcolumn in samfile.pileup("chr1", 100, 120):
    print ("\ncoverage at base %s = %s" %
           (pileupcolumn.pos, pileupcolumn.n))
for pileupread in pileupcolumn.pileups:
    if not pileupread.is_del and not pileupread.is_refskip:
# query position is None if is_del or is_refskip is set.
         print ('\tbase in read %s = %s' %
          (pileupread.alignment.query_name,
             pileupread.alignment.query_sequence[pileupread.query_position]))
samfile.close()

```

## API

### class pysam.AlignmentFile

```python
AlignmentFile(filepath_or_object, 
              mode=None, 
              template=None, 
              reference_names=None, 
              reference_lengths=None, 
              text=NULL, 
              header=None, 
              add_sq_text=False, 
              check_header=True, 
              check_sq=True,
              reference_filename=None, 
              filename=None, 
              index_filename=None, 
              filepath_index=None, 
              require_index=False,
              duplicate_filehandle=True, 
              ignore_truncation=False, 
              threads=1)
```

- mode：常见的有 `w` `rb` ...
- header: when writing, build header from a multilevel dictionary. The first level are the four types (‘HD’, ‘SQ’, . . . ). The second level are a list of lines, with each line being a list of tag-value pairs. The header is constructed first from all the defined fields, followed by user tags in alphabetical order. Alternatively, an AlignmentHeader object can be passed directly
- template: 根据指定的bam文件自动创建`header`, header格式如下

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210923102823681.png" alt="image-20210923102823681" style="zoom:80%;" />