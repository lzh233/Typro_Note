# fastqc结果查看

**html文件则为结果展示文件，**zip**文件储存了结果以及图片，下载查看**

![(Izh)  total  - rw-r  - rw-r  (Izh)  [root@VM-@-4-centos fastqcr]# 11  904  1 root root 571741 Jan 25 15  1 root root 348273 Jan 25 15  [root@VM-@-4-centos fastqcr]#  samplel 2 RI fastqc . html  samplel 2 RI fastqc .zip ](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image001.png)

2-1 fastqc结果（可以通过shell脚本批量比对，也可以把文件都跟在后面）

序列的基本信息

![Basic Statistics  Measure  Fil ename  File type  Encoding  Total Sequences  Sequences flagged as poor quality  Sequence length  NGC  Value  fastq. g z  Conventional base calls  Sanger / Illu_mina I. g  250000 ](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image002.png)

（2）每一个碱基的质量分布情况，一般分布在绿色区域（大于Q30）较好黄色（大于Q20），红色区域质量较差

![Per base sequence quality  Quality scores across all bases (Sanger / Illumina 1.9 encoding)  15-1 g  30-34  45-4g  60-64  Position in read (bp)  go-g4  105-109 120-124 135-139 150-151 ](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image003.png)

（3）每一条序列的质量情况，纵轴为reads数，横轴为质量，大多分布在Q30以上

![img](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image004.png)

（4）每条序列的GC含量、AT含量，好的结果为A=T以及G=C（如下图）

![I- 0 1 6 E I- E I Z I- 0 Z 1 6 01- 01 6 • 0 6  = SSO 」 02 u00  (dq) 92e 」 UOA!S0d  9 • 0 9  0 E  u.IöWO) Ed  61-gl 6 B 、 9 E Z I  0 01 ](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image005.png)

![1 23 4 56 7B g 10  Position in read (bp) ](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image006.png)

 

（4）每条序列的GC含量百分比，蓝色为经验值，红色为样品，符合经验值的分布（如下图）



（5）N碱基的数量

![Per base N content  100  15-1 g  30-34  N content across all bases  45-4g  60-64  Position in read (bp)  go-g4  105-109 120-124 135-139 150-15: ](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image008.png)

（6）序列长度的分布

![.판 ](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image009.png)

（6）重复序列的占比

![صصصصحححكد  صصصصصصم  صصصصصصصص  . صصصصحصصك  صىقتصسححه ](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image010.png)

（6）是否有大量重复的序列

![Overrepresented sequences  No overrepresented sequences ](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image011.png)

（7）测序接头情况，存在接头（下下图为不存在接头）

![I-gel  V11-011  (dq) pea' u! u0A!S0d  66-G6 6B-SB 6L-SL 69-gg ag-gg 6b-gv 6E-SE 6t-gz 61-GI  epv news anos  aouanbas asesodsue'l  Eldepv news  Eldepv E news  Eldepv les uamun  001  Eldepv % ](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image012.png)

 

![10 11 12 13 14 15 15 17 1B 19 20 21 22 23 24 25 25 27 2B 29  Position in read (bp) ](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/clip_image013.png)