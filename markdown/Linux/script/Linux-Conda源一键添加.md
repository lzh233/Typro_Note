# Conda源一键添加

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
#两个源任选其一
#名称：conda_tuna.sh、conda_sjtu.sh
##############
```

## 脚本（上海交大源）

```shell
#!/bin/bash
#上海交大源
url=https://mirrors.sjtug.sjtu.edu.cn/anaconda
#pkgs
conda config --add channels ${url}/pkgs/free/
conda config --add channels ${url}/pkgs/main/
conda config --add channels ${url}/pkgs/mro/
conda config --add channels ${url}/pkgs/msys2/
conda config --add channels ${url}/pkgs/pro/
conda config --add channels ${url}/pkgs/r/
#cloud
#conda config --add channels ${url}/cloud/atztogo/
conda config --add channels ${url}/cloud/bioconda/
conda config --add channels ${url}/cloud/conda-forge/
#conda config --add channels ${url}/cloud/matsci/
#conda config --add channels ${url}/cloud/menpo/
#conda config --add channels ${url}/cloud/pytorch/
#conda config --add channels ${url}/cloud/pytorch-test/
#conda config --add channels ${url}/cloud/soumith/
#conda config --add channels ${url}/cloud/viscid-hub/

```

## 脚本（清华源）

```shell
#!/bin/bash
# 添加常用下载频道
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

# 添加清华镜像加速下载
site=https://mirrors.tuna.tsinghua.edu.cn/anaconda
conda config --add channels ${site}/pkgs/free/
conda config --add channels ${site}/pkgs/main/
conda config --add channels ${site}/cloud/conda-forge/
conda config --add channels ${site}/pkgs/r/
conda config --add channels ${site}/cloud/bioconda/
conda config --add channels ${site}/cloud/msys2/
conda config --add channels ${site}/cloud/menpo/
conda config --add channels ${site}/cloud/pytorch/
```

