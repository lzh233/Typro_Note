# 离线创建conda环境

1、首先在可以联网的服务器或虚拟机中出创建好相应的环境

2、打包所创建的环境：

tar -tar -zcf [压缩文件名] [环境所在目录] &

3、打包创建该环境所需的依赖：

tar -zcf [压缩文件名] ~/miniconda2/pkgs &

4、将两压缩包上传服务器并解压

5（安装依赖包）

5-1 将 ~/miniconda2/pkgs备份：mv ~/miniconda2/pkgs ~/miniconda2/pkgsbak

5-2 将解压后的新pkgs复制到~/miniconda2/

6（创建环境）conda create -n lzh --clone [指定步骤2中解压后目录的位置] --offline

7、测试环境没问题后，删除原来的包，rm -rf ~/miniconda2/pkgsbak/