# openssl升级

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
#需要自己下载openssl-1.1.1h
##############
```

## 脚本

```shell
#安装telnent防止升级失败为连不上服务器
rpm -qa telnet-server
rpm -qa xinetd
rpm -qa telnet

#均为安装的话使用yum进行安装
yum -y install xinetd telnet telnet-server
#修改telnet配置文件
touch /etc/xinetd.d/telnet
#vim /etc/xinetd.d/telnet
#写入如下内容
echo "service telnet
{
	disable = no
	flags = REUSE
	socket_type = stream
	wait = no
	user = root
	server = /usr/sbin/in.telnetd
	log_on_failure += USERID
}" >>  /etc/xinetd.d/telnet

#开启telnet
systemctl start telnet.socket
systemctl start xinetd

#查看23端口是否开启
netstat -tlun | grep 23

#查看防火墙状态
systemctl status firewalld
#运行中则暂时关闭
systemctl stop firewalld

#xshell测试：默认不能root登录，普通用户测试登录，可以登录则进行下一步
#安装依赖
yum -y install gcc gcc-c++ kernel-devel
#安装zlib
cd /root/ssh/zlib-1.2.11
./configure --prefix=/usr/local/zlib
make -j 4 
make install

#安装openssl-1.1.1
cd /root/ssh/openssl-1.1.1h
./config --prefix=/usr/local/ssl -d shared
make -j 4 
make install
echo '/usr/local/ssl/lib' >> /etc/ld.so.conf
ldconfig -v

#安装openssh-8.4p1
cd /root/ssh/openssh-8.4p1
mv /etc/ssh /etc/ssh.bak
./configure --prefix=/usr/local/openssh --sysconfdir=/etc/ssh --with-ssl-dir=/usr/local/ssl --with-zlib=/usr/local/zlib 
make -j 4 
make install

#修改配置文件
echo "PermitRootLogin yes" >> /etc/ssh/sshd_config

#备份原有ssh文件
mv /usr/sbin/sshd /usr/sbin/sshd.bak
mv /usr/bin/ssh /usr/bin/ssh.bak
mv /usr/bin/ssh-keygen /usr/bin/ssh-keygen.bak

#更新ssh文件
cp -a /usr/local/openssh/sbin/sshd /usr/sbin/sshd
cp -a /usr/local/openssh/bin/ssh /usr/bin/ssh
cp -a /usr/local/openssh/bin/ssh-keygen /usr/bin/ssh-keygen

#测试
ssh -V

#关闭原有ssh服务
systemctl stop sshd.service
#删除原有服务
rm -rf /lib/systemd/system/sshd.service 
#重读服务列表
systemctl daemon-reload

#sshd加入systemctl管理
cp ~/ssh/openssh-8.4p1/contrib/redhat/sshd.init /etc/init.d/sshd

#启动服务
/etc/init.d/sshd restart

#systemctl查看状态
systemctl status sshd

#重启服务测试systemctl
systemctl restart sshd

#使用其他用户登录没问题即升级成功

#重新开启防火墙
systemctl start firewalld
systemctl status firewalld

#停止telnet服务
systemctl stop telnet.socket
systemctl stop xinetd

#查看23端口：没有则成功关闭
netstat -tlun | grep 23

#卸载telnet和xinetd
rpm -e telnet-server
rpm -e xinetd
rpm -e telnet
rm -rf /etc/xinetd.d/telnet 
```