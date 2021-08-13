# 学习Git的猫

## Git常用命令

### 常用命令总览

| 命令名称                             | 作用           |
| ------------------------------------ | -------------- |
| git config --global user.name 用户名 | 设置用户签名   |
| git config --global user.email 邮箱  | 设置用户签名   |
| git init                             | 初始化本地库   |
| git status                           | 查看本地库状态 |
| git add 文件名                       | 添加到暂存区   |
| git commit -m "日志信息" 文件名      | 提交到本地库   |
| git reflog                           | 查看历史记录   |
| git reset --hard 版本号              | 版本穿梭       |

### 常用命令用法

- 设置用户名和邮箱

  ```shell
  #设定用户名为aironi
  git config --global user.name aironi
  #设定邮箱
  git config --global user.email lzh360370285@qq.com
  ```

  配置文件位于家目录下**(windows位置)**`C:\Users\lzh233`：

  ```shell
  #可以看到.gitconfig文件，即为配置文件，里面为设定好的用户名和邮箱
  [user]
  	name = aironi
  	email = lzh360370285@qq.com
  ```

- 初始化工作空间**`git init`**

  ```shell
  gin init [工作空间]
  #创建工作空间，gitdemo
  cd /g
  mkdir git2
  #初始化工作空间
  git init ./git2/
  #初始化后的工作空间中有了.git目录
  lzh233@DESKTOP-S5MC5J5 MINGW64 /g/git2 (master)
  $ ll -a
  total 21
  drwxr-xr-x 1 lzh233 197609  0 Apr 28 21:32 ./
  drwxr-xr-x 1 lzh233 197609  0 Apr 28 20:54 ../
  drwxr-xr-x 1 lzh233 197609  0 Apr 28 21:32 .git/
  ```

- 创建一个项目

  ```shell
  #创建hello.sh
  cat hello.sh
  ```

  使用**`git status`**查看当前工作空间的状态

  ```shell
  #On branch master:当前位于主分支
  #No commits yet：还没有文件提交到本地库
  #Untracked files:检测到一个未追踪的文件，即刚创建的项目，因为还未提交暂存区和本本地库所以显示Untracked files
  $ git status
  On branch master
  No commits yet
  Untracked files:
    (use "git add <file>..." to include in what will be committed)
          hello.sh
  
  nothing added to commit but untracked files present (use "git add" to track)
  ```

  使用**`git add`**将指定的代码提交到暂存区

  ```shell
  git add hello.sh
  #warning信息为自动转换了换行符
  #warninD: LF will be replaced by CRLF in hello.sh.
  #The file will have its original line endings in your working directory
  
  #再次查看当前状态则不再显示Untracked files
  On branch master
  No commits yet
  Changes to be committed:
    (use "git rm --cached <file>..." to unstage)
          new file:   hello.sh
  ```

  使用**`git rm --cached <file>`**删除指定的缓存区文件

  ```shell
  #当前缓存区有两个文件，删除hello2.sh
  $ git status
  On branch master
  No commits yet
  Changes to be committed:
    (use "git rm --cached <file>..." to unstage)
          new file:   hello.sh
          new file:   hello2.sh
  #删除暂存区文件hello2.sh
  git rm --cached hello2.sh
  #查看状态，hello2.sh变成了Untracked files，而hello.sh依旧在暂存区
  git status
  On branch master
  No commits yet
  Changes to be committed:
    (use "git rm --cached <file>..." to unstage)
          new file:   hello.sh
  Untracked files:
    (use "git add <file>..." to include in what will be committed)
          hello2.sh
  ```

  使用**`git commit -m "<discrption>"`**提交到本地库

  ```shell
  git commit -m "This is v1.0" hello.sh
  #再次查看git状态
  lzh233@DESKTOP-S5MC5J5 MINGW64 /g/git2 (master)
  $ git status
  On branch master
  nothing to commit, working tree clean
  ```
  
  使用**`git reflog`**和**`git log`**查看本地库信息
  
  ```shell
  git reflog
  ###6a6bada：文件的版本号；(HEAD -> master)：当前指针的位置
  6a6bada (HEAD -> master) HEAD@{0}: commit (initial): This is v1.0
  ```
  
  提交新版本时，**可以看到文件指针默认指向最新版本，当前工作空间的版本默认也是最新版本**
  
  ```shell
  edcd6f0 (HEAD -> master) HEAD@{0}: commit: This is v2.0
  6a6bada HEAD@{1}: commit (initial): This is v1.04
  ```
  
  使用**`git log`**查看版本详细信息
  
  ```shell
  git log
  ####
  lzh233@DESKTOP-S5MC5J5 MINGW64 /g/git2 (master)
  $ git log
  commit edcd6f02b459bd980429dc619e824358994f1b27 (HEAD -> master)
  Author: aironi <lzh360370285@qq.com>
  Date:   Thu Apr 29 09:11:57 2021 +0800
  	This is v2.0
  
  commit 6a6badadc68f2da2d5389c76bd3be06a180bf451
  Author: aironi <lzh360370285@qq.com>
  Date:   Thu Apr 29 09:07:19 2021 +0800
  	 This is v1.0
  ```
  
  使用**`git reset --hard <index>`**实现版本穿梭
  
  ```shell
  git reset --hard 6a6bada
  ###
  HEAD is now at 6a6bada This is v1.0
  ```
  
  **常用命令关系总结**

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210429111455175.png" alt="image-20210429111455175" style="zoom:80%;" />



## Git分支操作

### 分支的基本定义

在版本控制过程中，同时推进多个任务，为每个任务，我们就可以创建每个任务的单独分支。使用分支意味着程序员可以把自己的工作从开发主线上分离开来，开发自己分支的时候，不会影响主线分支的运行。同时并行推进多个功能开发，提高开发效率。各个分支在开发过程中，如果某一个分支开发失败，不会对其他分支有任何影响。失败的分支删除重新开始即可。    

### 分支原理示意图

![d](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210429094151354.png)

### 分支常用命令

| 命令名称            | 作用                                               |
| ------------------- | -------------------------------------------------- |
| git branch 分支名   | 创建分支                                           |
| git branch -v       | 查看分支                                           |
| git checkout 分支名 | 切换分支                                           |
| git merge 分支名    | 把指定的分支合并到当前分支上（冲突合并和正常合并） |

- 使用**`git branch`**创建一个新的分支，**`-v`**选项可以查看分支，分支默认是克隆了一个master

  ```shell
  #创建一个hot-fix分支
  git branch hot-fix
  lzh233@DESKTOP-S5MC5J5 MINGW64 /g/git2 (master)
  git branch -v
    hot-fix 6a6bada This is v1.0
  * master  6a6bada This is v1.0
  ```

- 使用**`git checkout <branch name>`**切换分支

  ```shell
  git checkout hot-fix
  Switched to branch 'hot-fix'
  #可也看出当前已经切换到了hot-fix分支
  lzh233@DESKTOP-S5MC5J5 MINGW64 /g/git2 (hot-fix)
  #对hot-fix内的文件进行修改，增加一行echo "hello git!",后提交缓存区后提交本地库，以测试merge
  ```

- 使用**`git merge <branch name>`**进行合并, **合并到当前分支上**

  ```shell
  #正常合并，仅仅在分支上对代码进行了修改，但是master分支的代码并未修改
  #首先切换到master分支
  git checkout master
  #合并hot-fix分支
  git merge hot-fix
  ###提示信息，文件版本从 6a6bada迭代到6a761f2
  Updating 6a6bada..6a761f2
  Fast-forward
   hello.sh | 1 +
   1 file changed, 1 insertion(+)
  #查看合并后文件
  $ cat hello.sh
  #!/bin/bash
  echo "hello world!"
  echo "hello git!"
  ```

  **冲突合并**

  ```shell
  #在创建新的分支后，对master分支的文件同时进行了修改，git无法对其正常合并
  git merge hot-fix
  Auto-merging hello.sh
  CONFLICT (content): Merge conflict in hello.sh
  Automatic merge failed; fix conflicts and then commit the result.
  ##可以看到分支名称变成了 (master|MERGING)
  #因此使用vim手动修改后提交
  #!/bin/bash
  echo "hello world!"
  <<<<<<< HEAD
  echo "hello git!"aaaa
  =======
  echo "hello git!"
  echo "abcd"
  >>>>>>> hot-fix
  #手动修改
  #!/bin/bash
  echo "hello world!"
  echo "hello git!"aaaa
  echo "abcd"
  ###提交修改，提交暂存区，提交本地库即可完成合并
  git add hello.sh
  #提交本库的时候不能指定文件名，否则报错，可以看到(master|MERGING)变为(master)
  git commit -m "unorm test! fin"
  ```

  

## GitHub的使用

### Git团队合作的机制

### 创建远程仓库

注册github账号，并创建远程仓库`lzh233`

![image-20210429203646384](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210429203646384.png)

创建成功

**远程库地址：`https://github.com/lzh233/lzh233.git`**

![image-20210429203752569](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210429203752569.png)

### GitHub常用操作

| 命令名称                           | 作用                                                      |
| ---------------------------------- | --------------------------------------------------------- |
| git remote -v                      | 查看当前所有远程地址别名                                  |
| git remote add 别名 远程地址       | 给远程库起别名                                            |
| git push 别名 分支                 | 推送本地分支上的内容到远程仓库                            |
| git clone 远程地址                 | 将远程仓库的内容克隆到本地                                |
| git pull 远程库地址别名 远程分支名 | 将远程仓库对于分支最新内容拉下来后与 当前本地分支直接合并 |

使用**`git remote add <别名> <远程库地址>`**

```shell
#创建别名
git remote add test https://github.com/lzh233/lzh233.git
```

使用**`git remote -v`**查看当前存在的别名

```shell
git remote -v
##输出结果，fetch和push表示拉取和推送均使用该别名
lzh233@DESKTOP-S5MC5J5 MINGW64 /g/lzh233 (master)
test    https://github.com/lzh233/lzh233.git (fetch)
test    https://github.com/lzh233/lzh233.git (push)
```

使用**`git push <别名> <分支>`**将本地分支推送到托管中心

```shell
#将master分支推送到托管中心
git push test master
###提示信息
Enumerating objects: 3, done.
Counting objects: 100% (3/3), done.
Delta compression using up to 4 threads
Compressing objects: 100% (2/2), done.
Writing objects: 100% (3/3), 729 bytes | 729.00 KiB/s, done.
Total 3 (delta 0), reused 0 (delta 0), pack-reused 0
To https://github.com/lzh233/lzh233.git
 * [new branch]      master -> master
```

![image-20210429210118132](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210429210118132.png)

使用**`git pull <远程库别名> <远程库分支>`**对远程库进行拉取

```shell
#首选修改一下远程库，最后一行加点料
#拉取
git pull test master
###提示信息
remote: Enumerating objects: 5, done.
remote: Counting objects: 100% (5/5), done.
remote: Compressing objects: 100% (2/2), done.
remote: Total 3 (delta 1), reused 0 (delta 0), pack-reused 0
Unpacking objects: 100% (3/3), 675 bytes | 112.00 KiB/s, done.
From https://github.com/lzh233/lzh233
 * branch            master     -> FETCH_HEAD
   5fc103a..c8fadf1  master     -> test/master
Updating 5fc103a..c8fadf1
Fast-forward
 test.sh | 5 +++++
 1 file changed, 5 insertions(+)
#查看一下本地版本信息
git reflog
###当前指针已经指向最新版本
c8fadf1 (HEAD -> master, test/master) HEAD@{0}: pull test master: Fast-forward
5fc103a HEAD@{1}: commit (initial): This test v1.0
```

使用**`git clone <别名>/<远程库地址>`**将远程仓库的内容克隆到本地，使用自带WSL进行`git clone`，

1、拉取代码。 2、初始化本地仓库。 3、创建别名  

```shell
git clone https://github.com/lzh233/lzh233.git

###提示信息
Cloning into 'lzh233'...
remote: Enumerating objects: 6, done.
remote: Counting objects: 100% (6/6), done.
remote: Compressing objects: 100% (4/4), done.
remote: Total 6 (delta 1), reused 2 (delta 0), pack-reused 0
Unpacking objects: 100% (6/6), done.
##查看目录（所有目录都在）
lzh@DESKTOP-S5MC5J5:~$ cd lzh233/
lzh@DESKTOP-S5MC5J5:~/lzh233$ la
.git  test.sh
##git status /git reflog /git log等均能使用（本地库只有一个版本）
lzh@DESKTOP-S5MC5J5:~/lzh233$ git reflog
c8fadf1 (HEAD -> master, origin/master, origin/HEAD) HEAD@{0}: clone: from https://github.com/lzh233/lzh233.git
#git log
lzh@DESKTOP-S5MC5J5:~/lzh233$ git log
commit c8fadf1133b890ca7272046083fc69b0a6b3e02c (HEAD -> master, origin/master, origin/HEAD)
Author: aironi <83407205+lzh233@users.noreply.github.com>
Date:   Thu Apr 29 21:55:41 2021 +0800

    Update test.sh

commit 5fc103ad395de5549d7523751499d56be602e306
Author: aironi <lzh360370285@qq.com>
Date:   Thu Apr 29 20:40:56 2021 +0800

    This test v1.0
```

### Git团队合作

#### 团队内

#### 跨团队

#### SSH密钥登录

**和`linux`设置密钥登录一致**（输入 `ssh-keygen  -t rsa`），将公钥(用户目录下的`./ssh/[filename].pub`)写入到`github`即可

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210430111256832.png" alt="image-20210430111256832" style="zoom:67%;" />

