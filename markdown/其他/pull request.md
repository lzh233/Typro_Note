# 如何进行pull request

`Pull Request` 是一种通知机制。你修改了他人的代码，将你的修改通知原来的作者，希望他合并你的修改，这就是 `Pull Request`。"

[Pull Request 的命令行管理](https://www.ruanyifeng.com/blog/2017/07/pull_request.html)

## step1

`fork`别人的库，`github`对应库的主页的右上角

![image-20210901134823410](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210901134823410.png)

## step2

`fork`完成后，`clone`到本地，创建分支修改代码，修改后提交缓冲区提交本地库` push`到`fork`到自己那里的库

```shell
git clone git@github.com:lzh233/CeleScope-1.git
git checkout -b <new branch>
'....修改代码....'
git add <file>
git commit -m "..." 
git push git@github.com:lzh233/CeleScope-1.git <barnch>
```

## stpe3

点击`Pull Request` ,进入自己的`github`后进入自己`fork`的库

![image-20210901135312142](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210901135312142.png)

## step4

点击 `New pull Request`

![image-20210901135406188](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210901135406188.png)

## step5

选择你修改的分支，选择希望合并到哪条分支

![image-20210901135626876](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210901135626876.png)

## step5

选择完成后，下方会显示出做的所有修改，然后点击提交即可，(因为提交过所以按钮位置内容不一样)

![image-20210901135800631](https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210901135800631.png)