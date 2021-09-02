# glob模块

[glob](https://docs.python.org/zh-cn/3/library/glob.html#module-glob)

`glob` 模块可根据 Unix 终端所用规则找出所有匹配特定模式的路径名（支持`shell`的通配符），但会按不确定的顺序返回结果。波浪号扩展不会生效，但 `*`, `?` 以及表示为 `[]` 的字符范围将被正确地匹配。

## glob.glob()

返回匹配的目录、文件..，其中`recursive`表示是否递归查找，`True`代表递归查找，即查找所有目录，子目录以及文件，`pathname`支持`shell`风格的通配符，绝对路径和相对路径都支持

```python
usage: glob.glob(pathname, *, recursive=False)
#eg
print(glob.glob("./st*"))
print(glob.glob("/Personal/liuzihao/vscode_use/*.bam*",recursive=False))
#output
['./stat.txt']
['/Personal/liuzihao/vscode_use/test1_virus_Aligned.sortedByCoord.out.bam.bai', '/Personal/liuzihao/vscode_use/test1_virus_Aligned.sortedByCoord.out.bam']
```

## glob.iglob()

返回一个 `iterator`，它会产生与 [`glob()`](https://docs.python.org/zh-cn/3/library/glob.html#module-glob) 相同的结果，但不会实际地同时保存它们

```python
print(glob.iglob("./st*"))
print(glob.iglob("/Personal/liuzihao/vscode_use/*.bam*",recursive=False))
#output
<generator object _iglob at 0x7fbc35809410>
<generator object _iglob at 0x7fbc35809410>
#通过循环访问迭代器
for i in glob.iglob("/Personal/liuzihao/vscode_use/*.bam*",recursive=False):
    print(i)
#output
/Personal/liuzihao/vscode_use/test1_virus_Aligned.sortedByCoord.out.bam.bai
/Personal/liuzihao/vscode_use/test1_virus_Aligned.sortedByCoord.out.bam
```

## glob.escape

转义所有特殊字符 (`'?'`, `'*'` 和 `'['`)。 这适用于当你想要匹配可能带有特殊字符的任意字符串字面值的情况。

