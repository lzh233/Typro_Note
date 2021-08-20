# Rosalind

## Day4

孟德尔遗传定律计算

给定：一个群体中AA  Aa  aa 数量，假设两两之间均可交配，求后代是显性性状的概率，AA 2  Aa 2  aa2

输出：0.7833

可以计算1-隐性性状的概率

```python
def test(AA,Aa,aa):
    total = AA + Aa + aa
    num = 1-(aa/total*(aa-1)/(total - 1) + Aa/total * (Aa - 1)/(total -1)*0.25 + Aa/total * aa/(total -1)*0.5 + aa/total * Aa/(total - 1)*0.5)
    return round(num,4)

def main():
    print(test(AA = 2, Aa = 2, aa = 2))

if __name__=='__main__':
    main()
##-------------------output---------------------##
0.7833
```

