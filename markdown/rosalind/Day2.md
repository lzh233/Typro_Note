# Rosalind

## Day2

斐波那契数列（递归思想）

```python
def fab(n):
    if n == 1 or n == 2:
        return 1
    else:
        return fab(n-1) + fab(n-2)
print(fab(14))
##-------------------output---------------------##
55
```

