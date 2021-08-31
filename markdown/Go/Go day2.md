# Go Day2

## 基本数据类型

Go 是一种强类型的静态编译语言，类型是高级语言的基础，有了类型，高级语言才能对不同类型抽象出不同的运算，编程者才能在更高的抽象层次上操纵数据，而不用关注具体存储和运算细节。  

具体类型如下

![image-20210825234354854](https://aironi.oss-cn-beijing.aliyuncs.com/Typro_imgimage-20210825234354854.png)



### 布尔类型 bool

布尔类型的关键字是`bool`，只有`True`和`False`, 二者Go 内置的两个预声明标识符  ，声明的时候如果不指定值，默认`false`

```go
var isOk bool
isOk := True
```

**布尔型数据和整型数据不能进行相互转换**  

```go
var a bool
a = 1 // error
```

比较表达式和逻辑表达式的结果都是布尔类型数据  

```go
var a bool = (x > y) & (x < z)
```

**`if` 和 `for `语句的条件部分一定是布尔类型的值或表达式**  

![image-20210825235006573](https://aironi.oss-cn-beijing.aliyuncs.com/Typro_imgimage-20210825235006573.png)

### 整型

Go 语言内置了12种整数类型，分别是 `byte`、`int` 、 `int8` 、` int16` 、 `init32` 、 `int64` 、 `uint` 、` uint8` 、`uint64 `、 `uintptr`转换。其中`byte`是 `uint8`的别名，不同类型的整型必须进行强制类型, 整型支持算术运算和位操作，算术表达式和位操作表达式的结果还是整型。  

```go
var a int32 = 5
var b int64 = 6
a = b // error
```

### 浮点型

语言内置两种浮点数类型，分别是` float32 `和` float64`, 默认是`float64`

**计算机很难进行浮点数的精确表示和存储 ， 因此两个浮点数之间不应该使用`＝`或 `!＝`进行比较操作，高精度科学计算应该使用 math 标准库**  

### 复数类型

Go 语言 内 置的复数类型有两种， 分别是 `complex64 `和 `complex l28 `，复数在计算机里面使用两个浮点数表示 ， 一个表示实部， 一个表示虚部。   

### 字符串

Go 语言将字符串作为一种原生的基本数据类型， 字符串的初始化可以使用字符串字面量

```go
var a = "hello world!"
```

( 1 ) 字符串属于**常量**，可以通过类似于数组的方式去索引访问，但是**不能修改某个字节的值**

```go
a[1] = "A" //error
```

（ 2 ）宇符串转换为切片 ［]byte(s）要慎用，尤其是当数据量较大时（每转换一次都需复制内容）。例如：  

```go
var a = "hello world!"
b := []byte(a)
```

( 3 ）字符串尾部不包含 NULL 字符  

( 4 ）字符串类型底层实现是一个二元的数据结构，一个是指针指 向字节数组的起点，另 一个是长度 。  

![image-20210826220607652](https://aironi.oss-cn-beijing.aliyuncs.com/Typro_imgimage-20210826220607652.png)

( 5 ）基于字符串创建的切片和原字符串指向相同的底层字符数组 ， 一样不能修改 ， 对字符串的切片操作返回的子串仍然是`string`，而非 `slice`。   

( 6 ）字符串 和切片 的转换 ： 字符串可以转换为**字节数组 ，也可 以转换为 Unicode 的 字数组。例如 ：**  

```go
package main

import "fmt"

var a string = "你是猪"

func main() {
	b := []rune(a)
	c := []byte(a)
	fmt.Println(len(a))
	fmt.Printf("Unicode的字数组: %v\n", b)
	fmt.Printf("字节数组: %v\n", c)
}
9
Unicode的字数组: [20320 26159 29482]
字节数组: [228 189 160 230 152 175 231 140 170]
```





