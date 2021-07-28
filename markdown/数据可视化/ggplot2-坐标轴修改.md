# ggplot2-图表的修改

## R版本与运行环境信息

```R
Author:liuzihao
Date:2021-4-14
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)
```

## 载入相关包

```r
library("tidyverse")
#设定图片主题
plot_theme <- theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 13,color = "black"))
```

## 坐标轴的修改

```r
ggplot(mpg) + 
  geom_point(aes(cty,hwy,color=class,size = cyl)) + 
  scale_color_npg()
```

<img src="G:\Desktop\s_note\data\picture\image-20210720232139595.png" alt="image-20210720232139595" style="zoom:80%;" />

### 坐标轴标题

- 使用**`labs()`**对坐标轴和图例的标题进行修改

```r
ggplot(mpg) + 
  geom_point(aes(cty,hwy,color=class,size = cyl)) + 
  scale_color_npg() + 
  labs(x = "city drivibg(mpg)",y = "highway driving(mpg)",color = "type of\ncar",size = NULL) + 
  plot_theme
```

<img src="G:\Desktop\s_note\data\picture\image-20210720232556700.png" alt="image-20210720232556700" style="zoom:80%;" />

### 坐标轴刻度与标签

- 使用**`scale_x_continuous`**和**`scale_x_continuous`**可以对坐标轴的刻度进行重排，并设定标签其中

  **`names`**：指定坐标轴的名称

  **`breaks`**：指定刻度的显示

  **`labels`**：将刻度换成指定的标签

```r
ggplot(mpg) + 
  geom_point(aes(cty,hwy,color=class,size = cyl)) + 
  scale_color_npg() + 
  scale_x_continuous(name = "city drivibg(mpg)",breaks = c(10,15,20,25,30,35),labels = str_c("T_",c(1:6))) + 
  scale_y_continuous(name = "highway driving(mpg)",breaks = c(15,20,30,40)) + 
  plot_theme
```

<img src="G:\Desktop\s_note\data\picture\image-20210720234014830.png" alt="image-20210720234014830" style="zoom:80%;" />

