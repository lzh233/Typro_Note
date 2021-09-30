# cosmic脚本问题总结

#### 2021 - 9-18

问题1：

pyranges读取数据出错，不应该使用pyranges库去读取cosmic的数据库，应该去读取bed文件，因为应该是通过遍历数据库中的所有的位点信息，然后通过pyranges去判断数据库中的位点是否在包含在bed文件中，如果包含则保留这条数据库中的信息，所以应该是用pyranges去读取bed文件。

按照错误的方法去做就成了提取数据库中的位点信息，如`chr: 1	position_start: 114713961	position_end: 114716164`, 按照错误的办法来就成了 选出所有在上面范围的位点信息了，会造成数据的重复，同时无法得到除了位点以外的信息。原因如图，

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210922090056777.png" alt="image-20210922090056777" style="zoom:80%;" />

染色体7为例子，第3组的start是5518306 End是55181436 然后第四组的Start是55181312 End是55181443 两者之间有重叠 所以我找出来的数据里面有很多重复的

<img src="https://aironi.oss-cn-beijing.aliyuncs.com/typro_image/image-20210922090139256.png" alt="image-20210922090139256" style="zoom:80%;" />

```python
#错误的方法!
##错误的方法!
gr = pr.PyRanges(database_pos)
pos_effective = self.get_effective_position()
            
chr = pos_effective[0]
pos_start = pos_effective[1]
pos_end = pos_effective[2]
file_count = 0

for chr_s,start,end in zip(chr,pos_start,pos_end):
    print(f"chr: {chr_s}\tposition_start: {start}\tposition_end: {end}")
    file_count += 1
    
    chr_comsic = gr[str(chr_s),start:end]
    chr_comsic.to_csv(f"{self.result_dir}/{chr_s}_{file_count}_results.tsv",sep="\t")
```



问题2：

对Pandas了解不够多，如果要对数据框每行进行迭代，直接调用`iterrows()`方法就行，而且没有对数据框的header命名，导致各种索引混用，代码混乱



重写脚本：

```python
import pyranges as pr
import pandas as pd

class Cosmic:
    def __init__(self,db,bed):
        self.db = db
        self.bed = bed
    
    def get_cosmic_data(self):
        database = pd.read_table(self.db,
                                 header=None,
                                 sep = "\t",
                                 names = ['Chromosome','Start','End','Ref','Alt','Annotation'])
        
        bed = pd.read_table(self.bed,
                            usecols=[0,1,2],
                            names=['Chromosome', 'Start', 'End'])

        chr_use = set(bed.loc[:,"Chromosome"])
        database = database[database.loc[:,"Chromosome"].isin(chr_use)]

        gr = pr.PyRanges(bed)
        cosmic_result = []
        for _,row in database.iterrows():
            if gr[str(row['Chromosome']), row['Start']:row['Start']+1]:
                cosmic_result.append(row)

        df_cosmic_result = pd.DataFrame(cosmic_result)
        df_cosmic_result.to_csv("./cosmic_results.tsv",sep = "\t")
        

def main():
    a = Cosmic(db = "/SGRNJ/Database/script/database/annovar/humandb/hg38_cosmic70.txt",
              bed="./lung_COVER_Aligned.sortedByCoord.out.bed")
    a.get_cosmic_data()

main()
```



脚本位置

/SGRNJ03/randd/user/liuzihao/cosmic_find
