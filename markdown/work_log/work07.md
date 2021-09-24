# work-9.15

## 描述

输入数据：
--bam bam：/SGRNJ03/randd/RD20102301_DZH/Blood_panel/20210907/B0831_R1_V2FJ/06.target_metrics/B0831_R1_V2FJ_filtered_sorted.bam
--bed bed: /SGRNJ03/randd/test_rd/dxh/data/Blood/blood_Aligned.sortedByCoord.out.bed
--bp_min 50 (可选，默认为50bp)
--zl zl：/SGRNJ03/randd/RD20102301_DZH/Blood_panel/20210906/B0831_R1_V2ZL/06.analysis/B0831_R1_V2ZL_tsne_coord.tsv

需求：
  1，统计bam文件中在指定位置（根据bed文件）的reads和UMI的cover情况，输出统计表格
  position：bed文件中的位置对应的名称
  cell:在该位置cover的细胞数
  reads_number:在该位置的reads数（要求read 在指定位置bp至少大于bp_min）
  umi_number：在该位置的umi种类数
  mean_reads:对每个位置的每个cell的read数进行排序，reads的平均数
  medium_reads：对每个位置的每个cell的read数进行排序，reads的中位数
  mean_umi：对每个位置的每个cell的umi数进行排序，umi的平均数
  medium_umi:对每个位置的每个cell的umi数进行排序，umi的中位数

  \##example
  position  cell  reads_number  umi_number  mean_reads  medium_reads  mean_umi  medium_umi
  位置1  
  位置2
  位置3
  ...

2，根据position的cell barcode和zl文件，整合输出position ，cell barcode list，cluster 表格

\##example
position  cell barcode list  cluster
位置1    AAACATCGAAACATCGCCGAAGTA  1
位置2
位置3...