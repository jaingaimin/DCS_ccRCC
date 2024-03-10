

#save.image(file = "/t8a/allimage_data/Sdeath202302.rds")
#save.image(file = "/t8a/allimage_data/Sdeath202303.rds")
load("/t8a/allimage_data/Sdeath202302.rds")
load("/t8a/allimage_data/Sdeath202303.rds")


library(Seurat)
library(CellChat)
library(cellcall)
library(AUCell)
library(Scissor)


###拟使用的数据集 最好把大类已经标注好了
##直接做Sdeath score 高低评分的细胞 
#优先肿瘤细胞
#再者是做免疫和间质细胞

####singlecell addmodule score####


###AUC score#####


###cell chat####


####monocle2####


###immune diff and GSVA###