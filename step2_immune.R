

###immune cell algorithm1####
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

hmdat <- read.csv("/t8a/2022backup/xuzijun/vip39/huitu/FigureYa230immunelandscape/easy_input.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
type <- read.csv("/t8a/2022backup/xuzijun/vip39/huitu/FigureYa230immunelandscape/easy_input_type.csv", row.names = 1)
head(type)
table(rownames(sub)%in%rownames(hmdat))

comsam <- rownames(sub)
hmdat <- hmdat[comsam,]
risk <- sub[comsam,,drop = F]
dim(hmdat)

# 拆分不同算法结果，获得类的名字
#immMethod <- sapply(strsplit(colnames(hmdat),"_",fixed = T),"[",2) #用easy_input.csv列名里的算法信息
immMethod <- type$Methods # 用easy_input_type.csv的算法那一列

library(ComplexHeatmap) 
# 最新版ComplexHeatmap好像有一些bug，使用里面的pheatmap函数会报错
# 因此我们从脚本直接加载pheatmap函数
source("/t8a/2022backup/xuzijun/vip39/huitu/FigureYa230immunelandscape/pheatmap_translate.R") # 位于当前文件夹，出自ComplexHeatmap_2.7.9.tar.gz
ht_opt$message = FALSE
# 创建注释
annCol <- data.frame(
  Subtype = risk$Cluster,
  row.names = rownames(risk),
  stringsAsFactors = F)
annRow <- data.frame(row.names = colnames(hmdat),
                     Methods = factor(immMethod,levels = unique(immMethod)),
                     stringsAsFactors = F)

annColors <- list(
  "Subtype" = c("C1"="#BC3C29FF", "C2"="#0072B5FF","C3"="#E18727FF","C4"="#7876B1FF"))

# 数据标准化
indata <- t(hmdat)
indata <- indata[,colSums(indata) > 0] # 确保没有富集全为0的细胞
plotdata <- standarize.fun(indata,halfwidth = 2)

# 样本按risk score排序
samorder <- c(rownames(sub)[sub$Cluster=="C1"],rownames(sub)[sub$Cluster=="C2"],
              rownames(sub)[sub$Cluster=="C3"],rownames(sub)[sub$Cluster=="C4"])

# 拆分各算法的结果
plotdata1 <- plotdata[rownames(annRow[which(annRow$Methods == "TIMER"),,drop = F]),]
plotdata2 <- plotdata[rownames(annRow[which(annRow$Methods == "CIBERSORT"),,drop = F]),]
plotdata3 <- plotdata[rownames(annRow[which(annRow$Methods == "CIBERSORT-ABS"),,drop = F]),]
plotdata4 <- plotdata[rownames(annRow[which(annRow$Methods == "QUANTISEQ"),,drop = F]),]
plotdata5 <- plotdata[rownames(annRow[which(annRow$Methods == "MCPCOUNTER"),,drop = F]),]
plotdata6 <- plotdata[rownames(annRow[which(annRow$Methods == "XCELL"),,drop = F]),]
plotdata7 <- plotdata[rownames(annRow[which(annRow$Methods == "EPIC"),,drop = F]),]

# 分别画7次热图（参数基本同pheatmap里的pheatmap）
hm1 <- pheatmap(mat = as.matrix(plotdata1[,samorder]),
                border_color = NA,
                #color = bluered(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                annotation_col = annCol[samorder,,drop = F],
                annotation_colors = annColors,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = c(92,259,333,516),
                name = "TIMER") # 为子热图的图例命名

hm2 <- pheatmap(mat = as.matrix(plotdata2[,samorder]),
                border_color = NA,
                color = greenred(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = c(92,259,333,516),
                name = "CIBERSORT")

hm3 <- pheatmap(mat = as.matrix(plotdata3[,samorder]),
                border_color = NA,
                color = blueyellow(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = c(92,259,333,516),
                name = "CIBERSORT-ABS")

hm4 <- pheatmap(mat = as.matrix(plotdata4[,samorder]),
                border_color = NA,
                color = bluered(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = c(92,259,333,516),
                name = "QUANTISEQ")

hm5 <- pheatmap(mat = as.matrix(plotdata5[,samorder]),
                border_color = NA,
                color = inferno(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = c(92,259,333,516),
                name = "MCPCOUNTER")

hm6 <- pheatmap(mat = as.matrix(plotdata6[,samorder]),
                border_color = NA,
                color = viridis(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = c(92,259,333,516),
                name = "XCELL")

hm7 <- pheatmap(mat = as.matrix(plotdata7[,samorder]),
                border_color = NA,
                color = magma(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = c(92,259,333,516),
                name = "EPIC")

pdf("immune heatmap by ComplexHeatmap met.pdf", width = 10,height = 20) # 保存前请注意RGUI里不能有任何显示的图像，否则不会pdf打不开
draw(hm1 %v% hm2 %v% hm3 %v% hm4 %v% hm5 %v% hm6 %v% hm7, # 垂直连接子热图
     heatmap_legend_side = "bottom", # 热图注释放底部
     annotation_legend_side = "bottom") # 顶部注释放底部
invisible(dev.off())

save(sub,file="Sdeath_sub.rds")

##fig280
setwd("../")

#fig

####immune gene #####
type <- fread("/t8a/2022backup/xuzijun/vip39/data/immunegene5type.csv",header = T,data.table=F)
hmdat <- as.data.frame(t(KIRC_mRNA_fpkm[rownames(type),sub$Sample]))

head(type)
table(type$Gene%in%rownames(KIRC_mRNA_fpkm))
table(rownames(sub)%in%rownames(hmdat))
type$Gene <- rownames(type)

type <- type[type$Gene%in%rownames(KIRC_mRNA_fpkm),]
type <- type[!duplicated(type$Gene),]
rownames(type) <- type$Gene
type <- type[-1]
head(type)
head(hmdat)

comsam <- rownames(sub)
hmdat <- hmdat[comsam,]
risk <- sub[comsam,,drop = F]
dim(hmdat)
head(hmdat[,1:5])

# 拆分不同算法结果，获得类的名字
#immMethod <- sapply(strsplit(colnames(hmdat),"_",fixed = T),"[",2) #用easy_input.csv列名里的算法信息
immMethod <- type$Methods # 用easy_input_type.csv的算法那一列

library(ComplexHeatmap) 
# 最新版ComplexHeatmap好像有一些bug，使用里面的pheatmap函数会报错
# 因此我们从脚本直接加载pheatmap函数
source("/t8a/2022backup/xuzijun/vip39/huitu/FigureYa230immunelandscape/pheatmap_translate.R") # 位于当前文件夹，出自ComplexHeatmap_2.7.9.tar.gz

ht_opt$message = FALSE

# 创建注释
annCol <- data.frame(
  Subtype = risk$Cluster,
  row.names = rownames(risk),
  stringsAsFactors = F)
annRow <- data.frame(row.names = colnames(hmdat),
                     Methods = factor(immMethod,levels = unique(immMethod)),
                     stringsAsFactors = F)

annColors <- list(
  "Subtype" = c("C1"="#BC3C29FF", "C2"="#0072B5FF","C3"="#E18727FF","C4"="#7876B1FF"))

# 数据标准化
indata <- t(hmdat)
indata <- indata[,colSums(indata) > 0] # 确保没有富集全为0的细胞
plotdata <- standarize.fun(indata,halfwidth = 2)

# 样本按risk score排序
samorder <- c(rownames(sub)[sub$Cluster=="C1"],rownames(sub)[sub$Cluster=="C2"],
              rownames(sub)[sub$Cluster=="C3"],rownames(sub)[sub$Cluster=="C4"])

# 拆分各算法的结果
###进行差异分析 


## 用前面的自定义函数计算组间统计差异
setwd("/home/aim/Desktop/TCGA_work/vip33new/Sdeath_project/result/heatmap_gene")
head(Subtype)
Subtype <- data.frame(Subtype=sub$Cluster)
rownames(Subtype) <- sub$Sample
plotdata <- plotdata[,rownames(Subtype)]
comprTab <- cross_subtype_compr(expr = plotdata, # 或log2(mygene_data + 1)，如果用参数检验，请注意对数转化；若非参均可
                                subt = Subtype,
                                #two_sam_compr_method = "wilcox", # 两组"t.test", "wilcox"
                                multi_sam_compr_method = "kruskal", # 多组"anova", "kruskal"
                                res.path = ".")
# 用全部基因来画
n.show_top_gene <- nrow(plotdata)
# 或者取top 20个基因来画
#n.show_top_gene <- 20 
# 按分组排序
subt.order <- Subtype[order(Subtype$Subtype),,drop = F]
indata <- mygene_data[comprTab$gene[1:n.show_top_gene],rownames(subt.order)]

# 开始画图
# 数据标准化和边界设置
# plotdata <- t(scale(t(indata)))
# plotdata[plotdata > 2] <- 2
# plotdata[plotdata < -2] <- -2

# 调整行名
blank <- "    " # 行名和p值之间的间隔
p.value <- comprTab$adjusted.p.value[1:n.show_top_gene]
table(comprTab$gene==rownames(plotdata))
head(comprTab)
head(type)


sig.label <- ifelse(p.value < 0.001,"****",
                    ifelse(p.value < 0.005,"***",
                           ifelse(p.value < 0.01,"**",
                                  ifelse(p.value < 0.05,"*",""))))
rownames(plotdata) <- paste0(rownames(plotdata)," ",sig.label)

plotdata1 <- plotdata[rownames(annRow[which(annRow$Methods == "chemokine"),,drop = F]),]
plotdata2 <- plotdata[rownames(annRow[which(annRow$Methods == "chemokine_receptor"),,drop = F]),]
plotdata3 <- plotdata[rownames(annRow[which(annRow$Methods == "MHC"),,drop = F]),]
plotdata4 <- plotdata[rownames(annRow[which(annRow$Methods == "Immunoinhibitor"),,drop = F]),]
plotdata5 <- plotdata[rownames(annRow[which(annRow$Methods == "Immunostimulator"),,drop = F]),]


plotdata1 <- plotdata[rownames(annRow[which(annRow$Methods == "chemokine"),,drop = F]),]
rownames(plotdata1) <- paste0(rownames(plotdata1)," ",sig.label[which(annRow$Methods == "chemokine")])
plotdata2 <- plotdata[rownames(annRow[which(annRow$Methods == "chemokine_receptor"),,drop = F]),]
rownames(plotdata2) <- paste0(rownames(plotdata2)," ",sig.label[which(annRow$Methods == "chemokine_receptor")])
plotdata3 <- plotdata[rownames(annRow[which(annRow$Methods == "MHC"),,drop = F]),]
rownames(plotdata3) <- paste0(rownames(plotdata3)," ",sig.label[which(annRow$Methods == "MHC")])
plotdata4 <- plotdata[rownames(annRow[which(annRow$Methods == "Immunoinhibitor"),,drop = F]),]
rownames(plotdata4) <- paste0(rownames(plotdata4)," ",sig.label[which(annRow$Methods == "Immunoinhibitor")])
plotdata5 <- plotdata[rownames(annRow[which(annRow$Methods == "Immunostimulator"),,drop = F]),]
rownames(plotdata5) <- paste0(rownames(plotdata5)," ",sig.label[which(annRow$Methods == "Immunostimulator")])
# 分别画7次热图（参数基本同pheatmap里的pheatmap）
library(RColorBrewer)
library(circlize)
library(ggsci)
library(gplots)
library(viridis)
library(oompaBase)
hm1 <- pheatmap(mat = as.matrix(plotdata1[,samorder]),
                border_color = NA,
                #color = bluered(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                annotation_col = annCol[samorder,,drop = F],
                annotation_colors = annColors,
                cellwidth = 0.8,
                cellheight = 12,
                gaps_col = c(92,259,333,516),
                name = "chemokine") # 为子热图的图例命名

hm2 <- pheatmap(mat = as.matrix(plotdata2[,samorder]),
                border_color = NA,
                color = greenred(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 12,
                gaps_col = c(92,259,333,516),
                name = "chemokine_receptor")

hm3 <- pheatmap(mat = as.matrix(plotdata3[,samorder]),
                border_color = NA,
                color = blueyellow(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 12,
                gaps_col =c(92,259,333,516),
                name = "MHC")

hm4 <- pheatmap(mat = as.matrix(plotdata4[,samorder]),
                border_color = NA,
                color = bluered(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 12,
                gaps_col = c(92,259,333,516),
                name = "Immunoinhibitor")

hm5 <- pheatmap(mat = as.matrix(plotdata5[,samorder]),
                border_color = NA,
                color = inferno(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 12,
                gaps_col = c(92,259,333,516),
                name = "Immunostimulator")


pdf("immune heatmap by ComplexHeatmap gene RNAmodify.pdf", width = 10,height = 30) # 保存前请注意RGUI里不能有任何显示的图像，否则不会pdf打不开
draw(hm1 %v% hm2 %v% hm3 %v% hm4 %v% hm5, # 垂直连接子热图
     heatmap_legend_side = "bottom", # 热图注释放底部
     annotation_legend_side = "bottom") # 顶部注释放底部
invisible(dev.off())


###pericyte to endothail vessel normalization score######
rownames(plotdata)
sub$PtoEscore <- scale(as.numeric(plotdata["Cancer associated fibroblast_XCELL",sub$Sample]/plotdata["Endothelial cell_XCELL",sub$Sample]))

compaired <- list(c("C1", "C2"), 
                  c("C1", "C3"), 
                  c("C1", "C4"),
                  c("C2", "C3"),
                  c("C2", "C4"),
                  c("C3", "C4"))   #设置比较组别
library(ggsci)
#colnames(df.merge)
colnames(sub)
col_1<-pal_nejm()(8)
sub$PtoEscore[sub$PtoEscore > 1] <- 1
sub$PtoEscore[sub$PtoEscore < -1] <- -1
ggboxplot(sub, x="Cluster", y="PtoEscore",color = "Cluster",
          palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"), add = "jitter")+
  ### 添加多组比较的统计学结果
  stat_compare_means(label.y = max(sub$PtoEscore))+
  ### 添加每个两两比较的显著性标记位置信息
  stat_compare_means(comparisons=compaired, method = "wilcox.test",
                     label ="p.value")
ggsave(filename = "vessel_normalization_score.pdf",height = 5,width = 5)
###more signal from IOBR ####
library(IOBR)
names(sig_group)
signature_collection
length(signature_collection)
names(signature_collection)


###免疫 以及immune score####
#各种评分
##TIP
##CIBERSORT
##ssGSEA
##IOBR
##CYT
#fig163


###两个subtype cor heatmap####

#fig97 gene and immune cor
#fig277 immune gene subtype
library(ChAMPdata) # 用于提供甲基化注释文件
library(genefu) # 用于获取乳腺癌PAM50分型
#data("pam50.robust")
data("probe.features")
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
# 设置热图颜色
heatmap.BlWtRd <- c("#1F66AC", "#75AFD3", "grey90", "#FAB99B", "#B2192B")
# 设置感兴趣基因集
immunomodulator <- read.table("/t8a/2022backup/ssy088_202210/ssy088/tutulaile/FigureYa277Immunomodulator/immunomodulator.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# 数据处理 #
# 该部分是为产生最终用于绘图的文件
## 表达谱
expr <- read.table("TCGA.BRCA.sampleMap-HiSeqV2.gz",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


expr <- tpms
is.element(rownames(immunomodulator),rownames(tpms)) # 所有基因都在表达谱内

# pam50pred <- molecular.subtyping(sbt.model = "pam50",
#                                  data = t(expr),
#                                  annot = data.frame(Gene.Symbol = rownames(expr)),
#                                  do.mapping = FALSE)
# subtype <- as.character(pam50pred$subtype)
# 
# sinfo <- data.frame(row.names = colnames(expr),
#                     pam50 = subtype,
#                     subtype = sapply(subtype, switch,
#                                      "Basal"  = "C1",
#                                      "Her2"   = "C2",
#                                      "LumB"   = "C3",
#                                      "LumA"   = "C4",
#                                      "Normal" = "C5"))

sinfo <- data.frame(row.names=sub$Sample,
                    subtype =sub$Cluster)
head(sinfo)
rownames(sinfo)<- paste0(rownames(sinfo),"A")
table(sinfo$subtype)

expr <- tpms[rownames(immunomodulator),]
colnames(expr) <- paste0(colnames(expr),"A")
## 甲基化谱
meth <- fread("/t8a/vip39database/database/UCSC_TCGA/TCGA_methylation450/TCGA-KIRC.methylation450.tsv.gz", sep = "\t",check.names = F,stringsAsFactors = F,header = T,data.table = F)
meth[1:4,1:4]

rownames(meth) <- meth$`Composite Element REF`; meth <- meth[,-1]
meth <- as.data.frame(na.omit(meth))


probeOfInterest <- probe.features[which(probe.features$gene %in% rownames(immunomodulator)),]
probeOfInterest <- probeOfInterest[intersect(rownames(probeOfInterest), rownames(meth)),]
is.element(rownames(immunomodulator), probeOfInterest$gene) # 有些基因没有对应的甲基化探针

meth <- meth[rownames(probeOfInterest),]
meth$gene <- probeOfInterest$gene
meth <- as.data.frame(apply(meth[,setdiff(colnames(meth), "gene")], 2, function(x) tapply(x, INDEX=factor(meth$gene), FUN=median, na.rm=TRUE)))


## 拷贝数变异 (-2,-1,0,1,2: 2 copy del, 1 copy del, no change, amplification, high-amplification)
cna <- read.table("/home/aim/Desktop/urology/Urinary_data/kidney/KIRC_DATA/TCGA.KIRC.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
cna$gene <- sapply(strsplit(rownames(cna),"|",fixed = T),"[",1)
cna <- cna[!duplicated(cna$gene),]; cna <- cna[,setdiff(colnames(cna),"gene")]
is.element(rownames(immunomodulator),rownames(cna)) # 有些基因没有对应的拷贝数结果
cna <- cna[intersect(rownames(cna),rownames(immunomodulator)),]
cna[cna > 1] <- 1 # 统一扩增
cna[cna < -1] <- -1 # 统一缺失

## 提取共同样本
comsam <- intersect(colnames(expr), colnames(meth))
colnames(cna) <- paste0(colnames(cna),"A")
comsam <- intersect(comsam, colnames(cna))
comsam <- intersect(comsam, rownames(sinfo))
head(sinfo)
sinfo <- sinfo[comsam,,drop = F]
expr <- expr[,comsam]
meth <- meth[,comsam]
cna <- cna[,comsam]
setwd("./figya277/")
head(sinfo)

write.table(sinfo[,"subtype",drop = F], file = "easy_input_subtype.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(expr, file = "easy_input_expr.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(meth, file = "easy_input_meth.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(cna, file = "easy_input_cna.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 读取数据 #
sinfo <- read.table(file = "easy_input_subtype.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
expr <- read.table(file = "easy_input_expr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
meth <- read.table(file = "easy_input_meth.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
cna <- read.table(file = "easy_input_cna.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

sinfo <- sinfo[,"subtype",drop = F]
head(sinfo)

# 亚型数目（注意，亚型必须以C1，C2，C3...等命名）
(n.subt <- length(unique(sinfo$subtype))) # 获取亚型数目

subt <- unique(sinfo$subtype) # 获取亚型名

# 初始化绘图矩阵
expMat <- as.data.frame(t(expr[rownames(immunomodulator),]))
expMat$subtype <- sinfo[rownames(expMat), "subtype"]
expMat <- as.data.frame(t(apply(expMat[,setdiff(colnames(expMat), "subtype")], 2, 
                                function(x) 
                                  tapply(x, 
                                         INDEX = factor(expMat$subtype), 
                                         FUN = median, 
                                         na.rm = TRUE)))) # 对同一亚型内的样本取中位数
corExpMeth <- ampFreq <- delFreq <- 
  as.data.frame(matrix(NA,
                       nrow = nrow(immunomodulator),
                       ncol = n.subt, 
                       dimnames = list(rownames(immunomodulator), 
                                       unique(sinfo$subtype))))

## 表达谱与甲基化的相关性
for (i in rownames(immunomodulator)) {
  if(!is.element(i, rownames(expr)) | !is.element(i, rownames(meth))) { # 如果存在任意一方有缺失的基因
    corExpMeth[i,] <- NA # 则保持矩阵为NA
  } else { # 否则取出亚型样本，做表达和甲基化的相关性
    for (j in subt) {
      sam <- rownames(sinfo[which(sinfo$subtype == j),,drop = F])
      expr.subset <- as.numeric(expr[i, sam])
      meth.subset <- as.numeric(meth[i, sam])
      ct <- cor.test(expr.subset, meth.subset, method = "spearman") # 这里采用speaman相关性
      corExpMeth[i, j] <- ct$estimate
    }
  }
}

## 扩增/缺失频率
for (i in rownames(immunomodulator)) {
  if(!is.element(i, rownames(cna))) { # 同理，如果存在拷贝数中缺失某基因，则保持NA
    ampFreq[i,] <- NA 
    delFreq[i,] <- NA
  } else { # 否则
    # 计算i在总样本中的频率
    ampFreqInAll <- sum(as.numeric(cna[i,]) == 1)/ncol(cna) # 总样本中扩增的数目除以总样本数
    delFreqInAll <- sum(as.numeric(cna[i,]) == -1)/ncol(cna) # 总样本中缺失的数目除以总样本数
    for (j in subt) {
      # 计算i在亚型j中的频率
      sam <- rownames(sinfo[which(sinfo$subtype == j),,drop = F])
      cna.subset <- cna[, sam]
      ampFreqInSubt <- sum(as.numeric(cna.subset[i,]) == 1)/length(sam) # 该亚型中扩增的数目除以该亚型样本数
      delFreqInSubt <- sum(as.numeric(cna.subset[i,]) == -1)/length(sam) # 该亚型中缺失的数目除以该亚型样本数
      
      ampFreqInDiff <- ampFreqInSubt - ampFreqInAll # 根据原本，用亚型特异性扩增比例减去总扩增比例
      delFreqInDiff <- delFreqInSubt - delFreqInAll # 同理
      
      ampFreq[i, j] <- ampFreqInDiff
      delFreq[i, j] <- delFreqInDiff
    }
  }
}
write.table(expMat,"expMat.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(corExpMeth,"corExpMeth.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(ampFreq,"ampFreq.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(delFreq,"delFreq.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 创建列注释
annCol <- data.frame(subtype = subt,
                     row.names = subt)
annCol <- annCol[order(annCol$subtype),,drop = F] # 按照亚型排序
annColors <- list()
annColors[["subtype"]] <- c("C1" = "red",
                            "C2" = "yellow",
                            "C3" = "green",
                            "C4" = "cyan"
                            #"C5" = "blue"
                            )
top_anno <- HeatmapAnnotation(df                   = annCol,
                              col                  = annColors,
                              gp                   = gpar(col = "grey80"), # 每个单元格边框为灰色
                              simple_anno_size     = unit(3.5, "mm"), # 注释高3.5毫米
                              show_legend          = F, # 不显示亚型的图例，因为一目了然
                              show_annotation_name = F, # 不显示该注释的名称
                              border               = FALSE) # 不显示注释的外边框

# 创建行注释
annRow <- immunomodulator
annRow[which(annRow$Category == "Co-stimulator"),"Category"] <- "Co-stm" # 这里字符少一些，不会挤在一起，可以后期AI
annRow[which(annRow$Category == "Co-inhibitor"),"Category"] <- "Co-ihb"
annRow[which(annRow$Category == "Cell adhesion"),"Category"] <- "Cell\nadhesion" # 这里换行，不会挤在一起，可以后期AI
annRow[which(annRow$Category == "Antigen presentation"),"Category"] <- "Antigen\npresentation"
annRow$Category <- factor(annRow$Category, levels = c("Co-stm","Co-ihb","Ligand","Receptor","Cell\nadhesion","Antigen\npresentation","Other")) # 由于行需要按照类分割，所以需要定义因子顺序，否则按照字母表
annRow$ICI <- factor(annRow$ICI, levels = c("Inhibitory","N/A","Stimulatory"))
annRowColors <- list("ICI" = c("Inhibitory" = "black","N/A" = "#888888","Stimulatory" = "#E59E02"))
left_anno <- HeatmapAnnotation(df                   = annRow[,"ICI",drop = F],
                               which                = "row", # 这里是行注释（默认为列）
                               gp                   = gpar(col = "grey80"), # 每个单元格边框为灰色
                               col                  = annRowColors,
                               simple_anno_size     = unit(3.5, "mm"), # 注释宽3.5毫米
                               show_annotation_name = F,
                               border               = F)

## 绘制表达谱热图（参数下同）
col_expr <- colorRamp2(seq(min(na.omit(expMat)), max(na.omit(expMat)), length = 5), heatmap.BlWtRd) # 创建热图颜色（将热图输入矩阵的最大最小值取5个点，分配颜色红蓝色板；注意矩阵中可能存在的NA值）
hm.expr <- Heatmap(matrix             = as.matrix(expMat),
                   col                = col_expr,
                   border             = NA, # 无热图外边框
                   rect_gp = gpar(col = "grey80"), # 热图单元格边框为灰色
                   cluster_rows       = F, # 行不聚类
                   cluster_columns    = F, # 列不聚类
                   show_row_names     = T, # 显示行名
                   row_names_side     = "left", # 行名显示在左侧
                   row_names_gp       = gpar(fontsize = 10), # 行名字号为10
                   show_column_names  = F, # 不显示列名（可后期在颜色内AI使得亚型一目了然）
                   column_names_side  = "top", # 列名显示在顶部
                   row_split          = annRow$Category, # 行按照Category进行分割（因子顺序）
                   top_annotation     = top_anno, # 热图顶部注释
                   left_annotation    = left_anno, # 热图左侧注释
                   name               = "mRNA\nExpression", # 热图颜色图例的名称
                   width              = ncol(expMat) * unit(4, "mm"), # 热图单元格宽度（稍大于高度，因为所有注释都放在底部，平衡图形纵横比）
                   height             = nrow(expMat) * unit(3.5, "mm")) # 热图单元格高度

col_corExprMeth <- colorRamp2(seq(min(na.omit(corExpMeth)), max(na.omit(corExpMeth)), length = 5), heatmap.BlWtRd)
hm.corExprMeth <- Heatmap(matrix             = as.matrix(corExpMeth),
                          col                = col_corExprMeth,
                          border             = NA,
                          rect_gp = gpar(col = "grey80"),
                          cluster_rows       = F,
                          cluster_columns    = F,
                          show_row_names     = F,
                          row_names_side     = "left",
                          row_names_gp       = gpar(fontsize = 10),
                          show_column_names  = F,
                          column_names_side  = "top",
                          row_split          = annRow$Category,
                          row_title          = NULL,
                          top_annotation     = top_anno,
                          name               = "Expression\nvs. Methylation",
                          width              = ncol(expMat) * unit(4, "mm"),
                          height             = nrow(expMat) * unit(3.5, "mm"))

col_ampFreq <- colorRamp2(seq(min(na.omit(ampFreq)), max(na.omit(ampFreq)), length = 5), heatmap.BlWtRd)
hm.ampFreq <- Heatmap(matrix             = as.matrix(ampFreq),
                      col                = col_ampFreq,
                      border             = NA,
                      rect_gp = gpar(col = "grey80"),
                      cluster_rows       = F,
                      cluster_columns    = F,
                      show_row_names     = F,
                      row_names_side     = "left",
                      row_names_gp       = gpar(fontsize = 10),
                      show_column_names  = F,
                      column_names_side  = "top",
                      row_split          = annRow$Category,
                      row_title          = NULL,
                      top_annotation     = top_anno,
                      name               = "Amplification\nFrequency",
                      width              = ncol(expMat) * unit(4, "mm"),
                      height             = nrow(expMat) * unit(3.5, "mm"))

col_delFreq <- colorRamp2(seq(min(na.omit(delFreq)), max(na.omit(delFreq)), length = 5), heatmap.BlWtRd)
hm.delFreq <- Heatmap(matrix             = as.matrix(delFreq),
                      col                = col_delFreq,
                      border             = NA,
                      rect_gp = gpar(col = "grey70"),
                      cluster_rows       = F,
                      cluster_columns    = F,
                      show_row_names     = F,
                      row_names_side     = "left",
                      row_names_gp       = gpar(fontsize = 10),
                      show_column_names  = F,
                      column_names_side  = "top",
                      row_split          = annRow$Category,
                      row_title          = NULL,
                      top_annotation     = top_anno,
                      name               = "Deletion\nFrequency",
                      width              = ncol(expMat) * unit(4, "mm"),
                      height             = nrow(expMat) * unit(3.5, "mm"))

pdf(file = "complexheatmap of immunomodulator.pdf", width = 8,height = 12)
draw(hm.expr + hm.corExprMeth + hm.ampFreq + hm.delFreq, # 水平衔接各个子热图
     heatmap_legend_side = "bottom") # 热图颜色图例显示在下方
invisible(dev.off())


#fig248 突变高级热图
#fig136  fgsea
#fig76 gene corrlation
#fig282 CMAP_XSum
#fig238 corRiskMut
#fig161 stemness
###fig152 immune tow cor####
#tipcor
#C3 单独作为一组
setwd("./fig152/")

stepscore[1:4,1:4]

tcga_gsva1 <- stepscore[,sub$Sample[sub$Cluster=="C4"]]

tcga_gsva1[1:3,1:3] 
tcga_gsva1 <- as.data.frame(t(tcga_gsva1))##列为细胞或者通路分数
tcga_gsva2 <- stepscore[,sub$Sample[sub$Cluster!="C4"]]
tcga_gsva2[1:3,1:3] 
tcga_gsva2 <- as.data.frame(t(tcga_gsva2))

#这里将计算每列之间的相关性
#如果要计算行之间的相关性就运行下面这行转置（把行变成列，列变成行）
#tcga_gsva <- t(tcga_gsva)
# 计算相关系数
cor_r <- cor(tcga_gsva1)
cor_r[1:3,1:3]
write.csv(cor_r, "easy_input_R1.csv", quote = F)
cor_r <- cor(tcga_gsva2)
cor_r[1:3,1:3]
write.csv(cor_r, "easy_input_R2.csv", quote = F)

# 计算p value
tcga_gsva <- tcga_gsva1
cor_p <- matrix(0, nrow = ncol(tcga_gsva), ncol = ncol(tcga_gsva))
rownames(cor_p) <- colnames(tcga_gsva)
colnames(cor_p) <- colnames(tcga_gsva)
for (i in 1:ncol(tcga_gsva)){
  for (j in 1:ncol(tcga_gsva)){
    p <- cor.test(tcga_gsva[,i],tcga_gsva[,j])
    cor_p[i,j] <- p$p.value
  }
}
write.csv(cor_p, "easy_input_P1.csv", quote = F)

tcga_gsva <- tcga_gsva2
cor_p <- matrix(0, nrow = ncol(tcga_gsva), ncol = ncol(tcga_gsva))
rownames(cor_p) <- colnames(tcga_gsva)
colnames(cor_p) <- colnames(tcga_gsva)
for (i in 1:ncol(tcga_gsva)){
  for (j in 1:ncol(tcga_gsva)){
    p <- cor.test(tcga_gsva[,i],tcga_gsva[,j])
    cor_p[i,j] <- p$p.value
  }
}
write.csv(cor_p, "easy_input_P2.csv", quote = F)

#KICH是C4  KIRC是其他的
KICHR <- read.csv("easy_input_R1.csv", row.names = 1)
KICHR[1:3,1:3]
KICHP <- read.csv("easy_input_P1.csv", row.names = 1)
KICHP[1:3,1:3]
KIRCR <- read.csv("easy_input_R1.csv", row.names = 1)
KIRCP <- read.csv("easy_input_P1.csv", row.names = 1)

## 合并相关系数的数据
datR <- KICHR
for(i in 1:nrow(datR)){
  datR[i,1:i] <- KIRCR[i,1:i]
}
datR[1:3,1:3]

## 合并P值的数据
datP <- KICHP
for (i in 1:nrow(datP)) {
  datP[i,1:i] <- KIRCP[i,1:i]
}
datP[1:3,1:3]

# P>= 0.05时，把相关系数设为NA
datR[datP > 0.05] <- NA

## 定义左右两个不同相关系数图的颜色
# 定义右上部分图形的颜色
colCorRight <-  circlize::colorRamp2(c(-1, 0, 1), c("green", "white", "#ef3b2c"))
# 定义左上部分图形的颜色
colCorLeft <- circlize::colorRamp2(c(-1, 0, 1), c("yellow", "white", "#762a83"))

## 绘制基本图形
p1 <- Heatmap(datR, rect_gp = gpar(type = "none"), 
              show_heatmap_legend = F,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.rect(x = x, y = y, width = width, height = height,
                          gp = gpar(col = "grey", fill = NA))
                if(i == j) {
                  grid.circle(x = x, y = y, r = 0.5 * min(unit.c(width, height)), gp = gpar(fill = "grey", col = NA))
                }else if(i > j) {
                  grid.circle(x = x, y = y, r = abs(datR[i, j])/2 * min(unit.c(width, height)), 
                              gp = gpar(fill = colCorLeft(datR[i, j]), col = NA))
                } else {
                  grid.circle(x = x, y = y, r = abs(datR[i, j])/2 * min(unit.c(width, height)), 
                              gp = gpar(fill = colCorRight(datR[i, j]), col = NA))
                }
              },
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = T, show_column_names = T, 
              row_names_side = "right", 
              row_names_rot = 45,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8)
)

p1

## 绘制两个不同不同相关的图例
lgdRight <- Legend(col_fun = colCorRight, title = "KICH", 
                   direction = "horizontal")
lgdLeft <- Legend(col_fun = colCorLeft, title = "KIRC", 
                  direction = "horizontal")
pd = list(lgdRight, lgdLeft)

## 最后出图
pdf("DouleCorPlot.pdf", width = 5, height = 5.5)
draw(p1, annotation_legend_list = pd,
     annotation_legend_side = "top")
dev.off()

#fig158 MutationPattern
#fig227 boxdensity drug AUC
#fig292 HCCsubtype 和其他分型比较
#fig168 legoplot 突变位点


####immune cohort####
#进行聚类 或者NTP算法 进行再次分析
#Cancer cell 
#javalin
#GSE 

##免疫推测 TIDE submap
###immune violion cluster####

###gene and immune cell cor####
### 用ssGSEA来量化浸润水平
### 1.加载marker
load(file = "/t8a/2022backup/xuzijun/vip39/data/cellMarker_ssGSEA.Rdata")
### 2.加载表达量
#load(file = "tcga_panmRNA_expr.Rdata")
expr <- tpms[,sub$Sample]
expr[1:5,1:5]
expr <- as.matrix(expr)
library(GSVA)
### 挺耗时间的，调用了12个线程，17:12开始, 17:37结束
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")
gsva_data[1:4,1:4]
tcga_gsva <- as.data.frame(t(gsva_data))
tcga_gsva[1:4,1:4]

tcga_expr <- tpms[,sub$Sample]
gene <- gene_list
immuscore <- function(gene){
  y <- as.numeric(tcga_expr[gene,])
  colnames <- colnames(tcga_gsva)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_gsva[,x]), y , method="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}

#以FOXP3为例，测试一下函数
immuscore("OXSM")
genelist <- gene_list
data <- do.call(rbind,lapply(genelist,immuscore))
head(data)
data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")
data$pstar[1:20]
ggplot(data, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "black",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=pstar),col ="black",size = 5)+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(angle = 45, hjust = 1),# 调整x轴文字
        axis.text.y = element_text(size = 8))+#调整y轴文字
  #调整legen
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))


setwd("./ssgsea_cor/")
tcga_gsva[1:4,1:4]
tcga_gsva <- as.data.frame(t(tcga_gsva))

tcga_gsva1 <- tcga_gsva[,sub$Sample[sub$Cluster=="C4"]]
tcga_gsva1[1:3,1:3] 
tcga_gsva1 <- as.data.frame(t(tcga_gsva1))##列为细胞或者通路分数
tcga_gsva2 <- tcga_gsva[,sub$Sample[sub$Cluster!="C4"]]
tcga_gsva2[1:3,1:3] 
tcga_gsva2 <- as.data.frame(t(tcga_gsva2))

#这里将计算每列之间的相关性
#如果要计算行之间的相关性就运行下面这行转置（把行变成列，列变成行）
#tcga_gsva <- t(tcga_gsva)
# 计算相关系数
cor_r <- cor(tcga_gsva1)
cor_r[1:3,1:3]
write.csv(cor_r, "easy_input_R1.csv", quote = F)
cor_r <- cor(tcga_gsva2)
cor_r[1:3,1:3]
write.csv(cor_r, "easy_input_R2.csv", quote = F)

# 计算p value
tcga_gsva <- tcga_gsva1
cor_p <- matrix(0, nrow = ncol(tcga_gsva), ncol = ncol(tcga_gsva))
rownames(cor_p) <- colnames(tcga_gsva)
colnames(cor_p) <- colnames(tcga_gsva)
for (i in 1:ncol(tcga_gsva)){
  for (j in 1:ncol(tcga_gsva)){
    p <- cor.test(tcga_gsva[,i],tcga_gsva[,j])
    cor_p[i,j] <- p$p.value
  }
}
write.csv(cor_p, "easy_input_P1.csv", quote = F)

tcga_gsva <- tcga_gsva2
cor_p <- matrix(0, nrow = ncol(tcga_gsva), ncol = ncol(tcga_gsva))
rownames(cor_p) <- colnames(tcga_gsva)
colnames(cor_p) <- colnames(tcga_gsva)
for (i in 1:ncol(tcga_gsva)){
  for (j in 1:ncol(tcga_gsva)){
    p <- cor.test(tcga_gsva[,i],tcga_gsva[,j])
    cor_p[i,j] <- p$p.value
  }
}
write.csv(cor_p, "easy_input_P2.csv", quote = F)

#KICH是C4  KIRC是其他的
KICHR <- read.csv("easy_input_R1.csv", row.names = 1)
KICHR[1:3,1:3]
KICHP <- read.csv("easy_input_P1.csv", row.names = 1)
KICHP[1:3,1:3]
KIRCR <- read.csv("easy_input_R1.csv", row.names = 1)
KIRCP <- read.csv("easy_input_P1.csv", row.names = 1)

## 合并相关系数的数据
datR <- KICHR
for(i in 1:nrow(datR)){
  datR[i,1:i] <- KIRCR[i,1:i]
}
datR[1:3,1:3]

## 合并P值的数据
datP <- KICHP
for (i in 1:nrow(datP)) {
  datP[i,1:i] <- KIRCP[i,1:i]
}
datP[1:3,1:3]

# P>= 0.05时，把相关系数设为NA
datR[datP > 0.05] <- NA

## 定义左右两个不同相关系数图的颜色
# 定义右上部分图形的颜色
colCorRight <-  circlize::colorRamp2(c(-1, 0, 1), c("green", "white", "#ef3b2c"))
# 定义左上部分图形的颜色
colCorLeft <- circlize::colorRamp2(c(-1, 0, 1), c("yellow", "white", "#762a83"))

## 绘制基本图形
p1 <- Heatmap(datR, rect_gp = gpar(type = "none"), 
              show_heatmap_legend = F,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.rect(x = x, y = y, width = width, height = height,
                          gp = gpar(col = "grey", fill = NA))
                if(i == j) {
                  grid.circle(x = x, y = y, r = 0.5 * min(unit.c(width, height)), gp = gpar(fill = "grey", col = NA))
                }else if(i > j) {
                  grid.circle(x = x, y = y, r = abs(datR[i, j])/2 * min(unit.c(width, height)), 
                              gp = gpar(fill = colCorLeft(datR[i, j]), col = NA))
                } else {
                  grid.circle(x = x, y = y, r = abs(datR[i, j])/2 * min(unit.c(width, height)), 
                              gp = gpar(fill = colCorRight(datR[i, j]), col = NA))
                }
              },
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = T, show_column_names = T, 
              row_names_side = "right", 
              row_names_rot = 45,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8)
)

p1

## 绘制两个不同不同相关的图例
lgdRight <- Legend(col_fun = colCorRight, title = "KICH", 
                   direction = "horizontal")
lgdLeft <- Legend(col_fun = colCorLeft, title = "KIRC", 
                  direction = "horizontal")
pd = list(lgdRight, lgdLeft)

## 最后出图
pdf("DouleCorPlotssgea.pdf", width = 5, height = 5.5)
draw(p1, annotation_legend_list = pd,
     annotation_legend_side = "top")
dev.off()


###immune signature boxplot####
#多组间分数比较 标注多个p值!
RNA_stemness <- fread("/t8a/2022backup/xuzijun/vip39/data/pancancer/StemnessScores_RNAexp_20170127.2.tsv.gz",data.table = F)
#会有RNAss EREG.EXPss两个分数
RNA_stemness <- as.data.frame(t(RNA_stemness))
head(RNA_stemness)
colnames(RNA_stemness) <- RNA_stemness[1,]
RNA_stemness <- RNA_stemness[-1,]
head(RNA_stemness)
rownames(RNA_stemness) <- str_replace_all(rownames(RNA_stemness),"\\.","\\-")
table(rownames(sub)%in%rownames(RNA_stemness))
# FALSE  TRUE 
# 17   499
RNA_stemness$Sample <- rownames(RNA_stemness)

subdata <- merge(sub,RNA_stemness,by="Sample")
head(subdata)
#subdata$Cluster <- factor(subdata$Cluster,levels = c("C1","C2"))

DNA_stemness <- fread("/t8a/2022backup/xuzijun/vip39/data/pancancer/StemnessScores_DNAmeth_20170210.tsv.gz",data.table = F)
DNA_stemness <- as.data.frame(t(DNA_stemness))
head(DNA_stemness)
colnames(DNA_stemness) <- DNA_stemness[1,]
DNA_stemness <- DNA_stemness[-1,]
head(DNA_stemness)
library(stringr)
rownames(DNA_stemness) <- str_replace_all(rownames(DNA_stemness),"\\.","\\-")
table(rownames(sub)%in%rownames(DNA_stemness))

DNA_stemness$Sample <- rownames(DNA_stemness)
subdata2 <- merge(sub,DNA_stemness,by="Sample")
head(subdata2)

library(ggpubr)

my_comparisons <- list(c("C1", "C2"), c("C2", "C3"), c("C1", "C3"))
subdata$RNAss <- as.numeric(subdata$RNAss)
p <- ggboxplot(subdata, x="Cluster",
               y = "RNAss", color = "Cluster",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
ggviolin(subdata, x = "Cluster", y = "RNAss",
         fill = "Cluster", palette = c(pal_nejm()(4),'grey60'),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c(pal_nejm()(4),'grey60'))) + 
  stat_compare_means()
#+stat_compare_means(label.y = 0.5)

ggboxplot(subdata, x="Cluster", y="RNAss",color = "Cluster",
          palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"), add = "jitter")+
  ### 添加多组比较的统计学结果
  stat_compare_means(label.y = max(subdata$RNAss)+0.3)+
  ### 添加每个两两比较的显著性标记位置信息
  stat_compare_means(comparisons=compaired, method = "wilcox.test",
                     label ="p.value")
ggsave(filename = "RNAss.pdf",height = 5,width = 5)


subdata$EREG.EXPss <- as.numeric(subdata$EREG.EXPss)
subdata$Cluster <- factor(subdata$Cluster,levels = c("C1","C2","C3","C4"))
ggviolin(subdata, x = "Cluster", y = "EREG.EXPss",
         fill = "Cluster", palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"))) + 
  stat_compare_means(label = "p.signif")

ggboxplot(subdata, x="Cluster", y="EREG.EXPss",color = "Cluster",
          palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"), add = "jitter")+
  ### 添加多组比较的统计学结果
  stat_compare_means(label.y = max(subdata$EREG.EXPss)+0.5)+
  ### 添加每个两两比较的显著性标记位置信息
  stat_compare_means(comparisons=compaired, method = "wilcox.test",
                     label ="p.value")
ggsave(filename = "EREG.EXPss.pdf",height = 5,width = 5)

subdata2$DNAss <- as.numeric(subdata2$DNAss)
#subdata2$Cluster <- factor(subdata2$Cluster,levels = c("C1","C2","C3"))
ggviolin(subdata2, x = "Cluster", y = "DNAss",
         fill = "Cluster", palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"))) + 
  stat_compare_means( label = "p.signif")

ggboxplot(subdata2, x="Cluster", y="DNAss",color = "Cluster",
          palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"), add = "jitter")+
  ### 添加多组比较的统计学结果
  stat_compare_means(label.y = max(subdata2$DNAss)+0.1)+
  ### 添加每个两两比较的显著性标记位置信息
  stat_compare_means(comparisons=compaired, method = "wilcox.test",
                     label ="p.value")
ggsave(filename = "DNAss.pdf",height = 5,width = 5)

subdata2$DMPss<- as.numeric(subdata2$DMPss)
ggviolin(subdata2, x = "Cluster", y = "DMPss",
         fill = "Cluster", palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"))) + 
  stat_compare_means(label = "p.signif")

ggboxplot(subdata2, x="Cluster", y="DMPss",color = "Cluster",
          palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"), add = "jitter")+
  ### 添加多组比较的统计学结果
  stat_compare_means(label.y = max(subdata2$DMPss)+0.1)+
  ### 添加每个两两比较的显著性标记位置信息
  stat_compare_means(comparisons=compaired, method = "wilcox.test",
                     label ="p.value")
ggsave(filename = "DMPss.pdf",height = 5,width = 5)


subdata2$ENHss<- as.numeric(subdata2$ENHss)
subdata2$Cluster <- factor(subdata2$Cluster,levels = c("C1","C2","C3","C4"))
ggviolin(subdata2, x = "Cluster", y = "ENHss",
         fill = "Cluster", palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"))) + 
  stat_compare_means(label = "p.signif")

ggboxplot(subdata2, x="Cluster", y="ENHss",color = "Cluster",
          palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"), add = "jitter")+
  ### 添加多组比较的统计学结果
  stat_compare_means(label.y = max(subdata2$ENHss)+0.1)+
  ### 添加每个两两比较的显著性标记位置信息
  stat_compare_means(comparisons=compaired, method = "wilcox.test",
                     label ="p.value")
ggsave(filename = "ENHss.pdf",height = 5,width = 5)




HRD_score <- fread("/t8a/2022backup/xuzijun/vip39/data/pancancer/TCGA.HRD_withSampleID.txt.gz",data.table = F)
HRD_score <- as.data.frame(t(HRD_score))
head(HRD_score)
colnames(HRD_score) <- HRD_score[1,]
HRD_score <- HRD_score[-1,]
head(HRD_score)
#rownames(DNA_stemness) <- str_replace_all(rownames(DNA_stemness),"\\.","\\-")
table(rownames(sub)%in%rownames(HRD_score))

HRD_score$Sample <- rownames(HRD_score)
subdata4 <- merge(sub,HRD_score,by="Sample")
head(subdata4)

subdata4$HRD <- as.numeric(subdata4$HRD)
subdata4$Cluster <- factor(subdata4$Cluster,levels = c("C1","C2","C3","C4"))
ggviolin(subdata4, x = "Cluster", y = "HRD",
         fill = "Cluster", palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"))) + 
  stat_compare_means(label = "p.signif")

####美化版plot
library(ggpubr)
#load('df.merge.RData')
compaired <- list(c("C1", "C2"), 
                  c("C1", "C3"), 
                  c("C1", "C4"),
                  c("C2", "C3"),
                  c("C2", "C4"),
                  c("C3", "C4"))   #设置比较组别
library(ggsci)
#colnames(df.merge)
col_1<-pal_nejm()(8)
ggboxplot(subdata4, x="Cluster", y="HRD",color = "Cluster",
          palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"), add = "jitter")+
  ### 添加多组比较的统计学结果
  stat_compare_means(label.y = max(subdata4$HRD)+0.1)+
  ### 添加每个两两比较的显著性标记位置信息
  stat_compare_means(comparisons=compaired, method = "wilcox.test",
                     label ="p.value")
ggsave(filename = "HRD.pdf",height = 5,width = 5)

colnames(subdata4)
subdata4$`hrd-loh`
subdata4$hrd_loh <- as.numeric(subdata4$`hrd-loh`)
ggplot(subdata4,aes(x=Cluster,y=hrd_loh,fill=Cluster))+
  geom_bar(stat="summary",width = 0.5,fun.data = 'mean_sd')+
  scale_fill_manual(values = col_1)+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
               width = 0.25,size=0.8,position = position_dodge( .9))+
  geom_jitter( size =3,
               alpha = 0.1,
               shape = 21,stroke = 1) +
  theme_classic()+
  stat_compare_means(comparisons=compaired,method = "wilcox.test",
                     label = "p.forma",bracket.size = 0.8,size=4)+
  theme(legend.position = 'top')

ggboxplot(subdata4, x="Cluster", y="hrd_loh",color = "Cluster",
          palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"), add = "jitter")+
  ### 添加多组比较的统计学结果
  stat_compare_means(label.y = max(subdata4$hrd_loh)+5)+
  ### 添加每个两两比较的显著性标记位置信息
  stat_compare_means(comparisons=compaired, method = "wilcox.test",
                     label ="p.value")
ggsave(filename = "hrd_loh.pdf",height = 5,width = 5)



#####estimate IOBR related signature#####
#setwd("/home/data/vip39/TCGA/jiaowang_KIRC/met_model/met20211216/estimate")
dat=log2(KIRC_mRNA_fpkm[,sub$Sample]+1)
dat[1:4,1:4]
library(estimate)
estimate <- function(dat,pro){ 
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate) 
  filterCommonGenes(input.f=input.f, 
                    output.f=output.f ,
                    id="GeneSymbol") 
  estimateScore(input.ds = output.f,
                output.ds=output.ds, 
                platform="illumina")  
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro='Sdeath'
scores=estimate(dat,pro)
scores[1:4,]
gsva_es_imm_left <- data.frame(ID=sub$Sample,
                               TMEcluster=sub$Cluster)

scores <- as.data.frame(scores)
rownames(scores) <- str_replace_all(rownames(scores),"\\.","\\-")
table(gsva_es_imm_left$ID==rownames(scores))

scores <- scores[gsva_es_imm_left$ID,]
pdata <- cbind(gsva_es_imm_left,scores)

head(pdata)

ggviolin(pdata, x = "TMEcluster", y = "ESTIMATEScore",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means(label = "p.signif")
ggsave(filename = "ESTIMATEScore.pdf",height = 5,width = 6)

ggviolin(pdata, x = "TMEcluster", y = "ImmuneScore",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means(label = "p.signif")
ggsave(filename = "ImmuneScore.pdf",height = 5,width = 6)

ggviolin(pdata, x = "TMEcluster", y = "StromalScore",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means( label = "p.signif")
ggsave(filename = "StromalScore.pdf",height = 5,width = 6)


pdata[,3:ncol(pdata)] <- scale(pdata[,3:ncol(pdata)],scale = T,center = T)
#数据宽转长
pdata_melt <- melt(pdata,
                   id.vars = c("ID","TMEcluster"),
                   variable.name ="Signature",
                   value.name = "Signature_score")
head(pdata_melt)



# 使用ggplo2画图
c <- ggplot(pdata_melt,
            aes(x=Signature, y=Signature_score, 
                fill = TMEcluster, #按类填充颜色
                color = TMEcluster)) + #按类给边框着色
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", #outlier点用黑色
               outlier.size = 0.65) +
  #自定义配色
  scale_fill_manual(values= c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF")) +
  #ggtitle("Gene signature score", "stratified by TME-cluster") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top")

# 标注p value
c + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

# 标注*
# `****` = 1e-04, `***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1
p <- c + stat_compare_means(label = "p.signif")
p

ggsave(p, filename = "TME-relevant-estimate-boxplot.pdf", width = 8, height = 6)

####TIP####
stepscore <- fread("/t8a/2022backup/xuzijun/vip39/database/TIP/KIRC/ssGSEA.normalized.score.txt",header = T,data.table = F)

head(stepscore[,1:4])
stepscore <- stepscore%>%
  column_to_rownames("Steps")
colnames(stepscore) <- str_sub(colnames(stepscore),1,15)
rownames(stepscore)[1:3] <- c("Step1.Release of cancer cell antigens",
                              "Step2.Cancer antigen presentation",
                              "Step3.Priming and activation")
rownames(stepscore)[c(21,22,23)] <- c("Step5.Infiltration of cancer cells into tumors",
                                      "Step6.Recognitioin of cancer cells by T cells",
                                      "Step7.Killing of cancer cells")
#stepscore <- as.data.frame(t(stepscore))
gsva_es_imm <- stepscore[,sub$Sample]
gsva_es_imm <- as.data.frame(t(gsva_es_imm))
gsva_es_imm_left <- data.frame(ID=sub$Sample,
                               TMEcluster=sub$Cluster)
pdata <- cbind(gsva_es_imm_left,gsva_es_imm)

head(pdata)
pdata[,3:ncol(pdata)] <- scale(pdata[,3:ncol(pdata)],scale = T,center = T)
#数据宽转长
pdata_melt <- melt(pdata,
                   id.vars = c("ID","TMEcluster"),
                   variable.name ="Signature",
                   value.name = "Signature_score")
head(pdata_melt)


# 使用ggplo2画图
par(oma=c(1,1,1,1), mar=c(2,2,2,2))
c <- ggplot(pdata_melt,
            aes(x=Signature, y=Signature_score, 
                fill = TMEcluster, #按类填充颜色
                color = TMEcluster)) + #按类给边框着色
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", #outlier点用黑色
               outlier.size = 0.65) +
  #自定义配色
  scale_fill_manual(values= c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF")) +
  #ggtitle("Gene signature score", "stratified by TME-cluster") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top")

# 标注p value
c + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

# 标注*
# `****` = 1e-04, `***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1
p <- c + stat_compare_means(label = "p.signif")
p
ggsave(p, filename = "TME-relevant-TIP-boxplot.pdf", width = 20, height = 6)

#####fig146####
(load("/t8a/2022backup/xuzijun/vip39/huitu/FigureYa146TMEbox/signature.RData"))
head(signature)
gsym.expr[1:4,1:4]
gsva_es_imm <- gsva(as.matrix(gsym.expr), signature)
head(gsva_es_imm)

# 分组
gsva_es_imm <- gsva_es_imm[,sub$Sample]
gsva_es_imm <- as.data.frame(t(gsva_es_imm))
gsva_es_imm_left <- data.frame(ID=sub$Sample,
                               TMEcluster=sub$Cluster)
pdata <- cbind(gsva_es_imm_left,gsva_es_imm)

head(pdata)
pdata[,3:ncol(pdata)] <- scale(pdata[,3:ncol(pdata)],scale = T,center = T)
#数据宽转长
pdata_melt <- melt(pdata,
                   id.vars = c("ID","TMEcluster"),
                   variable.name ="Signature",
                   value.name = "Signature_score")
head(pdata_melt)


# 使用ggplo2画图
c <- ggplot(pdata_melt,
            aes(x=Signature, y=Signature_score, 
                fill = TMEcluster, #按类填充颜色
                color = TMEcluster)) + #按类给边框着色
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", #outlier点用黑色
               outlier.size = 0.65) +
  #自定义配色
  scale_fill_manual(values= c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF")) +
  #ggtitle("Gene signature score", "stratified by TME-cluster") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top")

# 标注p value
c + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

# 标注*
# `****` = 1e-04, `***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1
p <- c + stat_compare_means(label = "p.signif")
p

ggsave(p, filename = "TME-relevant-sigature-boxplot.pdf", width = 13, height = 6)

###TIDE####
#行是基因 列是样本 fpkm不需要log2处理
#setwd("~/TCGA/jiaowang_KIRC/met_model/met20211216/TIDE/")
dat_TIDE <- KIRC_mRNA_fpkm[,sub$Sample]
TIDE <- dat_TIDE
# 这里为了得到比较好的结果，采用two direction median centered
TIDE <- sweep(TIDE,2, apply(TIDE,2,median,na.rm=T))
TIDE <- sweep(TIDE,1, apply(TIDE,1,median,na.rm=T))
write.table(TIDE,"TIDE_input.self_subtract",sep = "\t",row.names = T,col.names = NA,quote = F)
TIDE.res <- read.csv("~/TCGA/jiaowang_KIRC/met_model/met20211216/TIDE/met_TIDE_result.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
TIDE.res[1:4,1:4]

###绘制不同分组的TIDE评分图
table(sub$Sample%in%rownames(TIDE.res))
TIDE.res <- TIDE.res[sub$Sample,]

gsva_es_imm_left <- data.frame(ID=sub$Sample,
                               TMEcluster=sub$Cluster)
pdata <- cbind(gsva_es_imm_left,TIDE.res)

head(pdata)
#save(pdata,file = "pdata_TIDE.rds")
# 默认方法为 method = "kruskal.test"
library(ggpubr)

ggboxplot(pdata, x="TMEcluster", y="TIDE",color = "TMEcluster",
          palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"), add = "jitter")+
  ### 添加多组比较的统计学结果
  stat_compare_means(label.y = max(pdata$TIDE)+5)+
  ### 添加每个两两比较的显著性标记位置信息
  stat_compare_means(comparisons=compaired, method = "wilcox.test",
                     label ="p.value")
ggsave(filename = "TIDE.pdf",height = 5,width = 5)

ggboxplot(pdata, x="TMEcluster", y="Dysfunction",color = "TMEcluster",
          palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"), add = "jitter")+
  ### 添加多组比较的统计学结果
  stat_compare_means(label.y = max(pdata$Dysfunction)+2)+
  ### 添加每个两两比较的显著性标记位置信息
  stat_compare_means(comparisons=compaired, method = "wilcox.test",
                     label ="p.value")
ggsave(filename = "Dysfunction.pdf",height = 5,width = 5)

ggboxplot(pdata, x="TMEcluster", y="MSI Expr Sig",color = "TMEcluster",
          palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"), add = "jitter")+
  ### 添加多组比较的统计学结果
  stat_compare_means(label.y = max(pdata$`MSI Expr Sig`)+1)+
  ### 添加每个两两比较的显著性标记位置信息
  stat_compare_means(comparisons=compaired, method = "wilcox.test",
                     label ="p.value")
ggsave(filename = "MSI.pdf",height = 5,width = 5)

ggboxplot(pdata, x="TMEcluster", y="Exclusion",color = "TMEcluster",
          palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"), add = "jitter")+
  ### 添加多组比较的统计学结果
  stat_compare_means(label.y = max(pdata$Exclusion)+3)+
  ### 添加每个两两比较的显著性标记位置信息
  stat_compare_means(comparisons=compaired, method = "wilcox.test",
                     label ="p.value")
ggsave(filename = "Exclusion.pdf",height = 5,width = 5)

ggboxplot(pdata, x="TMEcluster", y="CAF",color = "TMEcluster",
          palette = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"), add = "jitter")+
  ### 添加多组比较的统计学结果
  stat_compare_means(label.y = max(pdata$CAF)+0.5)+
  ### 添加每个两两比较的显著性标记位置信息
  stat_compare_means(comparisons=compaired, method = "wilcox.test",
                     label ="p.value")
ggsave(filename = "CAF.pdf",height = 5,width = 5)



ggviolin(pdata, x = "TMEcluster", y = "TIDE",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means(label = "p.signif")
ggsave(filename = "TIDE.pdf",height = 5,width = 6)
# ggviolin(pdata, x = "TMEcluster", y = "TIDE",
#          palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",
#          add = "boxplot",add.params = list(fill=c("#AFEEEE","#FFE4E1"))) + 
#   stat_compare_means()
ggviolin(pdata, x = "TMEcluster", y = "Dysfunction",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means( label = "p.signif")
ggsave(filename = "Dysfunction.pdf",height = 5,width = 6)

# ggviolin(pdata, x = "TMEcluster", y = "Dysfunction",
#          fill = "TMEcluster", palette = c("#48D1CC", "#E9967A"),alpha=.8,color = "white",add = "boxplot",
#          add.params = list(fill=c("#AFEEEE","#FFE4E1"))) + 
#   stat_compare_means()

ggviolin(pdata, x = "TMEcluster", y = "MSI Expr Sig",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means( label = "p.signif")
ggsave(filename = "MSI Expr Sig.pdf",height = 5,width = 6)


ggviolin(pdata, x = "TMEcluster", y = "Exclusion",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means(label = "p.signif")
ggsave(filename = "Exclusion.pdf",height = 5,width = 6)

ggviolin(pdata, x = "TMEcluster", y = "CAF",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means( label = "p.signif")
ggsave(filename = "CAF.pdf",height = 5,width = 6)


library(dplyr)
meta = data.frame(Sample=pdata$ID,
                  Response=pdata$Responder,
                  risk=pdata$TMEcluster)

#一句话搞定
library(ggstatsplot)
ggbarstats(data = meta,x=Response,y=risk)
ggsave(filename = "Immune_response.pdf",height = 6,width = 6)
###治疗反应###


##CYT评分 融细胞分数####
sub$GZMA <- as.numeric(tpms["GZMA",sub$Sample])
sub$PRF1 <- as.numeric(tpms["PRF1",sub$Sample])
sub$CYT=sqrt(sub$GZMA*sub$PRF1)##sqrt为开根号函数


sub$Cluster <- factor(sub$Cluster,levels = c("C1","C2","C3","C4"))##把MSI变成因子型，排序一下
library(ggplot2)
library(ggpubr)
pdf(file = 'TCGA-KIRC-CYT.pdf',width = 5,height = 5)
ggplot(sub, aes(x = Cluster, y = CYT, fill = Cluster)) +
  geom_violin(scale="width") +
  geom_boxplot(width=.05, 
               fill='grey40', 
               notch=F,
               outlier.color=NA, 
               col="grey40") +
  stat_summary(fun="median", 
               geom="point", shape=20, col="white")+ ##画中位数小白点
  labs(x='',y='CYT activity')+
  scale_fill_manual(values = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"))+
  theme(legend.position = 'none',
        axis.title.y = element_text(color = 'grey40'))+
  stat_compare_means(comparisons=compaired, method = "wilcox.test",
                     label ="p.value")#label.x调整P值的水平位置
dev.off()




####其他评分 pan cancer 项目#####



###和经典的KIRC分型关联 柱状图#####


#####fig255####
###fig255免疫热图
setwd("/home/aim/Desktop/TCGA_work/vip33new/Sdeath_project/result/immnue_heatmap")
library(utils)
library(GSVA)
library(ComplexHeatmap) # 用到里面的pheatmap画图然后拼图，需安装2.8以上版本的ComplexHeatmap
source("/t8a/2022backup/xuzijun/vip39/huitu/FigureYa255TIME/pheatmap_translate.R") # 如果不想放弃较老版本的R及其对应的老ComplexHeatmap，可以只加载新版ComplexHeatmap里的这个函数，该脚本出自2.9.4版ComplexHeatmap
library(circlize)
library(viridis)
library(gplots)
library(data.table)
library(estimate)
source("/t8a/2022backup/xuzijun/vip39/huitu/FigureYa255TIME/annTrackScale.R") # 加载函数，对数值进行标准化并对极端值做截断
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


# 加载自定义函数
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}


# 读入用MOVICS获取的muscle-invasive bladder cancer (MIBC)亚型

annCol.tcga <- data.frame(CMOIC=paste0("CS",str_sub(sub$Cluster,2,2)))
rownames(annCol.tcga) <- sub$Sample

# 读取450K甲基化数据

library(data.table)
library(impute)
library(ChAMP)
library(stringr)
library(tibble)
options(stringsAsFactors = F)
meth <- data.table::fread('/t8a/vip39database/database/UCSC_TCGA/TCGA_methylation450/TCGA-KIRC.methylation450.tsv.gz')
a <- meth
a = column_to_rownames(a,"Composite Element REF")
colnames(a)= str_sub(colnames(a),1,15)
#load("/home/data/vip39/TCGA/Immune_MOVICS_sub/KIRC_meth_metil.rds")
meth <- a

# 匹配亚型
colnames(meth) <- str_sub(colnames(meth),1,15)
#meth <- meth[,which(substr(colnames(meth), 14, 15) == "01")]
table(rownames(annCol.tcga)%in%colnames(meth)) #306
interid <- intersect(rownames(annCol.tcga),colnames(meth))
annCol.tcga <- as.data.frame(annCol.tcga)
annCol.tcga <-  as.data.frame(annCol.tcga[interid,])
rownames(annCol.tcga) <- interid
colnames(annCol.tcga) <- "CMOIC"
meth <- meth[,interid] 

MeTIL.marker <- c("cg20792833","cg20425130","cg23642747","cg12069309","cg21554552")
meth[1:4,1:4]
meth.metil <- meth[MeTIL.marker,]
MeTIL <- meth[MeTIL.marker,]
MeTIL <- t(scale(t(MeTIL)))

# 计算MeTIL得分
pca.MeTIL <- prcomp(MeTIL,center = F,scale. = F)
MeTIL <- pca.MeTIL$rotation[,1]
annCol.tcga$MeTIL <- as.numeric(MeTIL)

# 甲基化数据只要用到5个探针的marker就可以，所以这里我的输入数据是简化的，方便传输
#load("/t8a/2022backup/xuzijun/vip39/TCGA/Immune_MOVICS_sub/KIRC_meth_metil.rds")
#meth <- meth.metil

# 匹配亚型
#colnames(meth) <- substr(colnames(meth), start = 1,stop = 16)
#meth <- meth[,which(substr(colnames(meth), 14, 15) == "01")]
# meth <- meth[,rownames(annCol.tcga)] 
# 
# MeTIL.marker <- c("cg20792833","cg20425130","cg23642747","cg12069309","cg21554552")
# meth.metil <- meth[MeTIL.marker,]
# MeTIL <- meth[MeTIL.marker,]
# MeTIL <- t(scale(t(MeTIL)))
# 
# # 计算MeTIL得分
# pca.MeTIL <- prcomp(MeTIL,center = F,scale. = F)
# MeTIL <- pca.MeTIL$rotation[,1]
# annCol.tcga$MeTIL <- MeTIL



# 加载表达数据
tpm <- KIRC_mRNA_fpkm[,rownames(annCol.tcga)]
immune.signature <- read.table("/t8a/2022backup/xuzijun/vip39/huitu/FigureYa255TIME/Curated_Immune_Cell_Signature.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)

# 构建计算GSVA的列表
cell.type <- unique(immune.signature$CellType)
immune.sig.ccr <- list()
for (i in cell.type) {
  immune.sig.ccr[[i]] <- immune.signature[which(immune.signature$CellType == i),"Symbol"]
}

# 免疫检查点相关基因
imm.targets <- c("CD274","PDCD1","CD247","PDCD1LG2","CTLA4","TNFRSF9","TNFRSF4","TLR9") 

# 免疫细胞的排序
immune.sig.ccr.order <- c("T.cells.CD8", 
                          "T.cells.regulatory..Tregs.",
                          "T.cells.CD4.naive",
                          "T.cells.follicular.helper",
                          "B.cells.naive",
                          "B.cells.memory",
                          "T.cells.gamma.delta",
                          "Dendritic.cells.activated",
                          "Macrophages.M1",
                          "NK.cells.activated",
                          "Plasma.cells",
                          "T.cells.CD4.memory.resting",
                          "T.cells.CD4.memory.activated",
                          "Mast.cells.activated",
                          "NK.cells.resting",
                          "Macrophages.M0",
                          "Macrophages.M2",
                          "Eosinophils",
                          "Monocytes",
                          "Dendritic.cells.resting",
                          "Mast.cells.resting",
                          "Neutrophils",
                          "Endothelial cells",
                          "Fibroblasts")

# 计算immune/stromal得分
range(tpm)
indata <- tpm
#indata <- log2(tpm + 1)
# 保存到文件
write.table(indata,file = "TCGA_log2TPM_hugo.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

filterCommonGenes(input.f = "TCGA_log2TPM_hugo.txt" , output.f="TCGA_log2TPM_hugo_ESTIMATE.txt", id="GeneSymbol")

estimateScore("TCGA_log2TPM_hugo_ESTIMATE.txt","TCGA_log2TPM_hugo_estimate_score.txt", platform="affymetrix")

est.tcga <- read.table(file = "TCGA_log2TPM_hugo_estimate_score.txt",header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.tcga) <- est.tcga[,2]; colnames(est.tcga) <- est.tcga[1,]; est.tcga <- est.tcga[-1,c(-1,-2)];
est.tcga <- sapply(est.tcga, as.numeric); rownames(est.tcga) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.tcga.backup = as.data.frame(est.tcga); colnames(est.tcga.backup) <- colnames(indata)

# 对数值进行标准化并对极端值做截断
est.tcga <- annTrackScale(indata = est.tcga, halfwidth = 2, poolsd = F); est.tcga <- as.data.frame(t(est.tcga)) 
rownames(est.tcga) <- colnames(tpm)

tcga.immune.gsva <- gsva(as.matrix(log2(tpm + 1)),
                         immune.sig.ccr,
                         method = "gsva")

# 设置颜色
clust.col <- c("#DD492E","#40548A","#32A087","#EC7D21")
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")
blue <- "#5bc0eb"
gold <- "#ECE700"
cyan <- "#00B3D0"

annCol.tcga$ImmuneScore <- as.numeric(est.tcga[rownames(annCol.tcga),"ImmuneScore"])
annCol.tcga$StromalScore <- as.numeric(est.tcga[rownames(annCol.tcga),"StromalScore"])
annCol.tcga <- annCol.tcga[order(annCol.tcga$CMOIC),]
annColors.tcga <- list() # 构建热图的图例颜色
annColors.tcga[["CMOIC"]] <- c("CS1" = clust.col[1],
                               "CS2" = clust.col[2],
                               "CS3" = clust.col[3],
                               "CS4" = clust.col[4]
)
annColors.tcga[["ImmuneScore"]] <- inferno(64)
annColors.tcga[["StromalScore"]] <- viridis(64)

## 热图1：免疫检查点基因表达
indata <- log2(tpm[intersect(rownames(tpm),imm.targets),] + 1)
hm1 <- pheatmap(standarize.fun(indata[imm.targets,rownames(annCol.tcga)],halfwidth = 2), # 表达谱数据标准化
                border_color = NA, # 热图单元格无边框
                annotation_col = annCol.tcga[,c("CMOIC","StromalScore","ImmuneScore")],
                annotation_colors = annColors.tcga[c("CMOIC","StromalScore","ImmuneScore")],
                color = colorpanel(64,low=blue,mid = "black",high=gold),
                show_rownames = T, # 显示行名
                show_colnames = F, # 不显示列名
                cellheight = 12, # 热图高度固定
                cellwidth = 0.6, # 热图宽度固定
                name = "ICI", # 图例名字
                cluster_rows = F, # 行不聚类
                cluster_cols = F) # 列不聚类

#pdf("CheckPoint_heatmap.pdf",width = 10,height = 10)
hm1

dev.off()

## 热图2：肿瘤免疫微环境富集得分
hm2 <- pheatmap(standarize.fun(tcga.immune.gsva[immune.sig.ccr.order,rownames(annCol.tcga)],halfwidth = 1), # 富集得分标准化并排序
                border_color = NA,
                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                gaps_row = c(14,22), # 根据不同类别的免疫细胞分割
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 0.6,
                name = "TIME",
                cluster_rows = F,
                cluster_cols = F)

#pdf("TIMEheatmap.pdf",width = 10,height = 10)
hm2

#dev.off()

## 热图3：MeTIL得分
hm3 <- pheatmap(standarize.fun(t(annCol.tcga[,"MeTIL",drop = F]),halfwidth = 1),
                border_color = NA,
                color = NMF:::ccRamp(c(cyan,"black","#F12AFE"),64),
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 0.6,
                name = "MeTIL",
                cluster_rows = F,
                cluster_cols = F)

#pdf("MeTILheatmap.pdf",width = 10,height = 10)
hm3
#dev.off()

# 合并热图并输出
pdf("TIME.pdf",width = 10,height = 10)
draw(hm1 %v% hm2 %v% hm3, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right")
invisible(dev.off())



