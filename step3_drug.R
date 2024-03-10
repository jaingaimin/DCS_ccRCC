



###maftools+genomic variation####


###drug resistence MOVICS + new 小丫 ####
##fig105
setwd("/home/aim/Desktop/TCGA_work/vip33new/Sdeath_project/result/figureYa105")
library(pRRophetic)
library(ggplot2)
library(cowplot)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
#药物名字
GCP.drug <- read.table("/t8a/2022backup/ssy088_202210/ssy088/tutulaile/FigureYa105GDSC/drug.txt") #如果要例文的两种药物，就换成drug_eg.txt
GCP.drug <- GCP.drug$V1

###Sdeath作图

#load(file = "/home/data/ssy088/tutulaile/FigureYa105GDSC/Sdeath/Sdeath_sub.rds")
table(sub$Cluster)

load(file = "/t8a/2022backup/ssy088_202210/ssy088/database/KITC_tumor_exp.Rdata")

ann_KIRC <- sub
ann_KIRC$ImmClust <- ann_KIRC$Cluster
#ann_KIRC$ImmClust <- ifelse(ann_KIRC$Cluster=="C1","C1",ifelse(ann_KIRC$Cluster=="C2","C2","C3"))
#ann_KIRC$ImmClust <- ifelse(ann_KIRC$Cluster=="C3","C3","C12")

dat_KIRC <- tumor_exp[,rownames(ann_KIRC)]
dat_KIRC[1:4,1:4]
dat <- dat_KIRC
ann <- ann_KIRC
GCP.drug

setwd("/home/aim/Desktop/TCGA_work/vip33new/Sdeath_project/result/figureYa105")
# 自定义足够多的box的颜色，颜色数量至少等于分组数量
jco <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF")

### 药物敏感性预测 ###、
GCPinfo <- GCP.IC50 <- GCP.expr <- cvOut <- predictedPtype <- predictedBoxdat <- list() # 初始化列表

plotp <- list()

for (drug in GCP.drug) {
  set.seed(1248103) # 因为预测过程默认10-fold CV，所以设置种子以便结果可重复
  cat(drug," starts!\n") # 提示当前药物已开始分析
  
  # 预测IC50值，采用默认参数，详细可参考??pRRopheticPredict参数列表
  predictedPtype[[drug]] <- pRRopheticPredict(testMatrix = as.matrix(dat[,rownames(ann)]),
                                              drug = drug,
                                              tissueType = "allSolidTumors",
                                              selection = 1) # 1表示若有重复基因取均值处理
  
  if(!all(names(predictedPtype[[drug]])==rownames(ann))) {stop("Name mismatched!\n")} # 若名字不匹配则报错退出
  
  predictedBoxdat[[drug]] <- data.frame("est.ic50"=predictedPtype[[drug]],
                                        #"ImmClust"=ifelse(ann$ImmClust == "C3","C3","C12"), # 这里我修改了C1和C2的名
                                        "ImmClust"= ann$ImmClust, # 这里我修改了C1和C2的名
                                        row.names = names(predictedPtype[[drug]])) 
  #predictedBoxdat[[drug]]$ImmClust <- factor(predictedBoxdat[[drug]]$ImmClust,levels = c("C12","C3"),ordered = T) # 把类改成因子变量
  predictedBoxdat[[drug]]$ImmClust <- factor(predictedBoxdat[[drug]]$ImmClust,levels = c("C1","C2","C3","C4"),ordered = T) # 把类改成因子变量
  # 绘图
  p <- ggplot(data = predictedBoxdat[[drug]], aes(x=ImmClust, y=est.ic50))
  p <- p + geom_boxplot(aes(fill = ImmClust)) + 
    scale_fill_manual(values = jco[1:length(unique(ann$ImmClust))]) + #自定义box的配色
    theme(legend.position="none") + # 倾斜字体
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),plot.title = element_text(size = 12, hjust = 0.5)) + 
    xlab("") + ylab("Estimated IC50") + 
    ggtitle(drug) # 补上title
  
  plotp[[drug]] <- p # 保存在列表里供合并图片用
  cat(drug," has been finished!\n") # 提示当前药物已分析结束
}

# 合并图片
#适合展示两种药物
p1 <- plot_grid(plotp[[1]],plotp[[2]],labels = c("A","B"),nrow = 1) # title可以AI下拉到合适位置，就如例文所示
ggsave("boxplot of predicted IC50.pdf", width = 6, height = 5)

# 适合展示多种药物
p2 <- plot_grid(plotlist=plotp, ncol=11)
ggsave("./jiaowang/boxplot of predicted IC50_multiple.pdf", width = 25, height = 15)

## 检验组间差异
p <- vector()

for (drug in GCP.drug) {
  tmp <- wilcox.test(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C2"),"est.ic50"]),
                     as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C1"),"est.ic50"]),alternative = "less")$p.value
  p <- append(p,tmp) # 两组样本秩和检验p值
}

for (drug in GCP.drug) {
  tmp <- kruskal.test(list(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C1"),"est.ic50"]),
                           as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C2"),"est.ic50"]),
                           as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C3"),"est.ic50"]),
                           as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C4"),"est.ic50"])))$p.value
  # tmp <- wilcox.test(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C2"),"est.ic50"]),
  #                    as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C1"),"est.ic50"]),alternative = "less")$p.value
  p <- append(p,tmp) # 两组样本秩和检验p值
}
names(p) <- GCP.drug
print(p) #打印p值，因为这里有一个不显著，所以当时没有放在boxplot上，有需要的也可以加在ggplot的title里，或参考FigureYa12box直接画在图上。
# names(p)[p<0.000000000000000000000001]
# names(p)[p<0.001]

sort(p,decreasing = F)
names(sort(p,decreasing = F))[1:10]

#保存到文件
write.table(p,"./output_pvalue0711.txt", quote = F, sep = "\t")
#save.image(file = "./all.rds")
#load(file = "./all.rds")

###CS1 CS2作图

plot_grid(plotp[["PAC.1"]],plotp[["Vinorelbine"]],plotp[["CI.1040"]],
          plotp[["Embelin"]],plotp[["SL.0101.1"]],plotp[["PD.0325901"]],
          plotp[["Nutlin.3a"]],plotp[["LFM.A13"]],plotp[["GSK.650394"]],
          plotp[["Gefitinib"]],
          nrow = 2)

####fig227 只是用于box plot可视化 不是药敏分析的脚本####


#fig213  CMap绘图 基于差异gene 适合两组分析  C4单独一组 其他单独一组#####
setwd("/home/aim/Desktop/TCGA_work/vip33new/Sdeath_project/result/fig213/")
source("/t8a/2022backup/ssy088/tutulaile/FigureYa213customizeHeatmap/pheatmap_translate.R") #来自github版ComplexHeatmap
library(devtools)
library(ComplexHeatmap) # 用于绘制热图
library(circlize) # 用于配色
library(tidyverse) # 用于读取MAF文件
library(limma) # 用于差异表达

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
##计算差异gene 150up 150down
C4sam <- sub$Sample[sub$Cluster=="C4"]
C123sam <- sub$Sample[sub$Cluster!="C4"]

pd <- data.frame(Samples = c(C4sam,C123sam),
                 Group = rep(c("tumor","normal"),c(length(C123sam),length(C4sam))),
                 stringsAsFactors = FALSE)

design <-model.matrix(~ -1 + factor(pd$Group, levels = c("tumor","normal")))
colnames(design) <- c("tumor","normal")

gset <- tpms[,pd$Samples]
range(gset)

fit <- limma::lmFit(gset, design = design);
contrastsMatrix <- limma::makeContrasts(tumor - normal, levels = c("tumor", "normal"))
fit2 <- limma::contrasts.fit(fit, contrasts = contrastsMatrix)
fit2 <- limma::eBayes(fit2, 0.01)

resData <- limma::topTable(fit2, adjust = "fdr", sort.by = "B", number = 100000)
resData <- as.data.frame(subset(resData, select=c("logFC","t","B","P.Value","adj.P.Val")))
resData$id <- rownames(resData)
colnames(resData) <- c("log2fc","t","B","pvalue","padj","id")
resData$fc <- 2^resData$log2fc
resData <- resData[order(resData$padj),c("id","fc","log2fc","pvalue","padj")]

# 保存到文件
write.table(resData, file = "output_degs.txt", row.names = FALSE, sep = "\t", quote = FALSE)

ngene <- 150
degs <- na.omit(resData[order(resData$log2fc,decreasing = T),])
updegs <- rownames(degs)[1:ngene]
dndegs <- rownames(degs)[(nrow(degs)-ngene + 1):nrow(degs)]
cmap.input <- data.frame(up = updegs,
                         dn = dndegs,
                         stringsAsFactors = F)
write.table(cmap.input,"CMap_input.txt",sep = "\t",row.names = F,col.names = T,quote = F)

#开始绘图
cmap.res <- read.delim("/home/aim/Desktop/TCGA_work/vip33new/Sdeath_project/result/fig213/CMap_export.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

drug <- c("sunitinib","NU-1025","4,5-dianilinophthalimide","amisulpride","butein","arachidonamide")  ##需要自己选定药物 可以选Cmap score最高的五个drug
selres <- cmap.res[which(cmap.res$Name %in% intersect(cmap.res$Name,drug)),]
print(selres) 

## 左侧区块
# 背景颜色
dt1 <- matrix(c(0,1, # 1为深色，0为浅色
                0,1,
                1,1,
                0,0,
                0,0,
                1,1),
              ncol = 2,
              byrow = T,
              dimnames = list(c("sunitinib","NU-1025","4,5-dianilinophthalimide","amisulpride","butein","arachidonamide"),
                              c("Clinical status","Experimental evidence")))

# 文字标签
lb1 <- matrix(c("Phase2","Present",
                "Phase1","Present",
                "Launched","Present",
                "Preclinical","Absent",
                "Preclinical","Absent",
                "Launched","Present"),
              ncol = 2,
              byrow = T)

# 文字颜色
cl1 <- matrix(c("black","white",
                "black","white",
                "white","white",
                "black","black",
                "black","black",
                "white","white"),
              ncol = 2,
              byrow = T)

# 画图
# 如果这里报错，请返回开头，看“环境设置”
hm1 <- pheatmap(mat = dt1,
                color = c("#EFFAE8","#BAE5BC"),
                cluster_cols = F,
                cluster_rows = F,
                show_rownames = T,
                show_colnames = T,
                display_numbers = lb1,
                number_color = cl1,
                fontsize_number = 11,
                border_color = "white",
                cellwidth = 50,
                cellheight = 50,
                legend = F)


## 中间区块
# 背景颜色
# dt2 <- matrix(c(1,0,
#                 1,0,
#                 0,0,
#                 1,0,
#                 0,0,
#                 1,0),
#               ncol = 2,
#               byrow = T,
#               dimnames = list(c("BI-2536","Leptomycin B","Methotrexate","Narciclasine","SR-II-138A","Vincristine"),
#                               c("Log2FC of mRNA expression","Log2FC of protein expression")))
# 
# # 文字标签
# lb2 <- matrix(c("1.23","0.42",
#                 "1.21","0.43",
#                 "0.89","-0.15",
#                 "1.07","0.24",
#                 "0.59","0.06",
#                 "1.54","-0.01"),
#               ncol = 2,
#               byrow = T)
# # 文字颜色
# cl2 <- matrix(c("white","black",
#                 "white","black",
#                 "black","black",
#                 "white","black",
#                 "black","black",
#                 "white","black"),
#               ncol = 2,
#               byrow = T)
# # 画图
# hm2 <- pheatmap(mat = dt2,
#                 color = c("#B7E4ED","#019AC9"),
#                 cluster_cols = F,
#                 cluster_rows = F,
#                 show_rownames = T,
#                 show_colnames = T,
#                 display_numbers = lb2,
#                 number_color = cl2,
#                 fontsize_number = 11,
#                 border_color = "white",
#                 cellwidth = 50,
#                 cellheight = 50,
#                 legend = F)

## 右侧区块
# 背景颜色
dt3 <- matrix(c(0, # 0为中间色
                -1, # -1为最浅色
                0,
                1, # 1为深色
                -1,
                0),
              ncol = 1,
              byrow = T,
              dimnames = list(c("sunitinib","NU-1025","4,5-dianilinophthalimide","amisulpride","butein","arachidonamide"),
                              c("CMap score")))

# 文字标签
# 根据selres的结果填写
# 因为有两个药物没有CMap结果，所以是NA
lb3 <- matrix(c("84.99",
                "-1.34",
                "-2.04",
                "-21.09",
                "-0.89",
                "-16.69"),
              ncol = 1,
              byrow = T)

# 文字颜色
cl3 <- matrix(c("black",
                "black",
                "black",
                "white",
                "black",
                "black"),
              ncol = 1,
              byrow = T)

# 画图
hm3 <- pheatmap(mat = dt3,
                color = c("white","#FDCEB9","#F26E5F"),
                cluster_cols = F,
                cluster_rows = F,
                show_rownames = T,
                show_colnames = T,
                display_numbers = lb3,
                number_color = cl3,
                fontsize_number = 11,
                border_color = "white",
                cellwidth = 50,
                cellheight = 50,
                legend = F)

## 水平合并热图
hm <- hm1  + hm3
draw(hm) 
dev.copy2pdf(file = "heatmap.pdf",width = 6,height = 6)


#fig212  ###这个代码很复杂需要  需要一份类似riskscore的分数  500个样本需要一天时间######
library(tidyverse) # 用于读取MAF文件
library(ISOpureR) # 用于纯化表达谱
library(impute) # 用于KNN填补药敏数据
library(pRRophetic) # 用于药敏预测
library(SimDesign) # 用于禁止药敏预测过程输出的信息
library(ggplot2) # 绘图
library(cowplot) # 合并图像

####fig282 CMAP plus####
setwd("/home/aim/Desktop/TCGA_work/vip33new/Sdeath_project/result/fig282/")
library(PharmacoGx)
library(parallel)
library(dplyr)
library(stringr)
library(tidyverse)
library(tibble)
library(clusterProfiler)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(gridExtra)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# 提取药物处理矩阵
load("CMAP_gene_signatures.RData")
camp_sig <- CMAP.genePerturbations[,,c("tstat")] %>% data.frame()

# 基因名转换
camp_sig$ENTREZID <- do.call(rbind, strsplit(rownames(camp_sig),'\\.'))[,2]

SYMBOL <- bitr(camp_sig$ENTREZID, fromType = "ENTREZID",
               toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

camp_sig <- merge(SYMBOL, camp_sig, by = "ENTREZID")
camp_sig <- column_to_rownames(camp_sig, var = "SYMBOL"); camp_sig <- camp_sig[,-1]

# 保存数据
saveRDS(camp_sig, "camp_sig.rds")


load("")
# 以下生成signature基因的代码只适用于TCGA来源的数据
#data <- read.table("LIHC_TPM.txt", sep="\t", check.names=F, header=T, row.names=1)

data <- tpms[,sub$Sample]
# 得到癌与癌旁比较的logfc
tumor <- tpms[, C123sam]
normal <- tpms[, C4sam]

dis_sig <- data.frame(id=rownames(data),
                      fc=log2(rowMeans(tumor)/rowMeans(normal)))
dis_sig <- dis_sig[dis_sig$fc != "Inf" & dis_sig$fc != "-Inf",]
dis_sig <- na.omit(dis_sig)

# 提取变化倍数最大的基因
# 例文（Yang et al）提到疾病分子特征数量选择100时可以获得较好的预测性能。但这篇研究是基于LINCS数据，cmap数据的维度与LINCS差别明显。这里我们建议疾病分子特征数量可以稍微多一点，以下演示使用top300基因进行XSum分析
dis_sig <- rbind(top_n(dis_sig, 150, fc), top_n(dis_sig, -150, fc))

# 将logfc转成二分类变量（logfc>0的基因用1表示，logfc小于0的基因用-1表示）
# 使用XSum时不需要考虑差异基因的差异倍数，这步分析是为了让大家更好的理解，并不是并要的

dis_sig$fc[dis_sig$fc>0] <- 1; dis_sig$fc[dis_sig$fc<0] <- -1
rownames(dis_sig) <- NULL

# 保存结果
write.table(dis_sig, "dis_sig.csv", sep=",", quote=F, row.names=F, col.names=T)

# dis_sig只要整理成以下格式都可以在后面的分析中使用
head(dis_sig)


# 读入drug signature
camp_sig <- readRDS("/home/aim/Desktop/xiaoyahuitu/FigureYa282CMAP_XSum/camp_sig.rds")

# 读入disease signature
dis_sig <- read.csv('dis_sig.csv', sep=',', header=TRUE)

# 读入XSum函数
source("/home/aim/Desktop/xiaoyahuitu/FigureYa282CMAP_XSum/Core_function.R")

# 选择XSum的topN（Yang et al的研究提到topN选择200效果可能比较好，但这个结论可能不适用与cmap的数据，这里我们选择topN = 500）
XLogFC <- eXtremeLogFC(camp_sig, N = 500)

up_gene <- dis_sig$id[dis_sig$fc == 1]
dn_gene <- dis_sig$id[dis_sig$fc == -1]

xsum <- data.frame(score=XSum(XLogFC, up_gene, dn_gene))
xsum <- rownames_to_column(xsum, var = "id")

# 把结果标准化至-1到1（这步也可不做）
xsum_pos <- xsum[xsum$score>0,]
xsum_pos$score <- xsum_pos$score/max(xsum_pos$score)

xsum_neg <- xsum[xsum$score<0,]
xsum_neg$score <- xsum_neg$score/min(xsum_neg$score) * -1

xsum <- rbind(xsum_pos, xsum_neg)

# 将结果从低到高排序
xsum <- xsum[order(xsum$score),]
head(xsum)

xsum$number <- 1:nrow(xsum)

# 突出显示top5的药物，标出药物名
select <- xsum[1:5,]

# 开始画图
ggplot(xsum, aes(number,score))+
  geom_point(size=3, color="grey50") + 
  geom_point(data = select, alpha = 1, 
             size = 5, color = "#5ec7dd") + 
  geom_label_repel(data = select, aes(label=id), 
                   color = "white",
                   alpha = 1, point.padding = 1, 
                   size = 5, fill = "#009bc7",
                   segment.size = 1, nudge_x=-0.5, 
                   segment.color = "grey50",
                   direction = "x",
                   hjust = 1) + 
  theme_classic()

ggsave("CMAP_XSum.pdf", width = 5, height = 7)




###validation of MOVICS ####


##JAPAN clinical 热图#####
setwd("~/Desktop/TCGA_work/vip33new/Sdeath_project/result/JAPAN_cli/")
head(japan_clin)
table(japan_clin$metastases)

###ICGC 热图clincal#####
setwd("~/Desktop/TCGA_work/vip33new/Sdeath_project/result/ICGC_cli/")




####突变分析  fig19  fig #####



