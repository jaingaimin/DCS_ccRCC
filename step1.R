


load("~/backup/t8a/allimage_data/Sdeath202302.rds")

.libPaths(c("~/R/forMOVICS"))
.libPaths()
.libPaths(c(#"/home/data/aim/R/x86_64-pc-linux-gnu-library/aimold/Rlib",
  "/home/data/aim/R/x86_64-pc-linux-gnu-library/4.3",
  "/home/data/aim/R/x86_64-pc-linux-gnu-library/4.2",
  "/usr/local/lib/R/library"
))
###新的分析流程 
#更新绘图方式和配色
#颜色协调一致

###安排一套最完美代码
# .libPaths(c("/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.1",
#             "/home/data/refdir/Rlib/",
#             "/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.0",
#             #"/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.2/",
#             "/usr/local/lib/R/library"))
# .libPaths(c("/home/data/refdir/Rlib/"))
library(RColorBrewer)
library(circlize)
library(gplots)
library(oompaBase)
library(viridis)
library(ggplot2)
library(clusterProfiler)
library(survival)
library(reshape2)
library(corrplot)
library(plyr)
library(igraph)
library(CMScaller)
#library(MOVICS,lib.loc = "~/R/forMOVICS/")
library(tidyverse)
library(MOVICS)
library(ggpubr)
library(maftools)
library(survminer)
library(pROC)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(tidyverse)
library(stringr)
library(EnhancedVolcano)
library(data.table)
library(monocle)
library(DDRTree)
library(ggsci)
library(ggstatsplot)
library(viridis)
library(scales)
library(RTN)
library(gplots)
library(grid)
###数据读入
load("~/backup/t8a/2022backup/xuzijun/vip39/TCGA/Copper_death/Copper_all.rds")

##保存需要的变量即可
gene_list <- c("NUBPL","NDUFA11","LRPPRC","OXSM","NDUFS1","GYS1","SLC7A11","SLC3A2","RPN1","NCKAP1")

save()
load("~/backup/t8a/vip39database/database/tcga_counts_fpkm_tpm/TCGA-KIRC_tpm_gene_symbol.Rdata")
tpms <- as.data.frame(tpms)
tpms[1:4,1:4]
tpms <- as.data.frame(log2(tpms+1))
range(tpms)
colnames(tpms) <- str_sub(colnames(tpms),1,15)
##HR 单因素多因素 多个KM in RCC###
Coxoutput <- data.frame(OS=sub$EVENT,
                        OS.time=sub$OS.time)

rownames(Coxoutput) <- sub$Sample
Coxoutput2 <- as.data.frame(t(tpms[unique(gene_list),sub$Sample]))
Coxoutput <- cbind(Coxoutput,Coxoutput2)
Coxoutput[1:4,1:4]

realdata <- Coxoutput
realdata[1:3,1:6]

setwd("~/Desktop/TCGA_work/vip33new/Sdeath_project/")
dir.create("COX")
setwd("./COX/")
Coxoutput=data.frame()
for(i in colnames(realdata[,3:ncol(realdata)])){
  cox <- coxph(Surv(OS.time, OS) ~ realdata[,i], data = realdata)
  coxSummary = summary(cox)
  Coxoutput=rbind(Coxoutput,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                  z=coxSummary$coefficients[,"z"],
                                  pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                  lower=coxSummary$conf.int[,3],
                                  upper=coxSummary$conf.int[,4]))
}
for(i in c(2:6)){
  Coxoutput[,i] <- as.numeric(as.vector(Coxoutput[,i]))
}

# Coxoutput <- arrange(Coxoutput,pvalue)  %>% #按照p值排序
#   filter(pvalue < 0.01) 
#Coxoutput <- arrange(Coxoutput,pvalue)
#Coxoutput <- Coxoutput[Coxoutput$pvalue<0.01,]
Coxoutput[1:4,1:4]


input_data <- Coxoutput2
input_data[1:4,1:4]
corr <- cor(input_data, method = "spearman")
corrplot(corr,title = "", 
         method = "pie", #或"circle" (default), "square", "ellipse", "number", "pie", "shade" and "color"
         outline = T, addgrid.col = "darkgray", 
         order="hclust", addrect = 4, #hclust聚为4类，根据数据的具体情况调整
         mar = c(4,0,4,0), #撑大画布，让细胞名显示完全
         rect.col = "black", rect.lwd = 5, cl.pos = "b", 
         tl.col = "black", tl.cex = 1.08, cl.cex = 1.5, tl.srt=60)
cor.mtest <- function(corr, ...) {
  corr <- as.matrix(corr)
  n <- ncol(corr)
  p.corr <- matrix(NA, n, n)
  diag(p.corr) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(corr[, i],method = "spearman", corr[, j], ...)
      p.corr[i, j] <- p.corr[j, i] <- tmp$p.value
    }
  }
  colnames(p.corr) <- rownames(p.corr) <- colnames(corr)
  p.corr
}

p.corr <- cor.mtest(input_data) 
head(p.corr[, 1:5])

#合并相关系数和P值
rr <- as.data.frame(corr)
rr$ID <- rownames(rr)
cor <- melt(rr,"ID",value.name = "cor"); #head(cor)

pp <- as.data.frame(p.corr)
pp$ID <- rownames(pp)
pvalue <- melt(pp,"ID",value.name = "pvalue"); #head(pvalue)
colnames(pvalue) <- c("from","to","pvalue")
corpvlue <- cbind(pvalue, cor)
head(corpvlue)

corpvlue <- corpvlue[, -c(4:5)]
head(corpvlue)
dim(corpvlue)

#去掉相关性较弱的连接
corpvlue <- corpvlue[corpvlue$pvalue < 0.0001,] #只保留pvalue < 0.0001的
dim(corpvlue)
corpvlue$weight <- corpvlue$pvalue
corpvlue$weight <- -log10(corpvlue$weight)
head(corpvlue)
#去掉相关系数为1，也就是两个相同变量之间的连接
corpvlue <- corpvlue[!corpvlue$cor==1,]
dim(corpvlue)
#去掉相关系数一样的连接--也就是重复计算的连接
summary(duplicated(corpvlue$weight))
corpvlue <- corpvlue[!duplicated(corpvlue$weight),]
dim(corpvlue)
#相关系数的正负用不同颜色表示
corpvlue$color <- ifelse(corpvlue$cor<0, negcol, poscol)

#保存到文件，便于查看
write.csv(corpvlue, "output_links.csv")
##进行聚类 细胞聚类 或者gene 聚类  可不运行
# cellcluster <- as.data.frame(t(input_data))
# #cellcluster[1:5,1:5]
# hc <- hclust(dist((cellcluster)))
# hcd <- as.dendrogram(hc)
# (clus4 <- cutree(hc, 4)) #分4类
# A <- as.character(rownames(as.data.frame(subset(clus4,clus4==1))))
# B <- as.character(rownames(as.data.frame(subset(clus4,clus4==2))))
# C <- as.character(rownames(as.data.frame(subset(clus4,clus4==3))))
# D <- as.character(rownames(as.data.frame(subset(clus4,clus4==4))))
# cls <- list(A,B,C,D)
# 
# nodes <- as.data.frame(unlist(cls))
# nodes$type <- c(rep("B",9),rep("A",4),rep("C",5),rep("D",5))
# names(nodes) <- c("media","type.label")
# 
# #以hclust的结果为基础，调整部分细胞所属的cluster
# nodes$type.label[nodes$media=="T cells follicular helper"] <- "B"
# nodes$type.label[nodes$media=="B cells naive"] <- "A"
# nodes$type.label[nodes$media=="T cells CD4 naive"] <- "A"
# nodes$type.label[nodes$media=="Plasma cells"] <- "A"
# nodes$type.label[nodes$media=="Dendritic cells resting"] <- "C"
# nodes$type.label[nodes$media=="Eosinophils"] <- "C"
# nodes$type.label[nodes$media=="Mast cells resting"] <- "A"
# 
# nodes <- as.data.frame(nodes)
# nodes$media <- as.character(nodes$media)
# nodes

nodes <- fread("/home/aim/Desktop/TCGA_work/vip33new/input/media_cor.csv",header = T,data.table = F)
nodes

# 合并生存分析的数据和细胞分类的数据
bb <- Coxoutput
bb$ID <- Coxoutput$gene

#bb$Cell.types <- as.character(bb$Cell.types) 
#colnames(bb)[1] <- c("ID")
#用pvalue控制节点圆的大小
bb$weight <- abs(log10(bb$pvalue))
#用HR标圆心点的颜色
bb$weight_HR <- (as.numeric(bb$HR)-1)*100
bb$colr <- ifelse(bb$weight_HR<0, "green", "black")
head(bb)

summary(nodes$media %in% Coxoutput$gene)

summary(nodes$media %in% bb$ID) #检查细胞名是否一致

nodes <- merge(nodes, bb, by.x = "media", "ID", all.x = T, all.y = T) #按细胞名merge

nodes$Fraction <- abs(nodes$weight_HR)

nodes$id <- paste("S", 01:10, sep = "")

nodes <- nodes[order(nodes$type.label),]
nodes <- nodes[,c(ncol(nodes),1:ncol(nodes)-1)]
nodes <- nodes[order(nodes$type.label),]
nodes

#建立nodes和links的连接id，把细胞名换成ID
paste0("'",nodes$media,"'","=","'",nodes$id,"'",collapse = ",")
corpvlue$from <- revalue(corpvlue$from,c('NCKAP1'='S3','RPN1'='S8','SLC3A2'='S9',
                                         'SLC7A11'='S10','GYS1'='S1','LRPPRC'='S2',
                                         'NDUFA11'='S4','NDUFS1'='S5','NUBPL'='S6','OXSM'='S7'))

corpvlue$to <- revalue(corpvlue$to,c('NCKAP1'='S3','RPN1'='S8','SLC3A2'='S9',
                                     'SLC7A11'='S10','GYS1'='S1','LRPPRC'='S2',
                                     'NDUFA11'='S4','NDUFS1'='S5','NUBPL'='S6','OXSM'='S7'))
(links <- corpvlue)
View(links)
links[24,6]
dim(links)
#利用nodes和links构建网络的input文件
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
#table(links$from%in%nodes$media)
# Generate colors based on cell clusters:
V(net)$color <- revalue(nodes$type.label,c("A"=mycol[1],"B"=mycol[2]))
# Compute node degrees (#links) and use that to set node size:
# Set edge width based on weight-log10(p_value):
V(net)$size <- (1 + V(net)$weight)*3 #节点圆的大小，可根据自己的数据再做调整
V(net)$label <- V(net)$media #设置标签
E(net)$arrow.mode <- 0 #不需要箭头
E(net)$edge.color <- "tomato" # tomato gray80
E(net)$width <- 1+E(net)$weight/6  #连接之间权重

pdf("Immune_network.pdf", width = 9.75, height = 8.78 )
plot(net,
     layout=layout_in_circle, #按圆圈布局
     edge.curved=.2, #画弯曲的连线
     vertex.label.color=V(net)$color, #细胞名的颜色
     vertex.label.dist= -2, #标签和节点的位置错开，后期还是要用AI调整
     edge.color=links$color)

#cluster的图例
legend("topright", #图例的位置
       c("Cell cluster-A", "Cell cluster-B"),
       pch=21, col="black", pt.bg=mycol, pt.cex=3,
       cex=1.3, bty="n", ncol=1)

#节点圆大小的图例，参考了FigureYa75base_volcano
f <- c(0.05, 0.001, 0.00001, 0.00000001)
s <- sqrt(abs(log10(f)))*3
legend("bottomright", 
       inset=c(0,-.1), #向下移
       legend=f, text.width = .2, 
       title = "logrank test, P value", title.adj = -.5,
       pch=21, pt.cex=s, bty='n',
       horiz = TRUE, #横向排列
       col = "black")

#连线的图例
legend("bottomright",
       c("Positive correlation with P < 0.0001", 
         "Negative correlation with P < 0.0001"),
       col = c(poscol, negcol), bty="n", 
       cex = 1, lty = 1, lwd = 5)

dev.off()


###cluster####
setwd("~/Desktop/TCGA_work/vip33new/Sdeath_project/")
##提取出T1到T3
table(sub$pstage)
# T1  T2  T3  T4 
# 264  68 173  11 
expr <- tpms[gene_list,sub$Sample]
dim(expr)
setwd("./result/")
dir.create('ConsensusCluster/')
results = ConsensusClusterPlus(as.matrix(expr),
                               maxK=9,
                               reps=100,
                               pItem=0.8,
                               pFeature=1,
                               title='ConsensusCluster/',
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123456,
                               plot="png")
Kvec = 2:9
x1 = 0.1; x2 = 0.9 
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") 
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK
#[1] 4
PAC <- as.data.frame(PAC)
PAC$K <- 2:9
library(ggplot2)
ggplot(PAC,aes(factor(K),PAC,group=1))+
  geom_line()+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_point(size=4,shape=21,color='darkred',fill='skyblue')+
  ylab('Proportion of ambiguous clustering')+
  xlab('Cluster number K')
library(export)
library(eoffice)
ggsave(filename = "ConsensusCluster/PAC.pdf",width = 6,height = 5)

icl <- calcICL(results,title = 'ConsensusCluster/',plot = 'png')
clusterNum=4    
cluster=results[[clusterNum]][["consensusClass"]]

sub <- data.frame(Sample=names(cluster),Cluster=cluster)
sub$Cluster <- paste0('C',sub$Cluster)
table(sub$Cluster)
# C1  C2  C3  C4 
# 92 167  74 183 
head(sub)
rownames(sub) <- sub$Sample
sub$OS <- as.vector(allclin[rownames(sub),"OS"])
sub$OStime <- as.vector(allclin[rownames(sub),"OS.time"])
sub <- sub[sub$OS>30,]
class(sub)

meta <- sub[,1:4]
head(sub)
colnames(meta) <- c("ID","cluster","event","time")
meta <- meta[,-1]
meta$time <- meta$time/365
head(meta)
sfit <- survfit(Surv(time, event)~cluster, data=meta)
ggsurvplot(sfit, conf.int=F, pval=TRUE)

sfit <- survfit(Surv(OS,EVENT) ~ Cluster,data = sub)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
sfit <- survfit(Surv(OS,EVENT)~Cluster,data = sub)
sfit <- survfit(Surv(OS.time,EVENT)~Cluster,data = sub)

sub$PFI <- allclin[rownames(sub),"PFI"]
sub$PFI.time <- allclin[rownames(sub),"PFI.time"]
sfit <- survfit(Surv(PFI.time, PFI)~Cluster, data=sub)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
ggsurvplot(sfit,
           #data = sub,
           palette= c(pal_nejm()(4),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("C1","C2","C3","C4"), 
           legend.title="cluster",
           xlab="Time (years)",
           #ylab='Overall survival',
           ylab='Overall survival',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)


###PFI
meta$PFI <- allclin[rownames(meta),"PFI"]
meta$PFI.time <- allclin[rownames(meta),"PFI.time"]/365
sfit2 <- survfit(Surv(PFI.time, PFI)~cluster, data=meta)

ggsurvplot(sfit2,
           palette= c(pal_nejm()(4),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("C1","C2","C3","C4"), 
           legend.title="cluster",
           xlab="Time (years)",
           ylab='Progression Free Survival',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)

ggsave(filename = "PFI_twogroup.pdf",height = 6,width = 6)


#加上正常组织三组绘图
dir.create("complexheatmap")
setwd("/home/aim/Desktop/TCGA_work/vip33new/Sdeath_project/result/complexheatmap/")

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
###准备两个数据 mygene_data 和 Subtype 分组信息

normal_id <- colnames(KIRC_mRNA_fpkm)[str_sub(colnames(KIRC_mRNA_fpkm),14,15)=="11"]
c(normal_id,sub$Sample)
Subtype <- data.frame(Subtype=c(rep("N",length(normal_id)),sub$Cluster),
                      id=c(normal_id,rownames(sub)))
rownames(Subtype) <- Subtype$id
Subtype <- Subtype[-2]

tpms[1:4,1:4]
table(rownames(Subtype)%in%colnames(tpms))
table(gene_list%in%rownames(tpms))
mygene_data <- tpms[gene_list,rownames(Subtype)]
mygene_data[1:4,1:4]
#Subtype <- Subtype[com_sam,,drop = F]
head(Subtype)
table(Subtype$Subtype)
## 用前面的自定义函数计算组间统计差异
comprTab <- cross_subtype_compr(expr = mygene_data, # 或log2(mygene_data + 1)，如果用参数检验，请注意对数转化；若非参均可
                                subt = Subtype,
                                #two_sam_compr_method = "wilcox", # 两组"t.test", "wilcox"
                                multi_sam_compr_method = "kruskal", # 多组"anova", "kruskal"
                                res.path = ".")


# 用全部基因来画
n.show_top_gene <- nrow(mygene_data)

# 按分组排序
subt.order <- Subtype[order(Subtype$Subtype),,drop = F]
indata <- mygene_data[comprTab$gene[1:n.show_top_gene],rownames(subt.order)]

# 数据标准化和边界设置
plotdata <- t(scale(t(indata)))
plotdata[plotdata > 2] <- 2
plotdata[plotdata < -2] <- -2

# 调整行名
blank <- "    " # 行名和p值之间的间隔
p.value <- comprTab$adjusted.p.value[1:n.show_top_gene]
sig.label <- ifelse(p.value < 0.001,"****",
                    ifelse(p.value < 0.005,"***",
                           ifelse(p.value < 0.01,"**",
                                  ifelse(p.value < 0.05,"*",""))))
p.label <- formatC(p.value, # 将p值变成保留两位小数的科学计数法
                   format = "e",
                   digits = 2)

add.label <- str_pad(paste0(rownames(plotdata),sig.label), # 固定行名宽度并再右侧补齐" "
                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                     side = "right")

annCol <- subt.order # 获得排序后的亚型注释信息，这里只有一个变量需要注释
colnames(annCol)[1] <- paste(str_pad(colnames(annCol)[1], # 注释列名补上"P-value"，宽度和刚才一致
                                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                                     side = "right"),
                             "P-value",
                             sep = blank)

annColors <- list(c( "C1"="#BC3C29FF", "C2"="#0072B5FF","C3"="#E18727FF","C4"="#7876B1FF","N"="#33A02CFF")) # 如果有多个变量要注释颜色请补充c()
names(annColors) <- colnames(annCol)[1] # 如果有多个变量要注释颜色请补充每张list的name

# 绘制热图
table(Subtype$Subtype)
pheatmap(cellheight = 10, cellwidth = 1,
         mat = plotdata, # 输入数据
         scale = "none", # 不标准化因为数据已经被标准化
         annotation_col = annCol, # 列注释信息
         annotation_colors = annColors, # 列注释对应的颜色
         cluster_cols = F, # 列不聚类
         cluster_rows = F, # 行不聚类
         show_colnames = F, # 不显示列名
         show_rownames = T, # 显示基因名
         annotation_legend = F, # 不显示图例
         gaps_col = c(92,259,333,516),
         color = paletteer_c("scico::berlin", n = 100),
         labels_row = paste(add.label, p.label, sep=blank),
)
table(Subtype$Subtype)

####clinical heatmap+ table +subtype HR####
#View(rt)
rt <- sub[,c("Sample","Cluster")]
head(rt)
rt <- rt[-1]
colnames(rt) <- "risk"

rt=rt[order(rt$risk),,drop=F]
View(allclin)
rt$age <- allclin[rownames(rt),4]
rt$age <-ifelse(rt$age<=62,'<=65','>65')
rt$age=factor(rt$age,labels = c('<=65','>65'))
rt$gender <- allclin[rownames(rt),5]
table(rt$gender)
rt$gender <- factor(rt$gender,labels=c("MALE", "FEMALE"))

rt$T <- str_sub(pd[rownames(rt),"pathologic_T"],1,2)
table(rt$T)
rt$T <- factor(rt$T,labels=c("T1", "T2", "T3", "T4"))
rt$N <- str_sub(pd[rownames(rt),"pathologic_N"],1,2)
table(rt$N)
rt$N <- gsub("NX","unknown",rt$N)
rt$N <- factor(rt$N,labels=c("N0", "N1", "unknown"))
rt$M <- str_sub(pd[rownames(rt),"pathologic_M"],1,2)
table(rt$M)
rt$M <- gsub("MX","unknown",rt$M)
rt$M <- factor(rt$M,labels=c("M0", "M1","unknown"))

rt$stage <- allclin[rownames(rt),7]
table(rt$stage)
rt$stage <- gsub("\\[Discrepancy]","unknown",rt$stage)
table(rt$stage)
rt$stage <- factor(rt$stage,labels=c("Stage I", "Stage II", "Stage III", "Stage IV","unknown"))


myCol=c(
  '#98D5F4', '#3B7EF4', '#FEEDDE','#E6550D', '#B5B4B4', '#727270', '#444746', '#070707', '#00DBA0',
  '#15788E', '#C45132', '#EDF8E9', '#31A354' ,'#B5B4B4', '#727270', '#444746', '#070707')

ColorList=list(risk=c("C1"="#BC3C29FF", "C2"="#0072B5FF","C3"="#E18727FF","C4"="#7876B1FF"))

i=0
for(cli in colnames(rt[,2:ncol(rt)])){
  cliLength=length(levels(factor(rt[,cli])))
  cliCol=myCol[(i+1):(i+cliLength)]
  i=i+cliLength
  names(cliCol)=levels(factor(rt[,cli]))
  cliCol["unknown"]="grey75"
  ColorList[[cli]]=cliCol
}

library(ComplexHeatmap)
ha=HeatmapAnnotation(df=rt, col=ColorList)
zero_row_mat=matrix(nrow=0, ncol=nrow(rt))
Hm=Heatmap(zero_row_mat, top_annotation=ha,)

#输出热图
pdf(file="heatmap_clincal.pdf", width=7, height=3)
draw(Hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

###差异分析适合MOVICS分析 四组####

load("~/database/kirc.tcga5_omics.rds")
sub$pstage <- pd[sub$Sample,"pathologic_T"]
sub$pstage <- str_sub(sub$pstage,1,2)
head(sub)
newsub <- sub
newsub <- sub[c(1,2,3,6,7)]
newsub <- sub[c(1,2,3,4,7)]
head(newsub)
colnames(newsub) <- c("samID","cluster","fustat","futime","pstage")
head(newsub)
newsub$clust <-str_sub(newsub$cluster,2,2)
head(newsub)

id258 <- colnames(kirc.tcga5_omics$mRNA.expr)
table(id258%in%rownames(newsub))##只有247个
inter_sam <- id258[id258%in%rownames(newsub)]
inter_sub <- newsub[inter_sam,]

pseudo.moic.res<- list("clust.res" = inter_sub,"mo.method" = "PAM50")

# survival comparison
pam50.brca <- compSurv(moic.res         = pseudo.moic.res,
                       surv.info        = surv.info,
                       convt.time       = "y", # convert day unit to year
                       surv.median.line = "h", # draw horizontal line at median survival
                       fig.name         = "KAPLAN-MEIER CURVE OF PAM50 BY PSEUDO")

cmoic.brca <- kirc.tcga5_omics[1:5]
tmb.brca <- compTMB(moic.res     = pseudo.moic.res,
                    maf          = kirc.tcga5_omics$maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "DISTRIBUTION OF TMB AND TITV")

colnames(segment) <- c("sample","chrom","start","end","value")
head(segment)
segment <- kirc.tcga5_omics$segment
segment <- as.data.frame(segment)
segment <- segment[segment$sample%in%inter_sam,]
fga.brca <- compFGA(moic.res     = pseudo.moic.res,
                    segment      = segment,
                    iscopynumber = F, # this is a segmented copy number file
                    cnathreshold = 0.2, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "BARPLOT OF FGA")

maf <- kirc.tcga5_omics$mut.status[,colnames(kirc.tcga5_omics$mut.status)%in%inter_sam]
mut.brca <- compMut(moic.res     = pseudo.moic.res,
                    mut.matrix   = kirc.tcga5_omics$mut.status, # binary somatic mutation matrix
                    doWord       = TRUE, # generate table in .docx format
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
                    p.adj.cutoff = 0.05, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                    innerclust   = TRUE, # perform clustering within each subtype
                    #annCol       = annCol, # same annotation for heatmap
                    #annColors    = annColors, # same annotation color for heatmap
                    width        = 6, 
                    height       = 5,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")


drug.brca <- compDrugsen(moic.res    = pseudo.moic.res,
                         norm.expr   = KIRC_mRNA_fpkm[,pseudo.moic.res$clust.res$samID], # double guarantee sample order
                         drugs       = c("Sunitinib","Afatinib","Erlotinib","Gefitinib"), # a vector of names of drug in GDSC
                         tissueType  = "urogenital_system", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "BOXVIOLIN OF ESTIMATED IC50")


drug.brca <- compDrugsen(moic.res    = pseudo.moic.res,
                         norm.expr   = KIRC_mRNA_fpkm[,pseudo.moic.res$clust.res$samID], # double guarantee sample order
                         drugs       = c("Imatinib","Crizotinib","Saracatinib","Dasatinib","Lisitinib"), # a vector of names of drug in GDSC
                         tissueType  = "urogenital_system", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "BOXVIOLIN OF ESTIMATED  IC50 ")
###完善版本
drug.brca <- compDrugsen(moic.res    = pseudo.moic.res,
                         norm.expr   = KIRC_mRNA_fpkm[,pseudo.moic.res$clust.res$samID], # double guarantee sample order
                         drugs       = c("Saracatinib","Crizotinib","Axitinib","Erlotinib","Pazopanib","Temsirolimus"), # a vector of names of drug in GDSC
                         tissueType  = "urogenital_system", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "BOXVIOLIN OF ESTIMATED  IC50 ")

dir.create("MOVICS")
setwd("./MOVICS")
runDEA(dea.method = "edger",
       expr       = KIRC_mRNA_count, # raw count data
       moic.res   = pseudo.moic.res,
       prefix     = "TCGA-KIRC") # prefix of figure name
marker.up <- runMarker(moic.res      = pseudo.moic.res,
                       dea.method    = "edger", # name of DEA method
                       prefix        = "TCGA-KIRC", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = KIRC_mRNA_fpkm, # use normalized expression as heatmap input
                       #annCol        = annCol, # sample annotation in heatmap
                       #annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")
head(marker.up$templates)
marker.down <- runMarker(moic.res      = pseudo.moic.res,
                         dea.method    = "edger", # name of DEA method
                         prefix        = "TCGA-KIRC", # MUST be the same of argument in runDEA()
                         dat.path      = getwd(), # path of DEA files
                         res.path      = getwd(), # path to save marker files
                         p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                         p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                         dirct         = "down", # direction of dysregulation in expression
                         n.marker      = 100, # number of biomarkers for each subtype
                         doplot        = TRUE, # generate diagonal heatmap
                         norm.expr     = KIRC_mRNA_fpkm, # use normalized expression as heatmap input
                         #annCol        = annCol, # sample annotation in heatmap
                         #annColors     = annColors, # colors for sample annotation
                         show_rownames = FALSE, # show no rownames (biomarker name)
                         fig.name      = "DOWNREGULATED BIOMARKER HEATMAP")
head(marker.down$templates)

MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)

gsea.dn <- runGSEA(moic.res     = pseudo.moic.res,
                   dea.method   = "edger",
                   prefix       = "TCGA-KIRC",
                   msigdb.path  = MSIGDB.FILE,
                   norm.expr    = KIRC_mRNA_fpkm,
                   dirct        = "down",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method  = "ssgsea", # switch to ssgsea
                   norm.method  = "median", # switch to median
                   fig.name     = "DOWNREGULATED PATHWAY HEATMAP") 

gsea.up <- runGSEA(moic.res     = pseudo.moic.res,
                   dea.method   = "edger",
                   prefix       = "TCGA-KIRC",
                   msigdb.path  = MSIGDB.FILE,
                   norm.expr    = KIRC_mRNA_fpkm,
                   dirct        = "up",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method  = "ssgsea", # switch to ssgsea
                   norm.method  = "median", # switch to median
                   fig.name     = "UPREGULATED PATHWAY HEATMAP") 

CC <- ""
gsea.dn <- runGSEA(moic.res     = pseudo.moic.res,
                   dea.method   = "edger",
                   prefix       = "TCGA-KIRC",
                   msigdb.path  = CC,
                   norm.expr    = KIRC_mRNA_fpkm,
                   dirct        = "down",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method  = "ssgsea", # switch to ssgsea
                   norm.method  = "median", # switch to median
                   fig.name     = "DOWNREGULATED PATHWAY HEATMAP") 

gsea.up <- runGSEA(moic.res     = pseudo.moic.res,
                   dea.method   = "edger",
                   prefix       = "TCGA-KIRC",
                   msigdb.path  = CC,
                   norm.expr    = KIRC_mRNA_fpkm,
                   dirct        = "up",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method  = "ssgsea", # switch to ssgsea
                   norm.method  = "median", # switch to median
                   fig.name     = "UPREGULATED PATHWAY HEATMAP")

###Dasatinib 白血病
###Linsitinib  肾上腺肿瘤
##Saracatinib renal cancer
###Crizotinib 乳头状肾癌 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)00152-5/fulltext
##Imatinib 白血病
##Axitinib 肾癌
##Sunitinib  肾癌
##Afatinib 肺癌
###Erlotinib 可用于特殊肾癌
## Gefitinib  局部晚期或转移性非小细胞肺癌(NSCLC)
##Pazopanib 晚期肾癌
##Temsirolimus  肾癌

# convert beta value to M value for stronger signal
indata <- kirc.tcga5_omics[c(1,2,4,5)]
indata$meth.beta <- log2(indata$meth.beta / (1 - indata$meth.beta))

# data normalization for heatmap
feat   <- iClusterBayes.res$feat.res
feat1  <- feat[which(feat$dataset == "mRNA.expr"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA.expr"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "meth.beta"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "mut.status"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4)
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,F)) # no scale for mutation


mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, meth.col, mut.col)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = pseudo.moic.res$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")



##GSVA_免疫
GSET.FILE <- 
  system.file("extdata", "gene sets of interest.gmt", package = "MOVICS", mustWork = TRUE)
# run GSVA to estimate single sample enrichment score based on given gene set of interest


gsva.res <- 
  runGSVA(moic.res      = pseudo.moic.res,
          norm.expr     = KIRC_mRNA_fpkm,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva",  # method to calculate single sample enrichment score
          # annCol        = annCol,
          #annColors     = annColors,
          #color =magma(64),
          fig.path      = getwd(),
          show_rownames = T,
          fig.name      = "GENE SETS OF INTEREST HEATMAP",
          height        = 5,
          width         = 8)

###加上其他的数据集IOBR
library(IOBR)
IOBR::signature_metabolism
metabolism <- IOBR::signature_metabolism
GSET.FILE2 <- "/t8a/2022backup/xuzijun/vip39/TCGA/jiaowang_KIRC/met_model/data/try.gmt"
gsva.res <- 
  runGSVA(moic.res      = pseudo.moic.res,
          norm.expr     = KIRC_mRNA_fpkm,
          gset.gmt.path = GSET.FILE2, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          # annCol        = annCol,
          #annColors     = annColors,
          fig.path      = getwd(),
          show_rownames = T,
          fig.name      = "GENE SETS OF daixie HEATMAP",
          height        = 20,
          width         = 8)

get_gmt <- function(gmtinput,filename){
  output <- file(filename, open="wt")
  lapply(gmtinput[["name"]],function(name){
    outlines = paste0(c(name, gmtinput[["description"]][[name]],
                        gmtinput[["genes"]][[name]]),collapse='\t')
    writeLines(outlines, con=output)
  })
  close(output)
}
write.table(metabolism,file = "./data/Hallmark_my_geneset.gmt",sep = "\t",row.names = F,col.names = F,quote = F)
unlist(metabolism)

lst <- metabolism
library(GSEABase)
gsc <- GeneSetCollection(mapply(function(geneIds, keggId) {
  GeneSet(geneIds, geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
}, lst, names(lst)))
gsc
write.table(gsc,file = "./data/Hallmark_my_geneset.gmt")


sink("./data/try.gmt")
yourlist <- metabolism
for (i in 1:length(yourlist)){
  cat(names(yourlist)[i])
  cat('\tNA\t')
  cat(paste(yourlist[[i]], collapse = '\t'))
  cat('\n')
}
sink()

tme <- IOBR::signature_tme
sink("./data/tme.gmt")
yourlist <- tme
for (i in 1:length(yourlist)){
  cat(names(yourlist)[i])
  cat('\tNA\t')
  cat(paste(yourlist[[i]], collapse = '\t'))
  cat('\n')
}
sink()
GSET.FILE3 <- "/t8a/2022backup/xuzijun/vip39/TCGA/jiaowang_KIRC/met_model/data/tme.gmt"
gsva.res <- 
  runGSVA(moic.res      = pseudo.moic.res,
          norm.expr     = KIRC_mRNA_fpkm,
          gset.gmt.path = GSET.FILE3, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          # annCol        = annCol,
          #annColors     = annColors,
          fig.path      = getwd(),
          show_rownames = T,
          fig.name      = "GENE SETS OF tme HEATMAP",
          height        = 20,
          width         = 8)


tumor <- IOBR::signature_tumor
sink("./data/tumor.gmt")
yourlist <- tumor
for (i in 1:length(yourlist)){
  cat(names(yourlist)[i])
  cat('\tNA\t')
  cat(paste(yourlist[[i]], collapse = '\t'))
  cat('\n')
}
sink()
GSET.FILE4 <- "/t8a/2022backup/xuzijun/vip39/TCGA/jiaowang_KIRC/met_model/data/tumor.gmt"
gsva.res <- 
  runGSVA(moic.res      = pseudo.moic.res,
          norm.expr     = KIRC_mRNA_fpkm,
          gset.gmt.path = GSET.FILE4, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          # annCol        = annCol,
          #annColors     = annColors,
          fig.path      = getwd(),
          show_rownames = T,
          fig.name      = "GENE SETS OF tumor HEATMAP",
          height        = 6,
          width         = 8)


GSET.FILE5 <- "/t8a/2022backup/xuzijun/vip39/database/h.all.v7.5.symbols.gmt"
gsva.res <- 
  runGSVA(moic.res      = pseudo.moic.res,
          norm.expr     = KIRC_mRNA_fpkm,
          gset.gmt.path = GSET.FILE5, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          # annCol        = annCol,
          #annColors     = annColors,
          fig.path      = getwd(),
          show_rownames = T,
          fig.name      = "GENE SETS OF gsva HEATMAP",
          height        = 8,
          width         = 12)


rouliugenelist <- fread("/home/data/vip39/database/rouliu_genelist.csv",header = T,data.table = F)

rouliugenelist <- lapply(rouliugenelist, function(x) {
  unique(na.omit(x)) 
})
length(rouliugenelist)

sink("./data/rouliu.gmt")
yourlist <- rouliugenelist
for (i in 1:length(yourlist)){
  cat(names(yourlist)[i])
  cat('\tNA\t')
  cat(paste(yourlist[[i]], collapse = '\t'))
  cat('\n')
}
sink()
GSET.FILE6 <- "/home/data/vip39/TCGA/WWOX/data/rouliu.gmt"

gsva.res <- 
  runGSVA(moic.res      = pseudo.moic.res,
          norm.expr     = KIRC_mRNA_fpkm,
          gset.gmt.path = GSET.FILE6, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          # annCol        = annCol,
          #annColors     = annColors,
          fig.path      = getwd(),
          show_rownames = T,
          fig.name      = "GENE SETS OF rouliu HEATMAP",
          height        = 10,
          width         = 12)

####ceRNA####
