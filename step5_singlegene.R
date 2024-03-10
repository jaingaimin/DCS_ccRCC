


###单gene
####forest选出gene####
#使用
setwd("/home/aim/Desktop/TCGA_work/vip33new/Sdeath_project/result/prognosis_model")
library(survival)
library(randomForestSRC)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN) == 0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}
# 去除无表达的基因，log变换，z-score
head(expr)
dim(expr)
expr <- tpms[unique(gene_list),sub$Sample]
##
# expr <- as.data.frame(round(t(scale(t(log2(expr + 1)))),3))
# expr <- expr[,rownames(sub)]
dim(expr)
colnames(sub)
str(sub)
head(sub)
surv <- sub[,c("OStime","OS")]
colnames(surv) <- c("OS.time","OS")
head(surv)
head(sub)
cox.pcutoff <- 0.05 # cox的p阈值
Coxoutput.OS <- NULL
for (i in 1:nrow(expr)) {
  display.progress(index = i,totalN = nrow(expr)) # 显示进度
  
  # 产生临时变量存储生存以及变量表达值
  tmp <- data.frame(gene = as.numeric(expr[i,]),
                    OS.time = surv[,"OS.time"],
                    OS = surv[,"OS"],
                    stringsAsFactors = F)
  
  # 单变量cox比例风险模型
  cox <- coxph(Surv(OS.time, OS) ~ gene, data = tmp)
  coxSummary = summary(cox)
  
  # 生成cox结果数据框，包括基因名，风险比，z值，waldtest p值，以及HR置信区间
  Coxoutput.OS=rbind.data.frame(Coxoutput.OS,data.frame(gene=rownames(expr)[i],
                                                        HR=as.numeric(coxSummary$coefficients[,"exp(coef)"]),
                                                        z=as.numeric(coxSummary$coefficients[,"z"]),
                                                        pvalue=as.numeric(coxSummary$coefficients[,"Pr(>|z|)"]),
                                                        lower=as.numeric(coxSummary$conf.int[,3]),
                                                        upper=as.numeric(coxSummary$conf.int[,4]),
                                                        stringsAsFactors = F),
                                stringsAsFactors = F)
}
head(Coxoutput.OS)
write.csv(Coxoutput.OS,"univariate cox regression for gene filtering.csv",row.names = F,quote = F)

##随机森林进一步降维
gene.sel <- Coxoutput.OS[which(Coxoutput.OS$pvalue < cox.pcutoff),"gene"]
tmp <- expr[gene.sel,]

rownames(tmp) <- gsub("-","_",rownames(tmp)) # 防止出现“-”导致程序报错
dt.rf <- cbind.data.frame(surv[,c("OS","OS.time")],t(tmp))

ntree <- 1000
surv.rf <- rfsrc(Surv(OS.time, OS) ~ ., 
                 data = dt.rf, 
                 ntree = ntree,
                 importance = TRUE,
                 seed = 12345678)

#排列组化确定最优签名
num.imp <- 5
rel.imp <- sort(surv.rf$importance, decreasing = T)
rel.imp.sel <- rel.imp[1:num.imp] # 取出一定数量的基因
names(rel.imp.sel) <- gsub("_","-",names(rel.imp.sel)) # 还原基因名

outTab <- NULL
n.sum <- 0
for (i in 1:num.imp) {
  cat(paste0("combination using ",i," genes...\n"))
  tmp <- utils::combn(names(rel.imp.sel), m=i) # 获取当前基因个数下的排列组合
  n <- ncol(tmp)
  for (j in 1:n) {
    combgene <- tmp[,j] # 取出每一次组合的基因名
    combexpr <- cbind.data.frame(t(expr[combgene,]), # 构建数据库做多变量cox
                                 OS.time = surv[,"OS.time"],
                                 OS = surv[,"OS"],
                                 stringsAsFactors = F)
    cox <- coxph(Surv(OS.time, OS) ~ ., data = combexpr)
    coxSummary <- summary(cox)
    coeff <- coxSummary$coefficients[,1] # 取出系数
    riskscore <- as.matrix(combexpr[,combgene]) %*% coeff # 计算riskscore
    riskscore <- data.frame(riskscore = as.numeric(riskscore[,1]),
                            group = ifelse(riskscore[,1] > median(riskscore[,1]),"HRisk","LRisk"), # 根据中位数分组
                            row.names = rownames(riskscore),
                            OS.time = combexpr$OS.time,
                            OS = combexpr$OS,
                            stringsAsFactors = F)
    fitd <- survdiff(Surv(OS.time, OS) ~ group,
                     data = riskscore,
                     na.action = na.exclude)
    p.val <- 1-pchisq(fitd$chisq, length(fitd$n) - 1) # log-rank检验
    
    outTab <- rbind.data.frame(outTab,
                               data.frame(num.gene = ifelse(i == 1, paste0(i," gene"), paste0(i," genes")), # 当前基因数目
                                          km.pvalue = p.val, # KM曲线p值
                                          core.gene = paste(combgene,collapse = " | "), # 该组合下的基因
                                          stringsAsFactors = F),
                               stringsAsFactors = F)
  }
  n.sum <- n + n.sum # 校验排列组合的总数
}
if(n.sum == 2^num.imp-1) { # 如果总和不等则报错
  write.csv(outTab,"combination of important genes with KM pvalues.csv",row.names = F,quote = F)
} else (message("Wrong combination!!!"))


sigpoints <- Coxoutput.OS[which(Coxoutput.OS$pvalue < cox.pcutoff),]
unsigpoints <- Coxoutput.OS[which(Coxoutput.OS$pvalue >= cox.pcutoff),]

pdf("volcano.pdf",width = 5,height = 5)
par(bty = "o", mgp = c(2,.6,0), mar = c(3,3,1,1), las = 1, font.axis = 1) # 基础参数
plot(log(Coxoutput.OS$HR),
     -log10(Coxoutput.OS$pvalue),
     xlab = "Univariate Cox coefficient",
     ylab = bquote("-log"[10]~"(P value)"),
     xlim = c(-2,2))
points(log(sigpoints$HR),
       -log10(sigpoints$pvalue),
       col = ggplot2::alpha("#E53435",0.8),
       pch = 19)
points(log(unsigpoints$HR),
       -log10(unsigpoints$pvalue),
       col = ggplot2::alpha("#21498D",0.8),
       pch = 19)
abline(h = -log10(cox.pcutoff), lty = 2, col = "grey60")
invisible(dev.off())


xrange <- range(pretty(range(rel.imp.sel))) # 根据重要性区间确定x轴范围
yrange <- c(1,length(rel.imp.sel))  # 根据重要变量个数确定y轴范围

pdf("variable importance.pdf",width = 5,height = 5)
par(bty = "o", mgp = c(1.5,.33,0), mar = c(3,7,1,2),las = 1, tcl = -.25)
plot(NULL,NULL,
     xlim = xrange,
     ylim = yrange,
     xlab = "Variable Importance",
     ylab = "",
     yaxt = "n",
     las = 1)
axis(side = 2,at = 1:length(rel.imp.sel),rev(names(rel.imp.sel))) # 补齐y轴
for (i in 1:length(rel.imp.sel)) { # 循环添加线
  lines(c(xrange[1],rev(rel.imp.sel)[i]),
        c(i,i),
        lwd = 2.5,
        col = "#33A02CFF")
}
invisible(dev.off())


## 绘制排列组合p值


num.comb <- 20
outTab2 <- outTab[order(outTab$km.pvalue),][1:num.comb,]
head(outTab2)
#outTab2$km.pvalue[1] <- 1.100223e-18
pdf("combination barplot.pdf",width = 5,height = 5)
par(bty = "o", mgp = c(1.5,.33,0), mar = c(1,4,3,1),las = 1, tcl = -.25)
barplot(rev(-log10(outTab2$km.pvalue)),
        horiz = T, # 柱状图横向
        names.arg = rev(outTab2$num.gene), # 添加y轴名称
        xaxt = "n", # 取消下方x轴
        col = "#FF7F00FF")
axis(side = 3) # 在上方添加x轴
mtext(side = 3, bquote("-log"[10]~"(P value)"), line = 1) # 添加x轴名称
invisible(dev.off())
outTab2$core.gene[1]



####ntp外部验证#####
setwd("~/Desktop/TCGA_work/vip33new/Sdeath_project/result/complexheatmap/MOVICS/")


head(sub)
sub$pstage <- pd[sub$Sample,"pathologic_T"]
sub$pstage <- str_sub(sub$pstage,1,2)
head(sub)
newsub <- sub[,c(1,2,3,9,7)]
colnames(newsub) <- c("samID","cluster","fustat","futime","pstage")

newsub$clust <- ifelse(newsub$cluster=="C1","1","2")
head(newsub)

id258 <- colnames(kirc.tcga5_omics$mRNA.expr)
table(id258%in%rownames(newsub))##只有247个
inter_sam <- id258[id258%in%rownames(newsub)]
inter_sub <- newsub[inter_sam,]

pseudo.moic.res<- list("clust.res" = inter_sub,"mo.method" = "PAM50")
pseudo.moic.res<- list("clust.res" = newsub,"mo.method" = "PAM50")
dir.create("NTP")
setwd("./NTP/")
valid91_id <- rownames(valid91)

count <- round(KIRC_mRNA_count)
runDEA(dea.method = "deseq2",
       expr       = count,
       moic.res   = pseudo.moic.res,
       prefix     = "TCGA-KIRC")
# choose edgeR result to identify subtype-specific up-regulated biomarkers
marker.up <- runMarker(moic.res      = pseudo.moic.res,
                       dea.method    = "deseq2", # name of DEA method
                       prefix        = "TCGA-KIRC", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = KIRC_mRNA_fpkm, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")
marker.down <- runMarker(moic.res      = pseudo.moic.res,
                         dea.method    = "deseq2", # name of DEA method
                         prefix        = "TCGA-KIRC", # MUST be the same of argument in runDEA()
                         dat.path      = getwd(), # path of DEA files
                         res.path      = getwd(), # path to save marker files
                         p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                         p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                         dirct         = "down", # direction of dysregulation in expression
                         n.marker      = 100, # number of biomarkers for each subtype
                         doplot        = TRUE, # generate diagonal heatmap
                         norm.expr     = KIRC_mRNA_fpkm, # use normalized expression as heatmap input
                         annCol        = annCol, # sample annotation in heatmap
                         annColors     = annColors, # colors for sample annotation
                         show_rownames = FALSE, # show no rownames (biomarker name)
                         fig.name      = "UPREGULATED BIOMARKER HEATMAP")
#> --all samples matched.
#> --log2 transformation done for expression data.
###运行maker基因
load(system.file("extdata", "brca.yau.RData",  package = "MOVICS", mustWork = TRUE))

clin.info=clin_JW[valid91_id,c("EVENT","OS")]
colnames(clin.info) <- c("fustat","futime")
clin.pfi=allclin[valid91_id,c("PFI","PFI.time")]
colnames(clin.pfi) <- c("fustat","futime")
japan_clin_need <- sub5
colnames(japan_clin_need) <- c("ID","fustat","futime","lasso","risk_group")
head(japan_clin_need)
japan_clin_need$futime <- as.numeric(japan_clin_need$futime*30)
kirc.japan <- list(mRNA.exp=japan_mRNA_fpkm[,rownames(japan_clin_need)],
                   clin.info=japan_clin_need)




yau.ntp.pred <- runNTP(expr       = kirc.japan$mRNA.exp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR japan") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.japan$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR japan") 




####cancercell验证####
load("~/backup/t8a/2022backup/xuzijun/vip39/database/KIRC_immune_therapy/cancercellpaper/cancercellRCC.rds")

load("/home/data/aim/TCGA/RCCdata/cancercellRCC.rds")
Cancercell_cli$futime <- as.numeric(Cancercell_cli$PFS_MONTHS)*30
Cancercell_cli$fustat <- ifelse(Cancercell_cli$PFS_CENSOR=="TRUE",1,0)

kirc.Cancercell <- list(mRNA.exp=Cancercell_tpm,
                        clin.info=Cancercell_cli)


yau.ntp.pred <- runNTP(expr       = kirc.Cancercell$mRNA.exp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR Cancercell") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.Cancercell$clin.info,
                     #clust.col = c(pal_nejm()(2),'grey60'),
                     
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR Cancercell2") 

Cancercell_cluster <- yau.ntp.pred$clust.res
Cancercell_cluster$clust <- paste0("DCS",Cancercell_cluster$clust)
write_csv(Cancercell_cluster,file = "~/TCGA/Sdeath_project/supp2024/Cancercell_cluster.csv")

##复旦队列验证####
load("/t8a/2022backup/xuzijun/vip39/database/KIRC_immune_therapy/Fudan_RCC/FudanRCCall.rds")
load("/home/data/aim/TCGA/RCCdata/FudanRCCall.rds")

group <- str_detect(colnames(FDRCC_exp2),"TA")

FDRCC_exptumor <- FDRCC_exp2[,!str_detect(colnames(FDRCC_exp2),"TA")]
rownames(FDRCC_clin) <- paste0(rownames(FDRCC_clin),"_T")
table(paste0(rownames(FDRCC_clin),"_T")%in%colnames(FDRCC_exptumor))
interid_FD <- intersect(rownames(FDRCC_clin),colnames(FDRCC_exptumor))
FDRCC_exp <- FDRCC_exptumor[,interid_FD]
FDRCC_clin <- FDRCC_clin[interid_FD,]

FDRCC_clin$futime <- as.numeric(FDRCC_clin$`OS(month)`)*30
FDRCC_clin$fustat <- ifelse(FDRCC_clin$`Live Status`=="living",0,1)

kirc.FDRCC <- list(mRNA.exp=FDRCC_exp,
                   clin.info=FDRCC_clin)


yau.ntp.pred <- runNTP(expr       = kirc.FDRCC$mRNA.exp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR FDRCC") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.FDRCC$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR FDRCC") 

FD_cluster <- yau.ntp.pred$clust.res
FD_cluster$clust <- paste0("DCS",FD_cluster$clust)
write_csv(FD_cluster,file = "~/TCGA/Sdeath_project/supp2024/FD_cluster.csv")



###checkmate验证####
load("/t8a/2022backup/xuzijun/vip39/database/KIRC_immune_therapy/PD1/CheckmateRCC.rds")
#os
Checkmate_clin$futime <- as.numeric(Checkmate_clin$OS)*30
Checkmate_clin$fustat <- ifelse(Checkmate_clin$OS_CNSR=="0",0,1)

kirc.Checkmate <- list(mRNA.exp=Checkmate_exp,
                       clin.info=Checkmate_clin)


yau.ntp.pred <- runNTP(expr       = kirc.Checkmate$mRNA.exp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR Checkmate OS") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.Checkmate$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR Checkmate OS") 


#PFS
Checkmate_clin$futime <- as.numeric(Checkmate_clin$PFS)*30
Checkmate_clin$fustat <- ifelse(Checkmate_clin$PFS_CNSR=="0",0,1)

kirc.Checkmate <- list(mRNA.exp=Checkmate_exp,
                       clin.info=Checkmate_clin)


yau.ntp.pred <- runNTP(expr       = kirc.Checkmate$mRNA.exp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR Checkmate PFS") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.Checkmate$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR Checkmate PFS") 



####JAVELIN验证####
load("/t8a/2022backup/xuzijun/vip39/database/KIRC_immune_therapy/JAVELIN_Renal_101_trial/JAVELIN_RCCall.rds")
JAVELIN_cli$futime <- as.numeric(JAVELIN_cli$PFS_P)*30
JAVELIN_cli$fustat <- JAVELIN_cli$PFS_P_CNSR

kirc.JAVELIN <- list(mRNA.exp=JAVELIN_exp,
                     clin.info=JAVELIN_cli)


yau.ntp.pred <- runNTP(expr       = kirc.JAVELIN$mRNA.exp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR JAVELIN") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.JAVELIN$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR JAVELIN") 


####GSE22541验证#####
load("/t8a/2022backup/xuzijun/vip39/database/KIRC_immune_therapy/GSE22541/GSE22541_RCCall.rds")

GSE22541_cli$futime <- as.numeric(GSE22541_cli$futime)
GSE22541_cli$fustat <- GSE22541_cli$fustat
GSE22541_exp <- GSE22541_exp[,rownames(GSE22541_cli)]


kirc.GSE22541 <- list(mRNA.exp=GSE22541_exp,
                      clin.info=GSE22541_cli)


yau.ntp.pred <- runNTP(expr       = kirc.GSE22541$mRNA.exp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR GSE22541") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.GSE22541$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR GSE22541")



####GSE29609####
load("/t8a/2022backup/xuzijun/vip39/database/KIRC_immune_therapy/GSE29609/GSE29609_RCCall.rds")

GSE29609_cli$futime <- as.numeric(GSE29609_cli$futime)
GSE29609_cli$fustat <- GSE29609_cli$fustat
GSE29609_exp <- GSE29609_exp[,rownames(GSE29609_cli)]


kirc.GSE29609 <- list(mRNA.exp=GSE29609_exp,
                      clin.info=GSE29609_cli)


yau.ntp.pred <- runNTP(expr       = kirc.GSE29609$mRNA.exp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR GSE29609") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.GSE29609$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR GSE29609")

###IMvigor210#####
load(file = "/t8a/2022backup/xuzijun/vip39/database/KIRC_immune_therapy/IMvigor210/IMvigor210all.rds")


IMvigor210_cli$futime <- as.numeric(IMvigor210_cli$os)*30
IMvigor210_cli$fustat <- IMvigor210_cli$censOS
IMvigor210_exp <- IMvigor210_exp[,rownames(IMvigor210_cli)]


kirc.IMvigor210 <- list(mRNA.exp=IMvigor210_exp,
                        clin.info=IMvigor210_cli)


yau.ntp.pred <- runNTP(expr       = kirc.IMvigor210$mRNA.exp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR IMvigor210") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.IMvigor210$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR IMvigor210")



####ETAB3267####
load("/t8a/2022backup/xuzijun/vip39/database/KIRC_immune_therapy/E-MTAB-3267/MTAB3267_RCCall.rds")

MTAB3267_cli$futime <- as.numeric(MTAB3267_cli$`Characteristics[progression free survival]`)*30
MTAB3267_cli$fustat <- MTAB3267_cli$`Characteristics[progression]`
MTAB3267_exp <- MTAB3267_exp[,rownames(MTAB3267_cli)]


kirc.MTAB3267 <- list(mRNA.exp=MTAB3267_exp,
                      clin.info=MTAB3267_cli)


yau.ntp.pred <- runNTP(expr       = kirc.MTAB3267$mRNA.exp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR MTAB3267") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.MTAB3267$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR MTAB3267")


####CPTAC队列验证####
setwd("~/TCGA/Sdeath_project/supp2024/")
load("~/TCGA/cptac_CCRCC/CTPACneed.rds")
dim(CTPAC_mRNA_exp);dim(CTPAC_cli)
kirc.CTPAC <- list(mRNA.exp=CTPAC_mRNA_exp,
                   clin.info=CTPAC_cli)

kirc.CTPAC$mRNA.exp[1:4,1:4]
table(marker.up$templates$probe%in%rownames(kirc.CTPAC$mRNA.exp))

yau.ntp.pred <- runNTP(expr       = kirc.CTPAC$mRNA.exp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR CTPAC") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.CTPAC$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR CTPAC") 


####蛋白层面

load("~/TCGA/cptac_CCRCC/CTPACneed_pro.rds")
setwd("~/TCGA/Sdeath_project/supp2024/")

kirc.CTPAC <- list(mRNA.exp=CTPAC_pro_exp,
                   clin.info=CTPAC_cli)

yau.ntp.pred <- runNTP(expr       = kirc.CTPAC$mRNA.exp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR CTPAC based on Pro") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.CTPAC$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR CTPAC based on Pro") 

###gene in ccRCC   SLC7A11####


###实验图#####
