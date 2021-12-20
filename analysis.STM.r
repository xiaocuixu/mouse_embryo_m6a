library(DESeq2)
require("ggrepel")

###PCA
mycounts<- read.table("/home/xuxc/projects/m6a/new_pipe/validate.M3.STM/featureCounts.RNA.M3.STM.txt",header=T)
rownames(mycounts)=mycounts[,1]
mycounts<-mycounts[,-1:-6]
colnames(mycounts) <- gsub("bam.RNA.M3.","",colnames(mycounts))
colnames(mycounts) <- gsub(".bam","",colnames(mycounts))
sample.vec <- gsub(".[12345]$","",colnames(mycounts))
sample.level <- unique(sample.vec)
condition <- factor(sample.vec,levels=sample.level)
colData <- data.frame(row.names=colnames(mycounts), condition)
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
vsd <- vst(dds,blind=FALSE)
pcaData <- plotPCA(vsd, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2,color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_text_repel(aes(label = condition),size =4,max.overlaps=50)+
  theme_bw(18)+theme(legend.position="none")
ggsave("PCA.RNA.M3.STM.pdf")


###zga gene boxplot
zga_nopeak <- read.table("/home/xuxc/projects/m6a/new_pipe/global/fpkm_cutoff/genelist.zga_nopeak.txt")
zga_zpeak <- read.table("/home/xuxc/projects/m6a/new_pipe/global/fpkm_cutoff/genelist.zga_zpeak.txt")
zga_nopeak$stat <- "zga.m6a-"
zga_zpeak$stat <- "zga.m6a+"
zga <- rbind(zga_nopeak,zga_zpeak)
colnames(zga)[1] <- c("gene")
fp <- myfun_f_cbind(c("RNA.M3.2C-DMSO.uniq.ave.fpkm.txt",
	"RNA.M3.2C-STM.uniq.ave.fpkm.txt",
	"RNA.M3.4C-DMSO.uniq.ave.fpkm.txt",
	"RNA.M3.4C-STM.uniq.ave.fpkm.txt",
	"RNA.M3.Morula-DMSO.uniq.ave.fpkm.txt",
	"RNA.M3.Morula-STM.uniq.ave.fpkm.txt"),
    "/home/xuxc/projects/m6a/new_pipe/validate.M3.STM/fpkm.mean",
    "uniq.ave.fpkm.txt")
colnames(fp) <- gsub("RNA.M3.","",colnames(fp))
colnames(fp) <- gsub("-","_",colnames(fp))
colnames(fp) <- gsub("2C","C2",colnames(fp))
colnames(fp) <- gsub("4C","C4",colnames(fp))
zga_fc <- inner_join(zga,fp) %>% 
mutate(logfc_2C=log2((0.0001+C2_STM)/(0.0001+C2_DMSO)),
	logfc_4C=log2((0.0001+C4_STM)/(0.0001+C4_DMSO)),
	logfc_Morula=log2((0.0001+Morula_STM)/(0.0001+Morula_DMSO))) %>%
select(gene,stat,logfc_2C:logfc_Morula) %>%
gather(key="sample",value="log2fc_fpkm",-gene,-stat)
stat.test <- zga_fc  %>% group_by(sample) %>% wilcox_test(log2fc_fpkm~stat)
stat.test <- stat.test %>%   add_xy_position(x = "sample", dodge = 0.8)
stat.test$y.position <- 3
ggplot(zga_fc,aes(sample,log2fc_fpkm))+geom_boxplot(aes(fill=stat),outlier.color=NA)+
coord_cartesian(ylim=c(-3,3))+
scale_fill_manual(values=c("#6D9EC1","#E46726"))+
geom_hline(aes(yintercept=0),col="#990000",linetype="dashed")+
theme_classic()+theme(axis.text.y=element_text(size=15),axis.title.y = element_text(size = 15),
	axis.text.x= element_text(size = 15),legend.text = element_text(size = 15),
	legend.title = element_text(size = 0),plot.title=element_text(size=25))+xlab("")+
    ylab("log2FC(STM/DMSO)")+
stat_pvalue_manual(
    stat.test, label = "p", tip.length = 0.002)
ggsave("bxp.FC.zga_gene.RNA.M3.STM.pdf")


##volcano plot
samples <- c("2C.STM.2C.DMSO","4C.STM.4C.DMSO","Morula.STM.Morula.DMSO")
for(i in 1:3)
{
	s <- samples[i]
df <- read.table(paste0("DEG/",s,".DEG.res.tsv"),header=T)
df2 <- df %>% dplyr::filter(!is.na(log2FoldChange) & !is.na(padj)) 
df2$threshold = factor(ifelse(df2$padj < 0.05 & abs(df2$log2FoldChange) >= 1, ifelse(df2$log2FoldChange>= 1 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
df2 %<>% mutate(log10padj=-log10(padj)) %>% mutate(log10padj=ifelse(log10padj>10,10,log10padj))
up_num <- length(df2[df2$threshold == "Up",]$gene)
down_num <- length(df2[df2$threshold == "Down",]$gene)
ggplot(df2,aes(x=log2FoldChange,y=log10padj,color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+
  theme_classic()+
  theme(
    legend.title = element_blank()
  )+labs(title=paste(s,"\n","up",up_num,"down",down_num))+
  geom_vline(xintercept=c(-1,1),lty="dashed",col="grey",lwd=1) +
  geom_hline(yintercept = -log10(0.05),lty="dashed",col="grey",lwd=1) +
  theme(axis.text.y=element_text(size=15),axis.title.y = element_text(size = 15),
	axis.title.x = element_text(size = 15),
	axis.text.x= element_text(size = 15),legend.text = element_text(size = 15),
	plot.title=element_text(size=25))+xlab("log2FC(STM/DMSO)")+ylab("-log10(padj)")
	ggsave(file=paste("volcano.DEG",s,"pdf",sep="."))
}



###MERVL boxplot
mervl <- myfun_f_rbind(
	c("RNA.M3.2C-DMSO.tab","RNA.M3.2C-STM.tab","RNA.M3.4C-DMSO.tab",
		"RNA.M3.4C-STM.tab","RNA.M3.Morula-DMSO.tab","RNA.M3.Morula-STM.tab"),
	"/home/xuxc/projects/m6a/new_pipe/validate.M3.STM",
	"tab")
mervl %<>% separate(type,c("stage","treatment"),sep="-")
colnames(mervl)[7] <- "RPKM"
colnames(mervl)[3] <- "id"
mervl %<>% select(id,stage,treatment,RPKM)
mervl$stage <- gsub("RNA.M3.","",mervl$stage)
stat.test <- mervl  %>% group_by(stage) %>% wilcox_test(RPKM ~ treatment)
stat.test <- stat.test %>%   add_xy_position(x = "stage", y.trans = function(x){log2(x+1)})
p1 <- ggplot(mervl,aes(stage,log2(RPKM+1))) +
geom_boxplot(aes(fill=treatment),outlier.colour=NA)+ylab("log2(RPKM+1)") +
scale_fill_manual(values=c("black","grey"))+
theme_classic()+theme(axis.text.x = element_text(size = 15),axis.text.y=element_text(size=15),
axis.title.y = element_text(size = 15),legend.title=element_blank(), legend.text = element_text(size = rel(1.5)))+
xlab("")+labs(title="M3.STM MERVL")+
stat_pvalue_manual(
    stat.test, label = "p", tip.length = 0.002)

exp.mean <- mervl %>% group_by(stage,treatment) %>% summarise(m=mean(RPKM)) %>%
	mutate(rel=ifelse(treatment=="DMSO",
	m/exp.mean[exp.mean$stage=="2C" & exp.mean$treatment=="DMSO",]$m,
	m/exp.mean[exp.mean$stage=="2C" & exp.mean$treatment=="STM",]$m))

p2 <- ggplot(exp.mean,aes(x=stage,y=rel,col=treatment))+
geom_line(aes(group=treatment),size=1)+
geom_point(size=2)+  theme_classic(15) +
  scale_color_manual(values=c("black","grey"))+
  ylab("relative to 2C")+ggtitle("M3.STM MERVL")

ggarrange(p1,p2,labels = c("A","B"),
          nrow = 2, ncol = 2)
ggsave("MERVL.M3.STM.pdf")


##zga decay gene plot
zga_fp <- inner_join(zga,fp)
zga_decay2 <- zga_fp %>% 
	mutate(decay=ifelse((C2_DMSO+0.1)/(Morula_DMSO+0.1)>3,"yes","no")) %>%  #define transient or  sustained gene
	select(gene,stat,decay,C2_DMSO,C4_DMSO,Morula_DMSO,C2_STM,C4_STM,Morula_STM) %>% 
	filter(!grepl("Mir",.$gene),rowMeans(.[4:9])>0) %>% na.omit()
out <- zga_decay2 %>% 
	mutate(type=ifelse(stat=="zga.m6a-",paste0("zga.m6a.minus.decay.",decay),paste0("zga.m6a.plus.decay.",decay))) %>%
	select(gene,type) %>%
  group_nest(type)
out %>%
  pwalk(
    ~ write.table2(..2, paste0( ..1, ".genelist.txt"))
  )

summary <- zga_decay2 %>% group_by(stat,decay) %>% summarise(n=n())
p <- fisher.test(matrix(summary$n,ncol=2))$p
p8 <- zga_decay2 %>% ggplot(aes(decay,fill=stat))+geom_bar(stat="count",position='fill')+
scale_fill_manual(name="",values=c("#6D9EC1","#E46726"))+
geom_text(stat='count',aes(label=..count..), color="white", size=3.5,position=position_fill(0.5))+
theme_minimal(15)+ylab("proportion")+xlab("")+theme(axis.text.x= element_text(angle=45,hjust=1))+
labs(title=paste0("zga decay\nfisher.p:",p))

zga_decay2_fc <-  zga_decay2 %>%
mutate(logfc_2C=log2((0.1+C2_STM)/(0.1+C2_DMSO)),
	logfc_4C=log2((0.1+C4_STM)/(0.1+C4_DMSO)),
	logfc_Morula=log2((0.1+Morula_STM)/(0.1+Morula_DMSO))) %>%
select(gene,stat,decay,logfc_2C:logfc_Morula) %>%
gather(key="sample",value="log2fc_fpkm",-gene,-stat,-decay)
zga_decay2_fc$sample <- gsub("logfc_","",zga_decay2_fc$sample)
stat.test <- zga_decay2_fc  %>% group_by(decay,sample) %>% wilcox_test(log2fc_fpkm~stat)
stat.test <- stat.test %>%   add_xy_position(x = "sample", dodge = 0.8)
stat.test$y.position <- 3.5
p9 <- ggplot(zga_decay2_fc,aes(sample,log2fc_fpkm))+geom_boxplot(aes(fill=stat),outlier.color=NA)+
coord_cartesian(ylim=c(-2,4))+
scale_fill_manual(values=c("#6D9EC1","#E46726"))+
geom_hline(aes(yintercept=0),col="#990000",linetype="dashed")+
theme_classic()+theme(axis.text.y=element_text(size=15),axis.title.y = element_text(size = 15),
	axis.text.x= element_text(size = 15),legend.text = element_text(size = 15),
	legend.title = element_text(size = 0),plot.title=element_text(size=25))+xlab("")+
    ylab("log2FC(STM/DMSO)")+
facet_wrap(~decay)+
stat_pvalue_manual(
    stat.test, label = "p", tip.length = 0.002)

zga_decay2.long <- zga_decay2 %>% select(gene,decay,C2_DMSO:Morula_STM) %>%
gather(key="sample",value="FPKM",-gene,-decay) %>% 
		separate(sample,c("stage","treatment"),sep="_") 

p11 <- zga_decay2.long %>% filter(treatment=="DMSO") %>%
				ggplot()+geom_boxplot(aes(fill=stage,y=log2(FPKM+1),x=decay),outlier.colour=NA)+
				theme_classic(15)+
				scale_fill_brewer(palette = "Accent")+
				coord_cartesian(ylim=c(0,7.5))
				

		ggarrange(p11,p8,p9)
ggsave("zga.decay.V2.pdf")




###MERVL MTA IP boxplot
mervl.m6a <- myfun_f_rbind(
	c("MERVL.M3.l2c.D.tab","MERVL.M3.l2c.IP-D.tab",
		"MERVL.M3.l2c.IP-S.scale.tab","MERVL.M3.l2c.S.scale.tab"),
	"/home/xuxc/projects/m6a/new_pipe/validate.M3.STM",
	"tab")
colnames(mervl.m6a)[6] <- "RPKM"
colnames(mervl.m6a)[2] <- "id"
mervl.m6a %<>% select(id,type,RPKM)
mervl.m6a$type <- gsub("MERVL.M3.l2c.","",mervl.m6a$type)
mervl.m6a$type <- gsub(".scale","",mervl.m6a$type)
mervl.m6a %<>% mutate(treatment=ifelse(type=="D" | type=="S","input","IP"))
mervl.m6a$type <- gsub("IP-","",mervl.m6a$type)
stat.test <- mervl.m6a %>% group_by(treatment) %>% wilcox_test(RPKM ~ type)
stat.test <- stat.test %>%   add_xy_position(x = "treatment", y.trans = function(x){log2(x+1)})
p12 <- ggplot(mervl.m6a,aes(treatment,log2(RPKM+1))) +
geom_boxplot(aes(fill=type),outlier.colour=NA)+ylab("log2(RPKM+1)") +
scale_fill_manual(values=c("black","grey"))+
theme_classic()+theme(axis.text.x = element_text(size = 15),axis.text.y=element_text(size=15),
axis.title.y = element_text(size = 15),legend.title=element_blank(), legend.text = element_text(size = rel(1.5)))+
xlab("")+labs(title="M3.STM MERVL")+
stat_pvalue_manual(
    stat.test, label = "p", tip.length = 0.002)

mervl.m6a.ratio  <- mervl.m6a %>% spread(key="treatment",value="RPKM") %>% mutate(log2fc=log2((IP+0.1)/(input+0.1)))
stat.test <- mervl.m6a.ratio %>% wilcox_test(log2fc ~ type)
stat.test <- stat.test %>%   add_xy_position(x = "type")
p13 <- ggplot(mervl.m6a.ratio,aes(type,log2fc)) +
geom_boxplot(aes(fill=type),outlier.colour=NA)+ylab("log2fc(IP/input)") +
scale_fill_manual(values=c("black","grey"))+
theme_classic()+theme(axis.text.x = element_text(size = 15),axis.text.y=element_text(size=15),
axis.title.y = element_text(size = 15),legend.title=element_blank(), legend.text = element_text(size = rel(1.5)))+
xlab("")+labs(title="M3.STM MERVL")+
stat_pvalue_manual(
    stat.test, label = "p", tip.length = 0.002)
mta.m6a <- myfun_f_rbind(
	c("MTA.M3.l2c.D.tab","MTA.M3.l2c.IP-D.tab",
		"MTA.M3.l2c.IP-S.scale.tab","MTA.M3.l2c.S.scale.tab"),
	"/home/xuxc/projects/m6a/new_pipe/validate.M3.STM",
	"tab")
colnames(mta.m6a)[6] <- "RPKM"
colnames(mta.m6a)[2] <- "id"
mta.m6a %<>% select(id,type,RPKM)
mta.m6a$type <- gsub("MTA.M3.l2c.","",mta.m6a$type)
mta.m6a$type <- gsub(".scale","",mta.m6a$type)
mta.m6a %<>% mutate(treatment=ifelse(type=="D" | type=="S","input","IP"))
mta.m6a$type <- gsub("IP-","",mta.m6a$type)

stat.test <- mta.m6a %>% group_by(treatment) %>% wilcox_test(RPKM ~ type)
stat.test <- stat.test %>%   add_xy_position(x = "treatment", y.trans = function(x){log2(x+1)})
p14 <- ggplot(mta.m6a,aes(treatment,log2(RPKM+1))) +
geom_boxplot(aes(fill=type),outlier.colour=NA)+ylab("log2(RPKM+1)") +
scale_fill_manual(values=c("black","grey"))+
theme_classic()+theme(axis.text.x = element_text(size = 15),axis.text.y=element_text(size=15),
axis.title.y = element_text(size = 15),legend.title=element_blank(), legend.text = element_text(size = rel(1.5)))+
xlab("")+labs(title="M3.STM mta")+
stat_pvalue_manual(
    stat.test, label = "p", tip.length = 0.002)

mta.m6a.ratio  <- mta.m6a %>% spread(key="treatment",value="RPKM") %>% mutate(log2fc=log2((IP+0.1)/(input+0.1)))
stat.test <- mta.m6a.ratio %>% wilcox_test(log2fc ~ type)
stat.test <- stat.test %>%   add_xy_position(x = "type")
stat.test$y.position <-3
p15 <- ggplot(mta.m6a.ratio,aes(type,log2fc)) +
geom_boxplot(aes(fill=type),outlier.colour=NA)+ylab("log2fc(IP/input)") +
scale_fill_manual(values=c("black","grey"))+
theme_classic()+theme(axis.text.x = element_text(size = 15),axis.text.y=element_text(size=15),
axis.title.y = element_text(size = 15),legend.title=element_blank(), legend.text = element_text(size = rel(1.5)))+
xlab("")+labs(title="M3.STM mta")+
stat_pvalue_manual(
    stat.test, label = "p", tip.length = 0.002)
ggarrange(p12,p13,p14,p15)
ggsave("M3.STM.l2c.IP.input.box.MERVL.pdf")


###gene IP boxplot
df.gl <- data.frame()
for (gl in c('decay_mpeak','zga_zpeak','zga.m6a.plus.decay.yes','zga.m6a.plus.decay.no')) 
{
gl.m6a <- myfun_f_rbind(
	paste0(gl,".M3.l2c.",c("D.tab","IP-D.tab",
		"IP-S.scale.tab","S.scale.tab")),
	"/home/xuxc/projects/m6a/new_pipe/validate.M3.STM",
	"tab")
colnames(gl.m6a)[6] <- "RPKM"
colnames(gl.m6a)[2] <- "id"
gl.m6a %<>% select(id,type,RPKM)
gl.m6a$type <- gsub(paste0(gl,".M3.l2c."),"",gl.m6a$type)
gl.m6a$type <- gsub(".scale","",gl.m6a$type)
gl.m6a %<>% mutate(treatment=ifelse(type=="D" | type=="S","input","IP"))
gl.m6a$type <- gsub("IP-","",gl.m6a$type)
gl.m6a$gl <- gl
df.gl <- rbind(df.gl,gl.m6a)
}
stat.test <- df.gl %>% group_by(gl,treatment) %>% wilcox_test(RPKM ~ type)
stat.test <- stat.test %>%   add_xy_position(x = "treatment",y.trans = function(x){log2(x+1)})
p16 <- ggplot(df.gl,aes(treatment,log2(RPKM+1))) +
geom_boxplot(aes(fill=type),outlier.colour=NA)+ylab("log2(RPKM+1)\nstop_codon+-200bp") +
scale_fill_manual(values=c("black","grey"))+
theme_classic()+theme(axis.text.x = element_text(size = 15),axis.text.y=element_text(size=15),
axis.title.y = element_text(size = 15),legend.title=element_blank(), legend.text = element_text(size = rel(1.5)))+
xlab("")+facet_wrap(~gl)+
stat_pvalue_manual(
    stat.test, label = "p", tip.length = 0.002)
p16
ggsave("M3.STM.l2c.IP.input.box.gene.stop_codon.pdf")




