#TE heatmap
en <- read.table("peakEnrich.noGenic.TE.txt",header=T)
exp <- read.table("TE.noGenic.rpkm.average.name.txt")
colnames(exp) <- c("repName","chr","s","e","GV","MII","Zygote","2cell","4cell","ESC","morula","blast")
exp.l <- melt(exp,id.vars=c("repName","chr","s","e"))
colnames(exp.l) <- c("repName","chr","s","e","sample","rpkm")
exp.m <- exp.l %>% group_by(repName,sample) %>% summarise(rpkm.m=round(mean(rpkm),3))
exp.m.w <- dcast(exp.m,repName ~ sample)
en$enrich <- log2(en$enrich)
colnames(en)[1] <- "repName"
en_log2.w <- dcast(en,repName ~ sample)
comb <- merge(exp.m.w,en_log2.w,by="repName")
colnames(comb) <- gsub("\\.y",".e",colnames(comb))
colnames(comb) <- gsub("\\.x","",colnames(comb))
colnames(comb) <- gsub("2cell","c2",colnames(comb))
colnames(comb) <- gsub("4cell","c4",colnames(comb))
comb2 <- comb %>% select(repName,GV,MII,Zygote,c2,c4,morula,blast,ESC,GV.e,MII.e,Zygote.e,c2.e,c4.e,ESC.e)
comb_f <- comb2 %>% filter_at(vars(GV,MII,Zygote,c2,c4,morula,blast), any_vars(. >5))
scl <- myfun_scl(comb_f,c(2:9))
colnames(scl) <- gsub("\\.e","",colnames(scl))
colnames(scl) <- gsub("\\.z","",colnames(scl))
rownames(scl)=scl$repName
library(ComplexHeatmap)
library(circlize)
col_fpkm <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
col_fun = colorRamp2(c(-5, 0, 2,6), c("#6D9EC1", "white","#fff5eb", "#E46726"))
ht1 <- Heatmap(as.matrix(scl[,16:23]),name="z-score",col = col_fpkm,use_raster=FALSE,cluster_columns = FALSE,show_row_dend=FALSE,column_names_rot = 45)
ht2 <- Heatmap(as.matrix(scl[,10:15]),name="peak enrichment\nlog2(obs/exp)",col = col_fun,use_raster=FALSE,cluster_columns = FALSE,column_names_rot = 45)
ht1+ht2
dev.copy2pdf(file="TE.byRepName.rpkm5.zScore.enrich.ht.pdf")



##maternal and zga TE reads count plot
files <- c("maternal.TE.GV.count.txt","maternal.TE.MII.count.txt","zga.TE.2cell.count.txt","zga.TE.4cell.count.txt")
path <- "/home/xuxc/projects/m6a/new_pipe/transposon/byRepName"
df <- myfun_f_rbind(files,path,pattern="count.txt")
df$type <- gsub("TE.","",df$type)
df3 <- df %>% separate(type,c("class","stage"),sep="\\.")
df3.m <- df3 %>% filter(class=="maternal")
df3.z <- df3 %>% filter(class=="zga")
p1 <- ggplot(df3.m)+
  geom_bar(aes(x=factor(V1,levels=c("MTA_Mm-int","MTA_Mm","RLTR10-int","RLTR10","MT-int","RLTR4_Mm")),
      y=V2,fill=stage),
      stat="identity",position="dodge")+
  guides(fill=guide_legend(title=NULL))+
  theme_classic(18)+
  theme(axis.text.x= element_text(angle=45,hjust=1))+
  ylab("reads count")+xlab("")+
  scale_fill_manual(values=brewer.pal(8,"Dark2")[1:2])
p2 <- ggplot(df3.z)+
  geom_bar(aes(x=factor(V1,levels=c("MERVL-int","ORR1A0-int","ORR1A1-int","ORR1A0","MT2_Mm")),
      y=V2,fill=stage),
      stat="identity",position="dodge")+
  guides(fill=guide_legend(title=NULL))+
  theme_classic(18)+
  theme(axis.text.x= element_text(angle=45,hjust=1))+
  ylab("reads count")+xlab("")+
  scale_fill_manual(values=brewer.pal(8,"Dark2")[4:5])
ggarrange(p1,p2)
ggsave("TE.readsnum.bar.pdf")



##maternal and zga TE average RMPM and copy number plot
ann <- read.table("~/ann/mm9.rmsk.repName.TE.noGenic.srt.bed")
ann <- ann %>% mutate(l=V3-V2) %>% select(l,V4,V5)
ann.s <- ann %>% group_by(V4,V5) %>% summarise(L=sum(l),n=n())
colnames(ann.s) <- c("repName","repClass","L","N")
exp <- read.table("TE.noGenic.rpkm.average.name.txt")
colnames(exp) <- c("repName","chr","s","e","GV","MII","Zygote","2cell","4cell","ESC","morula","blast")
exp.l <- exp %>% gather(key="sample",value="rpkm",-repName,-chr,-s,-e)
exp.m <- exp.l %>% group_by(repName,sample) %>% summarise(rpkm.m=round(mean(rpkm),3))
exp.m.w <- spread(exp.m,key="sample",value="rpkm.m")
df <- full_join(exp.m.w,ann.s)
df[is.na(df)] <- 0
colnames(df)[2] <- "cell2"
colnames(df)[3] <- "cell4"
dat <- df %>% select(repName,GV,MII,cell2,cell4,repClass,L,N)
gv <- dat %>% dplyr::filter(GV>=5)
c2 <- dat %>% dplyr::filter(GV<5 & cell2>=5 & repName != "B2_Mm1a")
maternal <- gv %>% select(repName,GV,MII,N)
maternal.l <- maternal %>% select(-N) %>% gather(key="stage",value="meanRPKM",-repName)
p1 <- ggplot(maternal)+
  geom_bar(aes(x=factor(repName,levels=c("MTA_Mm-int","MTA_Mm","RLTR10-int","RLTR10","MT-int","RLTR4_Mm")),
      y=N),fill="#0099B4FF",
      stat="identity")+
  theme_classic(18)+
  theme(axis.text.x= element_text(angle=45,hjust=1))+
  ylab("copy number")+xlab("")+ggtitle("Maternal TE")
p2 <- ggplot(maternal.l)+
  geom_bar(aes(x=factor(repName,levels=c("MTA_Mm-int","MTA_Mm","RLTR10-int","RLTR10","MT-int","RLTR4_Mm")),
      y=meanRPKM,fill=stage),
      stat="identity",position="dodge")+
  guides(fill=guide_legend(title=NULL))+
  theme_classic(18)+
  theme(axis.text.x= element_text(angle=45,hjust=1))+
  xlab("")+
  scale_fill_manual(values=brewer.pal(8,"Dark2")[1:2])
zga <- c2 %>% select(repName,cell2,cell4,N)
zga.l <- zga %>% select(-N) %>% gather(key="stage",value="meanRPKM",-repName)
p3 <- ggplot(zga)+
  geom_bar(aes(x=factor(repName,levels=c("MERVL-int","ORR1A0-int","ORR1A1-int","ORR1A0","MT2_Mm")),
      y=N),fill="#0099B4FF",
      stat="identity")+
  theme_classic(18)+
  theme(axis.text.x= element_text(angle=45,hjust=1))+
  ylab("copy number")+xlab("")+ggtitle("ZGA TE")
p4 <- ggplot(zga.l)+
  geom_bar(aes(x=factor(repName,levels=c("MERVL-int","ORR1A0-int","ORR1A1-int","ORR1A0","MT2_Mm")),
      y=meanRPKM,fill=stage),
      stat="identity",position="dodge")+
  guides(fill=guide_legend(title=NULL))+
  theme_classic(18)+
  theme(axis.text.x= element_text(angle=45,hjust=1))+
  xlab("")+
  scale_fill_manual(values=brewer.pal(8,"Dark2")[4:5])
ggarrange(p1,p2,p3,p4)
ggsave("TE.zga.maternal.RPKM.copy.num.bar.pdf")


