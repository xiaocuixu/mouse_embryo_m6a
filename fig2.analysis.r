fp <- read.table(file="../../all.embryo.fpkm.txt")
colnames(fp)=c("gene","GV","MII","Zygote","cell2","cell4")

#define zga and decay gene
decay <- fp %>% dplyr::filter((GV>5 | MII>5) & GV/(cell2+0.0001) >4)
zga <- fp %>% dplyr::filter((GV<1 & MII<1 ) & (Zygote>1 | cell2>1) & cell2/(GV+0.0001) >4)
decay$type <- "decay"
zga$type <- "zga"
fp_fil <- rbind(decay,zga)


#combine m6a information
files = list.files("../../genesWithExonPeak", pattern = "em.", full.names = TRUE) ##gene with exon peak
peak <- map(files,read.table,header=T) %>% reduce(full_join,by="gene")
colnames(peak) <- c("gene","cell2.p","cell4.p","GV.p","MII.p","Zygote.p")
peak <- select(peak,"gene","GV.p","MII.p","Zygote.p","cell2.p","cell4.p")
myfun_mu <- function(a){b=ifelse(is.na(a),0,a);return(b)}
all <- left_join(fp_fil,peak) %>% mutate_at(vars("GV.p","MII.p","Zygote.p","cell2.p","cell4.p"),myfun_mu) %>% 
mutate(GV.p=ifelse(GV<1,0,GV.p),
	MII.p=ifelse(MII<1,0,MII.p),
	Zygote.p=ifelse(Zygote<1,0,Zygote.p),
	cell2.p=ifelse(cell2<1,0,cell2.p),
	cell4.p=ifelse(cell4<1,0,cell4.p))   #####filter out peak with fpkm<1


#classify decay and zga gene into m6a+ and m6a-
all2 <- all %>% mutate(m.p=ifelse(GV.p=="1" | MII.p=="1",1,0),z.p = ifelse(Zygote.p=="1" |cell2.p=="1" ,1,0))
all2$gr=paste(all2$type,all2$m.p,all2$z.p,sep=".")
all3 <- all2 %>% mutate(biggr=fct_collapse(as.factor(gr),
						zga_nopeak=c("zga.0.0"),
						zga_zpeak=c("zga.0.1"),
						decay_nopeak=c("decay.0.0","decay.0.1"),
						decay_mpeak=c("decay.1.0","decay.1.1")))

#plot heatmap
scl <- myfun_scl(all3,2:6)
ht5 <- Heatmap(as.matrix(scl[,17:21]),name="z-score",col = col_fpkm,use_raster=FALSE,cluster_rows=FALSE,cluster_columns = FALSE,
	row_split=factor(scl$biggr,levels=c("decay_nopeak","decay_mpeak","zga_nopeak","zga_zpeak")),cluster_row_slices = FALSE,
	row_title_rot = 0,column_names_rot=45,show_row_names=FALSE)
ht6 <- Heatmap(as.matrix(scl[,8:12]),name="m6a",col = col_m6a,use_raster=FALSE,cluster_columns = FALSE,column_names_rot=45,show_row_names=FALSE)
ht5+ht6
dev.copy2pdf(file="cluster.ht.pdf")

##out put genelist
write.table(file="genelist.decay_nopeak.txt",scl[scl$biggr=="decay_nopeak",1],quote=F,col.names=F,row.names=F,sep="\t")
write.table(file="genelist.decay_mpeak.txt",scl[scl$biggr=="decay_mpeak",1],quote=F,col.names=F,row.names=F,sep="\t")
write.table(file="genelist.zga_nopeak.txt",scl[scl$biggr=="zga_nopeak",1],quote=F,col.names=F,row.names=F,sep="\t")
write.table(file="genelist.zga_zpeak.txt",scl[scl$biggr=="zga_zpeak",1],quote=F,col.names=F,row.names=F,sep="\t")

#heatmap for m6A gained genes
denovo <- scl %>% filter(gr == "zga.0.1" | gr == "decay.0.1")
ht7 <- Heatmap(as.matrix(denovo[,17:21]),name="z-score",col = col_fpkm,use_raster=FALSE,cluster_rows=FALSE,cluster_columns = FALSE,
	row_split=factor(denovo$gr,levels=c("decay.0.1","zga.0.1")),cluster_row_slices = FALSE,
	row_title_rot = 0,column_names_rot=45,show_row_names=FALSE,
	column_title = "N(98,754)")
ht8 <- Heatmap(as.matrix(denovo[,8:12]),name="m6a",col = col_m6a,use_raster=FALSE,cluster_columns = FALSE,column_names_rot=45,show_row_names=FALSE)
dev.copy2pdf(file="denovo.m6A.ht.pdf")


###GO analysis
library("org.Mm.eg.db")
library(clusterProfiler)
gene.zga_nopeak  <- bitr(scl[scl$biggr=="zga_nopeak",]$gene, fromType = "SYMBOL",toType = c("ENSEMBL","ENTREZID"),OrgDb = org.Mm.eg.db)
gene.zga_zpeak  <- bitr(scl[scl$biggr=="zga_zpeak",]$gene, fromType = "SYMBOL",toType = c("ENSEMBL","ENTREZID"),OrgDb = org.Mm.eg.db)
gene.decay_mpeak  <- bitr(scl[scl$biggr=="decay_mpeak",]$gene, fromType = "SYMBOL",toType = c("ENSEMBL","ENTREZID"),OrgDb = org.Mm.eg.db)
gene.decay_nopeak  <- bitr(scl[scl$biggr=="decay_nopeak",]$gene, fromType = "SYMBOL",toType = c("ENSEMBL","ENTREZID"),OrgDb = org.Mm.eg.db)
cp <- list(decay_nopeak=gene.decay_nopeak$SYMBOL,decay_mpeak=gene.decay_mpeak$SYMBOL,
	zga_nopeak=gene.zga_nopeak$SYMBOL,zga_zpeak=gene.zga_zpeak$SYMBOL)
xx <- compareCluster(cp,fun="enrichGO",pAdjustMethod = "BH",ont= "BP",OrgDb= org.Mm.eg.db,keyType= 'SYMBOL',
	qvalueCutoff = 0.05)
dotplot(xx,showCategory=15,includeAll=TRUE)
dev.copy2pdf(file="multi_go.dotplot.pdf")
write.table(file="multigo.txt",as.data.frame(xx),quote=F,sep="\t",row.names=F)
####plot heatmap for GO results
library(data.table)
go_dt <- setDT(as.data.frame(xx))
go_sel <- read.table("selected.GOterm.csv",header=T,sep="\t")
go_final <- go_dt[ID %in% go_sel$ID,.(Cluster,Description,p.adjust)]
go_final$log10 <- -log10(go_final$p.adjust)
go_final <- go_final[,-3]
go_final.w <- dcast(go_final,Description ~ Cluster)
go_final.w[is.na(go_final.w)] <- 1
go_final.l <- melt(go_final.w,id.vars=c("Description"))
go_final.l <- go_final.l %>%  mutate(log10P_gr= case_when(value<=1 ~ "1",
												value >1 & value <=2 ~ "0.1",
												value >2 & value <=3 ~ "0.01",
												value >3 & value <=4 ~ "0.001",
												value >4  ~ "0.0001")) 
go_final.l <- go_final.l[,-3]
go_final.w2 <- reshape2::dcast(go_final.l,Description ~ variable)
go_final.s <- go_final.w2 %>% arrange(decay_nopeak,decay_mpeak,zga_nopeak,zga_zpeak)
rownames(go_final.s ) <- go_final.s$Description
# c(brewer.pal(9,"Blues")[3],brewer.pal(9,"Reds")[3:6])
col_fun <- c("1"="#C6DBEF", "0.1"="#FCBBA1","0.01"="#FC9272","0.001"="#FB6A4A","0.0001"="#EF3B2C")
ht <- Heatmap(as.matrix(go_final.s[,2:5]),
	name="P.adj",col = col_fun,use_raster=FALSE,
	cluster_rows=FALSE,cluster_columns = FALSE,rect_gp = gpar(col = "white", lwd = 2),
	row_names_max_width = max_text_width(
        rownames(go_final.s), 
        gp = gpar(fontsize = 12)))
dev.copy2pdf(file="multi_go.ht.pdf")


##M-decay z-decay definition
m.log2 <- fp %>% dplyr::filter(GV>2) %>% mutate_if(is.numeric,myfun_log2)
m.log2 <- m.log2 %>% mutate(cluster= case_when(GV>Zygote+1 & Zygote < cell2+1 ~ "M-decay",
GV<Zygote+1 & GV>Zygote-1 & Zygote >cell2+1 ~ "Z-decay1",
GV>Zygote+1 & Zygote >cell2+1 ~ "Z-decay2",
GV<Zygote+1 & GV>Zygote-1 & Zygote<cell2+1 & Zygote > cell2-1 ~ "stable"
))  %>% select(gene,cluster) 
all3.decay <- all3 %>% dplyr::filter(biggr=="decay_nopeak" || biggr == "decay_mpeak") %>% select(-type,-gr,-m.p,-z.p)
all3.decay.type <- inner_join(all3.decay,m.log2,by="gene")
all3.decay.type.scl <- myfun_scl(all3.decay.type,2:6)
#MD gene pie plot
all3.decay.type.sum <- all3.decay.type %>% group_by(cluster) %>% summarize(n=n()) %>% mutate(pro=n/sum(n))
all3.decay.type.sum[is.na(all3.decay.type.sum)] <- "notdetermine"
all3.decay.type.sum.2 <- all3.decay.type.sum %>% dplyr::select(cluster,n)
all3.decay.type.sum.2.l <- melt(all3.decay.type.sum.2)
all3.decay.type.sum.2.w <- dcast(all3.decay.type.sum.2.l,variable~cluster)
colnames(all3.decay.type.sum.2.w) <- c("variable","Mdecay","notdetermine","stable","Zdecay1","Zdecay2")
all3.decay.type.sum.2.w <- all3.decay.type.sum.2.w %>% dplyr::mutate(Zdecay=Zdecay1+Zdecay2,other=notdetermine+stable) %>% 
dplyr::select(Mdecay,Zdecay,other)
all3.decay.type.sum.2.l2 <- melt(all3.decay.type.sum.2.w)
all3.decay.type.sum.2.l2 <- all3.decay.type.sum.2.l2 %>% dplyr::mutate(p=value/sum(value))
pie.scales <- all3.decay.type.sum.2.l2$p
names(pie.scales) <- paste(all3.decay.type.sum.2.l2$variable,all3.decay.type.sum.2.l2$value,sep="\n")
pie(pie.scales, col = brewer.pal(3,"Set1"),main="MD genes")
dev.copy2pdf(file="MD.gene.classification.2.pie.pdf") 
#MD gene bar plot
all3.decay.fil <- all3.decay.type %>% dplyr::filter(cluster != "stable" &  !is.na(cluster)) %>% 
dplyr::mutate(type=ifelse(cluster=="M-decay","Mdecay","Zdecay"))
all3.decay.fil %>% ggplot(aes(type,fill=biggr))+geom_bar(stat="count",position='fill')+
scale_fill_manual(name="",values=c("#6D9EC1","#E46726"))+
geom_text(stat='count',aes(label=..count..), color="white", size=3.5,position=position_fill(0.5))+
theme_minimal()+ylab("proportion")+xlab("")+theme(axis.text.x= element_text(angle=45,hjust=1))+labs(title="fisher.test:3.94e-05")
dev.copy2pdf(file="MD.gene.classification.m6a.stat.bar.pdf") 
write.table3(all3.decay.fil[,c(1,12,ncol(all3.decay.fil))],"MD.ov.mdecay.zdecay.txt")
#MD gene box plot
all3.decay.fil.l <- melt(all3.decay.fil,measure.vars=c("GV","MII","Zygote","cell2","cell4"))
p<-ggplot(all3.decay.fil.l,aes(x=variable,y=log2(1+value),fill=biggr))+
geom_boxplot(outlier.color=NA)+
facet_grid(~type)+labs(title="MD gene")+xlab("")+ylab("log2(fpkm+1)")+
scale_fill_manual(name="m6a",values=c("#6D9EC1","#E46726"))+theme_bw()
p+theme(axis.text.y=element_text(size=15),axis.title.y = element_text(size = 20),
	axis.title.x = element_text(size = 20),axis.text.x= element_text(size = 20,angle=45,hjust=1))
dev.copy2pdf(file="MD.gene.classification.fpkm.box.pdf")
#pval mdecay GV-c4: 0.3139789,0.5127467,0.7574035,0.3054662,0.007083418
wilcox.test(all3.decay.fil[biggr=="decay_mpeak" & type=="Mdecay"]$GV,all3.decay.fil[biggr=="decay_nopeak" & type=="Mdecay"]$GV)$p.val
wilcox.test(all3.decay.fil[biggr=="decay_mpeak" & type=="Mdecay"]$MII,all3.decay.fil[biggr=="decay_nopeak" & type=="Mdecay"]$MII)$p.val
wilcox.test(all3.decay.fil[biggr=="decay_mpeak" & type=="Mdecay"]$Zygote,all3.decay.fil[biggr=="decay_nopeak" & type=="Mdecay"]$Zygote)$p.val
wilcox.test(all3.decay.fil[biggr=="decay_mpeak" & type=="Mdecay"]$cell2,all3.decay.fil[biggr=="decay_nopeak" & type=="Mdecay"]$cell2)$p.val
wilcox.test(all3.decay.fil[biggr=="decay_mpeak" & type=="Mdecay"]$cell4,all3.decay.fil[biggr=="decay_nopeak" & type=="Mdecay"]$cell4)$p.val

#pval Zdecay GV-c4:0.003649774,0.02058349,4.106473e-06,0.000233244,0.5693552
wilcox.test(all3.decay.fil[biggr=="decay_mpeak" & type=="Zdecay"]$GV,all3.decay.fil[biggr=="decay_nopeak" & type=="Zdecay"]$GV)$p.val
wilcox.test(all3.decay.fil[biggr=="decay_mpeak" & type=="Zdecay"]$MII,all3.decay.fil[biggr=="decay_nopeak" & type=="Zdecay"]$MII)$p.val
wilcox.test(all3.decay.fil[biggr=="decay_mpeak" & type=="Zdecay"]$Zygote,all3.decay.fil[biggr=="decay_nopeak" & type=="Zdecay"]$Zygote)$p.val
wilcox.test(all3.decay.fil[biggr=="decay_mpeak" & type=="Zdecay"]$cell2,all3.decay.fil[biggr=="decay_nopeak" & type=="Zdecay"]$cell2)$p.val
wilcox.test(all3.decay.fil[biggr=="decay_mpeak" & type=="Zdecay"]$cell4,all3.decay.fil[biggr=="decay_nopeak" & type=="Zdecay"]$cell4)$p.val
