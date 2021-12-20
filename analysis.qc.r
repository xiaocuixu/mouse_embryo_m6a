####PCC plot
library(ggcorrplot)
files.input <- list.files("pcc.ratio",pattern="input.stop_codon.txt",full.names=T)
input <- myfun_f_cbind(files.input,"pcc.ratio","input.stop_codon.txt")
input.l <- input %>% gather("sample","value",-gene)
files.IP <- list.files("pcc.ratio",pattern="IP.stop_codon.txt",full.names=T)
IP <- myfun_f_cbind(files.IP,"pcc.ratio","IP.stop_codon.txt")
IP.l <- IP %>% gather("sample","value",-gene)

j.w <- inner_join(IP.l,input.l,by = c("gene", "sample"),suffix = c(".IP", ".input")) %>%
       mutate(ratio=(value.IP+1)/(value.input+1)) %>% 
       select(-value.input,-value.IP) %>% 
       pivot_wider(names_from =sample , values_from = ratio)
j.w$CV <- apply(j.w[,-1],1,CoefVar)

j.w_fil <- j.w %>% slice_max(CV,n=floor(nrow(j.w)*(1-0.05))) %>% slice_min(CV,n=floor(nrow(j.w)*(1-0.1)))
cor <- round(cor(j.w_fil[,2:18]),3)
 colnames(cor) <- gsub("2cell","L2C",colnames(cor))
 colnames(cor) <- gsub("4cell","4C",colnames(cor))
 colnames(cor) <- gsub("Zygote","L1C",colnames(cor))
 rownames(cor) <- gsub("2cell","L2C",rownames(cor))
 rownames(cor) <- gsub("4cell","4C",rownames(cor))
 rownames(cor) <- gsub("Zygote","L1C",rownames(cor))
 colnames(cor) <- gsub("-"," rep",colnames(cor))
 rownames(cor) <- gsub("-"," rep",rownames(cor))
p <- ggcorrplot(cor,tl.cex = 13,hc.order=TRUE)
p + scale_fill_gradient2(limit = c(0.3,1), low = "#6D9EC1", high =  "#E46726", mid = "white", 
    midpoint = 0.65,name="PCC(IP/input)\nstop_codon\n+-200bp")+
theme(legend.title = element_text(size = 12))
ggsave(file="PCC.ratio.pdf")


#enrich radar plot
a=read.table("peakEnrich.txt",header=T)
log <- a %>% mutate(enrich=log2(enrich))
a2=dcast(log,sample ~ element, value = 'enrich')
colnames(a2)[7] <- "stop"
colnames(a2)[8] <- "TSS"
ggradar(a2,grid.min = -1.2,grid.mid = 0,grid.max = 6,values.radar = c("-1.2", "0", "6"),
	plot.title="peak enrich at gene",group.line.width = 1, group.point.size = 2)+
scale_color_brewer(palette="Dark2",breaks=c("s1", "s2", "s3","s4","s5","s6"),
labels=c("GV", "MII", "Zygote","2cell","4cell","ESC"))
ggsave(file="peakEnrich.radar.pdf",width=6,height=6)


library(tidyverse)
library("ggalluvial")
library(reshape2)
library(ggthemes)
# cols <- c("#F8766D","#00BFC4") ##ggplot default 2 cols
cols <- c("#E46726","#6D9EC1")
myfun_plus <- function(a){b=a+0.0001;return(b)}
fp <- read.table(file="../../all.embryo.fpkm.txt")
colnames(fp)=c("gene","GV","MII","Zygote","cell2","cell4")
fp1 <- fp %>% filter_if(is.numeric,any_vars(. >1)) %>% mutate_if(is.numeric,myfun_plus) #filter any fpkm>1

files = list.files("../../genesWithExonPeak", pattern = "em.*.genesWithExonPeak.id.txt", full.names = TRUE) ##gene with  exon peak
peak <- map(files,read.table,header=T) %>% reduce(full_join,by="gene")
colnames(peak) <- c("gene","cell2.p","cell4.p","GV.p","MII.p","Zygote.p")
peak <- select(peak,"gene","GV.p","MII.p","Zygote.p","cell2.p","cell4.p")
peak[is.na(peak)] <- 0
peak_filter <- merge(fp1,peak)[,c(1,7:11)]

##alluvial plot for gene
peak_m <- melt(peak_filter)
peak_m2 <- melt(peak_m,measure.vars=c("gene"))
allu <- peak_m2[,c(1,2,4)]
colnames(allu) <- c("stage","value","gene")
allu <- allu %>% mutate(value=ifelse(value=="0","unmarked","marked"))
allu$stage <- gsub(".p","",allu$stage)
p <- ggplot(allu,aes(x = factor(stage,levels=c("GV","MII","Zygote","cell2","cell4")), stratum = value, alluvium = gene,
           fill = value, label = value)) +
  # scale_fill_manual(values = c("#F8766D","#00BFC4")) +
  scale_fill_manual(values = cols) +
  geom_flow(stat = "alluvium", lode.guidance = "forward",alpha= .1) +
  geom_stratum(color = "white",size=2) +
  theme(legend.position = "bottom") + theme_classic() 
p+theme(axis.text.y=element_text(size=25),axis.title.y = element_text(size = 25),
	axis.title.x = element_text(size = 25),axis.text.x= element_text(size = 25,angle=45,hjust=1),
	legend.text = element_text(size = 25),legend.title = element_text(size = 25),
	plot.title=element_text(size=25))+labs(title="m6a marked genes(fpkm>1) N=6940")+xlab("")
  ggsave(file="m6a.alluvial.gene.fpkm1.pdf",width=10, height=5)
#####

 ##gene gain or loss barplor
 dy <- peak_filter %>% mutate(MII=paste(GV.p,MII.p,sep=""),
 	Zygote=paste(MII.p,Zygote.p,sep=""),
 	cell2=paste(Zygote.p,cell2.p,sep=""),
 	cell4=paste(cell2.p,cell4.p,sep="")) %>% select(gene,MII:cell4)
 dy.long <- melt(dy,id.vars=c("gene")) %>% filter(value=="10" | value=="01")
 p <- ggplot(dy.long,aes(x=variable,fill=value)) +
 geom_bar(position=position_dodge(0.6),width=0.5)+ 
 scale_fill_manual(values=cols,breaks=c("01", "10"),labels=c("Gain", "Loss"))+theme_classic()+
 geom_text(stat='count',aes(label=..count..), color="black", size=3,position=position_dodge(0.5),vjust=-0.5)
 p+theme(axis.text.y=element_text(size=12),axis.title.y = element_text(size = 15),
	axis.title.x = element_text(size = 15),axis.text.x= element_text(size = 15,angle=45,hjust=1),
	legend.text = element_text(size = 15),legend.title=element_blank(),
	plot.title=element_text(size=15))+xlab("")+ylim(0,3000)
   ggsave(file="m6a.gain.or.loss.gene.fpkm1.pdf", width=5, height=5, dpi=300)

##dynamic peak analysis
files <- list.files(path=".",pattern="ref.peak.stat")
peak <- map(files,read.table,header=T) %>% reduce(full_join,by="id")
df <- peak %>% mutate_at(vars(-id), as.character) %>% select("id","GV","MII","Zygote","X2cell","X4cell")
df[is.na(df)] <- "unmarked"
allu <- melt(df,id.vars=c("id"))
colnames(allu) <- c("id","stage","value")
p <- ggplot(allu,aes(x = factor(stage,levels=c("GV","MII","Zygote","X2cell","X4cell")), 
	stratum = value, alluvium = id,fill = value, label = value)) +
  scale_fill_manual(values = cols) +
  geom_flow(stat = "alluvium", lode.guidance = "forward",alpha= .1) +
  geom_stratum(color = "white",size=2) +
  theme(legend.position = "bottom") + theme_classic() 
p+theme(axis.text.y=element_text(size=25),axis.title.y = element_text(size = 25),
	axis.title.x = element_text(size = 25),axis.text.x= element_text(size = 25,angle=45,hjust=1),
	legend.text = element_text(size = 25),legend.title = element_text(size = 25),
	plot.title=element_text(size=25))+labs(title="peaks N=43980")+xlab("")
ggsave(file="m6a.alluvial.peak.pdf",width=10,height=5,dpi=100)

 dy <- df %>% mutate(MII.dy=paste(GV,MII,sep=""),
 	Zygote.dy=paste(MII,Zygote,sep=""),
 	cell2.dy=paste(Zygote,X2cell,sep=""),
 	cell4.dy=paste(X2cell,X4cell,sep="")) %>% select(id,MII.dy:cell4.dy)
 colnames(dy) <- gsub(".dy","",colnames(dy))
 dy.long <- melt(dy,id.vars=c("id")) %>% filter(value=="markedunmarked" | value=="unmarkedmarked")
 p <- ggplot(dy.long,aes(x=variable,fill=factor(value,levels=c("unmarkedmarked","markedunmarked")))) +
 geom_bar(position=position_dodge(0.6),width=0.5)+ 
 scale_fill_manual(values=cols,breaks=c("unmarkedmarked", "markedunmarked"),labels=c("Gain", "Loss"))+theme_classic()+
 geom_text(stat='count',aes(label=..count..), color="black", size=3,position=position_dodge(0.5),vjust=-0.5)
 p+theme(axis.text.y=element_text(size=12),axis.title.y = element_text(size = 15),
	axis.title.x = element_text(size = 15),axis.text.x= element_text(size = 15,angle=45,hjust=1),
	legend.text = element_text(size = 15),legend.title=element_blank(),
	plot.title=element_text(size=15))+xlab("")#+ylim(0,3000)
   ggsave(file="m6a.gain.or.loss.peak.pdf", width=5, height=5, dpi=300)




