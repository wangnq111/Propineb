#package loading####
library(LorMe)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(vegan)
library(dplyr)
library(tidyr)
library(magrittr)

#data loading####
bac=read.csv("./input/bac.csv",header = T)
metafile=read.table("./input/metafile.txt",header = T,sep="\t")
my_col=c("Control"="#E31A1C","Propineb"="#377DB8")

#data processing####
Bactax_summary=tax_summary(groupfile = metafile,inputtable = bac[,2:25],reads = T,taxonomytable = bac[,c(1,26)],outputtax = "standard")
Bac_plan=object_config(taxobj = Bactax_summary,treat_location = 2,rep_location = 3,treat_col = my_col)

#community analysis####
#pcoa
structure=structure_plot(taxobj = Bac_plan,taxlevel = "Genus",diagram = "stick")
pcoa=structure$PCoA_Plot+labs(title="")
structure$PERMANOVA_statistics

#composition
community=community_plot(taxobj=Bac_plan,taxlevel = "Phylum",n = 10,palette = "Paired")
com_plot=community$barplot
community$alluvialplot # just view
#alpha
alpha=Alpha_diversity_calculator(taxobj = Bac_plan,taxlevel = "Genus")
alpha_plot=(alpha$plotlist$Plotobj_Shannon$Boxplot|alpha$plotlist$Plotobj_ACE$Boxplot)*guides(fill="none")

mean(alpha$alphaframe$Indexvalue[alpha$alphaframe$Indexname=="ACE"&alpha$alphaframe$Treatment=="Propineb"])/
  mean(alpha$alphaframe$Indexvalue[alpha$alphaframe$Indexname=="ACE"&alpha$alphaframe$Treatment=="Control"])

mean(alpha$alphaframe$Indexvalue[alpha$alphaframe$Indexname=="Shannon"&alpha$alphaframe$Treatment=="Propineb"])/
  mean(alpha$alphaframe$Indexvalue[alpha$alphaframe$Indexname=="Shannon"&alpha$alphaframe$Treatment=="Control"])

(((alpha_plot|pcoa+guides(fill="none"))+plot_layout(widths = c(1,1,2)))/com_plot)*guides(fill="none")
ggsave("./Figure3/community.pdf",width = 6,height = 6)

##detailed phyla comparsion####
inputdata<- subset(community$Grouped_Top10Phylum,tax != "Others")
for(i in unique(inputdata[,4])){
  sub_facet<-inputdata[which((inputdata[,4])==i),] 
  results<-auto_signif_test(data =sub_facet,treatment_col =2,value_col =5,prior = T)
  
}
mean_frame<- aggregate(inputdata[,5],by=list(inputdata[,2],inputdata[,4]),FUN=mean)
Sd<- aggregate(inputdata[,5],by=list(inputdata[,2],inputdata[,4]),FUN=sd)%>% .[,'x']
Treatment_Name<- mean_frame$Group.1
N<- table(inputdata[,2]) %>% as.numeric()
Mean<- mean_frame[,'x']
SEM<- Sd/(N^0.5)
input_mean_frame=data.frame(Treatment_Name,N,Mean,Sd,SEM,mean_frame[,'Group.2'])
colnames(input_mean_frame)[6]=colnames(inputdata)[4]
ggplot(input_mean_frame,aes(x=as.factor(Treatment_Name),y=Mean))+
  theme_zg()+
  geom_bar(stat = 'identity',size=0.2,width=0.5,aes(fill=as.factor(Treatment_Name)),color='#000000',alpha=0.8)+
  geom_errorbar(aes(ymin=Mean-Sd,ymax=Mean+Sd),size=0.2,width=0.2)+
  scale_y_continuous(expand = c(0.00,0.01))+
  scale_color_manual(values=my_col)+
  scale_fill_manual(values=my_col)+
  labs(x='',fill='',y='Relative abundance')+
  facet_wrap(~get(colnames(inputdata)[4]),scales ='free_y',strip.position ='top',as.table =T,ncol=5)+
  stat_compare_means(data=inputdata,aes(x=as.factor(inputdata[,2]),y=inputdata[,5],label = after_stat(p.signif)),label.x.npc = 0.5,label.y.npc  = 0.8,size=4,method = "wilcox.test")+
  theme(legend.position = 'right',)+
  guides(color='none',fill='none')
ggsave("./Supp_figure/Figire_S3.pdf",width = 8,height = 4)

#indicator analysis####
set.seed(999)
indicator_anno=indicator_analysis(taxobj =Bac_plan,taxlevel = "Genus",func = "r.g")
library(indicspecies)
library(permute)

dat <- as.data.frame(Bac_plan$Genus)
rownames(dat) <- make.unique(as.character(dat$Genus))
dat$Genus <- NULL
counts <- as.matrix(dat)
storage.mode(counts) <- "numeric"
counts <- counts[rowSums(counts) > 0, , drop = FALSE]
group <- factor(
  ifelse(
    grepl("^Control", colnames(counts)),
    "Control",
    "Propineb"
  ),
  levels = c("Control", "Propineb")
)
counts_pc <- counts + 30
log_counts <- log10(counts_pc)
clr <- sweep(
  log_counts,
  2,
  colMeans(log_counts),
  FUN = "-"
)
clr <- t(clr)
set.seed(999)
ind <- multipatt(
  x = clr,
  cluster = group,
  func = "r.g",
  duleg = TRUE,
  control = how(nperm = 999)
)
# 提取结果
res <- ind$sign
res$Genus <- rownames(res)
# BH校正
res$padj <- p.adjust(res$p.value, method = "BH")
indicator_anno1=left_join(indicator_anno[,6:14],res)
indicator_anno1$tag="None"
indicator_anno1$tag[indicator_anno1$s.Control==1&indicator_anno1$padj<0.05]="Control"
indicator_anno1$tag[indicator_anno1$s.Propineb==1&indicator_anno1$padj<0.05]="Propineb"
indicator_anno1$stat[indicator_anno1$s.Control==1]=-indicator_anno1$stat[indicator_anno1$s.Control==1]

#volcano plot
vol_plot=volcano_plot(inputframe =indicator_anno1,cutoff = 1,aes_col = my_col )
vol=vol_plot$FC_FDR+guides(fill="none")
#mannhaton plot
man_plot=manhattan(inputframe = indicator_anno1,taxlevel = "Phylum",control_name = "Control",mode = "select",palette =community$filled_color[-11] ,select_tax = names(community$filled_color)[-11])
man=man_plot$manhattan

write.table(indicator_anno1,file="./Supp_Table//Table_S3.txt",row.names = F,sep="\t")
##beneficial microbes compasrison####
library(tidyr)
select_beneficial=c("g__Pseudomonas","g__Sphingomonas","g__Streptomyces","g__Nocardioides","g__Brevundimonas","g__Flavobacterium")
select_beneficial_frame=Bac_plan$Genus_percent[Bac_plan$Genus_percent$Genus %in% select_beneficial,]

grouped_select_beneficial_frame=data.frame(Bac_plan$groupfile,t(select_beneficial_frame[,-1]))
colnames(grouped_select_beneficial_frame)[-c(1:3)]=select_beneficial_frame$Genus
grouped_select_beneficial_frame_long=gather(grouped_select_beneficial_frame,"Tax","Rel",-c(SAMPLE.ID,Treatment,Rep))
grouped_select_beneficial_frame_long$Tax=factor(grouped_select_beneficial_frame_long$Tax,
                                                levels=select_beneficial[c(2,1,3,4,5,6)])

beneficial_box=ggplot(grouped_select_beneficial_frame_long,aes(x=Treatment,y=Rel))+
  scale_y_continuous(expand = c(0.00,0.001))+
  scale_fill_manual(values=my_col)+
  theme_zg()+
  facet_wrap(~Tax,scales ='fixed',strip.position ='top',as.table =T,nrow=1)+
  geom_boxplot(aes(fill=Treatment),alpha=0.8,width=0.5,outlier.color =NA,color='#000000',linewidth=0.2)+
  labs(x='',fill='',y="Relative abundance")+
  stat_compare_means(aes(label = after_stat(p.signif)),size=4,label.x.npc = 0.5,label.y.npc = .9)+
  theme(legend.position = 'right')+guides(fill="none")

up=((vol|man)+plot_layout(widths = c(2,4)))
up/beneficial_box
ggsave("./Figure4/figure4.pdf",width = 7,height = 6)


#Kegg
kegg=read.csv("./input/tax4fun.csv",header = T)
kegg_tax=kegg[,c(1,26:28)]
kegg_matrix_pct=data.frame(kegg[,2:25],row.names = kegg$layer3)

##indicator genes####
library(magrittr)
indicator_gene<-multipatt(as.data.frame(t(kegg_matrix_pct)),Bac_plan$groupfile$Treatment,func = "r.g",control=how(nperm=1000)) %$%
  as.data.frame(sign) #%>% subset(.,.$p.value<0.05) 
indicator_gene_anno=left_join(data.frame(indicator_gene,layer3=rownames(indicator_gene)),kegg_tax)
indicator_gene_anno$threshold[indicator_gene_anno$s.Control==1]="Control"
indicator_gene_anno$threshold[indicator_gene_anno$s.Propineb==1]="Propineb"
indicator_gene_anno$threshold[indicator_gene_anno$p.value>0.05]="None"
indicator_gene_anno=indicator_gene_anno[!is.na(indicator_gene_anno$stat),]
indicator_gene_anno$stat[indicator_gene_anno$s.Control==1]=-indicator_gene_anno$stat[indicator_gene_anno$s.Control==1]
write.table(indicator_gene_anno,file="./Supp_Table//Table_S4.txt",row.names = F,sep="\t")
vol_gene=ggplot(indicator_gene_anno, aes(x = stat, y = -log10(p.value))) + 
  geom_point(aes(fill = threshold), pch=21,alpha = 0.8,show.legend = F) +  
  labs(x = "Enriched factor", y = "-Log10(q-value)",color = "") + 
  geom_hline(yintercept = -log10(0.05), lty = 2, color = "grey") + 
  geom_vline(xintercept = c(0), lty = 2, color = "grey") + 
  scale_fill_manual(values = c(my_col,"None"="gray")) + 
  theme_zg()+guides(color="none")

indicator_gene_anno_select=indicator_gene_anno[indicator_gene_anno$L3=="Metabolism",]
select_layer2= table(indicator_gene_anno_select$L2[indicator_gene_anno_select$s.Control==1]) %>% sort(.,decreasing=T) %>% names() %>% .[1:7]
indicator_gene_anno_select=subset(indicator_gene_anno_select,L2 %in% select_layer2)

indicator_gene_anno_select$L2=factor(indicator_gene_anno_select$L2,levels = select_layer2)

man_gene=ggplot(indicator_gene_anno_select, aes(x=L2, y=-log10(p.value))) +  #
  geom_hline(yintercept=-log10(0.05), linetype=2, color="lightgrey") +# 添加显著阈值线
  geom_point(aes(color=L2,shape=threshold),alpha=.6,position=position_jitter(0.5)) + #,stroke=2#绘制散点, position调整离散程度，stroke参数调整点的路径的粗细
  scale_shape_manual(values=c("Control"=25, "None"=20, "Propineb"=17))+#富集、下调、无显著性分别为空心上三角、实心下三角、圆
  #scale_size_manual(values=c("IIII"=2,"III"=1.5,"II"=1,"I"=0.5))+
  scale_color_manual(values=color_scheme("Plan8"))+
  labs(x=NULL, y="-log10(P)")+#删掉x轴标签
  theme_zg()+guides(color="none")+
  theme(legend.position="top", #图例置顶
        panel.grid = element_blank(), #去掉背景网格
        axis.text.x = element_text(angle = -45, hjust = 0.1, vjust = 0.1))#x轴坐标标签倾斜

(vol_gene|man_gene)+plot_layout(widths = c(2,5))
ggsave("./Supp_figure/Figire_S7.pdf",width = 7,height = 4)

##correlation analysis####
library(Hmisc)
library(pheatmap)
kegg_select=kegg_matrix_pct[rownames(kegg_matrix_pct) %in% indicator_gene_anno$layer3[indicator_gene_anno$threshold!="None"&indicator_gene_anno$L3=="Metabolism"],]
genus_all_percent=Bac_plan$Genus_percent
indicator_anno_select=indicator_anno1[indicator_anno1$Phylum %in% names(community$filled_color)[-11],]
genues_select=genus_all_percent[genus_all_percent$Genus %in% indicator_anno_select$Genus[indicator_anno_select$tag!="None"],]
genues_select=data.frame(genues_select[,-1],row.names = genues_select$Genus)

cor_object=rcorr(genues_select %>% as.matrix()%>% t(),t(kegg_select), "spearman")
cor_r=cor_object$r 
cor_r[cor_object$P>0.05]=0
heat_cor=cor_r[-c(1:nrow(genues_select)),1:nrow(genues_select)]
#filter
temp=function(x){max(abs(x))}
heat_cor1=heat_cor[which(apply(heat_cor,1,FUN=temp)>0.7),]
heat_cor1=heat_cor1[,which(apply(heat_cor,2,FUN=temp)>0.7)]

#anno
genus_anno=data.frame(Genus=colnames(heat_cor1)) %>% 
  left_join(.,indicator_anno[,c("Phylum","Genus","tag")])
genus_anno=data.frame(genus_anno[,-1],row.names = genus_anno$Genus)
genus_anno$Phylum[!genus_anno$Phylum %in% (table(genus_anno$Phylum) %>% sort()%>% names() %>% tail(5))]="Others"
kegg_anno=data.frame(layer3=rownames(heat_cor1))%>% 
  left_join(.,indicator_gene_anno[,c("layer3","L2","threshold")])
kegg_anno=data.frame(kegg_anno[,-1],row.names = kegg_anno$layer3)
kegg_anno$L2[!kegg_anno$L2 %in%(table(kegg_anno$L2) %>% sort()%>% names() %>% tail(5))]="Others"

L2_col=color_scheme("Plan7",6)
names(L2_col)=unique(kegg_anno$L2)
ann_color=list(Phylum=community$filled_color[unique(genus_anno$Phylum)],
               tag=my_col,
               threshold=my_col,
               L2=L2_col)
#visualization
pheatmap(heat_cor1[order(kegg_anno$threshold,kegg_anno$L2),
                   order(genus_anno$tag,genus_anno$Phylum)],
         show_rownames = F,show_colnames = F,
         annotation_col =genus_anno,annotation_row =kegg_anno,annotation_colors = ann_color,
         cluster_rows = F,cluster_cols = F,
         gaps_row=37,gaps_col = 68) #save 10*4,S7 bottom


