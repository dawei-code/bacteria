#* sample_data ：一个data.frame，包含了所有样本的表型信息，行名必须匹配otu_table 中的样本名；

#* otu_table ：一个数字矩阵 matrix，包含了 OTU 在每个样本中的丰度信息；

#* tax_table ：一个字符矩阵 matrix，包含了 OTU 的物种信息，行名必须匹配otu_table 中的 OTU 名。
library(tidyverse)
library(readxl)
library(phyloseq)
library(xlsx)
library(ape)
library(vegan)


#sample_data
X80_TA_chaifen <- read_excel("X80_TA_chaifen.xlsx")
X160_HC_chaifen <- read_excel("X160_HC_chaifen.xlsx")

sample_data_raw <- rbind(X80_TA_chaifen[,1:6],X160_HC_chaifen)
sample_data_raw <- column_to_rownames(sample_data_raw,var = "ID")

#write.xlsx2(sample_data_raw, "sample_data_raw.xlsx")

#sample_shunxu <- rownames(sample_data_raw)

#otu_table
otu_table_raw <- read_excel("otu_table.xlsx")
otu_table_raw <- as.data.frame(otu_table_raw)
otu_table_raw <- column_to_rownames(otu_table_raw,var = "OTU ID")
#otu_shunxu <- rownames(otu_table_raw)

#otu_table_raw <- otu_table_raw[sample_shunxu]
otu_table_raw <- as.matrix(otu_table_raw)

#write.xlsx2(otu_table_raw, "otu_table_raw.xlsx")


#tax_table
tax_table_raw <- read.delim('percent.split_tax_rarefac.otu_taxa_table.xls', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE,na.strings = "")
dput(names(tax_table_raw))
tax_table_raw <- tax_table_raw[c("kingdom", 
                         "phylum", "class", "order", "family", "genus", "species")]
#tax_table_raw <- tax_table_raw[otu_shunxu,]


tax_table_raw <- as.matrix(tax_table_raw)


#write.xlsx2(tax_table_raw, "tax_table_raw.xlsx")

#接下来，我们需要将他们组合成一个 phyloseq 对象：

#进行数据转化
OTU <-  otu_table(otu_table_raw, taxa_are_rows=TRUE)
TAX <-  tax_table(tax_table_raw)
SAMPLE <- sample_data(sample_data_raw)
# 使用phyloseq函数创建一个phyloseq对象。
my_ph_raw <-  phyloseq(OTU, TAX)
my_ph_raw

#使用ape包建立 OTU 系统发育树并导入 phyloseq 对象：?rtree
#注意下面这个tree的生成是随机的，每次都不一样.rooted=true时生成的是有根树，根据有根树计算距离时，每次距离结果不会变。根据无根树计算距离时，每次距离会变。
#下面这个树不是系统发生树，系统发生树需要另外计算生成
TREE <-  rtree(ntaxa(my_ph_raw), rooted=TRUE, tip.label=taxa_names(my_ph_raw))

#我们有了otu_table，sample_data，tax_table，phy_tree这四类数据，
#可以使用merge_phyloseq函数在之前创建的phyloseq对象中加入sample_data和phy_tree数据；
my_ph <-  merge_phyloseq(my_ph_raw, SAMPLE, TREE)
my_ph

##############################################################################################
##############################################################################################

#  ?`phyloseq-package`
#  ?rarefy_even_depth
#测序深度进行统一,即抽平
nr <- nrow(otu_table_raw)
nc <- ncol(otu_table_raw)
rarefied_otu_table <- matrix(rep(0,nr*nc),nrow = nr)
for( i in 1:100){
rarefied_my_ph <- rarefy_even_depth(my_ph,rngseed = i,sample.size = min(sample_sums(my_ph)), trimOTUs = FALSE)#抽平数为拥有otu数量最少的样本的otu数量
rarefied_otu_table <- otu_table(rarefied_my_ph) + rarefied_otu_table
}

rarefied_otu_table <- round(rarefied_otu_table/100,0)
#write.xlsx2(rarefied_otu_table,"rarefied_otu_table.xlsx")
OTU_rarefied <-  otu_table(rarefied_otu_table, taxa_are_rows=TRUE)
my_ph_raw_rarefied <-  phyloseq(OTU_rarefied, TAX)
TREE_rarefied <-  rtree(ntaxa(my_ph_raw_rarefied), rooted=TRUE, tip.label=taxa_names(my_ph_raw_rarefied))
my_ph_rarefied <-  merge_phyloseq(my_ph_raw_rarefied, SAMPLE, TREE_rarefied)
my_ph_rarefied
sample_sums(my_ph_rarefied)



##############################################################################################
##############################################################################################
##基于phyloseq的微生物群落距离的计算
###############################################################
#my_ph_rarefied
#unifrac距离 ?distance
d_unifrac_rarefied <- distance(my_ph_rarefied,method = "unifrac")

# ordinate ?ordinate
##unifrac距离 PCoA
PCoA_unifrac <- ordinate(my_ph_rarefied,"PCoA",distance = d_unifrac_rarefied)
str(PCoA_unifrac)
plot_ordination(my_ph_rarefied, PCoA_unifrac,axes = c(1,2), color="group") +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 	#添加置信椭圆，注意不是聚类
  geom_vline(xintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  theme_bw()

##unifrac距离 NMDS
NMDS_unifrac <- ordinate(my_ph_rarefied,"NMDS",distance = d_unifrac_rarefied)
str(NMDS_unifrac)
NMDS_unifrac$stress
plot_ordination(my_ph_rarefied, NMDS_unifrac,axes = c(1,2), color="group") +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 	#添加置信椭圆，注意不是聚类
  geom_vline(xintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  theme_bw()

######
#wunifrac距离 ?distance
d_wunifrac_rarefied <- distance(my_ph_rarefied,method = "wunifrac")

# ordinate ?ordinate
##wunifrac距离 PCoA
PCoA_wunifrac <- ordinate(my_ph_rarefied,"PCoA",distance = d_wunifrac_rarefied)
str(PCoA_wunifrac)
plot_ordination(my_ph_rarefied, PCoA_wunifrac,axes = c(1,2), color="group") +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 	#添加置信椭圆，注意不是聚类
  geom_vline(xintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  theme_bw()

##wunifrac距离 NMDS
NMDS_wunifrac <- ordinate(my_ph_rarefied,"NMDS",distance = d_wunifrac_rarefied)
str(NMDS_wunifrac)
plot_ordination(my_ph_rarefied, NMDS_wunifrac,axes = c(1,2), color="group") +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 	#添加置信椭圆，注意不是聚类
  geom_vline(xintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  theme_bw()

#PERMANOVA检验整体差异
##PERMANOVA 分析（所有分组间比较，即整体差异）
#根据 group这一列分组进行 PERMANOVA 分析检验组间差异，基于 999 次置换，详情 ?adonis
#现有的距离矩阵，这里为 unifrac距离

adonis_result <- adonis(d_unifrac_rarefied~group, sample_data_raw, permutations = 9999)
adonis_result
str(adonis_result)
summary(adonis_result)
adonis_result$aov.tab



######PCA计算   
#使用rarefied_otu_table
#先进行hellinger转化，注意是对每个样本转化，margin=2, ?decostand
otu_table_hell <-  decostand(rarefied_otu_table,"hellinger",MARGIN = 2 )
OTU_hell <-  otu_table(otu_table_hell, taxa_are_rows=TRUE)
# 使用phyloseq函数创建一个phyloseq对象。
my_ph_hell <-  phyloseq(OTU_hell, TAX,SAMPLE)
my_ph_hell

##PCA,使用vegan中的rda函数,不用距离参数
PCA <- ordinate(my_ph_hell,"RDA",scale=FALSE)
str(PCA)
summary(PCA)
plot_ordination(my_ph_hell, PCA,axes = c(1,2), color="group") +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 	#添加置信椭圆，注意不是聚类
  geom_vline(xintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  theme_bw()

##############################################################################################
#计算alpha多样性 #使用原始OTu

alpha_output <-  estimate_richness(my_ph, measures=c("Shannon", "Simpson", "Chao1", "ACE", "Observed", "InvSimpson", "Fisher"))



#####################################
#绘制稀疏曲线  #使用原始OTu

library(amplicon)
library(ggalt)
?alpha_rare_all
result <-  alpha_rare_all(ps = my_ph, group = "group", method = "observed", start = 1000, step = 1000)
p <- result[[4]]
result$table
ggplot(result$table,aes(i,index,fill=ID))+
  geom_xspline(aes(color=Group))+
  guides(fill="none")


###绘制Venn图  #使用原始OTu
library(MiscMetabar)
venn_phyloseq(my_ph,fact="group",print_values = TRUE)

p <- venn_phyloseq(my_ph,fact="group",print_values = FALSE)

str(p)





