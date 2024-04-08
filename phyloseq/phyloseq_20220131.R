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

sample_shunxu <- rownames(sample_data_raw)

#otu_table
otu_table_raw <- read_excel("otu_table.xlsx")
otu_table_raw <- as.data.frame(otu_table_raw)
otu_table_raw <- column_to_rownames(otu_table_raw,var = "OTU ID")
otu_shunxu <- rownames(otu_table_raw)

otu_table_raw <- otu_table_raw[sample_shunxu]
otu_table_raw <- as.matrix(otu_table_raw)

#write.xlsx2(otu_table_raw, "otu_table_raw.xlsx")


#tax_table
tax_table_raw <- read.delim('percent.split_tax_rarefac.otu_taxa_table.xls', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE,na.strings = "")
dput(names(tax_table_raw))
tax_table_raw <- tax_table_raw[c("kingdom", 
                         "phylum", "class", "order", "family", "genus", "species")]
tax_table_raw <- tax_table_raw[otu_shunxu,]


tax_table_raw <- as.matrix(tax_table_raw)


#write.xlsx2(tax_table_raw, "tax_table_raw.xlsx")

#接下来，我们需要将他们组合成一个 phyloseq 对象：

#进行数据转化
OTU <-  otu_table(otu_table_raw, taxa_are_rows=TRUE)
TAX <-  tax_table(tax_table_raw)
SAMPLE <- sample_data(sample_data_raw)
# 使用phyloseq函数创建一个phyloseq对象。
my_ph <-  phyloseq(OTU, TAX,SAMPLE)
my_ph
view(otu_table(my_ph))
view(sample_data(my_ph))
view(tax_table(my_ph))

##############################################################################################
#计算alpha多样性 #使用原始OTU, ?estimate_richness

alpha_output <-  estimate_richness(my_ph, measures=c("Shannon", "Simpson", "Chao1", "ACE", "Observed", "InvSimpson", "Fisher"))



#####################################
#绘制稀疏曲线  #使用原始OTU

library(amplicon)
library(ggalt)
?alpha_rare_all
result <-  alpha_rare_all(ps = my_ph, group = "group", method = "observed", start = 0, step = 1000)#step减小到100会报错，超出画图范围了
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

#rarefied_otu_table <- round(rarefied_otu_table/100,0)#因为抽样100次，所以除以100，但公司并未除
#write.xlsx2(rarefied_otu_table,"rarefied_otu_table.xlsx")
OTU_rarefied <-  otu_table(rarefied_otu_table, taxa_are_rows=TRUE)
my_ph_rarefied <-  phyloseq(OTU_rarefied, TAX, SAMPLE)
my_ph_rarefied
sample_sums(my_ph_rarefied)


######PCA计算   
#使用rarefied_otu_table
#先进行hellinger转化，注意是对每个样本转化，margin=2, ?decostand
otu_table_hell <-  decostand(rarefied_otu_table,"hellinger",MARGIN = 2)
#使用转化后的物种数据执行 PCA，无需标准化,?rda 注意：主成分分析对数据的要求是：行是标本，列是参数
otu_table_hell <- t(otu_table_hell)
pca_otu <- rda(otu_table_hell, scale = FALSE)
pca_otu_sum <- summary(pca_otu)
str(pca_otu_sum)
pca_otu_sum$cont$importance["Proportion Explained",c("PC1","PC2")]#PC1,PC2所占比例

#提取对象排序坐标
site.scaling <- pca_otu_sum$sites[ ,c("PC1","PC2")]%>%as.data.frame %>% rownames_to_column(var = "sample")

site.scaling


#添加分组信息
site.fenzu <- sample_data_raw%>%rownames_to_column(var = "sample")
site.scaling <- site.scaling%>%left_join(site.fenzu,by="sample")


#ggplot2 作图

ggplot(site.scaling, aes(PC1, PC2)) +
  geom_point(aes(color = group)) +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 	#添加置信椭圆，注意不是聚类
  labs(x = "PC1 [27.1%]", y = "PC2 [12.4%]") +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  theme_bw()


##############################################################################################
##############################################################################################
##基于phyloseq的微生物群落距离的计算
####基于距离的差异检验方法
####注意group的顺序十分重要，必须和距离中样本的顺序一致

d_bray_matrix <-read.delim('bray_curtis_dm.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
d_bray_matrix <- d_bray_matrix[sample_shunxu,sample_shunxu]
d_bray <-as.dist(d_bray_matrix)

d_unifrac_matrix <-read.delim('unweighted_unifrac_dm.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
d_unifrac_matrix <- d_unifrac_matrix[sample_shunxu,sample_shunxu]
d_unifrac <-as.dist(d_unifrac_matrix)

d_wunifrac_matrix <-read.delim('weighted_unifrac_dm.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
d_wunifrac_matrix <- d_wunifrac_matrix[sample_shunxu,sample_shunxu]
d_wunifrac <-as.dist(d_wunifrac_matrix)

#d_bray
set.seed(1)
adonis_result <- adonis(d_bray~group, sample_data_raw, permutations = 999)
adonis_result
anosim_result <- anosim(d_bray, sample_data_raw$group, permutations = 999)
anosim_result
summary(anosim_result)
mrpp_result <- mrpp(d_bray, sample_data_raw$group, permutations = 999)
mrpp_result

#d_unifrac
set.seed(1)
adonis_result <- adonis(d_unifrac~group, sample_data_raw, permutations = 999)
adonis_result
anosim_result <- anosim(d_unifrac, sample_data_raw$group, permutations = 999)
anosim_result
summary(anosim_result)
mrpp_result <- mrpp(d_unifrac, sample_data_raw$group, permutations = 999)
mrpp_result

#d_wunifrac
set.seed(1)
adonis_result <- adonis(d_wunifrac~group, sample_data_raw, permutations = 999)
adonis_result
anosim_result <- anosim(d_wunifrac, sample_data_raw$group, permutations = 999)
anosim_result
summary(anosim_result)
mrpp_result <- mrpp(d_wunifrac, sample_data_raw$group, permutations = 999)
mrpp_result


####PCoA ，详情 ?cmdscale
PCoA_unifrac <- cmdscale(d_unifrac, k = (nrow(d_unifrac_matrix) - 1),eig = TRUE)
PCoA_unifrac$eig
PCoA_unifrac$points
##查看结果
#各 PCoA 轴的特征值
PCoA_unifrac_eig <- PCoA_unifrac$eig

PCoA_unifrac_exp <- PCoA_unifrac$eig/sum(PCoA_unifrac$eig)
#查看PC1，PC2所占比例
PCoA_unifrac_exp[1:2]#写到lab里
#先评估下负特征值（末尾几个轴）
barplot(PCoA_unifrac_eig)

#如果负特征值影响甚微，则继续
#如果负特征值非常明显，则应当先校正

##ggplot2 作图
#添加分组信息
site.fenzu <- sample_data_raw%>%rownames_to_column(var = "sample")
PCoA_unifrac_site <- PCoA_unifrac$point[,c(1,2)]%>%as.data.frame %>% rownames_to_column(var = "sample")%>%left_join(site.fenzu,by="sample")

#ggplot2 作图

ggplot(PCoA_unifrac_site, aes(V1, V2)) +
  geom_point(aes(color = group)) +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 	#添加置信椭圆，注意不是聚类
  labs(x = "PC1 [17.5%]", y = "PC2 [10.8%]") +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  theme_bw()



##unifrac距离 NMDS
#NMDS 排序，定义 2 个维度，详情 ?metaMDS
NMDS_unifrac <- metaMDS(d_unifrac, k = 2)
#应力函数值，一般不大于 0.2 为合理
NMDS_unifrac$stress

#样方得分
NMDS_unifrac_site <- data.frame(NMDS_unifrac$points)

#添加分组信息
site.fenzu <- sample_data_raw%>%rownames_to_column(var = "sample")

NMDS_unifrac_site <- NMDS_unifrac_site %>% rownames_to_column(var = "sample")%>%left_join(site.fenzu,by="sample")

#ggplot2 作图

ggplot(NMDS_unifrac_site, aes(MDS1, MDS2)) +
  geom_point(aes(color = group)) +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 	#添加置信椭圆，注意不是聚类
  geom_vline(xintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  theme_bw()


#########使用weighted_unifrac计算

####PCoA ，详情 ?cmdscale
PCoA_wunifrac <- cmdscale(d_wunifrac, k = (nrow(d_wunifrac_matrix) - 1),eig = TRUE)
PCoA_wunifrac$eig
PCoA_wunifrac$points
##查看结果
#各 PCoA 轴的特征值
PCoA_wunifrac_eig <- PCoA_wunifrac$eig

PCoA_wunifrac_exp <- PCoA_wunifrac$eig/sum(PCoA_wunifrac$eig)
#查看PC1，PC2所占比例
PCoA_wunifrac_exp[1:2]#写到lab里
#先评估下负特征值（末尾几个轴）
barplot(PCoA_wunifrac_eig)

#如果负特征值影响甚微，则继续
#如果负特征值非常明显，则应当先校正

##ggplot2 作图
#添加分组信息
site.fenzu <- sample_data_raw%>%rownames_to_column(var = "sample")
PCoA_wunifrac_site <- PCoA_wunifrac$point[,c(1,2)]%>%as.data.frame %>% rownames_to_column(var = "sample")%>%left_join(site.fenzu,by="sample")

#ggplot2 作图

ggplot(PCoA_wunifrac_site, aes(V1, V2)) +
  geom_point(aes(color = group)) +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 	#添加置信椭圆，注意不是聚类
  labs(x = "PC1 [48.5%]", y = "PC2 [21.9%]") +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  theme_bw()



##wunifrac距离 NMDS
#NMDS 排序，定义 2 个维度，详情 ?metaMDS
NMDS_wunifrac <- metaMDS(d_wunifrac, k = 2)
#应力函数值，一般不大于 0.2 为合理
NMDS_wunifrac$stress

#样方得分
NMDS_wunifrac_site <- data.frame(NMDS_wunifrac$points)

#添加分组信息
site.fenzu <- sample_data_raw%>%rownames_to_column(var = "sample")
NMDS_wunifrac_site <- NMDS_wunifrac_site %>% rownames_to_column(var = "sample")%>%left_join(site.fenzu,by="sample")

#ggplot2 作图

ggplot(NMDS_wunifrac_site, aes(MDS1, MDS2)) +
  geom_point(aes(color = group)) +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 	#添加置信椭圆，注意不是聚类
  geom_vline(xintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5,linetype=2) +
  theme_bw()



#转化为相对丰度
my_ph_relative <- transform_sample_counts(my_ph, function(x) x / sum(x) )

group_compare_table <-  tax_table(my_ph_relative) %>% data.frame%>% rownames_to_column(var = "OTU")%>%
left_join((otu_table(my_ph_relative) %>% data.frame%>% rownames_to_column(var = "OTU")),by="OTU")


group_compare_table_phylum <- group_compare_table %>% select(phylum,TA1:HC160)%>%group_by(phylum) %>%
  summarise(across(TA1:HC160, ~ sum(.x)))%>%
  column_to_rownames(var = "phylum")%>%t%>%data.frame%>% rownames_to_column(var = "sample")%>%
  left_join(sample_data(my_ph_relative)%>%data.frame%>% rownames_to_column(var="sample"),by="sample")%>%
  select(-c(sample,sex,age,BMI,dis_val))

#把相对丰度之和<0.1的合并为others
others_names <- names(which(group_compare_table_phylum%>%select(-group)%>%colSums<0.1))
others <- rowSums(group_compare_table_phylum[others_names])
group_compare_table_phylum <- group_compare_table_phylum%>%mutate(others=others)%>%select(-others_names)

###画图
group_compare_table_phylum <- group_compare_table_phylum%>%pivot_longer(cols= names(group_compare_table_phylum%>%select(-group)),names_to="phylum",values_to="value")

  ggplot(group_compare_table_phylum,aes(phylum,value,fill=group))+
  geom_boxplot()



##计算统计值
gvar <- colnames(group_compare_table_phylum%>%select(group))
convars <- colnames(group_compare_table_phylum%>%select(-group))
mydata <- group_compare_table_phylum
convartest_p <- c()#p值
convartest_s <- c()#s值
for(i in convars){
  k <- formula(paste(i,"~",gvar))
  suanfa_3 <- wilcox.test(k,data = mydata,exact=FALSE)
  convartest_p[i] <- suanfa_3$p.value 
  convartest_s[i] <- suanfa_3$statistic 
}

group_compare_table_phylum_output <- data.frame(convartest_p,convartest_s)
group_compare_table_phylum_output%>%arrange(-desc(convartest_p)) %>%mutate(across(convartest_p,~p.adjust(.x,method = "fdr"),.names = "q_value"))


#random_forest
otutable_rf <- otu_table(my_ph)%>%t%>%data.frame
sampledata_rf <- sample_data(my_ph)%>%data.frame
sampledata_rf <- sampledata_rf[rownames(otutable_rf),]
sampledata_rf$group <- sampledata_rf$group%>%factor(ordered = T,levels = c("HC","TA"))
otutable_rf$group <- sampledata_rf$group
sampledata_rf_dis <- subset(sampledata_rf,dis_val=="dis")
sampledata_rf_val <- subset(sampledata_rf,dis_val=="val")
otutable_rf_dis <- otutable_rf[rownames(sampledata_rf_dis),]
otutable_rf_val <- otutable_rf[rownames(sampledata_rf_val),]

library(randomForest)
set.seed(315)
rf_train <- randomForest(group~.,
                         data = otutable_rf_dis,
                         ntree=1000,
                         importance=TRUE,
                         proximity=TRUE)
rf_train
plot(rf_train)

imp <- data.frame(rf_train$importance)
imp = imp[order(imp$MeanDecreaseAccuracy,decreasing = T),]
head(imp)#此处“Mean Decrease Accuracy”和“Mean Decrease Gini”为随机森林模型中的两个重要指标。
#其中，“mean decrease accuracy”表示随机森林预测准确性的降低程度，该值越大表示该变量的重要性越大；
#“mean decrease gini”计算每个变量对分类树每个节点上观测值的异质性的影响，从而比较变量的重要性。该值越大表示该变量的重要性越大

varImpPlot(rf_train)

#ggplot 美化贡献度绘图, 筛选前20个用于绘图
imp_sub=imp[1:20,]
imp_sub$taxa<-rownames(imp_sub)
# 添加因子保持顺序
imp_sub$taxa=factor(imp_sub$taxa,order=T,levels = rev(imp_sub$taxa))
p=ggplot(data = imp_sub, mapping = aes(x=taxa,y=MeanDecreaseAccuracy)) +
  geom_bar(stat="identity")+coord_flip()+theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        panel.border=element_rect(colour = "black"),
        legend.key = element_blank(),plot.title = element_text(hjust = 0.5))
p

#交叉验证确定最佳预测微生物数量：
#英文：10-fold cross validation：
#用来测试算法准确性。是常用的测试方法。将数据集分成十份，
#轮流将其中9份作为训练数据，1份作为测试数据，进行试验。每次试验都会得出相应的正确率（或差错率）。
#10次的结果的正确率（或差错率）的平均值作为对算法精度的估计，
#一般还需要进行多次10折交叉验证（例如10次10折交叉验证），再求其均值，作为对算法准确性的估计。
#之所以选择将数据集分为10份，是因为通过利用大量数据集、使用不同学习技术进行的大量试验，
#表明10折是获得最好误差估计的恰当选择，而且也有一些理论根据可以证明这一点。但这并非最终诊断，争议仍然存在。而且似乎5折或者20折与10折所得出的结果也相差无几。

#这里采用10 fold，重复5次
set.seed(315)
result <-  replicate(5, rfcv(otutable_rf_dis[,-ncol(otutable_rf_dis)],otutable_rf_dis$group,cv.fold=10), simplify=FALSE)
error.cv <-  sapply(result, "[[", "error.cv")
matplot(result[[1]]$n.var, cbind(rowMeans(error.cv), error.cv), type="l",
        lwd=c(2, rep(1, ncol(error.cv))), col=1, lty=1, log="x",
        xlab="Number of variables", ylab="CV Error")
head(error.cv)
cv_result <- data.frame(error.cv)%>%rownames_to_column(var = "number_of_variables")%>% as_tibble()%>%mutate(mean_error_rate=rowMeans(error.cv))%>% rowwise()%>%mutate(sd = sd(across(starts_with("X"))))


#确定用于模型构建的数量后，
#再优化模型选择最佳mtry，也就是每课树中的预测变量数量
#这里选择前10个重要的变量，并优化选择最佳mtry

rf.formula <- as.formula(paste0("group ~",paste(row.names(imp)[1:10],collapse="+")))
rf.formula

#校准模型mtry
library(caret)
?train
train.res <- train(rf.formula,
                data = otutable_rf_dis, # Use the train data frame as the training data
                method = 'rf',# Use the 'random forest' algorithm
                trControl = trainControl(method='repeatedcv',
                                         number=10,
                                         repeats=5,
                                         search='grid')) # Use 10 folds for cross-validation 重复5次
train.res

?randomForest
model <- randomForest(rf.formula, # 新建的模型
                      data=otutable_rf_dis, # 使用训练集构建随机森林
                      ntree = 500, # 决策树的数量，默认的是500
                      mtry = 2, # 每个分组中随机的变量数，默认是变量数开根
                      importance = TRUE, # 是否评估预测变量的重要性
                      proximity = TRUE) # 是否计算行之间的接近度

model

#对dis组进行预测
pred1 <- predict(model, newdata = otutable_rf_dis,type="response")
confusionMatrix(pred1,otutable_rf_dis$group)

pred2 <- predict(model, newdata = otutable_rf_dis,type="vote")
library(pROC)
roc.info<-roc(otutable_rf_dis$group,
              pred2[,2], #提取随机森林模型中对应的预测指标
              plot=TRUE,
              legacy.axes=TRUE,
              percent=FALSE,
              xlab="False positive percentage",
              ylab="True postive percentage",
              col="#4daf4a",
              lwd=4,
              print.auc=TRUE)

#对val组进行预测
pred1 <- predict(model, newdata = otutable_rf_val,type="response")
confusionMatrix(pred1,otutable_rf_val$group)

pred2 <- predict(model, newdata = otutable_rf_val,type="vote")
library(pROC)
roc.info<-roc(otutable_rf_val$group,
              pred2[,2], #提取随机森林模型中对应的预测指标
              plot=TRUE,
              legacy.axes=TRUE,
              percent=FALSE,
              xlab="False positive percentage",
              ylab="True postive percentage",
              col="#4daf4a",
              lwd=4,
              print.auc=TRUE)

