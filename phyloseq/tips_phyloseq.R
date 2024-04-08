plot_bar(my_ph_raw, fill = "phylum")
plot_tree(my_ph, color="group", label.tips="taxa_names", ladderize="left", plot.margin=0.3)
plot_heatmap(my_ph)

#查看phyloseq的参数
my_ph
#查看OTU有多少种
ntaxa(my_ph)
#查看样品数量
nsamples(my_ph)

#查看样品名称
sample_names(my_ph)
#查看OTU名称
taxa_names(my_ph)
#查看tax的分类等级信息
rank_names(my_ph)
#查看样品参数名称
sample_variables(my_ph)
#查看每种OTU数量总和
taxa_sums(my_ph)
#查看每个样品的OTU总和
sample_sums(my_ph)
min(sample_sums(my_ph))

#
otu_table(my_ph)
tax_table(my_ph)
sample_data(my_ph)
phy_tree(my_ph)#进化树


#
get_taxa(my_ph,1:2)
get_sample (my_ph,1:2)
get_variable(my_ph,1:2)


#按照丰度提取前十个OTU，并可视化进化树
myTaxa <-  names(sort(taxa_sums(my_ph), decreasing = TRUE)[1:10])
ex1 <-  prune_taxa(myTaxa, my_ph)
plot(phy_tree(ex1), show.node.label = TRUE)
plot_tree(ex1, color = "group", label.tips = "phylum", ladderize = "left", justify = "left" , size = "Abundance")





#数据预处理

#转化OTUcount数为相对丰度
ph_r  <-  transform_sample_counts(my_ph, function(x) x / sum(x) )
otu_table(ph_r)[1:5][1:5]
#提取丰度大于十万份之一的OTU
ph_fr <-  filter_taxa(ph_r, function(x) mean(x) > 1e-5, TRUE)
ntaxa(ph_fr)
#提取指定分类的OTU
ph_ch1 <-  subset_taxa(my_ph, phylum=="p__Bacteroidota")
ntaxa(ph_ch1)
nsamples(ph_ch1)
min(sample_sums(ph_ch1))

#提取总OTU数量大于500的样品
ph_ch1 <-  prune_samples(sample_sums(ph_ch1)>=500, ph_ch1)
nsamples(ph_ch1)
#合并指定OTU：taxa_names(GP.chl)[1:5]为一个OTU
ph_ch1_merged <- merge_taxa(ph_ch1, taxa_names(ph_ch1)[1:5])


#选取【OTU数量>3的样品个数/总样品个数】大于20%的OTU,注意编写函数的形式
ph <-  filter_taxa(my_ph, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
ntaxa(ph)
#对样品分组文件mapping添加一列
sample_data(ph)$dis <- factor( get_variable(ph, "dis_val") %in% c("dis") )


#变异系数过滤OTU，减少一些意外OTU误差
ph_sf <- filter_taxa(ph_s, function(x) sd(x)/mean(x) > 3.0, TRUE)
#提取指定门类OTU
ph_sfb <- subset_taxa(ph_sf, phylum=="p__Bacteroidota")
title <- "plot_bar; Bacteroidetes-only"

plot_bar(ph_sfb, "group", "Abundance", title=title)
plot_bar(ph_sfb, "group", "Abundance","family", title=title)
plot_bar(ph_sfb, "group", "Abundance","family", title=title,facet_grid = "sex")


#我们进行查看发现有很多的距离算法
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

#需要tree的距离算法
dist_methods[1:3]


删除需要tree的距离算法
# Remove them from the vector
dist_methods <- dist_methods[-(1:3)]
# 删除需要用户自定义的距离算法
dist_methods["designdist"]
dist_methods <- dist_methods[-which(dist_methods=="ANY")]
?ordinate

plist <- vector("list", length(dist_methods))


names(plist) <- dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(my_ph, method=i)
  # Calculate ordination
  iPCoA  <- ordinate(my_ph, "PCoA", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(my_ph, iPCoA, color="group", shape="dis_val")
  # Add title to each plot
  p <- p + ggtitle(paste("PCoA using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] <- p
}

df <- ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
unique(df$distance)
p <- ggplot(df, aes(Axis.1, Axis.2, color=group, shape=dis_val))
p <- p + geom_point(size=3, alpha=0.5)
p <- p + facet_wrap(~distance, scales="free")
p <- p + ggtitle("PCoA on various distance metrics for my_ph dataset")
p

?ldply
#Bray-Curtis
print(plist[["bray"]])

##############################################################################################
#计算alpha多样性


p <- plot_richness(my_ph, x="group", measures=c("shannon"))
p$data

p <- plot_richness(my_ph, x="group", measures=c("simpson"))
p$data

p <- plot_richness(my_ph, x="group", measures=c("Chao1"))
p$data

p <- plot_richness(my_ph, x="group", measures=c("ACE"))
p$data

p <- plot_richness(my_ph, x="group", measures=c("Observed"))
p$data