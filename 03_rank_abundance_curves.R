#rank- abundance curve
library(tidyverse)
library(BiodiversityR)

####导入数据

otu_table_raw <-as.data.frame( read_excel("otu_table_raw.xlsx"))
otu_table_raw <- column_to_rownames(otu_table_raw,var = "...1")

sample_data_raw <-as.data.frame( read_excel("sample_data_raw.xlsx"))
names(sample_data_raw)[1] <- "sample"



plot_otu <- t(otu_table_raw)

##Rank-abundance 曲线
#统计（BiodiversityR 包 rankabundance() 实现 OTU 排序）
otu_relative <- plot_otu / rowSums(plot_otu)
rank_dat <- data.frame()
for (i in rownames(otu_relative)) {
  rank_dat_i <- data.frame(rankabundance(subset(otu_relative, rownames(otu_relative) == i), digits = 6))[1:2]
  rank_dat_i$sample <- i
  rank_dat <- rbind(rank_dat, rank_dat_i)
}
rank_dat <- subset(rank_dat, rank_dat$abundance != 0)

rank_dat <- left_join(rank_dat,sample_data_raw,by="sample")

#ggplot2 作图
ggplot(rank_dat, aes(rank, log(abundance, 10), group = sample,color=group)) +
  geom_line() +
  labs(x = 'OTUs rank', y = 'Relative adundance (%)', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  scale_y_continuous(breaks = 0:-5, labels = c('100', '10', '1', '0.1', '0.01', '0.001'), limits = c(-5, 0))
