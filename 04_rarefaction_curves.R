#rarefaction curve 稀释曲线（丰富度曲线）



# 同上，使用read.delim函数导入TSV文件
tsv_data <- read.delim("rarefaction_plot_group_data.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

dput(names(tsv_data))
#ggplot2 作图
ggplot(tsv_data, aes(num, count, group = X.SampleID,color=group)) +
  
  labs(x = 'Number of sequences', y = 'Number of OTUs', color = NULL) +
  geom_line()+
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
    scale_x_continuous(breaks = seq(0, 100000, 10000), labels = as.character(seq(0, 100000, 10000)))
