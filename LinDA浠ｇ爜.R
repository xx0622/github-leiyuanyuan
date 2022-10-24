# LinDA 计算差异菌 ------------------------------------
library(tidyverse)
library(magrittr)
library(LinDA)
library(readxl)

# meta2 ---------------------------------------------
# 读入分组文件
gp <- read_xlsx('meta2/2_sample_metadata.xlsx') %>% 
  rename(seqid = index)

# 读入丰度文件
df <-  read_tsv('meta2/2_all_mpa_to_abundance_metagenome.tsv') %>% 
  rename(name = index) %>% 
  select(name,all_of(gp[['seqid']])) 
# 转换为数据框
df %<>% data.frame(row.names = 1)
# 菌名替换为简短名
anno = data.frame(index = 1:length(rownames(df)),anno = rownames(df))
# 去除行名
rownames(df) = NULL
# 分组文件改为数据框
gp %<>% data.frame(row.names = 'seqid')

# 开始LinDA
meta <- gp %>% dplyr::select('group') %>% filter(!is.na(group)) %>% mutate(group = factor(group))
ab <- df %>% select(all_of(rownames(meta)))

linda.obj <- linda(otu.tab = ab, meta = meta, 
                   formula = '~group', 
                   alpha = 0.05,type = "proportion",
                   prev.cut = 0.0001, 
                   #lib.cut = 1,
                   winsor.quan =NULL)


picked = rownames(linda.obj$output[[1]])[which(linda.obj$output[[1]]$reject)]
anno[picked,] %>% View()

inner_join(
  y = linda.obj$output[[1]] %>% 
    #.[which(linda.obj$output[[1]]$reject),] %>% 
    as_tibble(rownames = 'index'),
  x = anno %>% as_tibble() %>% mutate(index = as.character(index)),
  by = 'index'
) %>%
  arrange(pvalue) %>% #View()
  write_tsv('meta2/results.tsv')

linda.plot(linda.obj,variables.plot =  'group',
           titles ='NASH_F vs NASH_F_DSF', 
           alpha = 0.05, lfc.cut = 1,
           legend = T, 
           directory = 'meta2', 
           width = 15, height = 12)

# meta3  ---------------------------------------------
# 读入分组文件
gp <- read_xlsx('meta3/1_sample_metadata.xlsx') %>% 
  rename(seqid = index)

# 读入丰度文件
df <-  read_tsv('meta3/3_all_mpa_to_abundance_metagenome.tsv') %>% 
  rename(name = index) %>% 
  select(name,all_of(gp[['seqid']])) 
# 转换为数据框
df %<>% data.frame(row.names = 1)
# 菌名替换为简短名
anno = data.frame(index = 1:length(rownames(df)),anno = rownames(df))
# 去除行名
rownames(df) = NULL
# 分组文件改为数据框
gp %<>% data.frame(row.names = 'seqid')

# 开始LinDA
meta <- gp %>% dplyr::select('group') %>% filter(!is.na(group)) %>% mutate(group = factor(group))
ab <- df %>% select(all_of(rownames(meta)))

linda.obj <- linda(otu.tab = ab, meta = meta, 
                   formula = '~group', 
                   alpha = 0.05,type = "proportion",
                   prev.cut = 0.0001, 
                   #lib.cut = 1,
                   winsor.quan =NULL)


picked = rownames(linda.obj$output[[1]])[which(linda.obj$output[[1]]$reject)]
anno[picked,] %>% View()

inner_join(
  y = linda.obj$output[[1]] %>% 
    #.[which(linda.obj$output[[1]]$reject),] %>% 
    as_tibble(rownames = 'index'),
  x = anno %>% as_tibble() %>% mutate(index = as.character(index)),
  by = 'index'
) %>%
  arrange(pvalue) %>% #View()
  write_tsv('meta3/results.tsv')

linda.plot(linda.obj,variables.plot =  'group',
           titles ='NASH_F vs NASH_F_DSF', 
           alpha = 0.05, 
           lfc.cut = 1,
           legend = T, 
           directory = 'meta3', 
           width = 15, height = 12)

# 16s  ---------------------------------------------
# 读入分组文件
gp <- read_xlsx('16s1/1-16s-sample-metadata-20221010.xlsx') %>% 
  rename(seqid = index)

# 读入丰度文件
df <-  read_csv('16s1/1-16s-level-6-20221010.csv') %>% 
  data.table::transpose(make.names = 'index',keep.names = 'name') %>% 
  as_tibble() %>% 
  select(name,all_of(gp[['seqid']])) 
# 转换为数据框
df %<>% data.frame(row.names = 1)
# 菌名替换为简短名
anno = data.frame(index = 1:length(rownames(df)),anno = rownames(df))
# 去除行名
rownames(df) = NULL
# 分组文件改为数据框
gp %<>% data.frame(row.names = 'seqid')

# 开始LinDA
meta <- gp %>% 
  dplyr::select('group') %>% 
  filter(!is.na(group)) 


ab <- df[,rownames(meta)]

# ab2 <- apply(X = ab,
#             MARGIN =2,
#             FUN =  function(x){x/sum(x)}) %>% 
#   data.frame() #%>% colSums()

linda.obj <- linda(otu.tab = ab, 
                   meta = meta, 
                   formula = '~group',
                   alpha = 0.05,
                   type = "count",
                   prev.cut = 0.0001 
                   #lib.cut = 1,
                   #winsor.quan =NULL
                   )


picked = rownames(linda.obj$output[[1]])[which(linda.obj$output[[1]]$reject)]
anno[picked,] %>% View()

inner_join(
  y = linda.obj$output[[1]] %>% 
    #.[which(linda.obj$output[[1]]$reject),] %>% 
    as_tibble(rownames = 'index'),
  x = anno %>% as_tibble() %>% mutate(index = as.character(index)),
  by = 'index'
) %>%
  arrange(pvalue) %>% #View()
  write_tsv('16s1/results.tsv')

linda.plot(linda.obj,variables.plot =  'group',
           titles ='group4', 
           alpha = 0.05, 
           lfc.cut = 1,
           legend = T, 
           directory = '16s1', 
           width = 15, height = 12)

