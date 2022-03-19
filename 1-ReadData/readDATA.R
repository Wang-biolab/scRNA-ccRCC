library(Seurat)
rm(list=ls())
dir <- dir("RawData/")
dir <- paste0("RawData/", dir)
#查看文件顺序
dir                         
#按文件顺序给样本命名，名称不要以数字开头，中间不能有空格 
group = c('p1','p2','p3','t1','t2','t3','rb1','rb2','rb3','rn1','rn2','rn3','rt1','rt2','rt3')

##使用循环命令批量创建seurat对象
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNA <- CreateSeuratObject(counts, project= group[i])
  save(scRNA, file= paste('RawData/scRNA_',group[i],'.Rdata',sep=''))}


#合并 seurat 对象
scRNA <- merge(x = seob_list[[1]], #第一个
              y = seob_list[-1], #其他的
              add.cell.ids = names(seob_list) #cell id 添加前缀
              )

  
