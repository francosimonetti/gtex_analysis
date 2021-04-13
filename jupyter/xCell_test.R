install.packages("devtools")
devtools::install_github('dviraran/xCell')

library(xCell)
# expr_file = "/cbscratch/franco/trans-eqtl/new_preprocess_feb2020_freeze/gtex_v8/expression/tpms/wb_tpms_qcfilter.txt.protein_coding_lncRNA_filtered.gene_names"
expr_file = "/cbscratch/franco/datasets/gtex_v8/expression/TPMs_phASER_GTEx_v8.matrix.txt.gene_names"
exprMatrix = read.table(expr_file,header=TRUE,row.names=1, as.is=TRUE)
# white_blood_cells = c('B-cells','CD4+ T-cells','CD4+ naive T-cells','CD4+ memory T-cells',
#                      'CD4+ Tcm','CD4+ Tem','CD8+ T-cells','CD8+ naive T-cells','CD8+ Tcm',
#                      'CD8+ Tem','Class-switched memory B-cells','Memory B-cells','Tgd cells',
#                      'Monocytes', 'Macrophages','Macrophages M1','Macrophages M2',
#                      'naive B-cells','NK cells','NKT','Neutrophils','Tregs')

res <- xCellAnalysis(exprMatrix, parallel.sz=16)
write.table(x=res, file="TPMs_GTEx_v8_phASER_all_xCell.txt")