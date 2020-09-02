# RECITublet R script v1.0 (last updated 02/09/2020)
# Author: Rachael Bashford-Rogers (Wellcome Centre for Human Genetics, University of Oxford, rbr1@well.ox.ac.uk)
# This script will identify scRNA-seq droplets that have features of doublets/multiplets based on VDJ-seq information


#### identifying VDJ doublets/multiplets
#### load requirements
library("Seurat")


concat = function(v) {
	res = ""
	for (i in 1:length(v)){res = paste(res,v[i],sep="")}
	res
}


############# load scRNA-seq data
pbmc  = readRDS(file="scRNAseq.file.rds") ### rds object of seurat object
type = "ALL"

cell_ids = rownames(pbmc@meta.data)
scale_data = pbmc@scale.data
raw_data = pbmc@raw.data
genes = rownames(scale_data)
orig.ident = pbmc@meta.data$orig.ident
cell_type = pbmc@meta.data$cell_type_overall
names(cell_type) = cell_ids

raw_data = raw_data[,colnames(scale_data)]
cluster = cell_type
clusters = sort(unique(cluster))
orig.idents = sort(unique(orig.ident))

pbmc.umap  = readRDS(file="Seurat_UMAP.umap") ### rds object of seurat UMAP (or other dimensionaltiy reduction) for the scRNA-seq data
xy = pbmc.umap$layout

######################## identify doublets from BCR/TCR information
Get_BCR_TCR_doublets<-function(m_VDJ, threshold_umis, cell_type){
	cell_types = sort(unique(cell_type))
	B_cell_types = cell_types[grep("B cell", cell_types)]
	T_cell_types = cell_types[grep("T cell", cell_types)]
	non_B_T_cell_types = setdiff(cell_types, c(B_cell_types, T_cell_types))
	wb1 = suppressWarnings (cell_ids [intersect(which(m_VDJ$BCR[,"constant_region1"]!='-'), which(as.numeric(m_VDJ$BCR[,"n_umis1"])>= threshold_umis))])
	wb2 = suppressWarnings (cell_ids [intersect(which(m_VDJ$BCR[,"constant_region2"]!='-'), which(as.numeric(m_VDJ$BCR[,"n_umis2"])>= threshold_umis))])
	wt1 = suppressWarnings (cell_ids [intersect(which(m_VDJ$TCR[,"constant_region1"]!='-'), which(as.numeric(m_VDJ$TCR[,"n_umis1"])>= threshold_umis))])
	wt2 = suppressWarnings (cell_ids [intersect(which(m_VDJ$TCR[,"constant_region2"]!='-'), which(as.numeric(m_VDJ$TCR[,"n_umis2"])>= threshold_umis))])
	wb12 = intersect(wb1,wb2)
	wt12 = intersect(wt1,wt2)
	wb= unique(c(wb1,wb2))
	wt= unique(c(wt1,wt2))
	high_confidence_clashes = intersect(intersect(wb1,wb2), intersect(wt1,wt2))### both chains present for TCR/BCR
	mid_confidence_clashes = unique(c(wb12 [which(wb12 %in% c(wt1,wt2))],wt12 [which(wt12 %in% c(wb1,wb2))]))
	lower_confidence_clashes = intersect(wb, wt)
	
	non_B_cell_BCRs = intersect(rownames(m_VDJ$BCR)[which(cell_type %in% B_cell_types==F)], c(wb1, wb2))
	non_B_cell_BCRs = non_B_cell_BCRs[which(as.numeric (m_VDJ$BCR[non_B_cell_BCRs,"n_umis1"])> threshold_umis)]
	non_T_cell_TCRs = intersect(rownames(m_VDJ$TCR)[which(cell_type %in% T_cell_types==F)], c(wt1, wt2))
	non_T_cell_TCRs = non_T_cell_TCRs[which(as.numeric(m_VDJ$TCR[non_T_cell_TCRs,"n_umis1"])> threshold_umis)]

	VDJ_BCR_TCR_doublets = c(list(high_confidence_clashes), list(mid_confidence_clashes), list(lower_confidence_clashes),list(non_B_cell_BCRs), list(non_T_cell_TCRs))
	names(VDJ_BCR_TCR_doublets) = c("IGH+IGK/L+TRA+TRB","3 of [IGH,IGK/L,TRA,TRB]", "IGH or IGK/L + TRA or TRB","non-B cell clustered + BCR", "non-T cell clustered + TCR" )
	return(VDJ_BCR_TCR_doublets)
}
Count_BCR_TCR_doublets_by_cell_type<-function(cell_type, cell_types, VDJ_BCR_TCR_doublets){
	headers = c("all","high","mid","low")
	counts = matrix(data = 0,nrow = length(cell_types),ncol = length(headers), dimnames = c(list(cell_types),list(headers)))
	t = table(cell_type)
	counts[names(t), "all"] = t
	t = table(cell_type[ VDJ_BCR_TCR_doublets $ higher_confidence])
	counts[names(t), "high"] = t
	t = table(cell_type[ VDJ_BCR_TCR_doublets $ mid_confidence])
	counts[names(t), "mid"] = t
	t = table(cell_type[ VDJ_BCR_TCR_doublets $ lower_confidence])
	counts[names(t), "low"] = t
	return(counts)}
Get_intra_BCR_TCR_doublets <-function(m_VDJ, threshold_umis){
	w = which(m_VDJ$BCR[,"mixed_contig_n_umis1"]!='-')
	Bchain1_doublets = rownames(m_VDJ$BCR[w,]) [which(as.numeric(m_VDJ$BCR[w,"mixed_contig_n_umis1"])> threshold_umis)]
	w = which(m_VDJ$BCR[,"mixed_contig_n_umis2"]!='-')
	Bchain2_doublets = rownames(m_VDJ$BCR[w,]) [which(as.numeric(m_VDJ$BCR[w,"mixed_contig_n_umis2"])> threshold_umis)]
	
	w = which(m_VDJ$TCR[,"mixed_contig_n_umis1"]!='-')
	Tchain1_doublets = rownames(m_VDJ$TCR[w,]) [which(as.numeric(m_VDJ$TCR[w,"mixed_contig_n_umis1"])> threshold_umis)]
	w = which(m_VDJ$TCR[,"mixed_contig_n_umis2"]!='-')
	Tchain2_doublets = rownames(m_VDJ$TCR[w,]) [which(as.numeric(m_VDJ$TCR[w,"mixed_contig_n_umis2"])> threshold_umis)]

	Tchain12_doublets = intersect(Tchain1_doublets, Tchain2_doublets)
	Bchain12_doublets = intersect(Bchain1_doublets, Bchain2_doublets)
	
	VDJ_intra_BCR_TCR_doublets = c(list(Bchain1_doublets),list(Bchain2_doublets),list(Tchain1_doublets), list(Tchain2_doublets),list(Bchain12_doublets), list(Tchain12_doublets))
	names(VDJ_intra_BCR_TCR_doublets) = c("2x IGHs","2x IGK/Ls","2x TRAs","2x TRBs","2x IGH and 2x IGK/L","2x TRAs and 2x TRBs")
	return(VDJ_intra_BCR_TCR_doublets)
}
Count_intra_BCR_TCR_doublets_by_cell_type <-function(cell_type, cell_types, VDJ_intra_BCR_TCR_doublets){
	headers = c("all",names(VDJ_intra_BCR_TCR_doublets))
	counts = matrix(data = 0,nrow = length(cell_types),ncol = length(headers), dimnames = c(list(cell_types),list(headers)))
	t = table(cell_type)
	counts[names(t), "all"] = t
	for(i in c(1:length(names(VDJ_intra_BCR_TCR_doublets)))){
		t = table(cell_type[ VDJ_intra_BCR_TCR_doublets [[i]]])
		counts [names(t),names(VDJ_intra_BCR_TCR_doublets)[i]] = t
	}
	return(counts)}
Boxplot_nUMI_ngene_mt_doublets<-function(doublets_list, pbmc, fileout1,metaD){
	headers = colnames(metaD)
	all_cells = rownames(metaD)
	library(RColorBrewer)
	cols = add.alpha (brewer.pal(8, "Dark2"), alpha = 0.5)
	w=2.5
	pdf(file=fileout1, height=w*1.1, width=w*4*0.9)
	par(mfrow= c(1,5), mar = c(12,4.5,3,1))
	for(h in c(1:length(headers))){
		groups = NULL
		names = NULL
		factors = c(names(doublets_list), "other")
		for(i in c(1:length(doublets_list))){
			x =metaD[doublets_list[[i]], headers[h]]
#			x = x[which(x<25)]
			if(length(x)>2){
				names = c(names, factors[i])
				groups = c(groups, list(x))}}
		x1 = unlist(groups)
		groups= c(groups, list(metaD[setdiff(all_cells, unlist(doublets_list)), headers[h]]))
		names = c(names, "other")
		x2 = metaD[setdiff(all_cells, unlist(doublets_list)), headers[h]]
		pval = wilcox.test(x1,y=x2,alternative = "two.sided")$p.value
		pval1 = signif(pval,digits =3)
		if(pval<1e-10){pval1 = "<1e-10"}
		# names = c(names, factors[i+1])
		boxplot(groups, ylab = "", col = cols,names= names, las = 2, main= concat(c(headers[h],"\np-value ",pval1)), outline =T, border = "grey",cex.lab = 0.9, lwd = 0.5, cex.main = 0.7, ylim = c(0,max(unlist(groups))))
	}
	dev.off()
}

#### get the VDJ information stored as table with the same rows as that in the Seurat object, and headers as follows for the BCR and TCR separately: 
# headers:  "X.cell","contig1","chain1","constant_region1","n_umis1","V_gene_10X1",
	"J_gene_10X1","cdr3_aa1","cdr3_nn1","V_gene1","J_gene1","V_mm1" , "J_mm1","mixed_contig_chain1","mixed_contig_n_umis1","contig2", "chain2", "constant_region2",
	"n_umis2","V_gene2","J_gene2", "cdr3_aa2","cdr3_nn2","V_gene2.1", "J_gene2.1","V_mm2","J_mm2","mixed_contig_chain2","mixed_contig_n_umis2"
# where 1 refers to IGH or TRA for B and T cells respectively, and 2 refers to IGK/L or TRB for B and T cells respectively

m_VDJ_BCR  = readRDS(file= m_VDJ_BCR.file.rds) ### rds object of BCR object
m_VDJ_TCR  = readRDS(file= m_VDJ_TCR.file.rds) ### rds object of TCR object

######################## normalise the variables that will be used in the predictor
matrix = pbmc@ raw.data
mitrib <- grep(pattern = "^MT-|^RP[SL]", x = rownames(x = matrix), value = TRUE, ignore.case = TRUE)
mito <- grep(pattern = "^MT-", x = rownames(x = matrix), value = TRUE, ignore.case = TRUE)
mitribcounts<- matrix[which(rownames(matrix) %in% mitrib), ]
mitoribo_ratio <- Matrix::colSums(mitribcounts[mito, , drop = FALSE])/Matrix::colSums(mitribcounts)
pbmc <- AddMetaData(object = pbmc, metadata = as.data.frame(mitoribo_ratio) ,col.name = "mito.ribo_ratio")

u1 = m_VDJ$BCR[,"n_umis1"]
u2 = m_VDJ$BCR[,"n_umis2"]
u3 = m_VDJ$TCR[,"n_umis1"]
u4 = m_VDJ$TCR[,"n_umis2"]
u1[which(u1=='-')]=0
u2[which(u2=='-')]=0
u3[which(u3=='-')]=0
u4[which(u4=='-')]=0
n_UMIs_VDJ = as.numeric(u1)+as.numeric(u2)+as.numeric(u3)+as.numeric(u4)
names(n_UMIs_VDJ) = cell_ids
pbmc <- AddMetaData(object = pbmc, metadata = as.data.frame(n_UMIs_VDJ) ,col.name = "n_UMIs_VDJ")

variables = pbmc@meta.data[,c("percent.mito", "nGene","nUMI","mitoribo_ratio","n_UMIs_VDJ")]
rownames(variables)= cell_ids
variables_norm = variables
library(compositions)
for(i in c(1:length(variables[1,]))){
	x = variables[,i]
	names(x) = cell_ids
	for(s in c(1:length(orig.idents))){
		w = cell_ids[which(orig.ident== orig.idents[s])]
		if(length(w)>0){
			# clr1 = clr(x[w])
			# clr1  = scale(as.numeric(clr1), center = TRUE, scale = TRUE)
			clr1 = x[w]/mean(x[w])
			variables_norm[w,i] = as.numeric(clr1)
			print(concat(c(i, " ", orig.idents[s] )))
		}}
}
######################## Indetify and count doublets/multiplets from the VDJ-seq data
VDJ_BCR_TCR_doublets = Get_BCR_TCR_doublets(m_VDJ, threshold_umis = 10, cell_type)
VDJ_intra_BCR_TCR_doublets = Get_intra_BCR_TCR_doublets(m_VDJ, threshold_umis = 10)
counts_BCR_TCR_doublets = Count_BCR_TCR_doublets_by_cell_type(cell_type, cell_types, VDJ_BCR_TCR_doublets)
counts_intra_BCR_TCR_doublets = Count_intra_BCR_TCR_doublets_by_cell_type(cell_type, cell_types, VDJ_intra_BCR_TCR_doublets)

doublets = c(VDJ_BCR_TCR_doublets, VDJ_intra_BCR_TCR_doublets)
metaD = variables_norm

######################## plot nGenes, nUMIs and mt-rb-ratio between these methods
fileout1 ="Seurat_CITE_nGene_VDJ_doublets.pdf"
Boxplot_nUMI_ngene_mt_doublets(doublets, pbmc, fileout1,metaD)

fileout1 ="Seurat_CITE_nGene_VDJ_doublets.pdf"
VDJ_doublets1 = c(list(sort(unique(unlist(doublets[c("IGH+IGK/L+TRA+TRB","3 of [IGH,IGK/L,TRA,TRB]", "IGH or IGK/L + TRA or TRB")])))),
list(sort(unique(unlist(doublets[c("2x IGHs")])))),
list(sort(unique(unlist(doublets[c("2x TRBs")])))),
list(doublets[["non-B cell clustered + BCR"]]),list(doublets[["non-T cell clustered + TCR"]]))
names(VDJ_doublets1) = c("BCR + TCR","multiple BCRs","multiple TCRs","non-B cell clustered + BCR","non-T cell clustered + TCR")
Boxplot_nUMI_ngene_mt_doublets(VDJ_doublets1, pbmc, fileout1,metaD)

######################## plot UMAP locations (BCR, TCR etc)

library(RColorBrewer)
fileout1="Seurat_CITE_nGene_VDJ_doublets_UMAP_overall.pdf"
w=2.5
pdf(file=fileout1, height=w*2, width=w*5)
par(mfrow= c(2,5), mar = c(4.5,4.5,2,2))

for(ind in c(1:length(doublets1))){
	cluster = rep("other", length(cell_ids))
	names(cluster) = cell_ids
	cluster[doublets1[[ind]]] = names(doublets1)[ind]
	clusters = sort(unique(cluster))
	if(length(clusters)>=2){
		cluster_num = match(cluster, clusters)
		cols=add.alpha(colorRampPalette(c('red',"grey"))(length(clusters)), alpha = 0.7)
		w = order(cluster_num)
		if(clusters[2]=="other"){w = order(cluster_num, decreasing = T)}
		cex1 = c(0.05,0.05)
		plot(xy[w,1], xy[w,2],pch = 21, col = cols[cluster_num[w]],bg = cols[cluster_num[w]] ,main = names(doublets1)[ind], cex = cex1[cluster_num[w]],xlab = "umap dim1",ylab = "umap dim2")
		print (ind)
}}
cluster = rep("other", length(cell_ids))
names(cluster) = cell_ids
cluster[sort(unique(unlist(doublets)))] = 'VDJ-doublet'
clusters = sort(unique(cluster))
cluster_num = match(cluster, clusters)
cols=add.alpha(colorRampPalette(c('red',"grey"))(length(clusters)), alpha = 0.7)
w = order(cluster_num)
if(clusters[2]=="other"){w = order(cluster_num, decreasing = T)}
cex1 = c(0.05,0.05)
plot(xy[w,1], xy[w,2],pch = 21, col = cols[cluster_num[w]],bg = cols[cluster_num[w]] ,main = "VDJ doublets", cex = cex1[cluster_num[w]],xlab = "umap dim1",ylab = "umap dim2")
dev.off()

######################## Write output
VDJ_doublets_summary = c(list(VDJ_BCR_TCR_doublets),list(VDJ_intra_BCR_TCR_doublets))
names(VDJ_doublets_summary) = c('VDJ_BCR_TCR_doublets', 'VDJ_intra_BCR_TCR_doublets')
saveRDS(file="Seurat_VDJ_doublets.rds", VDJ_doublets_summary)




