# RECITublet R script v1.0 (last updated 02/09/2020)
# Author: Rachael Bashford-Rogers (Wellcome Centre for Human Genetics, University of Oxford, rbr1@well.ox.ac.uk)
# This script will predict scRNA-seq doublets/multiplets based on the features of known doublets/multiplet features (from VDJ-seq and/or CITE-seq information)


#### identifying VDJ doublets/multiplets
#### load requirements
library("Seurat")

concat = function(v) {
	res = ""
	for (i in 1:length(v)){res = paste(res,v[i],sep="")}
	res
}

############# load scRNA-seq data
pbmc  = readRDS(file="scRNAseq.file.rds") ### rds object of seurat RNA-seq object

############# get features for learning doublet/multiplet features
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
		}}}

metaD = variables_norm

######################## load known doublets/multiplets (here we are combining both VDJ-seq and CITE-seq identified doublets/multiplets)
all_CITE_doublets  = readRDS(file="Seurat_CITE_doublets.rds") ### load CITE-seq doublets/multiplets (pre-computed)
VDJ_doublets_summary =readRDS(file="Seurat_VDJ_doublets.rds")  ### load VDJ-seq doublets/multiplets (pre-computed)

######################## combine known doublets and set up classifier
VDJ_doublets_total = sort(unique(c(unlist(VDJ_doublets_summary[[2]]))))
doublets_CITE_plus_VDJ = unique(sort(unlist(c(VDJ_doublets_total, all_CITE_doublets))))

variables = metaD
cell_ids = rownames(variables)
doublet_probabilities = rep(-1, length(cell_ids))
names(doublet_probabilities) = cell_ids

double_pos = sort(unique(c(unlist(doublets_CITE_plus_VDJ))))
single_pos = setdiff(cell_ids, double_pos)

### check if there is a difference between known doublets and the remainder of droplets
summary(variables_norm[double_pos,"nCount_RNA"])
summary(variables_norm[single_pos,"nCount_RNA"])

######################## exclude outliers from doublet detection (optional)
exclude_doublets = NULL
exclude_singlets = NULL
for(h in c(1:length(variables_norm[1,]))){ ### removes droplets that are outside of the median +/- sd of the distribution of features for the known doublets/multiplets and the remainder
	var = variables_norm[double_pos,h]
	exclude_doublets  = c(exclude_doublets , double_pos [c(which(var<median(var)-1*sd(var)), which(var>median(var)+1*sd(var)))])
	var = variables_norm[single_pos,h]
	exclude_singlets  = c(exclude_singlets , single_pos [c(which(var<median(var)-1*sd(var)), which(var>median(var)+1*sd(var)))])
}
exclude_doublets = unique(exclude_doublets)
exclude_singlets = unique(exclude_singlets)
double_pos1 = setdiff(double_pos, exclude_doublets)
single_pos1 = setdiff(single_pos, exclude_singlets)

### check if there is a difference between known doublets and the remainder of droplets
summary(variables_norm[double_pos1,"nCount_RNA"])
summary(variables_norm[single_pos1,"nCount_RNA"])

######################## set up classifier
cells_sub = c(double_pos1, single_pos1)
names(seurat_clusters) = cell_ids
library(MASS)
variables_norm1 = variables_norm
x = data.frame(variables_norm1[cells_sub,]) ### subset of data (if excluded outliers in step above)
x_all = data.frame(variables_norm1) ### all data
x1=x
class = c(rep("doublet", length(double_pos1)), rep("singlet", length(single_pos1)))
x1$class = factor(class)
colnames(x_all) = colnames(x1)[which(colnames(x1)!="class")]
x = x1[,which(colnames(x1)!="class")]
model <- glm(class ~ ., x1,family=binomial(link='logit'))

######################## check prediction accuracy of model based on the training data
p <- predict(model, newdata=x, type="response")
threshold = 0.5
glm.pred=rep("doublet", length(p))
glm.pred[p>= threshold]="singlet"
prediction_table0 = table(factor(factor(glm.pred)), factor(class))

######################## run model over all data to get doublet/multiplet predictions 
p <- predict(model, newdata=x_all, type="response")
glm.pred=rep("doublet", length(p))
glm.pred[p>=threshold]="singlet"
names(glm.pred) = cell_ids
table(glm.pred)*100/length(glm.pred)

print(concat(c(length(intersect(names(w), double_pos))*100/length(double_pos),"% doublets identified from original list")))
print(concat(c(length(intersect(names(w), VDJ_doublets_total))*100/length(VDJ_doublets_total),"% VDJ doublets identified from original list")))

### check if there is a difference between predicted doublets and singlets
summary(variables[which(glm.pred =="doublet"),"nCount_RNA"])
summary(variables[which(glm.pred =="singlet"),"nCount_RNA"])	

######################## summarise output 
prediction_table1 = table(glm.pred)*100/length(glm.pred)
doublet_probabilities_without_sample = 1-p
doublet_probabilities_without_sample_add_true_doublets = doublet_probabilities_without_sample
doublet_probabilities_without_sample_add_true_doublets[VDJ_doublets_total] = 1
pred_doublets = names(which(doublet_probabilities_without_sample_add_true_doublets>0.5))

doublet_predictions = c(list(prediction_table0), list(prediction_table1), list(doublet_probabilities_without_sample_add_true_doublets),list(pred_doublets), list(doublet_probabilities_without_sample))

names(doublet_predictions) = c("prediction_table_training", "prediction_table_total_data", "doublet_probabilities_including_true_doublets","predicted_doublets", "doublet_probabilities_raw")

doublets_pred = c(list(unlist(VDJ_doublets_total)), list(unlist(all_CITE_doublets)), list(pred_doublets))
names(doublets_pred) = c("VDJ doublets","CITE-seq doublets", "predicted doublets")

######################## plot nGenes, nUMIs and mt% between these methods
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
		boxplot(groups, ylab = "", col = cols,names= names, las = 2, main= concat(c(headers[h],"\np-value ",pval1)), outline =T, border = "grey",cex.lab = 0.9, lwd = 0.5, cex.main = 0.7, ylim = c(0,max(unlist(groups))))
	}
	dev.off()
}

fileout1 ="Seurat_CITE_nGene_predicted_doublets.pdf"))
Boxplot_nUMI_ngene_mt_doublets(doublets_pred, pbmc, fileout1,metaD)

######################## save doublet predictions
saveRDS(file="Probabilities_RECITublet.doublet_probabilities_CITE_VDJ_seq", doublet_predictions)
saveRDS(file="Probabilities_RECITublet.predicted_doublets_RECITlet", doublet_predictions) ## RDS object of doublets/multiplets from VDJ-seq,CITE-seq, and predicted doublets/multiplets using the VDJ-seq/CITE-seq inputs as training points




