# RECITublet R script v1.0 (last updated 02/09/2020)
# Author: Rachael Bashford-Rogers (Wellcome Centre for Human Genetics, University of Oxford, rbr1@well.ox.ac.uk)
# This script will identify scRNA-seq droplets that have features of doublets/multiplets based on CITE-seq information


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

pbmc.umap  = readRDS(file="Seurat_UMAP.umap") ### rds object of seurat UMAP (or other dimensionaltiy reduction) for the scRNA-seq data
xy = pbmc.umap$layout

############# load and scale CITE-seq data
Scale_CITE_data<-function(CITE_seq, CITE_data, orig.ident, orig.idents){
	library(compositions)
	CITE_data_scaled = CITE_data*NA
	for(i in c(1:length(CITE_seq))){
		x = CITE_data[i,]
		w1 = names(which(x!=0))
		for(s in c(1:length(orig.idents))){
			w = intersect(rownames(pbmc@meta.data)[which(orig.ident== orig.idents[s])], w1)
			if(length(w)>0){
				clr1 = clr(x[w])
				clr1  = scale(as.numeric(clr1), center = TRUE, scale = TRUE)
				CITE_data_scaled[i,w] = as.numeric(clr1)
				print(concat(c(CITE_seq[i], " ", orig.idents[s] )))
			}}
	}
	CITE_data_scaled = CITE_data_scaled[,rownames(pbmc@meta.data)]
	return(CITE_data_scaled)
}

CITE_data =readRDS(file="scCITEseq.file.rds") ### matrix of CITE-seq counts
CITE_data_scaled = Scale_CITE_data(CITE_seq, CITE_data, orig.ident, orig.idents)


############# CITE-seq - gene matching files
genes = rownames(pbmc@assays$ RNA@counts)
ref_file = "CITE_Seq_gene_name.txt" ### file containing the mapping of CITE-seq names with gene names
p <- as.matrix(read.csv(ref_file, head=T, sep="\t"))
gene_CITE =p[match(CITE_seq, p[,1]),2]
names(CITE_seq) = gene_CITE

############# functions: threshold CITE-seq data to determine. which cells are expressing each marker
Plot_CITE_seq_versus_gene_expression<-function(gene_CITE, gene_CITEs, CITE_data_all,orig.ident, orig.idents, CITE_seq,  scale_data,file_pre){
	for(i in c(1:length(CITE_seq))){
		fileout1=concat(c(file_pre, CITE_seq[i],".pdf"))
		library(RColorBrewer)
		cols =  add.alpha (brewer.pal(8, "Dark2"), alpha = 0.5)
		w=3.25
		pdf(file=fileout1, height=w*5, width=w*6)
		par(mfrow= c(5,6), mar = c(5,5,3,3))
		for(s in c(1:length(orig.idents))){
			w = which(orig.ident== orig.idents[s])
			x = CITE_data_all [CITE_seq[i], w]
			y = scale_data[gene_CITE[i], w]
			w1 = intersect(which(is.na(x)==F),which(is.na(y)==F))
			x = x[w1]
			y = y[w1]				
			if(length(x)>100){
				print (length(x))
				plot (y, x, col = cols[1], pch = 21, bg =  cols[1],main = concat(c(orig.idents[s],"\n", CITE_seq[i])), xlab = gene_CITE[i], ylab =  CITE_seq[i])
			}}
		print(concat(c("scp -p mfj169@rescomp1.well.ox.ac.uk:", fileout1," ./ " )))
		print (i)
		dev.off()
	}
}
Get_mean_gene_expression<-function(gene_CITE, clusters, cluster, gene_CITEs, scale_data){
	gene_CITEs = sort(unique(gene_CITE))
	m_mean_expression = t(matrix(data = 0,nrow = length(clusters),ncol = length(gene_CITEs), dimnames = c(list(clusters),list(gene_CITEs))))
	m_quantile_expression = t(matrix(data = 0,nrow = length(clusters),ncol = length(gene_CITEs), dimnames = c(list(clusters),list(gene_CITEs))))
	for(c in c(1:length(clusters))){
		w = which(cluster== clusters[c])
		m_mean_expression[,clusters[c]] = apply(scale_data[gene_CITEs,w], 1, mean)
		m_quantile_expression[,clusters[c]] = apply(scale_data[gene_CITEs,w], 1, function(x){quantile(x,prob = 0.9)})
		print(c)	}
	CITE_gene_expression = c(list(m_mean_expression), list(m_quantile_expression))
	names(CITE_gene_expression) = c("mean", "90th quantile")
	return(CITE_gene_expression)
}
CITE_seq_level_per_sample<-function(orig.ident, orig.idents, CITE_seq, clusters, cluster, CITE_data_scaled){
	m_CITE_expression_by_sample = t(matrix(data = 0,nrow = length(orig.idents),ncol = length(CITE_seq), dimnames = c(list(orig.idents),list(CITE_seq))))
	m_CITE_expression_by_cluster = t(matrix(data = 0,nrow = length(clusters),ncol = length(CITE_seq), dimnames = c(list(clusters),list(CITE_seq))))
	non_na = NULL
	for(i in c(1:length(CITE_seq))){
		non_na = c(non_na, list(which(is.na(CITE_data_scaled[CITE_seq[i],])==F)))}
	names(non_na) = CITE_seq
	for(c in c(1:length(orig.idents))){
		w = which(orig.ident== orig.idents[c])
		print (c)
		for(i in c(1:length(CITE_seq))){
			w1 = intersect(non_na[[i]],w )
			m_CITE_expression_by_sample[CITE_seq[i], orig.idents[c]] = mean(CITE_data_scaled[CITE_seq[i],w1])
		}}
	for(c in c(1:length(clusters))){
		print(c)
		w = which(cluster== clusters[c])
		for(i in c(1:length(CITE_seq))){
			w1 = intersect(non_na[[i]],w )
			m_CITE_expression_by_cluster[CITE_seq[i],clusters[c]] = mean(CITE_data_scaled[CITE_seq[i],w1])
		}}
	CITE_seq_expression = c(list(m_CITE_expression_by_sample), list(m_CITE_expression_by_cluster))
	names(CITE_seq_expression) = c("m_CITE_expression_by_sample", "m_CITE_expression_by_cluster")
	return(CITE_seq_expression)
}
Threshold_CITE_seq_positivity<-function(orig.idents, orig.ident, CITE_seq, sample_use, m_CITE_expression_by_sample, CITE_data, CITE_gene_expression, cell_ids){
	library(MASS)
	thresholds = t(matrix(data = 0,nrow = length(orig.idents),ncol = length(CITE_seq), dimnames = c(list(orig.idents),list(CITE_seq))))
	prediction_accuracy = t(matrix(data = 0,nrow = length(orig.idents),ncol = length(CITE_seq), dimnames = c(list(orig.idents),list(CITE_seq))))
	pos_accuracy = t(matrix(data = 0,nrow = length(orig.idents),ncol = length(CITE_seq), dimnames = c(list(orig.idents),list(CITE_seq))))
	neg_accuracy = t(matrix(data = 0,nrow = length(orig.idents),ncol = length(CITE_seq), dimnames = c(list(orig.idents),list(CITE_seq))))
	cell_ids_CITE_GEX = intersect(cell_ids, colnames(CITE_data))
	for(i in c(1:length(CITE_seq))){
		sample_use = names(which(CITE_seq_expression$m_CITE_expression_by_sample[CITE_seq[i],]!=0))
		quantile = CITE_gene_expression $ '90th quantile'[names(CITE_seq)[i],]
		quantile = quantile[which(is.na(quantile)==F)]
		mean = CITE_seq_expression$m_CITE_expression_by_cluster[CITE_seq[i],]
		mean = mean[which(is.na(mean)==F)]
		pos = intersect(names(which(quantile> quantile(quantile, 0.5))), names(which(mean > quantile(mean, 0.5))))
		neg = names(which(quantile <=  quantile(quantile, 0.1)))
		w_use_na = names(which(is.na(CITE_data[CITE_seq[i],])==F))
		# neg = setdiff(clusters, pos)
		for(s in c(1:length(sample_use))){
			w = cell_ids [which(orig.ident== sample_use[s])]
			w1 = intersect(w,  cell_ids [which(cluster %in% neg)])
			w2 = intersect(w,  cell_ids [which(cluster %in% pos)])
			w1 = intersect(intersect(w1, cell_ids_CITE_GEX), w_use_na)
			w2 = intersect(intersect(w2, cell_ids_CITE_GEX), w_use_na)
			if(min(c(length(w1), length(w2)))>5){
				x_pos = CITE_data_all[CITE_seq[i],w2]
				x_neg = CITE_data_all[CITE_seq[i],w1]
				x = c(x_neg, x_pos)
				if(length(x)>10){
					x1 = c(x_neg, x_pos)
					class = c(rep("neg", length(x_neg)), rep("pos", length(x_pos)))
					ids = c(w1,w2)
					class1 = data.frame(ids )
					class1$class = factor(class)
					class1$x = x
					model <- lda(class ~x, data = class1)
					predictions <- model %>% predict(class1)
					print(table(predictions $class))
					print(s)
					if(min(table(predictions $class))>5){
						x2 = data.frame(seq(from = max(x[which(predictions $class=="neg")]), to = min(x[which(predictions $class=="pos")]), length = 10000))
						colnames(x2) = 'x'
						pred = predict(model, newdata = x2)$class
						thr = mean(c(min(x2[which(pred=="pos"),]), max(x2[which(pred=="neg"),]))) 
						acc = length(which(predictions$class==class))*100/length(class)
						pacc = length(intersect(which(class=="pos"),which(x1>= thr)))*100/length(which(class=="pos"))
						nacc = length(intersect(which(class=="neg"),which(x1< thr)))*100/length(which(class=="neg"))
						thresholds[CITE_seq[i],sample_use[s]] =  thr
						prediction_accuracy [CITE_seq[i],sample_use[s]] = acc
						pos_accuracy [CITE_seq[i],sample_use[s]] =  pacc
						neg_accuracy [CITE_seq[i],sample_use[s]] =  nacc
						print(i)
			}}}
		}
	}
	CITE_seq_positivity_thresholds = c(list(thresholds), list(prediction_accuracy), list(pos_accuracy), list(neg_accuracy))
	names(CITE_seq_positivity_thresholds) = c(("thresholds"),("prediction_accuracy"),("pos_accuracy"),("neg_accuracy"))
	return(CITE_seq_positivity_thresholds)
}
Choose_CITE_seq_markers<-function(CITE_seq, CITE_seq_positivity_thresholds, Tpos_accuracy, Tneg_accuracy, Tprediction_accuracy){
	headers = c("thresholds", "prediction_accuracy", "pos_accuracy", "neg_accuracy")
	thrs = matrix(data = 0,nrow = length(CITE_seq),ncol = length(headers), dimnames = c(list(CITE_seq),list(headers)))
	for(i in c(1:length(CITE_seq))){
		w1 = intersect(which(CITE_seq_positivity_thresholds$thresholds[CITE_seq[i],]!=0), which(is.finite(CITE_seq_positivity_thresholds$thresholds[i,])))
		w2 = intersect(which(CITE_seq_positivity_thresholds$prediction_accuracy[i,]!=0), which(is.finite(CITE_seq_positivity_thresholds$prediction_accuracy[i,])))
		w3 = intersect(which(CITE_seq_positivity_thresholds$pos_accuracy[i,]!=0), which(is.finite(CITE_seq_positivity_thresholds$pos_accuracy[i,])))
		w4 = intersect(which(CITE_seq_positivity_thresholds$neg_accuracy[i,]!=0), which(is.finite(CITE_seq_positivity_thresholds$neg_accuracy[i,])))
		thrs[i,] = c(mean(CITE_seq_positivity_thresholds$thresholds[i,w1]), mean(CITE_seq_positivity_thresholds$prediction_accuracy[i,w2]), mean(CITE_seq_positivity_thresholds$pos_accuracy[i,w3]), mean(CITE_seq_positivity_thresholds$neg_accuracy[i,w4]))
	}
	CITE_seq_use = CITE_seq[intersect (intersect(which(thrs[,'pos_accuracy']> Tpos_accuracy), which(thrs[,'neg_accuracy']> Tneg_accuracy)), intersect(which(is.finite(thrs[,'thresholds'])==T), which(thrs[,'prediction_accuracy']> Tprediction_accuracy)))]
	info_CITE_seq_use = thrs [CITE_seq_use,]
	CITE_seq_choose = c(list(thrs), list(CITE_seq_use), list(info_CITE_seq_use))
	names(CITE_seq_choose) = c("thrs","CITE_seq_use","info_CITE_seq_use")
	return(CITE_seq_choose)
}
Plot_CITE_seq_thresholds<-function(CITE_seq_choose, sample_use, fileout1, CITE_seq_positivity_thresholds,CITE_data_all, scale_data){
	library(RColorBrewer)
	cols =  rep(add.alpha (brewer.pal(8, "Dark2"), alpha = 0.5), 100)
	w=4
	pdf(file=fileout1, height=w*4, width=w*4)
	par(mfrow= c(5,5), mar = c(5,5,3,3))
	sample_match = match(orig.ident, orig.idents)
	for(i in c(1:length(CITE_seq_choose$CITE_seq_use))){
		sample_use = names(which(CITE_seq_expression$m_CITE_expression_by_sample[CITE_seq_choose$CITE_seq_use[i],]!=0))
		if(length(sample_use)>0){
			for(s in c(1:length(sample_use))){
				w = which(orig.ident== sample_use[s])
				x = scale_data[names(CITE_seq_choose$CITE_seq_use)[i],w]
				y = CITE_data_all[(CITE_seq_choose$CITE_seq_use[i]),w]
				w1 = intersect(which(is.na(x)==F),which(is.na(y)==F))
				x = x[w1]
				y = y[w1]
				if(length(w1)>100){
					plot (x, y, col = cols[sample_match[w]], pch = 21, bg =  cols[sample_match[w]],main = concat(c(sample_use[s],"\n", CITE_seq_choose$CITE_seq_use[i])), 
					xlab = names(CITE_seq_choose$CITE_seq_use)[i], ylab = (CITE_seq_choose$CITE_seq_use[i]))
					segments(-100, CITE_seq_positivity_thresholds $thresholds[CITE_seq_choose$CITE_seq_use[i], sample_use[s]], 100000, CITE_seq_positivity_thresholds $thresholds[CITE_seq_choose$CITE_seq_use[i], sample_use[s]], lwd = 2, col = "red")
				}}
			}
		print (i)
		}
	dev.off()
}
Score_CITE_positivity<-function(CITE_seq_choose, sample_use, CITE_seq_positivity_thresholds, CITE_data_all, cell_ids){
	CITE_score = (CITE_data_all*NA)
	for(i in c(1:length(CITE_seq_choose$CITE_seq_use))){
		sample_use = names(which(CITE_seq_expression$m_CITE_expression_by_sample [CITE_seq_choose$CITE_seq_use[i],]!=0))
		if(length(sample_use)>0){
			for(s in c(1:length(sample_use))){
				w = cell_ids [which(orig.ident== sample_use[s])]
				threshold = CITE_seq_positivity_thresholds $thresholds[CITE_seq_choose$CITE_seq_use[i], sample_use[s]]
				w1 = intersect(w, names (which(CITE_data_all[CITE_seq_choose$CITE_seq_use[i],]> threshold)))
				if(length(w1)>0){
					x = CITE_data_all[CITE_seq_choose$CITE_seq_use[i],w1]
					#score = rank(x)/length(x)
					score = (x-min(x))
					score = (score/max(score))+1
					CITE_score[CITE_seq_choose$CITE_seq_use[i],w1] = score}}}
		print (i)
		}
	return(CITE_score)}
Plot_CITE_scores<-function(CITE_score, CITE_data_all, CITE_seq_choose, orig.ident, orig.idents, fileout1,CITE_seq_expression){
	library(RColorBrewer)
	cols =  add.alpha (brewer.pal(8, "Dark2"), alpha = 0.5)
	w=4
	pdf(file=fileout1, height=w*4, width=w*4)
	par(mfrow= c(5,5), mar = c(5,5,3,3))
	sample_match = match(orig.ident, orig.idents)
	for(i in c(1:length(CITE_seq_choose$CITE_seq_use))){
		sample_use = names(which(CITE_seq_expression$m_CITE_expression_by_sample[CITE_seq_choose$CITE_seq_use[i],]!=0))
		if(length(sample_use)>0){
			for(s in c(1:length(sample_use))){
				w = which(orig.ident== sample_use[s])
				scores = CITE_score[CITE_seq_choose$CITE_seq_use[i],w]
				scaled_expression = CITE_data_all[CITE_seq_choose$CITE_seq_use[i],w]
				scores[which(is.na(scores))] = 0
				plot (scaled_expression, scores, col = cols[1], pch = 21, bg =  cols[1],main = concat(c(sample_use[s],"\n", CITE_seq_choose$CITE_seq_use[i])))
				}}
		print (i)
		}
	dev.off()
}

############# run thresholding functions
file_pre = "CITE_expression_versus_gene_expression_"
Plot_CITE_seq_versus_gene_expression(gene_CITE, gene_CITEs, CITE_data_all,orig.ident, orig.idents, CITE_seq,  scale_data,file_pre)
CITE_seq_expression  = CITE_seq_level_per_sample(orig.ident, orig.idents, CITE_seq,clusters = cell_types, cluster= cell_type, CITE_data_all)
CITE_gene_expression = Get_mean_gene_expression(gene_CITE, clusters, cluster, gene_CITEs, scale_data)
which(gene_CITE %in% genes)

############# predict thresholds of positive CITE seq cells
CITE_seq_positivity_thresholds = Threshold_CITE_seq_positivity(orig.idents, orig.ident, CITE_seq, sample_use, CITE_seq_expression, CITE_data_all, CITE_gene_expression,cell_ids)

############# choose CITE seq markers that have good predictions
CITE_seq_choose = Choose_CITE_seq_markers(CITE_seq, CITE_seq_positivity_thresholds,Tpos_accuracy = 20, Tneg_accuracy = 50, Tprediction_accuracy = 70)

############# plot CITE-seq thresholds
fileout1="CITE_RNA_correlation_thresholds.pdf"
Plot_CITE_seq_thresholds(CITE_seq_choose, sample_use, fileout1, CITE_seq_positivity_thresholds, CITE_data_all, scale_data)

############# generate CITE-seq score per positive cell and save output
CITE_score = Score_CITE_positivity(CITE_seq_choose, sample_use, CITE_seq_positivity_thresholds, CITE_data_all, cell_ids)
saveRDS(file="CITE_seq_scale.CITE_score", CITE_score)

############# plot scores versus expression of CITE seq
fileout1="CITE_expression_scores_correlation.pdf"
Plot_CITE_scores(CITE_score, CITE_data_all, CITE_seq_choose, orig.ident, orig.idents, fileout1,CITE_seq_expression)

############# identify doublets/multiplets based on mutually exclusive markers 
CITE_score = readRDS(file="CITE_seq_scale.CITE_score")
CITE_score = t(CITE_score) [cell_ids,]
CITE_genes = colnames(CITE_score)
CITE_genes1 = strsplit(CITE_genes,"_",fixed = T)
for(i in c(1:length(CITE_genes1))){CITE_genes1[i] = CITE_genes1[[i]][1]}
CITE_genes1 =unlist(CITE_genes1)
names(CITE_genes1) = CITE_genes

CITE_seq_short = strsplit(names(CITE_genes1),"_")
for(i in c(1:length(CITE_seq_short))){CITE_seq_short[i] = CITE_seq_short[[i]][1]}
CITE_seq_short = unlist(CITE_seq_short)
CITE_seq_short1 = names(CITE_genes1)
names(CITE_seq_short1) = CITE_seq_short

ref_file = "CITE_seq_mutually_exclusive_expression.txt"
p <- as.matrix(read.csv(ref_file, head=T, sep="\t"))
p=p[which(p[,3] %in% c("Mutually exclusive","Weak")),]
p=p[intersect(which(p[,1] !="CD56"),which(p[,2] !="CD56")), ]
############# check that all markers are in the CITE-seq list
p[,1][which(p[,1] %in% CITE_seq_short==F)]
p[,2][which(p[,2] %in% CITE_seq_short==F)]

mutual1 = CITE_seq_short1 [p[,1]]
mutual2 = CITE_seq_short1 [p[,2]]
mutually_exclusive_CITE_seq = NULL
mutually_exclusive_CITE_seq_names = NULL
for(i in c(1:length(mutual1))){
	if(concat(c(mutual2[i],":", mutual1[i])) %in% mutually_exclusive_CITE_seq_names==F){
		mutually_exclusive_CITE_seq = c(mutually_exclusive_CITE_seq, list(c(mutual1[i], mutual2[i])))
		mutually_exclusive_CITE_seq_names=c(mutually_exclusive_CITE_seq_names, concat(c(mutual1[i],":", mutual2[i])))
}}

############# normalise the variables that will be used in the predictor
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

############# functions: identify double positive cells
Boxplot_nUMI_ngene_mt_cells1<-function(doublets_list, pbmc, fileout1,metaD){
	headers = colnames(metaD)
	all_cells = rownames(metaD)
	library(RColorBrewer)
	cols = add.alpha (brewer.pal(8, "Dark2"), alpha = 0.5)
	w=2.5
	pdf(file=fileout1, height=w*1.3, width=w*4*1.8)
	par(mfrow= c(1,5), mar = c(12,4.5,3,1))
	for(h in c(1:length(headers))){
		groups = NULL
		for(i in c(1:length(doublets_list))){
			x = metaD[doublets_list[[i]], headers[h]]
			groups = c(groups, list(metaD[doublets_list[[i]], headers[h]]))}
		factors = names(doublets_list)
		boxplot(groups, ylab = "", col = cols,names= factors, las = 2, main= headers[h], outline =T, border = "grey",cex.lab = 0.9, lwd = 0.5, cex.main = 0.7)
	}
	dev.off()
}
Identify_CITE_seq_doublets <-function(fileout1, orig.idents, orig.ident, mutually_exclusive_CITE_seq_names,pbmc, CITE_genes1, CITE_score, metaD){
	library(MASS)
	library(RColorBrewer)
	cols =  add.alpha (brewer.pal(8, "Dark2"), alpha = 0.5)
	variables = metaD
	variables_names = colnames(variables)
	doublets_list_all = list()
	for(s in c(1:length(orig.idents))){
		w = cell_ids [which(orig.ident== orig.idents[s])]
		double_pos = NULL
		single_pos = NULL
		double_pos_classification = NULL
		sample_positive_CITE = names(which(apply(CITE_score[w,],2,function(x){length(which(is.na(x)==F))})>0))
		for(cc in c(1:length(mutually_exclusive_CITE_seq))){
			cc1 = mutually_exclusive_CITE_seq[[cc]][1]
			cc2 = mutually_exclusive_CITE_seq[[cc]][2]
			cite_sub = CITE_genes1[which(CITE_genes1 %in% names(c(cc1,cc2)))]
			# plot_CITE = names(which(rowSums(raw_data[names(cite_sub),w])>0))
			plot_CITE = intersect(sample_positive_CITE, names(cite_sub))
			if(length(plot_CITE)==2){
				x = CITE_score[w,plot_CITE[1]]
				y = CITE_score[w,plot_CITE[2]]
				w1 = w[intersect(which(is.finite(x)),which(is.finite(y)))]
				w1 = w1[intersect(which(x[w1]>0.01), which(y[w1]>0.01))]
				w2 = w[setdiff(which(is.finite(x)),which(is.finite(y)))]
				w3 = w[setdiff(which(is.finite(y)),which(is.finite(x)))]
				w23 = sort(unique(c(w2,w3)))
				if(min(c(length(w1),length(w23)))>2){
					double_pos = c(double_pos, w1)
					single_pos = c(single_pos, w23)
					double_pos_classification = c(double_pos_classification, rep(mutually_exclusive_CITE_seq_names[cc], length(w1)))
				}}
				}
		single_pos = setdiff(w, double_pos)
		doublets_list = NULL
		names = mutually_exclusive_CITE_seq
		var = variables[,"nCount_RNA"]
		names(var) = rownames(variables)
		summary = NULL
		for(cc in c(1:length(mutually_exclusive_CITE_seq))){
			x = var [double_pos[which(double_pos_classification== mutually_exclusive_CITE_seq_names[cc])]]
			doublets_list = c(doublets_list, list(names(x)))
			if(length(summary)==0){summary = summary(x)
			}else{summary = rbind(summary, summary(x))}
		}
		x = var [single_pos]
		summary = rbind(summary, summary(x))
		doublets_list = c(doublets_list, list(single_pos))		
		rownames(summary) = c(mutually_exclusive_CITE_seq_names,"singlets")
		names(doublets_list)= c(mutually_exclusive_CITE_seq_names,"singlets")
		w = intersect(which(summary[,'Median']<1.5*summary['singlets','Median']), which(rownames(summary)!="singlets"))
		w1 = intersect(which(summary[,'Median']<1.5*summary['singlets','Mean']), which(rownames(summary)!="singlets"))
		w = sort(unique(c(w,w1)))
		exclude_CITE = rownames(summary)[w]
		doublets_list[w] = list(NULL) 
		######################## plot
		names(doublets_list) = gsub("_TotalSeqC","",names(doublets_list) )
		fileout1 =concat(c(out_dir,PLOTS,"/Seurat_CITE_nGene_CITEseq_doublets_", type,"_", batch,"_", orig.idents[s],".pdf"))
		Boxplot_nUMI_ngene_mt_cells1(doublets_list, pbmc, fileout1, variables)
		print(concat(c("scp -p mfj169@rescomp1.well.ox.ac.uk:", fileout1," ./ " )))
		
		w = which(double_pos_classification %in% exclude_CITE==F)
		doublet_info = c(list(cbind(double_pos[w], double_pos_classification[w])),list(single_pos), list(summary))
		names(doublet_info) = c("true_doublets", "singlets","summary_nUMIs")
		doublets_list_all = c(doublets_list_all, list(doublet_info))
	}
	names(doublets_list_all) = orig.idents
	return(doublets_list_all)
}
Combine_CITE_seq_and_VDJ_doublets<-function(orig.idents, orig.ident, doublets_list_all, VDJ_doublets_total, cell_ids, metaD){
	doublets_CITE_plus_VDJ = NULL
	for(s in c(1:length(orig.idents))){
		w = cell_ids [which(orig.ident== orig.idents[s])]
		d = sort(unique(c(doublets_list_all[[s]]$ true_doublets[,"double_pos"], intersect(VDJ_doublets_total, w))))
		doublets_CITE_plus_VDJ = c(doublets_CITE_plus_VDJ,list(d))
	}
	names(doublets_CITE_plus_VDJ) = orig.idents
	return(doublets_CITE_plus_VDJ)
}

############# identify double positive cells

doublets_list_all =Identify_CITE_seq_doublets(fileout1, orig.idents, orig.ident, mutually_exclusive_CITE_seq_names,pbmc, CITE_genes1, CITE_score, metaD)

types = NULL
for(i in c(1:length(doublets_list_all))){types = sort(unique(c(types, doublets_list_all[[i]][[1]][,2])))}
list_doublets_CITE = rep(list(c()), length(types))
names(list_doublets_CITE) = types
all_CITE_doublets = NULL
for(i in c(1:length(doublets_list_all))){
	if(length(doublets_list_all[[i]][[1]])>0){
		all_CITE_doublets = c(all_CITE_doublets, doublets_list_all[[i]][[1]][,1])}
	for(t in types){
		w = which(doublets_list_all[[i]][[1]][,2]==t)
		list_doublets_CITE[[t]] = c(list_doublets_CITE[[t]], doublets_list_all[[i]][[1]][w,1])}}
all_CITE_doublets= sort(unique(all_CITE_doublets))
saveRDS(file="CITE_doublets.rds")), all_CITE_doublets)

############# plot UMAP locations of CITE-seq doublets/multiplets
library(RColorBrewer)
fileout1="CITE_nGene_CITE_doublets_UMAP_overall.pdf"
w=2.5
pdf(file=fileout1, height=w*4, width=w*5)
par(mfrow= c(4,5), mar = c(4.5,4.5,1.5,1.5))
all_cells = cell_ids
for(ind in c(1:length(list_doublets_CITE))){
	cluster = rep("other", length(all_cells))
	names(cluster) = all_cells
	cluster[list_doublets_CITE[[ind]]] = names(list_doublets_CITE)[ind]
	clusters = sort(unique(cluster))
	if(length(clusters)>=2){
		cluster_num = match(cluster, clusters)
		cols=add.alpha(colorRampPalette(c('red',"grey"))(length(clusters)), alpha = 0.7)
		w = order(cluster_num)
		if(clusters[2]=="other"){w = order(cluster_num, decreasing = T)}
		cex1 = c(0.01,0.01)
		main = names(list_doublets_CITE)[ind]
		main = gsub("_TotalSeqC","", main)
		plot(xy[w,1], xy[w,2],pch = 21, col = cols[cluster_num[w]],bg = cols[cluster_num[w]] ,main =main, cex = cex1[cluster_num[w]],xlab = "umap dim1",ylab = "umap dim2")
		print (ind)
}}
dev.off()


############# plot nGenes, nUMIs and mt% between these methods
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

fileout1 ="CITE_nGene_CITE_doublets.pdf"
metaD = variables_norm
names(list_doublets_CITE) = gsub("_TotalSeqC","",names(list_doublets_CITE))
Boxplot_nUMI_ngene_mt_doublets(list_doublets_CITE, pbmc, fileout1,metaD)

######################## Write output
saveRDS(file="Seurat_CITE_doublets.rds", list_doublets_CITE)


