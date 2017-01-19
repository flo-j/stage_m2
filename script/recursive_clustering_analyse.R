#!/usr/bin/Rscript
require(optparse)
require(ape)
require(MASS)
#require(fields)
#================================================================================
# Ce script permets de faire un clustering récursif sur des données en table (table de fréquence de kmer)

# Exemple d'appel:
# Rscript --no-save --no-restore --verbose recursive_clustering_analyse.R --kmer_file ../data/nt011630.rev.monomers.fst_149.length_166_177.kmers5 --matepair 0.96
#================================================================================
# ARGUMENTS PARSER
#================================================================================
option_list = list(
	make_option(c("--kmer_file"), action="store", default=NULL, type='character',
							help="INPUT kmer table file"),
	make_option(c("--matepair"), action="store", default=NULL, type='character',
							help="Set mate pair threshold"))

opt = parse_args(OptionParser(option_list=option_list))
#================================================================================
# BUILDING-BLOCKS FUNCTION
#================================================================================
mysplit=function(str1){
  split_1 = gsub(pattern = ".No.*found", replacement = '' , x = str1)
  return(split_1)
}

dist2dataframe <- function(inDist) {
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  data.frame(
    row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    col = rep(B[-length(B)], (length(B)-1):1),
    value = as.vector(inDist))
}

get_kmer_table = function(input_file){
	kmer_table = read.csv(input_file, sep="", stringsAsFactors=FALSE, row.names = 1)
	return (kmer_table)
}

normalize_kmer = function (col_kmer){
  return ( ((col_kmer)-mean(col_kmer)) / sd(col_kmer))
}

input_species = function(input_file){
	sub_str = strsplit(input_file,split = "/")
	sub_str = unlist(sub_str)
	sub_str = sub_str[length(sub_str)]
	species = strsplit(sub_str,split = "\\.fa\\.|\\.fasta\\.")
	species = unlist(species)
	species = species[1]
	return (species)
}

clustering_1 = function(data){
	kmeans_1 = kmeans(data,centers = 2,iter.max = 500,nstart=100)
	return(kmeans_1$cluster)
}

clustering_2 = function(data,size){
	print("Taille de la table de kmers dans clustering_2")
	print(dim(data))
	print(data)
	if (nrow(data) < size){
		sample_1_id = rownames(data)
	}
	else {
		#		set.seed(1)
		sample_1_id = sample(rownames(data),size)
	}
	# 	sample_1_dist_matrix = dist(data[sample_1_id,1:100],method = "euclidean")
	print("Distance computation")
	#	sample_1_dist_matrix = rdist(data[sample_1_id,1:100])
	#	rownames(sample_1_dist_matrix) = sample_1_id
	#	colnames(sample_1_dist_matrix) = sample_1_id
	#	print("dist object")
	#	sample_1_dist= dist(sample_1_dist_matrix)
	if(nrow(data)<=100){
	sample_1_dist = dist(data[sample_1_id,1:100],method = "euclidean") # avant data[sample_1_id,1:100] sauf que probleme quand moins de 100 colonnes..
	print("HCLUST")
  	hc_1 = hclust(d = sample_1_dist, method = "ward.D2")
	hc_subtree_n = cutree(hc_1, 2)
}else{
	sample_1_dist = dist(data[sample_1_id,1:100],method = "euclidean") # avant data[sample_1_id,1:100] sauf que probleme quand moins de 100 colonnes..
	print("HCLUST")
		hc_1 = hclust(d = sample_1_dist, method = "ward.D2")
	hc_subtree_n = cutree(hc_1, 2)
}

	print("LDA")

	if (nrow(data) < size){
		return(list("classes" = hc_subtree_n, "subtree"  = hc_subtree_n, "DIST" = sample_1_dist))
	}
	else{
		lda_1 = lda(x = data[sample_1_id,1:100],grouping = hc_subtree_n)
		pred_lda = predict(lda_1, data[,1:100])
		pred_lda$class = as.numeric(pred_lda$class)
		names(pred_lda$class) = rownames(data)

		return(list("classes" = pred_lda$class, "subtree"  = hc_subtree_n, "DIST" = sample_1_dist))
	}
}

write_sample_id = function(data,wkfile){
	filename = paste(wkfile,".sample.txt", sep = '')
	write.table(x = data , file = filename, quote = F,row.names = F,col.names = F)
}

write_sample_id2 = function(data,wkfile,output){
	filename = paste(wkfile,".sample_",output,".txt", sep = '')
	write.table(x = data , file = filename, quote = F,row.names = F,col.names = F)
}

fastacmd = function(wkfile){
	cmd = paste("-i ",wkfile,".sample.txt -d ",opt$formatdb_file,"> ",wkfile,".sample.fa",sep = '')
	print(cmd)
	system2(command = "fastacmd", args = cmd)
}

run_muscle = function(wkfile){
	file_IN = paste(wkfile,".sample.fa",sep = '')
	file_OUT= paste(wkfile,".sample.aligned.fa",sep = '')
	cmd = paste("muscle -in ",file_IN," -out ",file_OUT, sep = '')
	print(cmd)
	system(cmd, ignore.stdout=T, ignore.stderr=T)
}

get_distance_K2P = function(wkfile){
	file_IN = paste(wkfile,".sample.aligned.fa",sep = '')
	dna_aligned = read.dna(file = file_IN,format = "fasta")
	dna_dist    = dist.dna(x = dna_aligned,model = "K80")
	print("DEBUG")
	print(sum(is.na(dna_dist)))
	dna_dist.df = dist2dataframe(dna_dist)
	dna_dist.df[,1] = as.character(dna_dist.df[,1])
	dna_dist.df[,2] = as.character(dna_dist.df[,2])
	write.table(x=dna_dist.df,file = "debug_getdistance_dist_file")
	return(dna_dist.df)
}


#================================================================================
#
#================================================================================
#evaluate_distance = function(sample_distance,wkfile,clusters.data.frame){
#	dist_file = sample_distance
#	clust_file = clusters.data.frame
#	dist_file$row = sapply(X = dist_file$row,FUN = mysplit)
#	dist_file$col = sapply(X = dist_file$col,FUN = mysplit)
#	tag2 = c()
#	for (i in 1:nrow(dist_file)){
#		name_A = dist_file[i,1]
#		name_B = dist_file[i,2]
#		clu_A = clust_file[name_A,]
#		clu_B = clust_file[name_B,]
#		if (clu_A < clu_B){
#			cur_tag = paste(as.character(clu_A),as.character(clu_B),sep = '_')
#		}
#		else{
#			cur_tag = paste(as.character(clu_B),as.character(clu_A),sep = '_')
#		}
#		tag2 = c(tag2,cur_tag)
#	}
#	dist_file_all = cbind(dist_file,tag2)
#	dist_1_1 = dist_file_all[dist_file_all$tag == "1_1",]
#	dist_1_2 = dist_file_all[dist_file_all$tag == "1_2",]
#	dist_2_2 = dist_file_all[dist_file_all$tag == "2_2",]
#	print(paste("INTRA_1 = ",intra1,"; INTRA_2 = ",intra2," INTER = ",inter,sep=''))
#	if ( (intra1 < inter)&(intra2 < inter) ) {return(1)}
#	else {return(0)}
#sum_dd = tapply(dist_file_all$value, dist_file_all$tag2, median)
#	print(sum_dd)
#	ii = "1"
#	jj = "2"
#	intra1=sum_dd[paste(ii, ii,sep ="_")]
#	intra2=sum_dd[paste(jj, jj,sep ="_")]
#	inter=sum_dd[paste(ii,jj,sep ="_")]
#	print(paste("INTRA_1 = ",intra1,"; INTRA_2 = ",intra2," INTER = ",inter,sep=''))
#	if (intra1 < inter){
#		if (intra2 > inter){return(0)}
#		else(return(1))
#	}
##	if (intra2 < inter){
#		if (intra1 > inter){return(0)}
#		else(return(1))
#	}
#	else{return(0)}
#}
#================================================================================
#
#================================================================================
test_selfpairmate = function(h_clust, mp){
	#Prepare Data
	print("Test self pair mate : preparation")
	seuil = mp
	f.mat = as.matrix(h_clust$DIST)
	#	rownames(f.mat) = sample_id
	#	colnames(f.mat) = sample_id
	f.mat[f.mat == 0] = 999
	#	diag(f.mat) = 999
	result <- t(sapply(seq(nrow(f.mat)), function(i) {
 		 j <- which.min(f.mat[i,])
  		c(paste(rownames(f.mat)[i], colnames(f.mat)[j], sep='/'), f.mat[i,j])
		}))
	result = as.data.frame(result)
	result$V1 = as.character(result$V1)
	fhcsub.df = as.data.frame(h_clust$subtree,row.names = names(h_clust$subtree))
	cname = strsplit(result$V1,split = "/")
	cname = as.data.frame(cname)
	cname_A = as.data.frame(t(cname[1,]))
	cname_A = fhcsub.df[as.character(cname_A$`1`),]
	cname_B = as.data.frame(t(cname[2,]))
	cname_B = fhcsub.df[as.character(cname_B$`2`),]
	cname_f = paste(cname_A,cname_B,sep = '')
	#	sum_inter1 = sum(cname_f=="12") / (sum(cname_f=="12")+sum(cname_f=="11"))
	#	sum_inter2 = sum(cname_f=="21") / (sum(cname_f=="21")+sum(cname_f=="22"))
	#	(sum_inter1+sum_inter2)/2
	intra_AA = sum(cname_f=="11") /(sum(cname_f=="12")+sum(cname_f=="11"))
	intra_BB = sum(cname_f=="22") /(sum(cname_f=="21")+sum(cname_f=="22"))
	print(paste("AA",intra_AA,"BB",intra_BB, "threshold", seuil))
	if ( (intra_AA > seuil)&(intra_BB > seuil)){
		return(list("test" = 1, "intra_AA" = intra_AA, "intra_BB" =  intra_BB, "pairmate" = result$V1))
	}
	else{return(list("test" = 0, "intra_AA" = intra_AA, "intra_BB" =  intra_BB,"pairmate" = result$V1))}
}

simulation_clustering = function(hc_subtree, pairmate){
  randomcluster = as.character(names(hc_subtree))
  p1 = sum(hc_subtree == 1)
  p2 = length(hc_subtree)-p1
  clu1r = rep(x = 1,p1)
  clu2r = rep(x = 2,p2)
  clu_random = c(clu2r,clu1r)
  clu_r = sample(clu_random)
  print("==")
  clu_random = cbind(randomcluster,as.character(clu_r))
  clu_random = as.data.frame(clu_random)
  rownames(clu_random) = clu_random$V1
  cname = strsplit(as.character(pairmate),split = "/")
  cname = as.data.frame(cname)
  cname_A = as.data.frame(t(cname[1,]))
  cname_A = clu_random[cname_A$`1`,]
  cname_B = as.data.frame(t(cname[2,]))
  cname_B = clu_random[cname_B$`2`,]
  cname_f = paste(cname_A$V2,cname_B$V2,sep = '')
  intra1 = sum(cname_f=="11") /(sum(cname_f=="12")+sum(cname_f=="11"))
  intra2 = sum(cname_f=="22") /(sum(cname_f=="21")+sum(cname_f=="22"))
  return(list("intra1" = intra1, "intra2" = intra2, "clur" = clu_random, "cnamef" = cname_f))
}

simulation_clustering2 = function(hc_subtree,sample_dist){
  randomcluster = as.character(names(hc_subtree))
  p1 = sum(hc_subtree == 1)
  p2 = length(hc_subtree)-p1
  clu1r = rep(x = 1,p1)
  clu2r = rep(x = 2,p2)
  clu_random = c(clu2r,clu1r)
  clu_r = sample(clu_random)
  names(clu_r) = randomcluster
  sil_obs = silhouette(x=hc_subtree, dist = sample_dist)
  sil_ran = silhouette(x=clu_r, dist = sample_dist)
  return(list("sil_obs"= sil_obs,"sil_ran"=sil_ran))
}

run_test = function(sample,wkfile,clusters.data.frame){
	print("WRITE FASTA FILE from sample sequences")
	fastacmd(wkfile)
	print("RUN MUSCLE for sample sequences")
	run_muscle(wkfile)
	print("Use ape package for sample sequences")
	sample_distance = get_distance_K2P(wkfile)
	print("TEST continue or not ?")
	test = evaluate_distance(h_clust)
	return(test)
}



plot_pca_noclu = function(pca_data,titlename,n){
  base_col = "#002b36"
  adj_col  = adjustcolor(col = base_col, alpha.f = 4/10)
  file_name = paste("./results/",titlename,"_pca_",n,".png",sep='')
  png(filename = file_name,width = 1000, height = 800)
  par(mfrow=c(2,2),mar=c(5,5,4,2)+0.1)
  plot(pca_data$sdev[1:100],main = "",type="h",xlab="PCn",ylab="Variances",cex.axis = 1.5, cex.lab = 1.5)
	#  title(main = paste(titlename))
  plot(pca_data$x[,c(1,2)], pch=16,cex=4/10,col = adj_col,cex.axis = 1.5, cex.lab = 1.5)
  plot(pca_data$x[,c(1,3)], pch=16,cex=4/10,col = adj_col,cex.axis = 1.5, cex.lab = 1.5)
  plot(pca_data$x[,c(1,4)], pch=16,cex=4/10,col = adj_col,cex.axis = 1.5, cex.lab = 1.5)
  dev.off()
}
#================================================================================
# MAIN FUNCTION / RECURSIVE CLUSTERING
#================================================================================
clusterize_me=function(data, n=0 , wkfile,filename, lim ,matepair) {
	print("MATEPAIR")
	print(matepair)
	mycol = c("#002b36","#dc322f")
	adj_col  = adjustcolor(col = mycol, alpha.f = 4/10)
	palette(adj_col)
	bcol = "#586e75"
	b_col  = adjustcolor(col = bcol, alpha.f = 4/10)
	size = lim
	check_kmer = colSums(data)
	#	threshold = size/100
	#	threshold = 1024
	threshold = sum(check_kmer != 0)
	number_row=dim(data)[1]
	clusters_alt=rep(1, number_row)
	names(clusters_alt)=row.names(data)
	fileid = filename
	#if (number_row < threshold){
	#	print(paste("Too small", threshold))
		#		size_A = sum(clusters == 1)
		#		size_B = sum(clusters == 2)
		#		maintitle = paste(size_A,size_B)
		#		filename2 = paste(fileid,"E",sep='_')
		#		png(paste(filename2, ".png", sep=""))
		#		par(mar=c(5,5,4,2)+0.1)
		#plot(pca_1$x[, c(1,2)],main =maintitle,cex.axis = 1.25, cex.lab = 1.25, col=clust_1$classes, cex = 0.5 , pch = 16 , xlab= labA, ylab = labB)
		#		plot(pca_1$scores[, c(1,2)],main =maintitle,cex.axis = 1.25, cex.lab = 1.25, col=clust_1$classes, cex = 0.5 , pch = 16 , xlab= labA, ylab = labB)
		#		dev.off()
	#	return(clusters_alt)
	#}

	sub_data = data[,check_kmer != 0]
	print("dim de sub_data")
	print(dim(sub_data))
	print(paste("START DEPTH =",n))
	#1-Make Clusters
	#	print("dim of data:")
	#	print(dim(data))
	#	write.table(x = data,file = paste("./results/",species,".data_pcabug.txt",sep=''))
	#	pca_1 = prcomp(data)
	#	pca_1 = princomp(data)

	pca_1 = prcomp(sub_data)
	print(dim(pca_1$x))
	print("Done")
	print("Clustering")
	#	clust_1 = clustering_2(pca_1$x, size = 20000)
	clust_1 = clustering_2(pca_1$x, size = 20000)# size utile pour la lda
	print("length pca : ")
	print(dim(pca_1$x))
	print("Done")
	#	clusters = clustering_1(pca_1[,1:100])
	print("Test Clusters")
	clusters = clust_1$classes
	print("Pair Mate")
	test = test_selfpairmate(clust_1, mp = matepair)
 	print ("Test Done")
	maintitle = paste(test$intra_AA,test$intra_BB)
    labA = paste("PC1 (",round(pca_1$sdev[1],5)*100,"% )")
    labB = paste("PC2 (",round(pca_1$sdev[2],5)*100,"% )")
	print(fileid)
	png(paste(fileid, ".%02d.png", sep=""))
	par(mar=c(5,5,4,2)+0.1)
	#	plot(pca_1$x[, c(1,2)],main = maintitle,cex.axis = 1.25, cex.lab = 1.25, col=clust_1$classes,cex = 0.5 , pch = 16, xlab= labA, ylab = labB)
	plot(pca_1$x[, c(1,2)],main = maintitle,cex.axis = 1.25, cex.lab = 1.25, col = clust_1$classes,cex = 0.5 , pch = 16, xlab= labA, ylab = labB)
	plot(pca_1$x[, c(1,2)],main = maintitle,cex.axis = 1.25, cex.lab = 1.25,col = b_col,cex = 0.5 , pch = 16, xlab= labA, ylab = labB)
	dev.off()

	#	print("simul")
	#	aX = c()
	#	bX = c()
	#	cluR = c()
	#	for (i in 1:100){
	#  	tets = simulation_clustering(hc_subtree = clusters,pairmate = test$pairmate)
	#  		aX = c(aX,tets$intra1)
	# 		bX = c(bX,tets$intra2)
	#	}
	#essai1 = cbind(aX,bX)

	print("###########################################")
	print(fileid)
	print("dim of data:")
	print(dim(data))
	print("proportions:")
	print(paste(sum(clusters == 1),sum(clusters == 2)))
	print("Pair Mate frequency:")
	print(paste(test$intra_AA,test$intra_BB))
	#	print("Simulation:")
	#	print(summary(essai1))
	print("###########################################")


	#  	file_name = paste("./",wkfile,"_pca_",n,"clust.png",sep='')
	# 	png(filename = file_name,width = 1000, height = 1000)
	#  	par(mfrow=c(2,2),mar=c(5,5,4,2)+0.1)
	#  	plot(pca_1$sdev[1:100],main = "",type="h",xlab="PCn",ylab="Variances",cex.axis = 1.5, cex.lab = 1.5)
	#	sample_id = names(clust_1$subtree)
	#  	plot(pca_1$x[sample_id,c(1,2)], pch=16,cex=4/10,col = clust_1$classes,cex.axis = 1.5, cex.lab = 1.5)
	#  	plot(pca_1$x[sample_id,c(1,3)], pch=16,cex=4/10,col = clust_1$classes,cex.axis = 1.5, cex.lab = 1.5)
	#	 	plot(pca_1$x[sample_id,c(1,4)], pch=16,cex=4/10,col = clust_1$classes,cex.axis = 1.5, cex.lab = 1.5)
	#  	dev.off()



	if (test$test) {
		print("Test Pairmate OK")
		print("###########################################")
		print(length(clusters[clusters==1]))
		print("###########################################")
		clusters1=clusterize_me(data[names(clusters[clusters==1]),], n=n+1,wkfile, filename = paste(fileid,n, "A", sep="_"),lim = size,matepair=matepair)
		clusters1bis=paste(clusters1, n, "A", sep="_")
		names(clusters1bis)=names(clusters1)

		print("###########################################")
		print(length(clusters[clusters==2]))
		print("###########################################")
		clusters2=clusterize_me(data[names(clusters[clusters==2]),], n=n+1,wkfile, filename=paste(fileid,n, "B", sep="_"),lim = size,matepair=matepair)
		clusters2bis=paste(clusters2, n, "B", sep="_")
		names(clusters2bis) = names(clusters2)
		return(c(clusters1bis, clusters2bis))
	}

	else {
		print("Test Pairmate fail")
		return(clusters_alt)
	}
}
#================================================================================
# MAIN
#================================================================================
print("debut du programme")
if ( !(is.null(opt$kmer_file) ) ) {
	species = input_species(opt$kmer_file)
	dir.create("results")
	dir.create("results/plot")
	file1 = paste("results/plot/",species,"_plot_pca",sep = '')
	print("START")
	print(species)
	options(warn=1)
	print("1-Reading File")
	kmer_table = get_kmer_table(opt$kmer_file)
	size = nrow(kmer_table)
	print("Make PCA")
#	pca_1 = prcomp(kmer_table)
#	kmer_norm = apply(X = kmer_table,MARGIN = 2,FUN = normalize_kmer)
	print(nrow(kmer_table))
	print("Done")
	kmer_id = kmer_table$id
	write.table(x = kmer_id,file = paste("results/",species,".kmer_table_id.txt",sep=''))
	print("2-RUN Clustering")
	all_in = clusterize_me(data = kmer_table , n=0 , wkfile = species, filename = file1 , lim =size, matepair = as.numeric(opt$matepair))
	filename = paste("results/",species,"mp",opt$matepair,".clustering_done.txt",sep='')
	write.table(x = all_in, file = filename)
	print("END")
}
