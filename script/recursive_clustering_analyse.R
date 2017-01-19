#!/usr/bin/Rscript
require(optparse)
require(ape)
require(MASS)
#================================================================================
# Ce script permets de faire un clustering récursif sur des données en table (table de fréquence de kmer)

# Exemple d'appel:
# Rscript --no-save --no-restore --verbose recursive_clustering_analyse.R --kmer_file ../data/nt011630.rev.monomers.fst_149.length_166_177.kmers5 --matepair 0.96

option_list = list(
	make_option(c("--kmer_file"), action="store", default=NULL, type='character',
							help="INPUT kmer table file"),
	make_option(c("--matepair"), action="store", default=NULL, type='character',
							help="Set mate pair threshold"))

opt = parse_args(OptionParser(option_list=option_list))

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

clustering_2 = function(data,size,nb_comp){
	print("dim:data:clustering_2")
	print(dim(data))
	if (nrow(data) < size){
		sample_1_id = rownames(data)
	}
	else {
		sample_1_id = sample(rownames(data),size)
	}
	# prendre les 10% des valeurs propres
	sample_1_dist = dist(data[sample_1_id,1:nb_comp],method = "euclidean") # avant data[sample_1_id,1:100] sauf que probleme quand moins de 100 colonnes..
  hc_1 = hclust(d = sample_1_dist, method = "ward.D2")
	hc_subtree_n = cutree(hc_1, 2)

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

test_selfpairmate = function(h_clust, mp){
	#Prepare Data
	seuil = mp
	f.mat = as.matrix(h_clust$DIST)

	f.mat[f.mat == 0] = 999
	#print(f.mat)
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
	intra_AA = sum(cname_f=="11") /(sum(cname_f=="12")+sum(cname_f=="11"))
	intra_BB = sum(cname_f=="22") /(sum(cname_f=="21")+sum(cname_f=="22"))
	print(paste("AA",intra_AA,"BB",intra_BB, "threshold", seuil))
	if ( (intra_AA > seuil)&(intra_BB > seuil)){
		return(list("test" = 1, "intra_AA" = intra_AA, "intra_BB" =  intra_BB, "pairmate" = result$V1))
	}
	else{return(list("test" = 0, "intra_AA" = intra_AA, "intra_BB" =  intra_BB,"pairmate" = result$V1))}
}




clusterize_me=function(data, n=0 , wkfile,filename, lim ,matepair) {
	print("______________________________________________________")
	print("dim:data:clusterize_me")
	print(dim(data))
	mycol = c("#002b36","#dc322f")
	adj_col  = adjustcolor(col = mycol, alpha.f = 4/10)
	palette(adj_col)
	bcol = "#586e75"
	b_col  = adjustcolor(col = bcol, alpha.f = 4/10)
	size = lim
	check_kmer = colSums(data)
	threshold = sum(check_kmer != 0)
	number_row=dim(data)[1]
	clusters_alt=rep(1, number_row)
	names(clusters_alt)=row.names(data)
	fileid = filename
	sub_data = data[,check_kmer != 0]
	print("dim de sub_data")
	print(dim(sub_data))
	print(paste("START DEPTH =",n))
	pca_1 = prcomp(sub_data)

#	nb_comp=ncol(pca_1$x)/10
	if(ncol(pca_1$x)>100){
		nb_comp=100
	}else{
		nb_comp=ncol(pca_1$x)
	}
	clust_1 = clustering_2(pca_1$x, size = 20000,nb_comp)# size utile pour la lda
	print("Done")
	print("Clustering")
		clusters = clust_1$classes
	test = test_selfpairmate(clust_1, mp = matepair)
	maintitle = paste(test$intra_AA,test$intra_BB)
    labA = paste("PC1 (",round(pca_1$sdev[1],5)*100,"% )")
    labB = paste("PC2 (",round(pca_1$sdev[2],5)*100,"% )")

	par(mar=c(5,5,4,2)+0.1)
	#	plot(pca_1$x[, c(1,2)],main = maintitle,cex.axis = 1.25, cex.lab = 1.25, col=clust_1$classes,cex = 0.5 , pch = 16, xlab= labA, ylab = labB)
	plot(pca_1$x[, c(1,2)],main = maintitle,cex.axis = 1.25, cex.lab = 1.25, col = clust_1$classes,cex = 1 , pch = 16, xlab= labA, ylab = labB)
	plot(pca_1$x[, c(1,2)],main = maintitle,cex.axis = 1.25, cex.lab = 1.25,col = b_col,cex = 1 , pch = 16, xlab= labA, ylab = labB)
	dev.off()
print("###########################################")
	print(fileid)
	print("dim of data:")
	print(dim(data))
	print("proportions:")
	print(paste(sum(clusters == 1),sum(clusters == 2)))
	print("Pair Mate frequency:")
	print(paste(test$intra_AA,test$intra_BB))

	if (test$test) {
		print("test ok")
		clusters1=clusterize_me(data[names(clusters[clusters==1]),], n=n+1,wkfile, filename = paste(fileid,n, "A", sep="_"),lim = size,matepair=matepair)
		clusters1bis=paste(clusters1, n, "A", sep="_")
		names(clusters1bis)=names(clusters1)
		clusters2=clusterize_me(data[names(clusters[clusters==2]),], n=n+1,wkfile, filename=paste(fileid,n, "B", sep="_"),lim = size,matepair=matepair)
		clusters2bis=paste(clusters2, n, "B", sep="_")
		names(clusters2bis) = names(clusters2)
		return(c(clusters1bis, clusters2bis))
	}

	else {
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
	options(warn=1)
	kmer_table = get_kmer_table(opt$kmer_file)
	size = nrow(kmer_table)
	kmer_id = kmer_table$id
	write.table(x = kmer_id,file = paste("results/",species,".kmer_table_id.txt",sep=''))
	all_in = clusterize_me(data = kmer_table , n=0 , wkfile = species, filename = file1 , lim =size, matepair = as.numeric(opt$matepair))
	filename = paste("results/",species,"mp",opt$matepair,".clustering_done.txt",sep='')
	write.table(x = all_in, file = filename)
}
