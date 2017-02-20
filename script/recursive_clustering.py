import os
import numpy as np
import matplotlib.pyplot as plt
import time
from sklearn.decomposition import PCA, IncrementalPCA
from sklearn.cluster import AgglomerativeClustering
import scipy
#import fastcluster
from rpy2 import *
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
import rpy2.robjects as robjects
def get_kmer_table(kmer_file):
    return np.genfromtxt(kmer_file,skip_header=1)

def compute_pca(kmer_table):
    pca=PCA()
    return pca.fit_transform(kmer_table)

def compute_dist(res_pca, method):
    temp=scipy.spatial.distance.pdist(res_pca, method)
    return scipy.spatial.distance.squareform(temp)

def plot_pca(res_pca):
    for i in range(len(res_pca)):
        plt.scatter(res_pca[i,0],res_pca[i,2])
    plt.show()

def compute_clustering(res_pca,linkage):
    clustering = AgglomerativeClustering(linkage=linkage,n_clusters=2)
    return clustering.fit(res_pca)

def compute_clustering_fast(distance,linkage):
    return fastcluster.ward(distance)

def plot_clustering(res_clustering,res_pca):
    plt.scatter(res_pca[:,0],res_pca[:,1],c=res_clustering.labels_,cmap="ocean")
    plt.show()

def compute_clustering_R(distance, linkage):
    return r.hclust(d = distance, method = "ward.D")

time1=time.time()
rpy2.robjects.numpy2ri.activate()
kmer_file=open("Homo_sapiens_HGSC.fa.regions.fst_np_L30.monomers.noN.fst_149.161_182.fst.kmers5.withoutfirstcol2","r")
#kmer_file=open("../data/nt011630.rev.monomers.fst_149.length_166_177.kmers5.withoutfirstcol","r")
#time1=time.time()
kmer_table=get_kmer_table(kmer_file)
time2=time.time()
print("obtention kmers : "),
print(time2-time1)
print("pca "),
res_pca=compute_pca(kmer_table)
time3=time.time()
print(time3-time2)
print("distance "),
distance=compute_dist(res_pca,"euclidean")
time4=time.time()
print(time4-time3)
print("clustering "),
#t1=time.time()
#res_clustering=compute_clustering_fast(distance, "ward")
#t2=time.time()
#print res_clustering

#print(t2-t1)

#t6=time.time()
#res_clust=compute_clustering(res_pca,"ward")
#t7=time.time()
#print res_clust
#print t7-t6
t8=time.time()
#res_clustR=compute_clustering_R(distance,"linkage")
r=robjects.r
stats=importr("stats")
dis=stats.as_dist(distance)
res_clusteringR=compute_clustering_R(dis,"ward.D")
t9=time.time()
print t9-t8
print("temps total ")
print t9-time1
# dt=np.genfromtxt(kmerfile,skip_header=1)
# pca=PCA()
# X_pca=pca.fit_transform(dt)
# for i in range(len(dt)):
#     plt.scatter(X_pca[i,0],X_pca[i,1])
#
# t1=time.time()
# kmerfile2=open("../data/Homo_sapiens_HGSC.fa.regions.fst_np_L30.monomers.noN.fst_149.161_182.fst.kmers5.withoutfirstcol","r")
# t2=time.time()
# print "temps ouverture fichier ",
# print t2-t1
# t3=time.time()
# dt2=np.genfromtxt(kmerfile2,skip_header=1)
# t4=time.time()
# print "temps creation genfromtxt ",
# print t4-t3
# pca=PCA()
# t5=time.time()
# X_pca2=pca.fit_transform(dt2)
# t6=time.time()
# print "temps creation pca ",
# print t6-t5
#
# # t7=time.time()
# # for i in range(len(dt2)):
# #     plt.scatter(X_pca2[i,0],X_pca2[i,1])
# # t8=time.time()
# # print t8-t7
#
# t9=time.time()
# distance=scipy.spatial.distance.pdist(X_pca2,'euclidean')
# t10=time.time()
# print "temps distance euclidean",
# print t10-t9
