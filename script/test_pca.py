import os
import numpy as np
import matplotlib.pyplot as plt
import time
from sklearn.datasets import load_iris
from sklearn.decomposition import PCA, IncrementalPCA
import scipy

def get_kmer_table(kmer_file):
    return np.genfromtxt(kmerfile,skip_header=1)

def compute_pca(kmer_table):
    pca=PCA()
    return pca.fit_transform(kmer_table)

def compute_dist(res_pca, method):
    return scipy.spatial.distance.pdist(res_pca, method)

def plot_pca(res_pca):
    for i in range(len(res_pca)):
        plt.scatter(res_pca[i,0],res_pca[i,2])
    plt.show()


kmer_file=open("../data/nt011630.rev.monomers.fst_149.length_166_177.kmers5.withoutfirstcol","r")
kmer_table=get_kmer_table(kmer_file)
res_pca=compute_pca(kmer_table)
distance=compute_dist(res_pca,"euclidean")

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
