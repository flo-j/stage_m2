import os
import numpy as np
import matplotlib.pyplot as plt
import time
from sklearn.decomposition import PCA, IncrementalPCA
from sklearn.cluster import AgglomerativeClustering
import scipy
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
#kmer_file=open("Homo_sapiens_HGSC.fa.regions.fst_np_L30.monomers.noN.fst_149.161_182.fst.kmers5.withoutfirstcol2","r")
kmer_file=open("../data/nt011630.rev.monomers.fst_149.length_166_177.kmers5.withoutfirstcol","r")
kmer_table=get_kmer_table(kmer_file)
time2=time.time()
print("obtention kmers : "),
print(time2-time1)
print("pca "),
res_pca=compute_pca(kmer_table)
time3=time.time()
print(time3-time2)
print("distance "),
distance=compute_dist(res_pca[:,range(0,100)],"euclidean") # Pour avoir les 100 premieres composante de l'acp
time4=time.time()
print(time4-time3)
print("clustering "),

t8=time.time()
r=robjects.r
stats=importr("stats")
dis=stats.as_dist(distance)
res_clusteringR=compute_clustering_R(dis,"ward.D")
t9=time.time()
print(t9-t8)
print("temps total ")
print(t9-time1)

#cut=r.cutree(res_clusteringR,2)
#vec=[]
#for i in range(len(cut)):
    # if cut[i]==1:
    #     vec.append(i)
# sub_data=kmer_table[vec,:]
