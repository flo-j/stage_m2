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


def get_kmer_index(kmer_file):
    return np.genfromtxt(kmer_file, skip_header=1, usecols=0, dtype=str) # on retourne que la 1ere colonne

def get_kmer_table(kmer_file, nb_col):
    return np.genfromtxt(kmer_file, skip_header=1, usecols = range(1,nb_col+1))

def compute_pca(kmer_table):
    pca = PCA()
    return pca.fit_transform(kmer_table)

def compute_dist(res_pca, method):
    temp = scipy.spatial.distance.pdist(res_pca, method)
    return scipy.spatial.distance.squareform(temp)

def plot_pca(res_pca):
    for i in range(len(res_pca)):
        plt.scatter(res_pca[i,0],res_pca[i,2])
    plt.show()

def compute_clustering(res_pca, linkage):
    clustering = AgglomerativeClustering(linkage=linkage,n_clusters=2)
    return clustering.fit(res_pca)

def compute_clustering_fast(distance, linkage):
    return fastcluster.ward(distance)

def plot_clustering(cut, res_pca):
    plt.scatter(res_pca[:,0],res_pca[:,1],c=cut)
    plt.show()

def compute_clustering_R(distance, linkage):
    return r.hclust(d = distance, method = "ward.D")

def compute_pairmate(cut,distance):
    distance[distance==0.0] = 999
    res = []
    for i in range(len(distance)):
        j=np.argmin(distance[i,])
        res.append(str(cut[i])+str(cut[j]))
        #print i, " ",j, " ",distance[i,j]
    mp1=float(res.count('11'))/(res.count('11')+res.count('12'))
    mp2 = float(res.count('22'))/(res.count('21')+res.count('22'))
    return mp1,mp2

time1=time.time()
rpy2.robjects.numpy2ri.activate()
#kmer_file=open("Homo_sapiens_HGSC.fa.regions.fst_np_L30.monomers.noN.fst_149.161_182.fst.kmers5.withoutfirstcol2","r")
#kmer_file=open("../data/nt011630.rev.monomers.fst_149.length_166_177.kmers5","r")
kmer_file="../data/nt011630.rev.monomers.fst_149.length_166_177.kmers5"
kmer_name=get_kmer_index(kmer_file)
kmer_table=get_kmer_table(kmer_file,1024)
fil=[[]]
fil[0].append(kmer_name)
fil[0].append(kmer_table)
seuil=0.98
r=robjects.r
stats=importr("stats")
i=0
while(len(fil)!=0):
    temp=[]
    kmer=fil.pop()
    res=[[]]
    res_pca=compute_pca(kmer[1])
    if len(res_pca) > 100:
        distance=compute_dist(res_pca[:,range(0,100)],"euclidean") # Pour avoir les 100 premieres composante de l'acp
    else :
        distance=compute_dist(res_pca,"euclidean") # Pour avoir les 100 premieres composante de l'acp
    dis=stats.as_dist(distance)
    res_clusteringR=compute_clustering_R(dis,"ward.D")
    cut=r.cutree(res_clusteringR,2)
    mp1,mp2=compute_pairmate(cut,distance)
    group1,group2=[],[]
    print(mp1,mp2)
    if mp1>=seuil and mp2 >=seuil:
        for i in range(len(cut)):
            if cut[i]==1:
                group1.append(i)
            elif cut[i]==2:
                group2.append(i)
            # else ERROR
        print(len(group1)," ", len(group2))
        print(kmer[0])
        temp.append(kmer[0][group1])
        temp.append(kmer[1][group1,:])
        fil.append(temp)
        temp.append(kmer[0][group2])
        temp.append(kmer[1][group2,:])
        fil.append(temp)
    else:
        temp2=[]
        temp2.append(i)
        temp2.append(kmer[0])
        temp2.append(kmer[1])
        res.append(temp2)
        i=i+1

# kmer_table=get_kmer_table(kmer_file,1024)
# time2=time.time()
# print("obtention kmers : "),
# print(time2-time1)
# print("pca "),
# res_pca=compute_pca(kmer_table)
# time3=time.time()
# print(time3-time2)
# print("distance "),
# distance=compute_dist(res_pca[:,range(0,100)],"euclidean") # Pour avoir les 100 premieres composante de l'acp
# time4=time.time()
# print(time4-time3)
# print("clustering "),
#
# t8=time.time()
# r=robjects.r
# stats=importr("stats")
# dis=stats.as_dist(distance)
# res_clusteringR=compute_clustering_R(dis,"ward.D")
# t9=time.time()
# print(t9-t8)
# print("temps total ")
# print(t9-time1)
#
# cut=r.cutree(res_clusteringR,2)
# vec=[]
# for i in range(len(cut)):
#      if cut[i]==1:
#          vec.append(i)
#
#
# mp1,mp2=compute_pairmate(cut,distance)
