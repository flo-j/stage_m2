import os, sys
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
import argparse
import logging
import fastcluster
import scipy.cluster.hierarchy

def get_kmer_table(kmer_file):
    ''' get the names of each sequence which correspond to the kmer_file except
        the first column. Column names are lost.
        ARG    : - kmer_file : the filename with the kmer table

        RETURN : - a 2D numpy.array with the kmer frequencies
    '''
    kfile=open(kmer_file,"r")
    line=kfile.readline()
    nb_col=len(line.split())-1
    kfile.close()
    if verbose > 2 :
        logging.info("get kmer table")
    return np.genfromtxt(kmer_file, skip_header=1, usecols = range(1,nb_col+1))

def get_sequences_name(kmer_file):
    ''' get the names of each sequence which correspond to the first column of
        the kmer_file
        ARG    : - kmer_file : the file with the kmer table
        RETURN : - a 1D numpy.array with the sequences' name
    '''
    if verbose > 2:
        logging.info("get sequences name")
    return np.genfromtxt(kmer_file, skip_header=1, usecols=0, dtype=str)

def compute_pca(kmer_table):
    ''' compute a Principal component analysis of the kmer frequencies
        ARG    : - kmer_table : an 2D array with the kmer frequencies
        RETURN : - An 2D numpy.array with the eigenvalues
    '''
    pca = PCA()
    print("compute pca")
    if verbose > 2:
        logging.info("compute pca")
    return pca.fit_transform(kmer_table)


def compute_dist2(res_pca, method = "euclidean"):
    ''' compute the distance between each sequence using the eigenvalues of the
        PCA. By default the method used is "euclidean"
        ARG    : - res_pca : an 2D array with the eigenvalues
                 - method  : the method used to compute the distance
        RETURN : - An 2D array with pairwise distances between sequences
    '''
    t1=time.clock()
    if verbose > 2:
        logging.info("compute distances with "+ method + " method")
    d =  scipy.spatial.distance.pdist(res_pca, method)
    t2=time.clock()
    print(t2-t1)
    return d

def compute_clustering_fast(distance):
    t1=time.clock()
    c=fastcluster.ward(distance)
    t2=time.clock()
    print(t2-t1)
    return c

def plot_clustering( res_pca, group1, group2, output):
    ''' plot the 2 first dimension of the pca colored by cluster.
        ARG    : - res_pca        : an 2D array with the eigenvalues
                 - group1, group2 : vector of indice
                 - output         : the name of the output
        NO RETURN
    '''
    output+="clustering.png"
    if verbose > 2:
        logging.info("plot clustering")
    plt.figure()
    for ind in range(len(res_pca)):
        if ind in group1:
            plt.scatter(res_pca[ind,0],res_pca[ind,1],c="r",alpha=0.2)
        elif ind in group2:
            plt.scatter(res_pca[ind,0],res_pca[ind,1],c="b",alpha=0.2)
    plt.savefig(output)
    plt.close()

def create_2_groups(cut):
    ''' Create two groups following clustering results.
        ARG    : - cut : a rpy2 IntVector with the repartition of the 2 groups
        RETURN : 2 lists one for each group
    '''
    plot_output="plot_"+str(i)
    if verbose > 1:
        logging.info("matepair test ok")
    g1,g2 = [],[]
    for x in range(len(cut)):
        (g1,g2)[cut[x] == 1].append(x)
    if verbose > 0:
        logging.info("length "+str(len(g1))+" "+str(len(g2)))
    return g1, g2


def compute_pairmate(cut, distance):
    ''' compute the pairmate for the two groups.Pairmate corresponds to the
        proportion of sequences which has their nearest neighbor in a same group

        ARG    : - cut      : a rpy2IntVector with the repartition of the groups
                 - distance : An 2D array with distances between each sequence
        RETURN : - a tupple with the matepair of each group
    '''
    if verbose > 2:
        logging.info("compute pairmate")
    distance[distance<=0.00000000001] = 999
    res = []
    for i in range(len(distance)):
        j = np.argmin(distance[i,])
        res.append(str(cut[i])+str(cut[j]))
    mp1 = float(res.count('11'))/(res.count('11')+res.count('12'))
    mp2 = float(res.count('22'))/(res.count('21')+res.count('22'))
    return mp1,mp2

def get_args():
    ''' get the argument from the command line
        NO ARG
        RETURN : a tupple with the different args used in the programm
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairmate", help = "Set mate pair threshold.\
                          default = 0.95", type = float, default = 0.95)
    parser.add_argument("input", help = "INPUT kmer table file")
    parser.add_argument("--nbcomppca", help = "Number of composante of pca\
                         used to compute distance. Default = 100",type=int,
                         default = 100)
    parser.add_argument("-o", "--output", help = "name of the outputfile. \
                        defaut results/results.txt",
                        default="results/results.txt")
    parser.add_argument("--sorting",help = "you can sort the output by sequence\
                         or by cluster. default = cluster",
                        choices = ["cluster","sequence"], default ="cluster")
    parser.add_argument("--verbose","-v", help ="degree of verbose you want",
                        type=int, default=0)
    parser.add_argument("--pca", help ="yes if you want to plot the \
                        pca else : no ", default="no")
    args = parser.parse_args()
    return args.pairmate, args.input, args.nbcomppca, args.output, \
           args.sorting, args.pca, args.verbose


# get args
pairmate, kmer_file, nbcomp_pca, outputfile, sort, plot, verbose = get_args()
res = {}
queue = []

# initialize some variables
cluster_number = 0
i = 0

kmer_table=get_kmer_table(kmer_file)
queue.append(range(len(kmer_table)))

while(queue):
    kmer_indice = queue.pop(0) # pop the first value of the queue
    pca = compute_pca(kmer_table[kmer_indice])
    if verbose > 0:
        logging.info(str(len(pca))+" "+str(i))
    # to reduce computation time, we use just nbcomp_pca components if possible
    if len(pca) > nbcomp_pca:
        distance = compute_dist2(pca[:,range(0, nbcomp_pca)])
    else :
        distance=compute_dist2(pca)
    clust_fast = compute_clustering_fast(distance)
    clusters = scipy.cluster.hierarchy.fcluster(clust_fast, 2, criterion="maxclust")
    mp1,mp2 = compute_pairmate(clusters, scipy.spatial.distance.squareform(distance))
    print("mp1 "+str(mp1)+' mp2 '+str(mp2))
    group1,group2 = [],[]
    if verbose > 0:
        logging.info("mp "+str(mp1)+" "+str(mp2))
    if mp1 >= pairmate and mp2 >= pairmate :
        group1, group2 =create_2_groups(clusters)
        if plot=="yes":
            plot_pca(pca, "plot/"+str(i))
            plot_clustering(pca, group1, group2, "plot/"+str(i))
        queue.append(group1)
        queue.append(group2)
    else:
        cluster_number = cluster_number+1
        if verbose > 1:
            logging.info("matepair test non ok")
        res[cluster_number] = kmer_indice
    i+=1
if sort == 'sequence':
    save_results_sorted_by_sequence(outputfile, res, sequence_name)
else:
    save_results_sorted_by_cluster(outputfile, res, sequence_name)
