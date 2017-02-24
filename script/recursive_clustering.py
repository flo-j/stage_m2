# "A computer program does what you tell it to do, not what you want it to do"

#AUTHOR : Sebastien Bridel, Florence Jornod
#CONTACT ; <florence@jornod.com>

################################################################################
#This work is licensed under the Creative Commons Attribution-ShareAlike 3.0
#Unported License. To view a copy of this license,
#visit http://creativecommons.org/licenses/by-sa/3.0/
#or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA
################################################################################

#python3 recursive_clustering.py --pairmate 0.98 ../data/nt011630.rev.monomers.fst_149.length_166_177.kmers5 --output TESTRES4 --sorting sequence

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
import argparse

def get_sequences_name(kmer_file):
    ''' get the names of each sequences who correspond of the first column of
        the kmer_file
        ARG    : - kmer_file : the file with the kmer table
        RETURN : - a 1D numpy.array with the name of the sequences
    '''
    return np.genfromtxt(kmer_file, skip_header=1, usecols=0, dtype=str)

def get_kmer_table(kmer_file, nb_col=1024): # XXX a terme lire la premiere ligne pour savoir le nombre de colonnes
    ''' get the names of each sequences who correspond of the kmer_file except
        the first column. The name of each col is lost. By default the number
        of column is 1024, the number of 5-mers possible.
        ARG    : - kmer_file : the filename with the kmer table
                 - nb_col    : the number of columns in the file
        RETURN : - a 2D numpy.array with the kmer frequencies
    '''
    return np.genfromtxt(kmer_file, skip_header=1, usecols = range(1,nb_col+1))

def compute_pca(kmer_table):
    ''' compute a Principal composant analysis of the kmer frequencies
        ARG    : - kmer_table : an 2D array with the kmer frequencies
        RETURN : - An 2D numpy.array with the eigenvalues
    '''
    pca = PCA()
    return pca.fit_transform(kmer_table)

def plot_pca(res_pca):
    ''' plot the 2 first dimension of the pca.
        ARG    : - res_pca : an 2D array with the eigenvalues
        NO RETURN
    '''
    for i in range(len(res_pca)):
        plt.scatter(res_pca[i,0],res_pca[i,1])
    plt.show()

def compute_dist(res_pca, method = "euclidean"):
    ''' compute the distance between each sequence using the eigenvalues of the
        PCA. By default the method used is "euclidean"
        ARG    : - res_pca : an 2D array with the eigenvalues
                 - method  : the method used to compute the distance
        RETURN : - An 2D numpy.array with the distance between each sequence
    '''
    temp = scipy.spatial.distance.pdist(res_pca, method)
    return scipy.spatial.distance.squareform(temp)

def compute_clustering_R(distance, linkage = "ward.D"):
    ''' compute clustering using the distance between the sequence using ward.D
        method. Cut this clustering in 2 groups.
        ARG    : - distance : An 2D array with distances between each sequence
                 - linkage  : The method used for the clustering
                 - cutoff   : the number of cluster wanted
        RETURN : - a rpy2 IntVector with the repartition of the 2 groups
    '''
    distance = stats.as_dist(distance)
    clustering = r.hclust(d = distance, method = linkage)
    return r.cutree(clustering, 2)

def compute_pairmate(cut, distance):
    ''' compute the pairmate for each 2 group. pairmate correspond to the
        proportion of each sequence who have the nearest neighbor in the same
        group.
        ARG    : - cut      : a rpy2IntVector with the repartition of the groups
                 - distance : An 2D array with distances between each sequence
        RETURN : - a tupple with the matepair for each group
    '''
    distance[distance==0.0] = 999
    res = []
    for i in range(len(distance)):
        j = np.argmin(distance[i,])
        res.append(str(cut[i])+str(cut[j]))
    mp1 = float(res.count('11'))/(res.count('11')+res.count('12'))
    mp2 = float(res.count('22'))/(res.count('21')+res.count('22'))
    return mp1,mp2

def plot_clustering( res_pca, group1, group2, output):
    ''' plot the 2 first dimension of the pca colored by cluster.
        ARG    : - res_pca :  an 2D array with the eigenvalues
                 - group1, group2 : vector of indice
                 - output : the name of the output
        NO RETURN
    '''
    for i in range(len(pca)):
        if i in group1:
            plt.scatter(pca[i,0],pca[i,1],c="r",alpha=0.2)
        elif i in group2:
            plt.scatter(pca[i,0],pca[i,1],c="b",alpha=0.2)
    plt.savefig(output)

def save_results_sorted_by_cluster(filename, res, seq_name):
    ''' save the result in a file. Each sequence is associated to its group.
        the result is order by cluster. This function is faster than sorting by
        sequence.
        ARG    : - filename : the result filename
                 - res      : a dict with cluster as keys and the seq as value
                 - seq_name : a 1D numpy.array with the name of the sequences
        NO RETURN
    '''
    filename = open(filename,"w")
    for key in res.keys():
        for indice in res[key]:
            towrite = str(seq_name[indice])+"  "+str(key)+"\n"
            filename.write(towrite)

def save_results_sorted_by_sequence(filename, res, kmer_name):
    ''' save the result in a file. Each sequence is associated to its group.
        the result is order by sequence. Attention! This function is slower than
        sorting by cluster.
        ARG    : - filename : the result filename
                 - res      : a dict with cluster as keys and the seq as value
                 - seq_name : a 1D numpy.array with the name of the sequences
        NO RETURN
    '''
    filename = open(filename,"w")
    for i in range(len(kmer_name)):
        for key in res.keys():
            if i in res[key]:
                filename.write(kmer_name[i]+"   "+str(key)+"\n")

def get_args():
    ''' get the argument from the command line
        NO ARG
        RETURN : a tupple with the different arg used in the programm
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairmate", help = "Set mate pair threshold.\
                          default = 0.95", type = float, default = 0.95)
    parser.add_argument("input", help = "INPUT kmer table file")
    parser.add_argument("--nbcomppca", help = "Number of composante of pca\
                         used to compute distance. Default = 100",type=int,
                         default = 100)
    parser.add_argument("-o", "--output", help = "name of the outputfile. \
                        defaut results.txt",default="results.txt")
    parser.add_argument("--sorting",help = "you can sort the output by sequence\
                         or by cluster. default = cluster",
                        choices = ["cluster","sequence"], default ="cluster")
    args = parser.parse_args()
    return args.pairmate, args.input, args.nbcomppca, args.output, args.sorting


if __name__ == '__main__':

    pairmate_threshold, kmer_file, nb_comp_pca, outputfile, sorting = get_args()

    rpy2.robjects.numpy2ri.activate()

    kmer_name = get_sequences_name(kmer_file)
    kmer_table = get_kmer_table(kmer_file)
    res = {}
    fil = []
    fil.append(range(len(kmer_table)))
    r = robjects.r
    stats = importr("stats")
    cluster = 0
    while(fil):
        kmer_indice = fil.pop()
        pca = compute_pca(kmer_table[kmer_indice])
        if len(pca) > nb_comp_pca:
            distance = compute_dist(pca[:,range(0, nb_comp_pca)])
        else :
            distance=compute_dist(pca)
        cut = compute_clustering_R(distance)
        mp1,mp2 = compute_pairmate(cut,distance)
        group1,group2 = [],[]
        print(mp1,mp2)
        if mp1 >= pairmate_threshold and mp2 >= pairmate_threshold:
            for i in range(len(cut)):
                if cut[i] == 1:
                    group1.append(i)
                elif cut[i] == 2:
                    group2.append(i)
            plot_clustering(pca,group1,group2,"test.png")
                # else ERROR
            print(len(group1)," ", len(group2))
            fil.append(group1)
            fil.append(group2)
        else:
            res[cluster] = kmer_indice
            cluster = cluster+1
    if sorting == 'sequence':
        save_results_sorted_by_sequence(outputfile,res,kmer_name)
    else:
        save_results_sorted_by_cluster(outputfile,res,kmer_name)
