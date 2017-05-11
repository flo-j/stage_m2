# Florence Jornod <florence@jornod.com>

#################################################################################

#This work is licensed under the Creative Commons Attribution-ShareAlike 3.0
#Unported License. To view a copy of this license,
#visit http://creativecommons.org/licenses/by-sa/3.0/
# or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#################################################################################

# A computer program does what you tell it to do, not what you want it to do
#-- Greer's Law

# example : python3 script_python.py ../data/nt011630.rev.monomers.fst_149.length_166_177.kmers5 -o res.txt -- --pairmate 0.90

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import time
import random
from sklearn.decomposition import PCA, IncrementalPCA
import scipy
import argparse
import fastcluster
import scipy.cluster.hierarchy
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

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
    if verbose > 2:
        logging.info("compute pca")
    pca.fit(kmer_table)
    expl_var_ratio = sum(pca.explained_variance_ratio_[0:nbcomp_pca])
    return pca.fit_transform(kmer_table), expl_var_ratio


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
    dist_time = t2 - t1
    return d, dist_time

def compute_clustering_fast(distance):
    t1=time.clock()
    c=fastcluster.ward(distance)
    t2=time.clock()
    return scipy.cluster.hierarchy.fcluster(c, 2, criterion="maxclust")

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

def create_2_groups(cut,kmer_indice):
    ''' Create two groups following clustering results.
        ARG    : - cut : a rpy2 IntVector with the repartition of the 2 groups
        RETURN : 2 lists one for each group
    '''
    plot_output="plot_"+str(i)
    if verbose > 1:
        logging.info("matepair test ok")
    g1,g2 = [],[]
    for x in range(len(cut)):
        (g1,g2)[cut[x] == 1].append(kmer_indice[x])
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
    print(max(cut),min(cut))
    if(max(cut)==min(cut)):
        print(cut)
        return(0,0)
    for i in range(len(distance)):
        j = np.argmin(distance[i,])
        res.append(str(cut[i])+str(cut[j]))
        if i == 1:
            print(min(distance[i,]))
    #    print(j)
    mp1 = float(res.count('11'))/(res.count('11')+res.count('12'))
    mp2 = float(res.count('22'))/(res.count('21')+res.count('22'))
    return mp1,mp2,res

def compute_pairmate2(cut, distance):
    ''' compute the pairmate for the two groups.Pairmate corresponds to the
        proportion of sequences which has their nearest neighbor in a same group

        ARG    : - cut      : a rpy2IntVector with the repartition of the groups
                 - distance : An 2D array with distances between each sequence
        RETURN : - a tupple with the matepair of each group
    '''
    if verbose > 2:
        logging.info("compute pairmate")
    res = []
    n = len(cut)
    if(max(cut)==min(cut)):
        return 0,0
    distance[distance<=0.00000000001] = 999
    line = distance[0:n-1]
    line = line.tolist()
    line.insert(0,999)
    j = np.argmin(np.asarray(line)) ###
    res.append(str(cut[0])+str(cut[j]))
    line = []
    for i in range(1,n):
        line = []
        pos = i - 1
        line.append(distance[pos])
        k = 1
        while k < i:
            pos = pos + n - k -1
            line.append(distance[pos])
            k=k+1
        line.append(999)
        pos = pos +n-k
        if i != n-1 :
            line2 = distance[pos:pos+n-k-1]
            line = line + line2.tolist()
        j = np.argmin(np.asarray(line)) # douze
        res.append(str(cut[i])+str(cut[j]))
    mp1 = float(res.count('11'))/(res.count('11')+res.count('12'))
    mp2 = float(res.count('22'))/(res.count('21')+res.count('22'))
    return mp1,mp2

def save_results_sorted_by_cluster(filename, res, seq_name):
    ''' save the result in a file. Each sequence is associated to its group.
        the result is order by cluster. This function is faster than sorting by
        sequence.
        ARG    : - filename : the result filename
                 - res      : a dict with cluster as keys and the seq as value
                 - seq_name : a 1D numpy.array with the name of the sequences
        NO RETURN
    '''
    if verbose > 2:
        logging.info("save file")
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
    if verbose > 2:
        logging.info("save file")
    filename = open(filename,"w")
    for i in range(len(kmer_name)):
        for key in res.keys():
            if i in res[key]:
                filename.write(kmer_name[i]+"   "+str(key)+"\n")

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
                         used to compute distance. Default = 1024", type=int,
                         default = 1024)
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
    parser.add_argument("--lda", help ="size of lda ", default="75000")
    parser.add_argument("--maxsizefamily", help ="size max of a family ",
                        type= float, default=0)
    args = parser.parse_args()
    return args.pairmate, args.input, args.nbcomppca, args.output, \
           args.sorting, args.pca, args.verbose, args.lda, args.maxsizefamily


# get args
t1=time.clock()
pairmate, kmer_file, nbcomp_pca, outputfile, \
sort, plot, verbose, lda, maxfam = get_args()
res = {}
queue = []

# initialize some variables
cluster_number = 0
i = 0
#size = len(kmer_indice)
sequence_name = get_sequences_name(kmer_file)
kmer_table=get_kmer_table(kmer_file)
queue.append(range(len(kmer_table)))
size = int(lda)
familysize = int(maxfam*len(kmer_table))
#print("size : ",size)
filename = open(outputfile,"w")
filename.write("input filename :" +kmer_file+" pairmate : " + str(pairmate) + \
 " nb comp pca : " + str(nbcomp_pca) + " size min family : " + str(familysize)\
 + " lda :"+ lda + "\n")

while(queue):
    kmer_indice = queue.pop(0) # pop the first value of the queue
    pca, expl_var_ratio = compute_pca(kmer_table[kmer_indice])
    nb_comp_dispo = len(pca[0])
    if verbose > 0:
        logging.info(str(len(pca))+" "+str(i))
    # to reduce computation time, we use just nbcomp_pca components if possible
    if(len(kmer_indice) > size):
        sample=range(len(kmer_indice))
        kmer_sample = random.sample(sample, size)
    else:
        kmer_sample=range(len(kmer_indice))
    if len(pca) > nbcomp_pca:
        pca=pca[:,range(nbcomp_pca)]
    distance, dist_time = compute_dist2(pca[kmer_sample,:])
    print(expl_var_ratio,'\t',dist_time,'\t',nb_comp_dispo,'\t',len(pca))
    clusters = compute_clustering_fast(distance)
    if(len(kmer_indice) > size):
        diff = list( set(kmer_indice) - set(kmer_sample) )
        clf = LinearDiscriminantAnalysis()
        clf.fit(pca[kmer_sample,:], clusters)
        clusters2 = clf.predict(pca[range(len(kmer_indice)),:])
    else:
        clusters2=clusters
    mp1, mp2 = compute_pairmate2(clusters, distance)
    group1, group2 = [],[]
    print(mp1,mp2)
    if verbose > 0:
        logging.info("mp "+str(mp1)+" "+str(mp2))
    if mp1 >= pairmate and mp2 >= pairmate :
        if max(clusters2)==min(clusters2):
            exit("LDA problem")
        group1, group2 =create_2_groups(clusters2, kmer_indice)
        if(len(group1)>familysize):
            queue.append(group1)
        else:
            cluster_number = cluster_number+1
            res[cluster_number] = group1
        if(len(group2)>familysize):
            queue.append(group2)
        else:
            cluster_number = cluster_number+1
            res[cluster_number] = group2
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
t2=time.clock()
print(t2-t1)
