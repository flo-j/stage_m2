

# ./gen_res_and_comparaison_nbcomp.sh recursive_clustering_analyse.R kmerfile
percent=100
nb_comp=-1
datafilename=$(basename $2)
name=$(echo $datafilename |sed "s/\..*//")
echo "total 3 10 50 100 500 " > res_comparaison_nbcomp.txt
res=""
for mp in $(seq 0.85 0.01 0.99);
do
  nb_comp=-1
  res=""
  matepair=$(echo $mp | sed "s/,/\./")
  echo "Rscript clustering analysis"
  #Rscript --no-save --no-restore --verbose $1 --kmer_file $2 --matepair $matepair --plotpca no
  echo "Rscript clustering analysis DONE"
  resfilename=$name"mp"$matepair'nb_comp_'$nb_comp'percent_'$percent
 temp=$(sed "s/\"//g" 'results/'$name'/mp_'$matepair'nbcomp_'$nb_comp'percent_'$percent'/'$resfilename.clustering_done.txt |sort | cut -d' ' -f2 |sort |uniq | wc -l)
  res+=$temp
  res+="  "
  for nb_comp in 3 10 50 100 500;
  do
    echo $nb_comp
    echo "Rscript clustering analysis"
   # Rscript --no-save --no-restore --verbose $1 --kmer_file $2 --matepair $matepair --nb_comp $nb_comp --plotpca no
    echo "Rscript clustering analysis DONE"
    resfilename=$name"mp"$matepair'nb_comp_'$nb_comp'percent_'$percent
    temp=$(sed "s/\"//g" 'results/'$name'/mp_'$matepair'nbcomp_'$nb_comp'percent_'$percent'/'$resfilename.clustering_done.txt |sort | cut -d' ' -f2 |sort | uniq | wc -l)
    res+=$temp
    res+="  "
  done
  echo $matepair $res >> res_comparaison_nbcomp.txt
done
