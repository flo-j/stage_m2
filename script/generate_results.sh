# A computer program does what you tell it to do, not what you want it to do -- Greer's Law

#author : Florence Jornod <florence@jornod.com>

# $1 nom du script
#$2 data à utiliser
for mp in $(seq 0.80 0.01 0.90);
do
  echo "matepair"
  matepair=$(echo $mp | sed "s/,/\./")
  echo "Rscript clustering analysis"
  Rscript --no-save --no-restore --verbose $1 --kmer_file $2 --matepair $matepair
  echo "Rscript clustering analysis DONE"
  echo "Comparaison"

  datafilename=$(basename $2)
  resfilename=$datafilename"mp"$matepair
  resdirectory=results/$datafilename/comparaison/
  ./script/comparaison_article_res.sh results/$datafilename/mp_$matepair/$resfilename.clustering_done.txt data/nt_names_tags_by_position_clustering.dat2 $matepair $resdirectory
  echo "Comparaison DONE"
done
