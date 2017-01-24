# A computer program does what you tell it to do, not what you want it to do -- Greer's Law

#author : Florence Jornod <florence@jornod.com>

# $1 nom du script
#$2 data à utiliser
for mp in $(seq 0.85 0.01 0.85);
do
  echo "matepair AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
  matepair=$(echo $mp | sed "s/,/\./")
  for nb_comp in 3 10 50 100 500;
  do
    echo $nb_comp
    echo "Rscript clustering analysis"
    Rscript --no-save --no-restore --verbose $1 --kmer_file $2 --matepair $matepair --nb_comp $nb_comp
    echo "Rscript clustering analysis DONE"
    echo "Comparaison"
    datafilename=$(basename $2)
    resfilename=$datafilename"mp"$matepair
    resdirectory=results/$datafilename/comparaison/$nb_comp/
    ./script/comparaison_article_res.sh results/$datafilename/mp_$matepair/$resfilename.clustering_done.txt data/nt_names_tags_by_position_clustering.dat2 $matepair $resdirectory
    echo "Comparaison DONE"
  done
  echo "Rscript clustering analysis"
  Rscript --no-save --no-restore --verbose $1 --kmer_file $2 --matepair $matepair --nb_comp $nb_comp
  echo "Rscript clustering analysis DONE"
  echo "Comparaison"
  datafilename=$(basename $2)
  resfilename=$datafilename"mp"$matepair
  resdirectory=results/$datafilename/comparaison/
  ./script/comparaison_article_res.sh results/$datafilename/mp_$matepair/$resfilename.clustering_done.txt data/nt_names_tags_by_position_clustering.dat2 $matepair $resdirectory
  echo "Comparaison DONE"
done
