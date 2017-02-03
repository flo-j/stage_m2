# A computer program does what you tell it to do, not what you want it to do -- Greer's Law

#author : Florence Jornod <florence@jornod.com>

# $1 nom du script
#$2 data à utiliser
percent=100
nb_comp=-1
for mp in $(seq 0.85 1.01 0.99);
do
  nb_comp=-1
  matepair=$(echo $mp | sed "s/,/\./")
  echo "Rscript clustering analysis"
  Rscript --no-save --no-restore --verbose $1 --kmer_file $2 --matepair $matepair
  echo "Rscript clustering analysis DONE"
  echo "Comparaison"
  datafilename=$(basename $2)
  echo "-----------------------------------------------------"
  resfilename=$datafilename"mp"$matepair'nb_comp_'$nb_comp'percent_'$percent
  echo $resfilename
  resdirectory=results/$datafilename/comparaison/
  ./../stage/programme/script/comparaison_article_res.sh 'results/'$datafilename'/mp_'$matepair'nbcomp_'$nb_comp'percent_'$percent'/'$resfilename.clustering_done.txt ../stage/programme/data/nt_names_tags_by_position_clustering.dat2 $matepair $resdirectory
  echo "Comparaison DONE"
  for nb_comp in 3 10 50 100 500;
  do
    echo $nb_comp
    echo "Rscript clustering analysis"
    Rscript --no-save --no-restore --verbose $1 --kmer_file $2 --matepair $matepair --nb_comp $nb_comp
    echo "Rscript clustering analysis DONE"
    echo "Comparaison"
    datafilename=$(basename $2)
    resfilename=$datafilename"mp"$matepair'nb_comp_'$nb_comp'percent_'$percent
    resdirectory=results/$datafilename/comparaison/$nb_comp/
    ./../stage/programme/script/comparaison_article_res.sh 'results/'$datafilename'/mp_'$matepair'nbcomp_'$nb_comp'percent_'$percent'/'$resfilename'.clustering_done.txt' ../stage/programme/data/nt_names_tags_by_position_clustering.dat2 $matepair $resdirectory
    echo "Comparaison DONE"
  done

done
