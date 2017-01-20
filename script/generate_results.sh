# A computer program does what you tell it to do, not what you want it to do -- Greer's Law

#author : Florence Jornod <florence@jornod.com>

# $1 nom du script
#$2 data à utiliser
for mp in $(seq 0.90 0.01 0.99);
do
  matepair=$(echo $mp | sed "s/,/\./")
  Rscript --no-save --no-restore --verbose $1 --kmer_file $2 --matepair $matepair
done
