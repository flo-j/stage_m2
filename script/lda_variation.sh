
for x in 1 2 3 4 5 6 7 8 9 10;
do
  echo $x
  for lda in 10 20 50 100 200 300 500 700 1000 1200 1400 1430 15000;
  do
    output='ntlda'$lda'V2mp090run'$x'.txt'
    python3 ../script/script_python.py --pairmate 0.90 ../data/nt011630.rev.monomers.fst_149.length_166_177.kmers5 --output $output --sorting sequence --lda $lda
  done
done
