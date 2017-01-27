# A computer program does what you tell it to do, not what you want it to do -- Greer's Law
#author : Florence Jornod <florence@jornod.com>


echo "3 10 50 100 500 total" > res_comparaison_nbcomp.txt
for mpair in $(seq 0.85 0.01 0.99);
do
  mp=$(echo $mpair | sed "s/,/\./")
  res=""
  for i in 3 10 50 100 500 -1;
  do
    #echo $mp >> AAA
    if (($i >= 0))
    then
      file=$i/mp_$mp"_res_comparaison.txt"
    else
      file=mp_$mp"_res_comparaison.txt"
    fi
    #sed -e "s/ *//" $file
   temp=$(sed -e "s/ *//" $file | cut -d' ' -f2 |cut -f1| sed -e "s/^[^0-9]//" | sort|uniq |wc -l )
   res+=$temp
    res+="  "
  done
   echo $mp $res >> res_comparaison_nbcomp.txt
done
#
