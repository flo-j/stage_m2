# A computer program does what you tell it to do, not what you want it to do -- Greer's Law
#author : Florence Jornod <florence@jornod.com>

file="results/nt011630.rev.monomers.fst_149.length_166_177.kmers5/mp_0.95nbcomp_100percent_100/nt011630.rev.monomers.fst_149.length_166_177.kmers5mp0.95nb_comp_100percent_100.clustering_done.txt"
article="data/nt_names_tags_by_position_clustering.dat2"
sort $article > $article.sort
cut -d' ' -f2 $article.sort | cut -d'_' -f2 > $article.col2
sort $file > filename2
paste filename2 $article.col2 > filename3
filename=filename3
nb_tot=$(grep -c "_A_0_B" $filename)
sed "s/\"//g" $filename |grep "_A_0_B" > prov
nb_A=$(grep -c "_A_1_A_0_B" $filename)
nb_B=$(grep -c "_B_1_A_0_B" $filename)
#cat prov
cut -d" " -f2 prov| sed "s/1_4_A_3_A_2_B_1_A_0_B/B1/" | sed "s/1_4_A_3_B_2_A_1_A_0_B/B2/" | sed "s/1_4_B_3_A_2_A_1_A_0_B/B3/" | sed "s/1_4_B_3_A_2_B_1_A_0_B/B4/" |sed "s/1_4_B_3_B_2_B_1_A_0_B/B5/" |sed "s/1_5_A_4_A_3_B_2_B_1_A_0_B/B6/"|sed "s/1_5_A_4_B_3_B_2_A_1_A_0_B/B7/"|sed "s/1_5_B_4_A_3_A_2_A_1_A_0_B/B8/"|sed "s/1_5_B_4_A_3_B_2_B_1_A_0_B/B9/"|sed "s/1_5_B_4_B_3_B_2_A_1_A_0_B/B10/"|sed "s/1_7_A_6_A_5_A_4_A_3_A_2_A_1_A_0_B/B11/"|sed "s/1_7_A_6_B_5_A_4_A_3_A_2_A_1_A_0_B/B12/"|sed "s/1_7_B_6_B_5_A_4_A_3_A_2_A_1_A_0_B/B13/"|sed "s/1_8_A_7_B_6_A_5_A_4_A_3_A_2_A_1_A_0_B/B14/"|sed "s/1_8_B_7_B_6_A_5_A_4_A_3_A_2_A_1_A_0_B/B15/" > prov1

cut -d" " -f1 prov > prov2
cut  -f2 prov >prov4
#echo ''> prov4
paste prov2 prov1 > prov3


awk -F_ 'Begin{print"begin"}{if($5==1) {temp=$3+$6} else{temp=$4-$6} print temp}' prov3 > prov5
echo $nb_A
echo $nb_B
paste  prov5  prov3 > prov6
paste prov6 prov4| sort > res4
cat res4
