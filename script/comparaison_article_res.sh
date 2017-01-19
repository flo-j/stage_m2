# ce script permet de voir pour chaque groupe obtenu par le "programme" a quel groupe il correspond dans l'article d'alexandrov.
# appel : ./comparaison_article_res.sh file1 file2 mp
# où file1 est notre fichier de res et file2 le fichier correspondant à l'article et mp le matepair utilisé


# on supprime les "" dans le fichier qui sort de R

sed "s/\"//g" $1 > $1.without_quote

# dans un premier temps on trie les deux fichiers par ordre alphabetique -> on a donc les sequences dans le meme ordre

sort $1.without_quote > $1.without_quote_sort
sort $2 > $2.sort

# on a maintenant les deux fichiers dans le même ordre

# on recupere la 2e colonnes pour chacun des deux fichiers
cut -d' ' -f2 $1.without_quote_sort > $1.col2
cut -d' ' -f2 $2.sort | cut -d'_' -f2 > $2.col2
# on combine les deux deuxiemes colonnes dans un seul fichier
paste $1.col2 $2.col2 > fusion.txt
# on garde que les lignes uniques pour plus de lisibilité
sort fusion.txt |uniq -c > mp_$3_res_comparaison.txt

#suppression des fichiers temporaires
rm $1.without_quote
rm $1.without_quote_sort
rm $2.sort
rm $1.col2
rm $2.col2
