while read i; do 
    echo $i; 
    esearch -db gene -query "$i [TITLE] AND Actinopterygii [ORGN]" < /dev/null | esummary | grep -w "Name"; 
done < ./description_list.txt > gene_symbols.txt

