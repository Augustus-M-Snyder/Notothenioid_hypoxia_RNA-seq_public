#!/usr/bin/bash

while read i; do 
    echo "Querying: $i"
    esearch -db gene -query "$i [TITLE] AND Actinopterygii [ORGN]" < /dev/null | \
    esummary | \
    grep -w "Name" | \
    awk -F'[<>]' '
    {
        if ($3 ~ /^LOC[0-9]+$/) {
            loc = $3;  # Capture LOC###### IDs as a fallback
        } else if ($3 != "") {
            symbol = $3;  # Capture real gene symbols
            if (!(symbol in seen)) {
                seen[symbol] = 1
                print symbol;  # Print only unique symbols
            }
        }
    }
    END {
        if (!symbol && loc && !(loc in seen)) print loc;  # If no real symbol, output unique LOC
    }'
done < description_list.txt > gene_symbols_unique.txt
