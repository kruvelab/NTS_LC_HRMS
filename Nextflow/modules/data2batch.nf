#!/bin/bash nextflow 

process DATA2BATCH {
    container = "quay.io/ida_rahu/ie4generation:v3.2.1"
    
    input:
    tuple val(batch_size), path(data)

    output:
    path("${data.baseName}_*.tsv")

    shell:
    '''
    basename=$(basename "!{data}" .tsv)
    
    awk -F'\\t' '
        BEGIN {OFS="\\t"} 
        NR==1 {
            for (i=1; i<=NF; i++) {
                if (\$i == "SMILES") colnum=i
            }
            if (!colnum) exit
        }
        NR>1 {print \$colnum}
    ' "!{data}" | split -l !{batch_size} - "output_chunk_"

    counter=1
    for file in output_chunk_*
    do
        echo -e "SMILES" > "${basename}_\${counter}.tsv"
        cat "\$file" >> "${basename}_\${counter}.tsv"
        rm "\$file"
        ((counter++))
    done
    '''
}