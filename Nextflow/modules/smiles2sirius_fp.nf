#!/bin/bash nextflow 

process SMILES2SIRIUS_FP {
    container = "quay.io/ida_rahu/ie4generation:v3.2.1"
    
    publishDir "${params.results}/", mode: 'copy', pattern: "*.tsv"
    
    input:
    tuple path(SMILES), val(mode), path(openbabel), path(ecfp_fp_hashes), path(biosmarts_aka_custom_made_fps), path(ringsystem_fps), path(csi_fingerid), path(csi_fingerid_neg)

    output:
    path('*.tsv')
    
    shell:
    '''
    Rscript !{baseDir}/bin/SMILES2SIRIUS_fp.R \
     --SMILES !{SMILES} \
     --mode !{mode} \
     --openbabel !{openbabel} \
     --ecfp_fp_hashes !{ecfp_fp_hashes} \
     --biosmarts_aka_custom_made_fps !{biosmarts_aka_custom_made_fps} \
     --ringsystem_fps !{ringsystem_fps} \
     --csi_fingerid !{csi_fingerid} \
     --csi_fingerid_neg !{csi_fingerid_neg}
    '''
}