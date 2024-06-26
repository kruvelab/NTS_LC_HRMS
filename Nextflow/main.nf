#!/bin/bash nextflow
nextflow.enable.dsl=2

// Importing modules
include { DATA2BATCH } from './modules/data2batch'
include { SMILES2SIRIUS_FP } from './modules/smiles2sirius_fp'

workflow {
    data_ch = Channel.from(params.batch_size)
                     .combine(Channel.fromPath(params.data))
   
    input_ch = DATA2BATCH(data_ch).flatten()
                                  .combine(Channel.from(params.mode))
                                  .combine(Channel.fromPath(params.openbabel))
                                  .combine(Channel.fromPath(params.ecfp_fp_hashes))
                                  .combine(Channel.fromPath(params.biosmarts_aka_custom_made_fps))
                                  .combine(Channel.fromPath(params.ringsystem_fps))
                                  .combine(Channel.fromPath(params.csi_fingerid))
                                  .combine(Channel.fromPath(params.csi_fingerid_neg))

    SMILES2SIRIUS_FP(input_ch)
}