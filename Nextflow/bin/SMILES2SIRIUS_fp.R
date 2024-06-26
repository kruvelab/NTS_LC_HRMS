library(rcdk)
library(rJava)
library(dplyr)
library(fingerprint)
library(data.table)
library(tidyr)
library(optparse)


option_list <- list(
    make_option('--SMILES', type='character',
    help='List of SMILES'), 
    make_option('--mode', type='character',
    help='Positive mode, negative mode, or overlapping fingerprint features'),
    make_option('--openbabel', type='character',
    help='File containing SMARTS corresponding to OpenBabel fingerprint features used by SIRIUS+CSI:FingerID'),
    make_option('--ecfp_fp_hashes', type='character',
    help='File containing hash codes corresponding to ECFP fingerprint features used by SIRIUS+CSI:FingerID'),
    make_option('--biosmarts_aka_custom_made_fps', type='character',
    help='File containing SMARTS corresponding to custom made fingerprint features used by SIRIUS+CSI:FingerID'),
    make_option('--ringsystem_fps', type='character',
    help='File containing SMARTS corresponding to ringsystem fingerprint features used by SIRIUS+CSI:FingerID'),
    make_option('--csi_fingerid', type='character',
    help='File containing information about positive mode SIRIUS fingerprint'),
    make_option('--csi_fingerid_neg', type='character',
    help='File containing information about negative mode SIRIUS fingerprint'));

option_parser <- OptionParser(option_list=option_list)
options <- parse_args(option_parser)
index <- options$index
SMILES_list <- as.data.frame(fread(options$SMILES, header=T)$SMILES)

mode <- options$mode

write_raw_fp=T
options(java.parameters = "- Xmx1024m")
# Function for reading in the structural patterns (from relevant files) that
# correspond to the SIRIUS+CSI:FingerID fingerprint features 
pattern_file_reader <- function(file_name, split_pattern) {
  con <- file(file_name)
  lines <- readLines(con)
  close(con)
  slines <- strsplit(lines, split_pattern)
  colCount <- max(unlist(lapply(slines, length)))
    
  patterns <- data.frame(matrix(nrow=0, ncol=colCount))
    
  for (i in 1:length(slines)) {
    line <- slines[[i]]
    for (j in 1:length(line)) {
      patterns[i, j] <- line[j]
    }
  }
  return(patterns)
}
  
# Fingerprints that cover all SIRIUS fingerprints 
# (files downloaded from: https://github.com/boecker-lab/sirius)
OpenBabelFP3_names <- paste0('AbsIdx_', c(0:54))
OpenBabelFP3_SMARTS <- unlist(pattern_file_reader(options$openbabel, '\t')[1], use.names=F)
  
CDKsubstructure_names <- paste0('AbsIdx_', c(55:361))
MACCS_names <- paste0('AbsIdx_', c(362:527))
PubChem_names <- paste0('AbsIdx_', c(528:1408))
KlekotaRoth_names <- paste0('AbsIdx_', c(1409:6268))
  
ECFP6_names <- paste0('AbsIdx_', c(6269:8178))
ECFP6_hashes <- unlist(read.table(options$ecfp_fp_hashes, header=F), use.names=F)
  
custommadeSMARTS_names <- paste0('AbsIdx_', c(8179:8461))
custommade_SMARTS <- unlist(pattern_file_reader(options$biosmarts_aka_custom_made_fps, '\n'), use.names=F)
  
ringsystems_names <- paste0('AbsIdx_', c(8462:8924))
ringsystems_SMARTS <- unlist(pattern_file_reader(options$ringsystem_fps, '\n'), use.names=F)
  
  
no_columns = 1 + length(OpenBabelFP3_names) + length(CDKsubstructure_names) +
             length(MACCS_names) + length(PubChem_names) + 
             length(KlekotaRoth_names) + length(ECFP6_names) + 
             length(custommadeSMARTS_names) + length(ringsystems_names)
  
  
final_fp_data <- data.frame(matrix(nrow=nrow(SMILES_list), ncol=no_columns))
  
colnames(final_fp_data) <- c('SMILES', OpenBabelFP3_names, 
                              CDKsubstructure_names, MACCS_names, 
                              PubChem_names, KlekotaRoth_names, ECFP6_names, 
                              custommadeSMARTS_names, ringsystems_names) 
final_fp_data[, 1] <- SMILES_list
  
fp_index <- data.frame(matrix(nrow=8, ncol=3))
colnames(fp_index) <- c('fingerprint', 'start_index', 'end_index')
fp_index$fingerprint <- c('FP3', 'substructure', 'maccs', 'pubchem', 'kr', 'ecfp6', 'custom', 'ring')
fp_index$start_index <- c(2, length(OpenBabelFP3_names)+2, 
                          length(OpenBabelFP3_names)+length(CDKsubstructure_names)+2,
                          length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+2,
                          length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+2,
                          length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+length(KlekotaRoth_names)+2,
                          length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+length(KlekotaRoth_names)+length(ECFP6_names)+2,
                          length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+length(KlekotaRoth_names)+length(ECFP6_names)+length(custommadeSMARTS_names)+2)
  
fp_index$end_index <- c(1+length(OpenBabelFP3_names), 
                        1+length(OpenBabelFP3_names)+length(CDKsubstructure_names),
                        1+length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names),
                        1+length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names),
                        1+length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+length(KlekotaRoth_names),
                        1+length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+length(KlekotaRoth_names)+length(ECFP6_names),
                        1+length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+length(KlekotaRoth_names)+length(ECFP6_names)+length(custommadeSMARTS_names),
                        1+length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+length(KlekotaRoth_names)+length(ECFP6_names)+length(custommadeSMARTS_names)+length(ringsystems_names))
  
i = 1
for (SMILES in SMILES_list[c(1:nrow(SMILES_list)),1]) {
  skip_to_next <- F
  tryCatch({
    mol <- parse.smiles(SMILES)[[1]]
    if(!is.null(mol)){
      openbabel_fingerprints <- get.fingerprint(mol, type='substructure',
                                                substructure.pattern=OpenBabelFP3_SMARTS)
      final_fp_data[i, c(fp_index$start_index[1]:fp_index$end_index[1])] <- strsplit(as.character(openbabel_fingerprints), "")[[1]]
        
      substr_fingerprints <- get.fingerprint(mol, type='substructure')
      final_fp_data[i, c(fp_index$start_index[2]:fp_index$end_index[2])] <- strsplit(as.character(substr_fingerprints), "")[[1]]
        
      maccs_fingerprints <- get.fingerprint(mol, type='maccs')
      final_fp_data[i, c(fp_index$start_index[3]:fp_index$end_index[3])] <- strsplit(as.character(maccs_fingerprints), "")[[1]]
        
      pubchem_fingerprints <- get.fingerprint(mol, type='pubchem')
      final_fp_data[i, c(fp_index$start_index[4]:fp_index$end_index[4])] <- strsplit(as.character(pubchem_fingerprints), "")[[1]]
        
      kr_fingerprints <- get.fingerprint(mol, type='kr')
      final_fp_data[i, c(fp_index$start_index[5]:fp_index$end_index[5])] <- strsplit(as.character(kr_fingerprints), "")[[1]]
        
      ecfp_fingerprints <- get.fingerprint(mol, type='circular', circular.type='ECFP6', fp.mode='count')
      for (idx in 1:length(ecfp_fingerprints@features)) {
        hash <- strsplit(as.character(ecfp_fingerprints@features[[idx]]), ':')[[1]][1]
        if (hash %in% ECFP6_hashes) {
          right_hash <- which(hash == ECFP6_hashes)
          column_name <- ECFP6_names[right_hash]
          final_fp_data[i, column_name] <- 1
        }
      }
        
      custommade_fingerprints <- get.fingerprint(mol, type="substructure",
                                                 substructure.pattern=custommade_SMARTS)
      final_fp_data[i, c(fp_index$start_index[7]:fp_index$end_index[7])] <- strsplit(as.character(custommade_fingerprints), "")[[1]]
        
      ring_fingerprints <- get.fingerprint(mol, type="substructure",
                                           substructure.pattern=ringsystems_SMARTS)
      final_fp_data[i, c(fp_index$start_index[8]:fp_index$end_index[8])] <- strsplit(as.character(ring_fingerprints), "")[[1]]
        
    }
    else {
      wrong_smile <- as.data.frame(SMILES)
      error_comp <- error_comp %>% rbind(wrong_smile)
      print(SMILES)
    }
    gc()
    i = i + 1
      
  }, error = function(e) { skip_to_next <<- T}
      
  )
  if(skip_to_next) {next}
}

final_fp_data <- final_fp_data %>% mutate_at(ECFP6_names, ~replace_na(as.numeric(.), 0))
final_fp_data <- final_fp_data %>% mutate_at(c(3:ncol(final_fp_data)), as.numeric)
  
if (write_raw_fp) {
    write.table(final_fp_data, paste0('raw_SIRIUS_fp_', tools::file_path_sans_ext(options$SMILES), '.tsv'), row.names=F, col.names=T, sep='\t', quote=F)
}

# Generating the names for features so calculated fingerprint features would 
# match with SIRIUS+CSI:FingerID absolute index.
positive_idxs <- pattern_file_reader(options$csi_fingerid, '\t')
colnames(positive_idxs) <- positive_idxs[1,]
positive_idxs <- positive_idxs[-1, ] 
  
negative_idxs <- pattern_file_reader(options$csi_fingerid_neg, '\t')
colnames(negative_idxs) <- negative_idxs[1,]
negative_idxs <- negative_idxs[-1, ] 
  
together_idx <- merge(positive_idxs, negative_idxs, by='absoluteIndex')
together_idx <- together_idx[order(as.numeric(as.character(together_idx$absoluteIndex))), ]
positive_idxs$absoluteIndex <- sub('^', 'AbsIdx_', positive_idxs$absoluteIndex)
negative_idxs$absoluteIndex <- sub('^', 'AbsIdx_', negative_idxs$absoluteIndex)
together_idx$absoluteIndex <- sub('^', 'AbsIdx_', together_idx$absoluteIndex)
  
# Generating data frames containing fingerprint features overlapping with those 
# outputted by SIRIUS+CSI:FingerID, based on the selected ionization mode ('pos', 'neg', 'overlapping' or 'raw').
if (mode == 'pos') {
  positive_mode_data <- final_fp_data %>% dplyr::select(c('SMILES', positive_idxs$absoluteIndex))
  write.table(positive_mode_data, paste0('pos_SIRIUS_fp_', tools::file_path_sans_ext(options$SMILES), '.tsv'), row.names=F, col.names=T, sep='\t', quote=F)
} else if (mode == 'neg') {
  negative_mode_data <- final_fp_data %>% dplyr::select(c('SMILES', negative_idxs$absoluteIndex))
  write.table(negative_mode_data, paste0('neg_SIRIUS_fp_', tools::file_path_sans_ext(options$SMILES), '.tsv'), row.names=F, col.names=T, sep='\t', quote=F)
} else {
  general_mode_data <- final_fp_data %>% dplyr::select(c('SMILES', together_idx$absoluteIndex))
  write.table(general_mode_data, paste0('overlapping_SIRIUS_fp_', tools::file_path_sans_ext(options$SMILES), '.tsv'), row.names=F, col.names=T, sep='\t', quote=F)
} 