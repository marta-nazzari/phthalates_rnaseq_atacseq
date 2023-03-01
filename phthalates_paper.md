This notebook contains the workflow used for the processing and analysis of RNA-Seq and ATAC-Seq data underlying the manuscript "Investigation of the epigenetic effects of phthalates on _in vitro_ thyroid models with RNA-Seq and ATAC-Seq".


# RNA-Seq data analysis
## 1 - Raw sequencing data processing
Used the script available at https://github.com/marta-nazzari/CODA/blob/main/CODA_preprocess.sh.


## 2 - Differential expression analysis
```{r}
# load required packages
library(DESeq2)
library(magrittr)
library(dplyr)
library(data.table)
library(edgeR)
library(ggrepel)
library(biomaRt)
library(PCAtools)
library(apeglm)
library(RUVSeq)

# Define %notin% function
`%notin%` = Negate(`%in%`) 

# General parameters
Control = 'DMSO'
FDR = 0.01
k_ruv = 2
mart = read.table('D:/lab/mouse_ensgid_to_gene_name.txt', header = T, fill = T)


Drugs = c('DEHP', 'DIDP', 'DINP', 'DnOP')
  
# Working directories
MAIN.DIR = '/project_folder/RNA-Seq/DE_analysis/'

if (dir.exists(MAIN.DIR)) { 
setwd(MAIN.DIR) 
print('Main directory exists')
} else { 
dir.create(MAIN.DIR)
setwd(MAIN.DIR) 
print('Created working directory')
}

# Raw gene read count table
genes = read.table('project_folder/RNA-Seq/DE_analysis/genes.data.tsv', header = TRUE)
rownames(genes) %<>% sub('\\.(\\d){1,2}', '', .) # remove gene version

# Raw miRNA read count table
mirna = read.table('project_folder/RNA-Seq/DE_analysis/miR_counts.txt', header = T)

# read the metadata table
samplekey = read.table('project_folder/RNA-Seq/DE_analysis/metadata.txt', header = T) %>% 
tibble::column_to_rownames(var = 'SampleName') %>% 
mutate_all(as.factor)

# Need to exclude DMSO_1 because it's an outlier
samplekey %<>% filter(rownames(.) != 'DMSO_ctrl_1')

# update gene and mirna tables
genes %<>% dplyr::select(-DMSO_ctrl_1) 
mirna %<>% dplyr::select(-DMSO_ctrl_1) 

# colorblind palette for plotting mapped reads
cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")


# Where list of genes that pass or fail the filters go
filters_dir = paste0(filtered_dir, 'Filtered_genes/')

if (dir.exists(unfiltered_dir)) { print('Directory exists') } else { 
  dir.create(unfiltered_dir)
  print('Created directory') 
}

if (dir.exists(filtered_dir)) { print('Directory exists') } else { 
  dir.create(filtered_dir)
  print('Created directory') 
}

if (dir.exists(filters_dir)) { print('Directory exists') } else { 
  dir.create(filters_dir)
  print('Created directory') 
}  

# PREPARE GENES TABLE
# check which samples are in `samplekey` and NOT in `genes`
rownames(samplekey) %notin% colnames(genes) %>% rownames(samplekey)[.]
colnames(genes) %notin% rownames(samplekey) %>% colnames(genes)[.]

# show colnames in alphabetical order
colnames(genes) %>% order %>% colnames(genes)[.]

# keep only samples in both `genes` and `samplekey`
common_samples = intersect(colnames(genes), rownames(samplekey))
samplekey %<>% .[rownames(samplekey) %in% common_samples, ] 
genes %<>% .[, colnames(genes) %in% common_samples] 

# reorder samples in genes to the same order of rownames(samplekey)
genes %<>% data.table::setcolorder(., rownames(samplekey))

# Check. Needed for DESeq2
colnames(genes) == rownames(samplekey)

# PREPARE MIRNA TABLE
# check which samples are in `samplekey` and NOT in `mirna`
rownames(samplekey) %notin% colnames(mirna) %>% rownames(samplekey)[.]
colnames(mirna) %notin% rownames(samplekey) %>% colnames(mirna)[.]

# show colnames in alphabetical orer
colnames(mirna) %>% order %>% colnames(mirna)[.]

# keep only samples in both `mirna` and `samplekey`
common_samples = intersect(colnames(mirna), rownames(samplekey))
samplekey %<>% .[rownames(samplekey) %in% common_samples, ] 
mirna %<>% .[, colnames(mirna) %in% common_samples] 

# reorder samples in mirna to the same order of samplekey$SampleName
mirna %<>% data.table::setcolorder(., as.character(rownames(samplekey)))


# plot number of aligned reads per samples
samplekey %<>% mutate(read_count_genes = colSums(genes))
samplekey %<>% mutate(read_count_mirna = colSums(mirna))

ggplot(samplekey, aes(x = rownames(samplekey), y = read_count_genes, fill = Compound)) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(breaks = samplekey$SampleCode, labels = rownames(samplekey)) +
  scale_fill_manual(values = cbPalette) +
  ggtitle('Genes aligned reads') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
ggsave('genes_read_count.png', width = 8000, height = 6000, dpi = 600, units = 'px')

ggplot(samplekey, aes(x = rownames(samplekey), y = read_count_mirna, fill = Compound)) +
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = cbPalette) +
  ggtitle('miRNA aligned reads') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
ggsave('mirna_read_count.png', width = 8000, height = 6000, dpi = 600, units = 'px')


# Remove snoRNA genes. To analyze snoRNAs, change `!=` to `==`
if (remove_snorna == T) { genes %<>% .[rownames(.) %in% mart$ensembl_gene_id[mart$gene_biotype != 'snoRNA'], ] }



###
# DE analysis genes
###

# PART 1: DE ANALYSIS
DEanalysisCompound = function (COMPARISON, DRUG, CONTROL, RNA, CONTRAST.FACTOR, OUTPUT.DIR) {	
  
  # Select output directory
  if (dir.exists(OUTPUT.DIR)) { 
    setwd(OUTPUT.DIR) 
  } else { 
    dir.create(OUTPUT.DIR)
    setwd(OUTPUT.DIR) }
  
  # Define DESeq2 design 
  DE_Design = samplekey[c(grep(DRUG[COMPARISON], samplekey$Compound), grep(CONTROL, samplekey$Compound)), ]
  
  # Select samples for DESeq2 analysis
  samples = RNA[, colnames(RNA) %in% rownames(DE_Design)]
  
  ColData = samplekey[rownames(samplekey) %in% colnames(samples), ] %>% droplevels
  
  print(paste(DRUG[COMPARISON], "vs", CONTROL))		
  
  dds = DESeqDataSetFromMatrix(countData = round(samples), 
                               colData = ColData,
                               design = ~ Compound)
  
  dds %<>% DESeq()
  
  res = results(dds, alpha = FDR, contrast = c(CONTRAST.FACTOR, DRUG[COMPARISON], CONTROL), pAdjustMethod = 'fdr')
  res = res[order(res$padj), ]	
  
  DEgenes = subset(res, res$padj < FDR)	
  DEgenes = DEgenes[order(DEgenes$padj), ]
  
  # Counts (all genes)
  normCounts = counts(dds, normalized = TRUE)
  
  ## Export tables
  FileName = paste(DRUG[COMPARISON], 'vs', CONTROL, 'FDR', FDR, sep = '_')		
  
  # results (all genes)
  write.table(res, file = paste0(FileName, '_Results_unfiltered.txt'), sep = '\t', quote = FALSE)
  
  # Normalized read count	(all genes)	
  write.table(normCounts, file = paste0(FileName, '_Norm_count_unfiltered.txt'), sep = '\t', quote = FALSE)
  
  # Normalized read count (DE genes)
  write.table(normCounts[rownames(normCounts) %in% rownames(DEgenes), ], file = paste0(FileName, '_Norm_count_DEG.txt'), sep = '\t', quote = F)
  
  # DE gene names only
  write.table(rownames(DEgenes), file = paste0(FileName, '_DEG_names.txt'), sep = '\t', quote = F, row.names = F, col.names = F)
  
  # assign variables to global env
  assign('dds', dds, envir = .GlobalEnv)
}	

# FILTER 1 - LOW COUNT: PASS if at least 75% of samples in either group have CPM >=1 
CPMfilter = function(COMPARISON, DRUG, CONTROL, DDS, CONTRAST.FACTOR, OUTPUT.DIR, ...) { 
  
  # Select output directory
  if (dir.exists(OUTPUT.DIR)) { 
    setwd(OUTPUT.DIR) 
  } else { 
    dir.create(OUTPUT.DIR)
    setwd(OUTPUT.DIR) }
  
  # Convert dds counts to CPM
  CPMdds = cpm(counts(DDS, normalized = TRUE))
  
  # Count how many DRUG and CONTROL samples there are
  nDrug = grep(DRUG[COMPARISON], colnames(CPMdds)) %>% length
  nCont = grep(CONTROL, colnames(CPMdds)) %>% length
  
  # Create df with only DRUG or CONTROL samples
  CPM_drug = CPMdds %>% as.data.frame %>% dplyr::select(grep(DRUG[COMPARISON], colnames(CPMdds)))
  CPM_drug[, 'PassFail'] = NA
  
  CPM_cont = CPMdds %>% as.data.frame %>% dplyr::select(grep(CONTROL, colnames(CPMdds)))
  CPM_cont[, 'PassFail'] = NA
  
  for (x in 1:nrow(CPM_drug)) {
    # Vector for gene x with TRUE if CPM >= 1, FALSE if CPM < 1
    countCPMpass_drug = CPM_drug %>% .[x,-ncol(CPM_drug)] >= 1  # Exclude PassFail column
    countCPMpass_cont = CPM_cont %>% .[x,-ncol(CPM_cont)] >= 1  # Exclude PassFail column
    
    # Sum samples for which CPM >= 1
    sumCPMpass_drug = countCPMpass_drug[countCPMpass_drug == TRUE] %>% length
    sumCPMpass_cont = countCPMpass_cont[countCPMpass_cont == TRUE] %>% length
    # If this value is HIGHER than 75% of tot samples, the gene passes the filter
    if (sumCPMpass_drug >= 0.75*nDrug) { CPM_drug$PassFail[x] = 'PASS' }
    else { CPM_drug$PassFail[x] = 'FAIL' }
    if (sumCPMpass_cont >= 0.75*nCont) { CPM_cont$PassFail[x] = 'PASS' }
    else { CPM_cont$PassFail[x] = 'FAIL' }
  }
  
  # Select genes that pass the filter in DRUG or CONTROL
  CPMpass = data.frame(genes = rownames(DDS),
                       drug_filter = CPM_drug$PassFail,
                       ctrl_filter = CPM_cont$PassFail) %>%
    mutate(PassFail = if_else(drug_filter == 'PASS' | ctrl_filter == 'PASS', 'PASS', 'FAIL'))
  
  CPMpass = CPMpass$genes[CPMpass$PassFail == 'PASS'] %>% #droplevels %>% as.vector
    as.vector()
  
  # Extract dds results that pass FILTER 1
  res = results(DDS[rownames(DDS) %in% CPMpass, ], 
                alpha = 0.01, cooksCutoff = F, independentFiltering = F, 
                contrast = c(CONTRAST.FACTOR, DRUG[COMPARISON], CONTROL), 
                pAdjustMethod = 'fdr')
  
  # Apply RUVSeq
  set = newSeqExpressionSet(counts(DDS))
  png(paste0(DRUG[COMPARISON], '_before_RUVg.png'), width = 7000, height = 4000, units = 'px', res = 600)
  par(mar = c(8,5,2,2), mfrow = c(1, 2))
  plotRLE(set, outline = FALSE, col = ifelse(grepl(DRUG[COMPARISON], colnames(DDS)), 'red', 'blue'), las = 2,
          main = 'DESeq2 normalized data, no RUVSeq applied')
  plotPCA(set, col = ifelse(grepl(DRUG[COMPARISON], colnames(DDS)), 'red', 'blue'), 
          las = 2, pch = 20, labels = F,
          main = 'DESeq2 normalized data, no RUVSeq applied')
  dev.off()
  
  set = betweenLaneNormalization(set[rownames(set) %in% CPMpass, ], which = "upper")
  not.sig = rownames(res)[which(res$pvalue > .1)]
  empirical = rownames(set)[ rownames(set) %in% not.sig ]
  set = RUVg(set, empirical, k = k_ruv)
  pData(set)
  
  png(paste0(DRUG[COMPARISON], '_after_RUVg.png'), width = 7000, height = 4000, units = 'px', res = 600)
  par(mar = c(8,5,2,2), mfrow = c(1, 2))
  plotRLE(set, outline = FALSE, col = ifelse(grepl(DRUG[COMPARISON], colnames(dds)), 'red', 'blue'), las = 2,
          main = paste0('DESeq2 normalized data, RUVg (k = ', k_ruv, ')'))
  plotPCA(set, col = ifelse(grepl(DRUG[COMPARISON], colnames(dds)), 'red', 'blue'), 
          las = 2, pch = 20, labels = F,
          main = paste0('DESeq2 normalized data, RUVg (k = ', k_ruv, ')'))
  dev.off()
  
  # plot the factors estimated by RUV
  # par(mfrow = c(2, 1), mar = c(3,5,3,1))
  # for (i in 1:2) {
  #   stripchart(as.formula(paste0('pData(set)[, i] ~ DDS$', CONTRAST.FACTOR)), vertical = TRUE, main = paste0("W", i))
  #   abline(h = 0)
  # }
  
  # control for these factors by adding them to the DESeqDataSet and to the design
  DDSRUV = DDS
  DDSRUV$W1 = set$W_1
  DDSRUV$W2 = set$W_2
  design(DDSRUV) = as.formula(paste0('~ W1 + W2 + ', CONTRAST.FACTOR))
  
  DDSRUV %<>% DESeq()
  
  # re-estimate he DESeq2 parameters and results
  res = results(DDSRUV[rownames(DDSRUV) %in% CPMpass, ], 
                alpha = 0.01, cooksCutoff = F, independentFiltering = F, 
                contrast = c(CONTRAST.FACTOR, DRUG[COMPARISON], CONTROL), 
                pAdjustMethod = 'fdr')  
  
  # Normalized read count for genes that pass FILTER 1
  normData_all = counts(DDSRUV[rownames(DDSRUV) %in% CPMpass, ], normalized = T)
  
  # DE genes (choose FDR)
  DEgenes = subset(res, res$padj < FDR)
  
  # Normalized read count of DE genes
  normData_DEG = normData_all[rownames(normData_all) %in% rownames(DEgenes), , drop = F]
  
  # Print info on filtering
  print(paste0(nrow(CPM_drug), ' genes in total.'))
  print(paste0(length(CPMpass), ' genes passed the CPM filter.'))
  print(paste0(nrow(CPM_drug[rownames(CPM_drug) %notin% rownames(res), ]), ' genes failed the CPM filter.'))
  print(paste0(nrow(DEgenes), ' genes that passed the CPM filter are differentilly expressed (fdr <', FDR, ').'))
  print(paste0(nrow(res[rownames(res) %notin% rownames(DEgenes), ]), ' genes that passed the CPM filter are not differentially expressed.'))
  
  # Save list of genes that pass FILTER 1
  write.table(CPMpass, file = paste0('Filtered_genes/', DRUG[COMPARISON], ..., '_vs_', CONTROL, '_CPMpass.txt'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Save results of DE analysis (for the volcano plot)
  FileName = paste(DRUG[COMPARISON], ...,  'vs', CONTROL, 'FDR', FDR, sep = '_')
  write.table(res, file = paste0(FileName, '_Results_all.txt'), sep = '\t', quote = F)
  write.table(normData_all, file = paste0(FileName, '_Norm_count_all.txt'), sep = '\t', quote = F)
  write.table(normData_DEG, file = paste0(FileName, '_Norm_count_DEG_prefilters.txt'), sep = '\t', quote = F)
  
  # assign variables to global env 
  assign('res', DEgenes, envir = .GlobalEnv)
  assign('normData_DEG', normData_DEG, envir = .GlobalEnv)
  assign('CPMpass', CPMpass, envir = .GlobalEnv)	
}


# FILTER 2 - SPIKE: DE gene PASSES if (Max/Sum) < 1.4*Samples of a group^(-0.66)
spikeFilter = function(COMPARISON, DRUG, CONTROL, NORM.DEG, OUTPUT.DIR, ...) {
  
  if (dir.exists(OUTPUT.DIR)) { 
    setwd(OUTPUT.DIR) 
  } else { 
    dir.create(OUTPUT.DIR)
    setwd(OUTPUT.DIR) }
  
  Drug = grep(DRUG[COMPARISON], colnames(NORM.DEG))
  Cont = grep(CONTROL, colnames(NORM.DEG))
  
  # Check genes in DRUG group
  Check_drug = data.frame(Max = apply(NORM.DEG[, Drug, drop = F], 1, max),
                          Sum = apply(NORM.DEG[, Drug, drop = F], 1, sum),
                          Length = length(Drug)) %>% 
    mutate(SpikeCheck = if_else( (Max/Sum) < (1.4*Length^(-0.66)), 'PASS', 'FAIL')) 
  
  
  # Extract genes that PASS in DRUG
  SpikePass_drug = Check_drug[grep('PASS', Check_drug$SpikeCheck), ] %>% rownames
  
  # Check genes in CONTROL group
  Check_cont = data.frame(Max = apply(NORM.DEG[, Cont, drop = F], 1, max),
                          Sum = apply(NORM.DEG[, Cont, drop = F], 1, sum),
                          Length = length(Cont)) %>% 
    mutate(SpikeCheck = if_else( (Max/Sum) < (1.4*Length^(-0.66)) , 'PASS', 'FAIL')) 
  
  # If read count in all replicates is 0, SpikeCheck = <NA>
  # Change all <NA> to PASS
  Check_drug[is.na(Check_drug)] = 'PASS' 
  Check_cont[is.na(Check_cont)] = 'PASS' 
  
  # Extract genes that PASS in CONTROL
  SpikePass_cont = Check_cont[grep('PASS', Check_cont$SpikeCheck), ] %>% rownames
  
  # Intersect the common genes: list of genes that pass the spike filter
  SpikePass = intersect(SpikePass_drug, SpikePass_cont) 
  SpikeFail = NORM.DEG[rownames(NORM.DEG) %notin% SpikePass, , drop = F] %>% rownames	# drop = F is needed when length(SpikeFail)=1
  
  # print info on filtering
  print(paste0(length(SpikePass), ' genes passed the Spike filter.'))
  print(paste0(length(SpikeFail), ' genes failed the Spike filter.'))
  
  # Save list of genes that pass FILTER 2
  write.table(SpikePass, file = paste(DRUG[COMPARISON], ..., 'vs', CONTROL, 'SpikePass.txt', sep = '_'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(SpikeFail, file = paste(DRUG[COMPARISON], ..., 'vs', CONTROL, 'SpikeFail.txt', sep = '_'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # assign variable to global env
  assign('SpikePass', SpikePass, envir = .GlobalEnv)
}

# Function that keeps normalized data of DE genes that pass FILTER 2 
filterNormData = function(COMPARISON, DRUG, CONTROL, SPIKEPASS, RES.DEG, NORM.DEG, OUTPUT.DIR, ...) {
  
  # Select output directory
  if (dir.exists(OUTPUT.DIR)) { 
    setwd(OUTPUT.DIR) 
  } else { 
    dir.create(OUTPUT.DIR)
    setwd(OUTPUT.DIR) }
  
  # Extract genes that pass FILTER 2
  idx = SPIKEPASS
  
  FileName = paste(DRUG[COMPARISON], ..., 'vs', CONTROL, 'FDR', FDR, sep = '_')	
  
  # Results of ONLY DE genes that PASS all filters
  res_filtered = RES.DEG[rownames(RES.DEG) %in% idx, , drop = F]
  
  # Normalized read count of ONLY DE genes that PASS all filters
  normData_filtered = NORM.DEG[rownames(NORM.DEG) %in% idx, , drop = F]
  
  # Info on all filters applied
  print(paste0(nrow(res_filtered), ' out of ', nrow(RES.DEG), ' DE genes passed the Spike filter - ', nrow(RES.DEG[rownames(RES.DEG) %notin% idx, ]), ' DE genes excluded.'))
  
  # Save tables
  write.table(res_filtered, file = paste0(FileName, '_Results_DEG.txt'), sep = '\t', quote = FALSE)
  write.table(rownames(res_filtered), file = paste0(FileName, '_DEG_names.txt'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)	
  write.table(normData_filtered, file = paste0(FileName, '_Norm_count_DEG.txt'), sep = '\t', quote = FALSE)
  
  # assign variable to global env
  assign('res_filtered', res_filtered, envir = .GlobalEnv)
}


###
# DE ANALYSIS miRNAs
###

# PART 1: DE ANALYSIS
DEanalysisCompound_miRNA = function (COMPARISON, DRUG, CONTROL, RNA, CONTRAST.FACTOR, OUTPUT.DIR) {	
  
  # Select output directory
  if (dir.exists(OUTPUT.DIR)) { 
    setwd(OUTPUT.DIR) 
  } else { 
    dir.create(OUTPUT.DIR)
    setwd(OUTPUT.DIR) }
  
  # Define DESeq2 design 
  DE_Design = samplekey[c(grep(DRUG[COMPARISON], samplekey$Compound), grep(CONTROL, samplekey$Compound)), ]
  
  # Select samples for DESeq2 analysis
  samples = RNA[, rownames(DE_Design) ]
  ColData = samplekey[rownames(samplekey) %in% colnames(samples), ] %>% droplevels
  
  print(paste(DRUG[COMPARISON], "vs", CONTROL))		
  
  dds = DESeqDataSetFromMatrix(countData = round(samples), 
                               colData = ColData,
                               design = ~ Compound)
  
  dds %<>% DESeq()
  
  res = results(dds, alpha = 0.01, contrast = c(CONTRAST.FACTOR, DRUG[COMPARISON], CONTROL), pAdjustMethod = 'fdr')
  res = res[order(res$padj), ]	
  
  DEgenes = subset(res, res$padj < 0.01)	
  DEgenes = DEgenes[order(DEgenes$padj), ]
  
  # Counts (all genes)
  normCounts = counts(dds, normalized = TRUE)
  
  ## Export tables
  FileName = paste('miRNA', DRUG[COMPARISON], 'vs', CONTROL, 'FDR_0.01', sep = '_')		
  
  # results (all genes)
  write.table(res, file = paste0(FileName, '_Results_unfiltered.txt'), sep = '\t', quote = FALSE)
  
  # Normalized read count	(all genes)	
  write.table(normCounts, file = paste0(FileName, '_Norm_count_unfiltered.txt'), sep = '\t', quote = FALSE)
  
  # Normalized read count (DE genes)
  write.table(normCounts[rownames(normCounts) %in% rownames(DEgenes), ], file = paste0(FileName, '_Norm_count_DEG.txt'), sep = '\t', quote = F)
  
  # DE gene names only
  write.table(rownames(DEgenes), file = paste0(FileName, '_DEG_names.txt'), sep = '\t', quote = F, row.names = F, col.names = F)
  
  # assign variables to global env
  assign('dds', dds, envir = .GlobalEnv)
}	

# FILTER 1 - LOW COUNT: PASS if at least 75% of samples in either group have CPM >=1 
CPMfilter_miRNA = function(COMPARISON, DRUG, CONTROL, DDS, CONTRAST.FACTOR, OUTPUT.DIR, ...) { 
  
  # Select output directory
  if (dir.exists(OUTPUT.DIR)) { 
    setwd(OUTPUT.DIR) 
  } else { 
    dir.create(OUTPUT.DIR)
    setwd(OUTPUT.DIR) }
  
  # Convert dds counts to CPM
  CPMdds = cpm(counts(DDS, normalized = TRUE))
  
  # Count how many DRUG and CONTROL samples there are
  nDrug = grep(DRUG[COMPARISON], colnames(CPMdds)) %>% length
  nCont = grep(CONTROL, colnames(CPMdds)) %>% length
  
  # Create df with only DRUG or CONTROL samples
  CPM_drug = CPMdds %>% as.data.frame %>% dplyr::select(grep(DRUG[COMPARISON], colnames(CPMdds)))
  CPM_drug[, 'PassFail'] = NA
  
  CPM_cont = CPMdds %>% as.data.frame %>% dplyr::select(grep(CONTROL, colnames(CPMdds)))
  CPM_cont[, 'PassFail'] = NA
  
  for (x in 1:nrow(CPM_drug)) {
    # Vector for gene x with TRUE if CPM >= 1, FALSE if CPM < 1
    countCPMpass_drug = CPM_drug %>% .[x,-ncol(CPM_drug)] >= 1  # Exclude PassFail column
    countCPMpass_cont = CPM_cont %>% .[x,-ncol(CPM_cont)] >= 1  # Exclude PassFail column
    
    # Sum samples for which CPM >= 1
    sumCPMpass_drug = countCPMpass_drug[countCPMpass_drug == TRUE] %>% length
    sumCPMpass_cont = countCPMpass_cont[countCPMpass_cont == TRUE] %>% length
    # If this value is HIGHER than 75% of tot samples, the gene passes the filter
    if (sumCPMpass_drug >= 0.75*nDrug) { CPM_drug$PassFail[x] = 'PASS' }
    else { CPM_drug$PassFail[x] = 'FAIL' }
    if (sumCPMpass_cont >= 0.75*nCont) { CPM_cont$PassFail[x] = 'PASS' }
    else { CPM_cont$PassFail[x] = 'FAIL' }
  }
  
  # Select genes that pass the filter in DRUG or CONTROL
  CPMpass = data.frame(genes = rownames(DDS),
                       drug_filter = CPM_drug$PassFail,
                       ctrl_filter = CPM_cont$PassFail) %>%
    mutate(PassFail = if_else(drug_filter == 'PASS' | ctrl_filter == 'PASS', 'PASS', 'FAIL'))
  
  CPMpass = CPMpass$genes[CPMpass$PassFail == 'PASS'] %>% #droplevels %>% as.vector
    as.vector()
  
  # Extract dds results that pass FILTER 1
  res = results(DDS[rownames(DDS) %in% CPMpass, ], 
                alpha = 0.01, cooksCutoff = F, independentFiltering = F, 
                contrast = c(CONTRAST.FACTOR, DRUG[COMPARISON], CONTROL), 
                pAdjustMethod = 'fdr')
  
  # Apply RUVSeq
  set = newSeqExpressionSet(counts(DDS))
  png(paste0('miRNA_', DRUG[COMPARISON], '_before_RUVg.png'), width = 7000, height = 4000, units = 'px', res = 600)
  par(mar = c(8,5,2,2), mfrow = c(1, 2))
  plotRLE(set, outline = FALSE, col = ifelse(grepl(DRUG[COMPARISON], colnames(DDS)), 'red', 'blue'), las = 2,
          main = 'DESeq2 normalized data, no RUVSeq applied')
  plotPCA(set, col = ifelse(grepl(DRUG[COMPARISON], colnames(DDS)), 'red', 'blue'), 
          las = 2, pch = 20, labels = F,
          main = 'DESeq2 normalized data, no RUVSeq applied')
  dev.off()
  
  set = betweenLaneNormalization(set[rownames(set) %in% CPMpass, ], which = "upper")
  not.sig = rownames(res)[which(res$pvalue > .1)]
  empirical = rownames(set)[ rownames(set) %in% not.sig ]
  set = RUVg(set, empirical, k = k_ruv)
  pData(set)
  
  png(paste0('miRNA_', DRUG[COMPARISON], '_after_RUVg.png'), width = 7000, height = 4000, units = 'px', res = 600)
  par(mar = c(8,5,2,2), mfrow = c(1, 2))
  plotRLE(set, outline = FALSE, col = ifelse(grepl(DRUG[COMPARISON], colnames(dds)), 'red', 'blue'), las = 2,
          main = paste0('DESeq2 normalized data, RUVg (k = ', k_ruv, ')'))
  plotPCA(set, col = ifelse(grepl(DRUG[COMPARISON], colnames(dds)), 'red', 'blue'), 
          las = 2, pch = 20, labels = F,
          main = paste0('DESeq2 normalized data, RUVg (k = ', k_ruv, ')'))
  dev.off()
  
  # control for these factors by adding them to the DESeqDataSet and to the design
  DDSRUV = DDS
  DDSRUV$W1 = set$W_1
  DDSRUV$W2 = set$W_2
  design(DDSRUV) = as.formula(paste0('~ W1 + W2 + ', CONTRAST.FACTOR))
  
  DDSRUV %<>% DESeq()
  
  # re-estimate he DESeq2 parameters and results
  res = results(DDSRUV[rownames(DDSRUV) %in% CPMpass, ], 
                alpha = 0.01, cooksCutoff = F, independentFiltering = F, 
                contrast = c(CONTRAST.FACTOR, DRUG[COMPARISON], CONTROL), 
                pAdjustMethod = 'fdr')  
  
  # Normalized read count for genes that pass FILTER 1
  normData_all = counts(DDS[rownames(DDS) %in% CPMpass, ], normalized = T)
  
  # DE genes (fdr < 0.01)
  DEgenes = subset(res, res$padj < 0.01)
  
  # Normalized read count of DE genes
  normData_DEG = normData_all[rownames(normData_all) %in% rownames(DEgenes), , drop = F]
  
  # Print info on filtering
  print(paste0(nrow(CPM_drug), ' genes in total.'))
  print(paste0(length(CPMpass), ' genes passed the CPM filter.'))
  print(paste0(nrow(CPM_drug[rownames(CPM_drug) %notin% rownames(res), ]), ' genes failed the CPM filter.'))
  print(paste0(nrow(DEgenes), ' genes that passed the CPM filter are differentilly expressed (fdr < 0.01).'))
  print(paste0(nrow(res[rownames(res) %notin% rownames(DEgenes), ]), ' genes that passed the CPM filter are not differentially expressed.'))
  
  # Save list of genes that pass FILTER 1	
  write.table(CPMpass, file = paste('Filtered_genes/miRNA', DRUG[COMPARISON], ..., 'vs', CONTROL, 'CPMpass.txt', sep = '_'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Save results of DE analysis (for the volcano plot)
  FileName = paste('miRNA', DRUG[COMPARISON], ...,  'vs', CONTROL, 'FDR_0.01', sep = '_')
  write.table(res, file = paste0(FileName, '_Results_all.txt'), sep = '\t', quote = F)
  write.table(normData_all, file = paste0(FileName, '_Norm_count_all.txt'), sep = '\t', quote = F)
  write.table(normData_DEG, file = paste0(FileName, '_Norm_count_DEG_prefilters.txt'), sep = '\t', quote = F)
  
  # assign variables to global env 
  assign('res', DEgenes, envir = .GlobalEnv)
  assign('normData_DEG',normData_DEG, envir = .GlobalEnv)
  assign('CPMpass', CPMpass, envir = .GlobalEnv)
}


# FILTER 2 - SPIKE: DE gene PASSES if (Max/Sum) < 1.4*Samples of a group^(-0.66)
spikeFilter_miRNA = function(COMPARISON, DRUG, CONTROL, NORM.DEG, OUTPUT.DIR, ...) {
  
  if (dir.exists(OUTPUT.DIR)) { 
    setwd(OUTPUT.DIR) 
  } else { 
    dir.create(OUTPUT.DIR)
    setwd(OUTPUT.DIR) }
  
  Drug = grep(DRUG[COMPARISON], colnames(NORM.DEG))
  Cont = grep(CONTROL, colnames(NORM.DEG))
  
  # Check genes in DRUG group
  Check_drug = data.frame(Max = apply(NORM.DEG[, Drug, drop = F], 1, max),
                          Sum = apply(NORM.DEG[, Drug, drop = F], 1, sum),
                          Length = length(Drug)) %>% 
    mutate(SpikeCheck = if_else( (Max/Sum) < (1.4*Length^(-0.66)), 'PASS', 'FAIL')) 
  
  
  # Extract genes that PASS in DRUG
  SpikePass_drug = Check_drug[grep('PASS', Check_drug$SpikeCheck), ] %>% rownames
  
  # Check genes in CONTROL group
  Check_cont = data.frame(Max = apply(NORM.DEG[, Cont, drop = F], 1, max),
                          Sum = apply(NORM.DEG[, Cont, drop = F], 1, sum),
                          Length = length(Cont)) %>% 
    mutate(SpikeCheck = if_else( (Max/Sum) < (1.4*Length^(-0.66)) , 'PASS', 'FAIL')) 
  
  # If read count in all replicates is 0, SpikeCheck = <NA>
  # Change all <NA> to PASS
  Check_drug[is.na(Check_drug)] = 'PASS' 
  Check_cont[is.na(Check_cont)] = 'PASS' 
  
  # Extract genes that PASS in CONTROL
  SpikePass_cont = Check_cont[grep('PASS', Check_cont$SpikeCheck), ] %>% rownames
  
  # Intersect the common genes: list of genes that pass the spike filter
  SpikePass = intersect(SpikePass_drug, SpikePass_cont) 
  SpikeFail = NORM.DEG[rownames(NORM.DEG) %notin% SpikePass, , drop = F] %>% rownames	# drop = F is needed when length(SpikeFail)=1
  
  # print info on filtering
  print(paste0(length(SpikePass), ' genes passed the Spike filter.'))
  print(paste0(length(SpikeFail), ' genes failed the Spike filter.'))
  
  # Save list of genes that pass FILTER 3	
  write.table(SpikePass, file = paste('miRNA', DRUG[COMPARISON], ..., 'vs', CONTROL, 'SpikePass.txt', sep = '_'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(SpikeFail, file = paste('miRNA', DRUG[COMPARISON], ..., 'vs', CONTROL, 'SpikeFail.txt', sep = '_'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # assign variable to global env
  assign('SpikePass', SpikePass, envir = .GlobalEnv)
}

# Function that keeps normalized data of DE genes that pass FILTER 2
filterNormData_miRNA = function(COMPARISON, DRUG, CONTROL, SPIKEPASS, RES.DEG, NORM.DEG, OUTPUT.DIR, ...) {
  
  # Select output directory
  if (dir.exists(OUTPUT.DIR)) { 
    setwd(OUTPUT.DIR) 
  } else { 
    dir.create(OUTPUT.DIR)
    setwd(OUTPUT.DIR) }
  
  # Extract genes that pass FILTER 2
  idx = SPIKEPASS
  
  FileName = paste('miRNA', DRUG[COMPARISON], ..., 'vs', CONTROL, 'FDR_0.01', sep = '_')	
  
  # Results of ONLY DE genes that PASS all filters
  res_filtered = RES.DEG[rownames(RES.DEG) %in% idx, ]
  
  # Normalized read count of ONLY DE genes that PASS all filters
  normData_filtered = NORM.DEG[rownames(NORM.DEG) %in% idx, ]
  
  # Info on all filters applied
  print(paste0(nrow(res_filtered), ' out of ', nrow(RES.DEG), ' DE genes passed the Spike filter - ', nrow(RES.DEG[rownames(RES.DEG) %notin% idx, ]), ' DE genes excluded.'))
  
  # Save tables
  write.table(res_filtered, file = paste0(FileName, '_Results_DEG.txt'), sep = '\t', quote = FALSE)
  write.table(rownames(res_filtered), file = paste0(FileName, '_DEG_names.txt'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)	
  write.table(normData_filtered, file = paste0(FileName, '_Norm_count_DEG.txt'), sep = '\t', quote = FALSE)
  
  # assign variable to global env
  assign('res_filtered', res_filtered, envir = .GlobalEnv)
}


# Perform DE analysis on genes
for (x in seq(1, length(Drugs))) { 
  
  DEanalysisCompound(COMPARISON = x, Drugs, Control, genes, 'Compound', unfiltered_dir) 
  CPMfilter(COMPARISON = x, Drugs, Control, dds, 'Compound', filtered_dir)
  spikeFilter(COMPARISON = x, Drugs, Control, normData_DEG, filters_dir)
  filterNormData(COMPARISON = x, Drugs, Control, SpikePass, res, normData_DEG, filtered_dir)
  
} 

# Perform DE analysis on miRNAs
for (x in seq(1, length(Drugs))) { 
  
  DEanalysisCompound_miRNA(COMPARISON = x, Drugs, Control, mirna, 'Compound', unfiltered_dir) 
  CPMfilter_miRNA(COMPARISON = x, Drugs, Control, dds, 'Compound', filtered_dir)
  spikeFilter_miRNA(COMPARISON = x, Drugs, Control, normData_DEG, filters_dir)
  filterNormData_miRNA(COMPARISON = x, Drugs, Control, SpikePass, res, normData_DEG, filtered_dir)
  
} 


# Add gene names to all the output files in the `filtered_dir` to make files inspection easier
setwd(filtered_dir)


for (i in seq(1, length(Drugs))) {
  
  # Result_all
  a = read.table(paste0(Drugs[i], '_vs_', Control, '_FDR_0.01_Results_all.txt'), header = T, sep = '\t')
  a %<>% tibble::rownames_to_column(., 'ensembl_gene_id')
  
  a = left_join(a, mart, by = 'ensembl_gene_id')
  write.table(a, file = paste0(Drugs[i], '_vs_', Control, '_FDR_0.01_Results_all.txt'), quote = F, sep ='\t', row.names = F)
  
  # Results_DEG
  a = read.table(paste0(Drugs[i], '_vs_', Control, '_FDR_0.01_Results_DEG.txt'), header = T, sep = '\t')
  a %<>% tibble::rownames_to_column(., 'ensembl_gene_id')
  
  a = left_join(a, mart, by = 'ensembl_gene_id')
  write.table(a, file = paste0(Drugs[i], '_vs_', Control, '_FDR_0.01_Results_DEG.txt'), quote = F, sep ='\t', row.names = F)
  
  # Norm_count_all
  a = read.table(paste0(Drugs[i], '_vs_', Control, '_FDR_0.01_Norm_count_all.txt'), header = T, sep = '\t')
  a %<>% tibble::rownames_to_column(., 'ensembl_gene_id')
  
  a = left_join(a, mart, by = 'ensembl_gene_id')
  write.table(a, file = paste0(Drugs[i], '_vs_', Control, '_FDR_0.01_Norm_count_all.txt'), quote = F, sep ='\t', row.names = F)
  
  # Norm_count_DEG
  a = read.table(paste0(Drugs[i], '_vs_', Control, '_FDR_0.01_Norm_count_DEG.txt'), header = T, sep = '\t')
  a %<>% tibble::rownames_to_column(., 'ensembl_gene_id')
  
  a = left_join(a, mart, by = 'ensembl_gene_id')
  write.table(a, file = paste0(Drugs[i], '_vs_', Control, '_FDR_0.01_Norm_count_DEG.txt'), quote = F, sep ='\t', row.names = F)
  
  # DEG_names
  a = scan(paste0(Drugs[i], '_vs_DMSO_FDR_0.01_DEG_names.txt'), what = character()) %>% as.data.frame()
  colnames(a) = 'ensembl_gene_id'
  
  mart_names = mart %>% dplyr::select(ensembl_gene_id, external_gene_name)
  
  a = left_join(a, mart_names, by = 'ensembl_gene_id')
  write.table(a, file = paste0(Drugs[i], '_vs_', Control, '_FDR_0.01_DEG_names.txt'), quote = F, sep ='\t', row.names = F)
  
}
```


## 3 - Gene Set Enrichment Analysis on Reactome
How to interpret GSEA results: http://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?_Interpreting_GSEA_Results 

```{r}
library(clusterProfiler)
library(enrichplot)
library(ggridges)
library(pathview)
library(org.Mm.eg.db)
library(tidyverse)
library(ReactomePA)


# set seed
set.seed(1234)

drug = 'DnOP' # can be any of 'DEHP', 'DIDP', 'DINP', 'DnOP'

res = read.table(paste0('project_folder/RNA-Seq/DE_analysis/DE_filtered/', drug, '_vs_DMSO_FDR_0.01_Results_all.txt'), header = T, fill = T) %>% 
  dplyr::rename(ENSEMBL = ensembl_gene_id)


# convert gene names to Entrez ID
ENTREZ = clusterProfiler::bitr(res$ENSEMBL, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db')

# Not all ENSEMBL ID are mapped to ENTREZ. Select only the `stats` value of mapped IDs
res_ENTREZ = left_join(res, ENTREZ, by = 'ENSEMBL') %>% 
  filter(!is.na(ENTREZID)) 

# named vector for analysis (ENTREZ IDs)
gl_entrez = res_ENTREZ$stat
names(gl_entrez) = make.names(res_ENTREZ$ENTREZID, unique = T) %>% 
  gsub('^X', '', .) # somehow, `X` is added at the beginning of the name
gl_entrez = gl_entrez[order(gl_entrez, decreasing = T)]


Reactome = ReactomePA::gsePathway(gl_entrez, 
                                  organism = 'mouse',
                                  pvalueCutoff = 0.2,
                                  pAdjustMethod = "BH", 
                                  verbose = FALSE,
                                  eps = 1e-50,  # or try using eps = 0
                                  seed = T)
View(Reactome@result)

# save table
write.table(Reactome@result, file = paste0('project_folder/RNA-Seq/DE_analysis/reactome_GSEA_', drug, '.txt'), row.names = F, quote = F, sep = '\t')


# Example of leading edge plot ('Metabolism of lipids') 
pw = 'Metabolism of lipids'

p = gseaplot2(Reactome, geneSetID = Reactome@result[Reactome@result$Description == 'Metabolism of lipids',]$ID, 
            title = paste0(drug, ' vs DMSO\n',
                           'Metabolism of lipids',
                           '\n(NES = ', round(Reactome@result[Reactome@result$Description == 'Metabolism of lipids',]$NES, digits = 3),
                           ', q-value = ', round(Reactome@result[Reactome@result$Description == 'Metabolism of lipids',]$qvalue, digits = 3), 
                           ')')) 

ggsave(plot = p, file = paste0('project_folder/RNA-Seq/DE_analysis/plots/Metabolism_of_lipids_', drug, '.png'), 
         dpi = 600, height = 3000, width = 4000, units = 'px')
```


# ATAC-Seq data analysis

## 1 - Process raw data with PEPATAC
The PEPATAC paper is available at https://doi.org/10.1093/nargab/lqab101.

Script: `1_atac_preprocess_pepatac.sh`.
```{bash}
#!/bin/bash 

# PEPATAC version 0.10.4

source /home/m.nazzari/miniconda3/etc/profile.d/conda.sh
conda activate pepatac

DataDir='/project_folder/ATAC-Seq/raw_data/'
cd ${DataDir}

FileNames=( DEHP_1* DEHP_2* DEHP_3* DEHP_4* DEHP_5* DEHP_6* DMSO_1* DMSO_3* DMSO_4* DMSO_5* DMSO_6* DINP_1* DINP_2* DINP_3* DINP_4* DINP_5* DINP_6* )

# get only sample names (e.g. DEHP_1)
Samples=($(printf "%s\n" "${FileNames[@]}" | sed "s/_R/ /" | awk -F"[ ]" "{print $1}" | sort -u))

# Run sample-level pipeline
for S in ${!Samples[@]}; do

echo "Analyzing sample: ${Samples[S]}"

INPUT_R1=${DataDir}${Samples[S]}"_R1_001.fastq.gz"
INPUT_R2=${DataDir}${Samples[S]}"_R2_001.fastq.gz"
	
/project_folder/pepatac/pipelines/pepatac.py -R \
--input ${INPUT_R1} \
--input2 ${INPUT_R2} \
--output-parent /project_folder/ATAC-Seq/processed_data/ \
--sample-name ${Samples[S]:0:6} \
--single-or-paired paired \
--trimmer pyadapt \
--aligner bowtie2 \
--deduplicator samtools \
--peak-caller hmmratac \
--genome GRCh38 \
--genome-size 2.7e9 \
--peak-type variable \
--prealignment-names GRCh38_chrM \
--prealignment-index GRCh38_chrM=/project_folder/genomes/ATAC-Seq_GRCh38/GRCh38_chrM_bowtie2 \
--genome-index /project_folder/genomes/ATAC-Seq_GRCh38/GRCh38_bowtie2 \
--chrom-sizes /project_folder/genomes/ATAC-Seq_GRCh38/GRCh38_chr_size.txt \
--TSS-name /project_folder/genomes/ATAC-Seq_GRCh38/GRCh38_TSS.bed \
--blacklist /project_folder/genomes/hg38-blacklist.v2.bed \
--anno-name /project_folder/genomes/ATAC-Seq_GRCh38/GRCh38_feat_annotations.bed

done

# Run project-level pipeline
/project_folder/pepatac/pipelines/pepatac_collator.py \
--config /project_folder/ATAC-Seq/analysis_config.yaml \
--output-parent /project_folder/ATAC-Seq/processed_data/ \
--results /project_folder/ATAC-Seq/processed_data/ \
--new-start

conda deactivate 
```


## 2 - Shift BAM files to account for transposase duplication
Script `2_shift_bam.sh`.
```{bash}
#!/bin/bash 

# using deeptools to shift the aligned BAM files generated by PEPATAC

# deeptools v3.5.1

Samples=(DEHP_1 DEHP_2 DEHP_3 DEHP_4 DEHP_5 DEHP_6 DMSO_1 DMSO_2 DMSO_3 DMSO_4 DMSO_5 DMSO_6 DINP_1 DINP_2 DINP_3 DINP_4 DINP_5 DINP_6)


for S in ${!Samples[@]}; do

echo "Analyzing sample: ${Samples[S]}"

# Perform ATAC shifting of the BAM file
/home/m.nazzari/Tools/deeptools3.5.1/bin/alignmentSieve \
--numberOfProcessors 16 \
--bam /project_folder/ATAC-Seq/processed_data/${Samples[S]}/aligned_GRCh38/${Samples[S]}_fixed_header.bam \
--ATACshift \
--outFile /project_folder/ATAC-Seq/shifted_bam/${Samples[S]}_shifted_presort.bam

# Sort the new ATAC shifted BAM file
samtools sort -@ 16 /project_folder/ATAC-Seq/shifted_bam/${Samples[S]}_shifted_presort.bam \
-o /project_folder/ATAC-Seq/shifted_bam/${Samples[S]}_ATACshift.bam

# Index the new ATAC shifted BAM file
samtools index /project_folder/ATAC-Seq/shifted_bam/${Samples[S]}_ATACshift.bam \
/project_folder/ATAC-Seq/shifted_bam/${Samples[S]}_ATACshift.bam.bai

rm /project_folder/ATAC-Seq/shifted_bam/${Samples[S]}_shifted_presort.bam

done
```


## 3 - Accessibility analysis using CSAW and sliding window approach
Modified from Sheik and Blais (2022) on bioRvix (https://doi.org/10.1101/2022.03.16.484118).
```{r}
# Set working directory 
setwd('/project_folder/ATAC-Seq/DA_analysis/')

# Load libraries
library(BiocParallel) # Enables parallelization
library(csaw) # Used for window based analysis
library(edgeR) # Used for differential testing
library(rtracklayer) # Used to import BED files as GRanges or export GRanges as BED files
library(RUVSeq) # Used for normalization, remove unwanted variance
library(dplyr)
library(magrittr)

# Specify directory to export results to
read_type = '5-prime-reads' # 'pe-reads' or '5-prime-reads'
window_size = 'window50'
compound = 'DEHP'
exportDir_analysis = paste0('/project_folder/ATAC-Seq/DA_analysis/', compound, '_', read_type, '_', window_size)

# Create export directory in case it doesn't exist
dir.create(exportDir_analysis, showWarnings = TRUE, recursive = TRUE)

# Create a BiocParallel object that uses specified number of CPU threads with a progress bar
my_multicoreParam = MulticoreParam(workers = 12, progressbar = TRUE, exportglobals = FALSE)
register(my_multicoreParam)

# Load tab-delimited TXT file with information on ATACSeq samples (i.e, name, treatment, drug, and the path to appropriate BAM file)
# ATAC-shifted BAM files are used.
atac_samples = read.table('metadata.txt', sep = '\t', header = TRUE) %>%
	as.data.frame() %>%
	dplyr::filter(Group == 'DMSO' | Group == compound)
	
# Set sample names as row names for later downstream analysis
rownames(atac_samples) = atac_samples$Sample

# Specify path to bam files of ATAC samples
bam_files = atac_samples$BAM_Path

# Design to compare samples, matching order of BAM files
# Adjust as needed depending on the study design
treatment = factor(atac_samples$Group)
design = model.matrix(~0 + treatment)
colnames(design) = levels(treatment)
treat_vs_ctrl = makeContrasts(DEHP - DMSO, levels = design)

# Specify chromosomes to limit analysis to (autosomes and sex chrs)
chrLimit = paste0('chr', c(1:22, 'X', 'Y'))

# Set regions to discard by importing blacklist loci BED file
blacklistPath = '/project_folder/genomes/hg38-blacklist.v2.bed'


##### CSAW-specific parameters
# Set parameters for CSAW on regions to discard and chromosomes to restrict analysis to, indicate if data is paired-end (see also CSAW manual)
# To count each read end independently (Tn5 insertions) (i.e. single-end), set `pe = 'none'`
# To use reads as paired end, set `pe = 'both'`
if (read_type == 'pe-reads') { my_pe_mode = 'both' } else if (read_type == '5-prime-reads') { my_pe_mode = 'none' }

# indicate whether to remove duplicates (TRUE) or not (FALSE)
duplicate_removal = F

# Set Window Size
windowSize = 50

# Set spacing between windows. Use `windowSize` to have adjacent windows without intervening gaps
windowSpacing = windowSize

# Set the minimum number of counts for a window to be retained 
minCount_filterSize = 50

# Set bin size for calculating global enrichment
bin_size = 10000

# Set the fold change value to filter windows for enrichment over global enrichment in the binned dataset
globalBckg_foldChange = 3

##### No further parameters to adjust #####


##### CSAW Analysis #####
blacklist = import.bed(blacklistPath)
insertParam = readParam(discard = blacklist, 
                        dedup = duplicate_removal, 
                        minq = NA, 
                        pe = my_pe_mode, 
                        restrict = chrLimit)

# Count genome into bins for filtering via fold change over global background using the insertParam
binned_data = windowCounts(bam_files,
                           bin = TRUE,  
                           param = insertParam, 
                           width = bin_size, 
                           BPPARAM = my_multicoreParam)

# Count Tn5 integration sites with CSAW
window_data = windowCounts(bam_files, 
                           bin = FALSE,
                           ext = 1,
                           filter = minCount_filterSize,
                           param = insertParam,
                           shift = 0,
                           spacing = windowSpacing, 
                           width = windowSize,      
                           BPPARAM = my_multicoreParam)


# Distinguish signal to noise based on global or local background?
filter_criteria = 'global'

if (filter_criteria == 'global') {
	# Filter data by fold change
  
	# Create a filter of windows that have signal well above global background
	keep_fcFilter = filterWindowsGlobal(window_data, background = binned_data)$filter > log2(globalBckg_foldChange)
	summary(keep_fcFilter)
	
	# Filter out unwanted windows by fold change
	filtered_winData = window_data[keep_fcFilter, ]

	# Calculate normalization factors on filtered high-abundance windows with CSAWs implementation of the TMM method to remove efficiency biases. 
	filtered_winData = normFactors(filtered_winData)

} else if (filter_criteria = 'local') {
	
  # filter uninteresting features (windows) by local enrichment
	# local background estimator: 2kb neighborhood
	# change width parameter as desired
	neighbor = resize(rowRanges(window_data), width = 2000, fix = 'center') %>% suppressWarnings()
	wider_data = regionCounts(bam_files, regions = neighbor, param = param) # count reads in neighborhoods
	keep_fcFilter = filterWindowsLocal(counts, wider_data)
	
	# threshold of 3-fold increase in enrichment over 2 kb neighborhood abundance; change as desired 
	filtered_winData = window_data[keep_fcFilter$filter > log2(3), ] 
	filtered_winData = normOffsets(filtered_winData, se.out = TRUE)
}

## Differential accessibility analysis
# Create an edgeR object from the normalized data, set appropriate column and row names
y_TMM = asDGEList(filtered_winData)
colnames(y_TMM$counts) = atac_samples$Sample
rownames(y_TMM$samples) = atac_samples$Sample

# Estimate dispersions and fit a GLM to data
y_TMM = estimateDisp(y_TMM, design)
fit_TMM = glmQLFit(y_TMM, design, robust = TRUE)

# Generate a SeqExpressionSet from TMM data to use with RUVs
set_TMM = newSeqExpressionSet(counts = y_TMM$counts, phenoData = atac_samples)
cpm_TMM = cpm(y_TMM, normalized.lib.sizes = TRUE, log = FALSE)
multiplier_TMM = 1e-6*(mean(y_TMM$samples$lib.size))
normCounts(set_TMM) = round(multiplier_TMM * cpm_TMM, digits = 0)

# See normalized count matrix
set_TMM@assayData$normalizedCounts 

# Plot PCA of data before RUV-Seq 
png(paste0(exportDir_analysis, 'before_RUVs.png'), width = 7000, height = 4000, units = 'px', res = 600)
par(mar = c(8,5,2,2), mfrow = c(1, 2))
plotPCA(set_TMM, col = c(rep('red', 6), rep('blue', 6)), pch = 20, labels = F)
plotRLE(set_TMM, col= c(rep('red',6), rep('blue',6)), las = 2)
dev.off()

# Define the replicate groups
replicateGroups = makeGroups(atac_samples$Group)
replicateGroups

# Create a RUVs set
set_RUVs = RUVs(set_TMM, rownames(set_TMM), k = 5, replicateGroups)

# plot PCA and RLE of data with RUV-Seq
png(paste0(exportDir_analysis, 'after_RUVs.png'), width = 7000, height = 4000, units = 'px', res = 600)
par(mar = c(8,5,2,2), mfrow = c(1, 2))
plotPCA(set_RUVs, col = c(rep('red', 6), rep('blue', 6)), main = 'k = 5', pch = 20, labels = F)
plotRLE(set_RUVs, col = c(rep('red',6), rep('blue',6)), las = 2, main = 'k = 5')
dev.off()

# Create a design matrix for differential testing that incorporates W, the variable representing unwanted variance. 
design_RUVs = model.matrix(~0 + treatment + W_1 + W_2 + W_3 + W_4 + W_5, data = pData(set_RUVs))

# Make contrasts for RUVs data based on the new design matrix
treat_vs_ctrl_RUVs = makeContrasts(treatmentDEHP - treatmentDMSO, levels = design_RUVs)


## Estimate dispersion with one command
y_RUVs = estimateDisp(y_TMM, design_RUVs)

# Fit QL model to data
fit_RUVs = glmQLFit(y_RUVs, design_RUVs)

# Test for differential gene expression
results_RUVs = glmQLFTest(fit_RUVs, contrast = treat_vs_ctrl_RUVs)

# Apply FDR to results table
resultsTable_RUVs = topTags(results_RUVs, adjust.method = 'fdr', n = Inf, sort.by = 'none')$table

# Calculate scores using the FDR values
resultsTable_RUVs$score = -10*log10(resultsTable_RUVs$FDR)

# Add FDR corrected results data to our filtered dataset of windows
winResults_Data = filtered_winData
rowData(winResults_Data) = cbind(rowData(winResults_Data), resultsTable_RUVs)

# Use base R to explore results table
winResults_Data@rowRanges %>% .[order(.$FDR),] %>% head
set_TMM@assayData$normalizedCounts # windows normalized counts 



# Merge nearby windows up to 'tol' distance apart: 500 bp in this case max merged window width: 5000 bp
# Reasoning: FDR correction should be done on `regions` and not `windows` (https://bioconductor.org/books/release/csawBook/correction-for-multiple-testing.html)
merged_peaks = mergeWindows(rowRanges(winResults_Data), tol = 500L, max.width = 5000L)

summary(width(merged_peaks$regions))

# For every region, use the window with lowest p-val as being representative of the whole region (https://bioconductor.org/books/release/csawBook/correction-for-multiple-testing.html). 
tab_best = getBestTest(merged_peaks$id, rowData(winResults_Data))

# Boxplot of regions size distribution
boxplot(width(merged_peaks$regions), main = paste0(compound, ' vs DMSO\nregions size distribution'), ylab = 'Window width, bp')

# Report the start location of the best window within each cluster. For example, the sequence of the DB subinterval can be extracted for motif discovery.
tab_best$rep.start = start(rowRanges(winResults_Data))[tab_best$rep.test]

# combine merged peaks window range with statistics
final_merged_peaks = merged_peaks$region
final_merged_peaks@elementMetadata = cbind(final_merged_peaks@elementMetadata, tab_best)

# give score to regions
final_merged_peaks$score = -10*log10(final_merged_peaks$PValue)
names(final_merged_peaks) = paste0('region_', 1:length(final_merged_peaks))

# Export table of merged regions 
# with dplyr...
final_merged_peaks %>% 
	as.data.frame() %>% 
	write.table(., paste0(exportDir_analysis, 'all_regions.txt'), quote = F, row.names = T, col.names = T, sep = '\t')

# ... or GRanges
export(final_merged_peaks, file = paste0(exportDir_analysis, 'all_regions.bed'))


# isolate significant regions, and ones that go up/down
PVAL = 0.01
final_merged_peaks_sig = subset(final_merged_peaks, abs(rep.logFC) >= 1 & PValue < PVAL)
names(final_merged_peaks_sig) = paste(names(final_merged_peaks_sig), final_merged_peaks_sig$direction, sep = '.')
final_merged_peaks_sig_up = subset(final_merged_peaks_sig, rep.logFC >= 1)
final_merged_peaks_sig_down = subset(final_merged_peaks_sig, rep.logFC <= -1)

# Isolate GRanges objects 
GR_Sig_All = rowRanges(final_merged_peaks_sig)
GR_Sig_Up = rowRanges(final_merged_peaks_sig_up)
GR_Sig_Down = rowRanges(final_merged_peaks_sig_down)

# export significant windows in .bed files
export(GR_Sig_All, file = paste0(exportDir_analysis, compound, '_merged_sigAll_', PVAL, '.bed'))
export(GR_Sig_Up, file = paste0(exportDir_analysis, compound, '_merged_sigUp_', PVAL, '.bed'))
export(GR_Sig_Up, file = paste0(exportDir_analysis, compound, '_merged_sigDown_', PVAL, '.bed'))

# Save image of workspace
save.image(file = paste0(exportDir_analysis, '/', compound, '.RData'))
```


## 4 - Annotate peaks with HOMER 
HOMER main page http://homer.ucsd.edu/homer/index.html.
Script `4_annotate_peaks.sh`.
```{bash}
#!/bin/bash 

# some of these variables depend on the parameters decided during the previous DA analysis (step 3)
READ_TYPE="5-prime-reads" # pe-reads, 5-prime-reads
WINDOW="window50"
DIRECTION="All" # All, Up, Down

# DEHP vs DMSO 
COMPOUND="DEHP"
MAIN_DIR="/project_folder/ATAC-Seq/DA_analysis/${COMPOUND}_${READ_TYPE}_${WINDOW}/"

mkdir -p "${MAIN_DIR}HOMER/"

annotatePeaks.pl \
${MAIN_DIR}/${COMPOUND}_merged_sig${DIRECTION}.bed \
hg38 \
> ${MAIN_DIR}HOMER/${COMPOUND}_sig${DIRECTION}.txt \
-go ${MAIN_DIR}HOMER/${PVAL}_${DIRECTION}_GO \
-annStats ${MAIN_DIR}HOMER/annStats_${DIRECTION}.txt

# DINP vs DMSO
COMPOUND="DINP"
MAIN_DIR="/project_folder/ATAC-Seq/DA_analysis/${COMPOUND}_${READ_TYPE}_${WINDOW}/"

mkdir -p "${MAIN_DIR}HOMER/"

annotatePeaks.pl \
${MAIN_DIR}/${COMPOUND}_merged_sig${DIRECTION}.bed \
hg38 \
> ${MAIN_DIR}HOMER/${COMPOUND}_sig${DIRECTION}.txt \
-go ${MAIN_DIR}HOMER/${PVAL}_${DIRECTION}_GO \
-annStats ${MAIN_DIR}HOMER/annStats_${DIRECTION}.txt
```


## 5 - Take screenshots of tracks on IGV normalized with BeCorrect
BeCorrect is used to scale track files based on the normalization values. It's used for visualization purposes. Here I used it to obtain scaled track files to visualize on the Integrative Genomics Viewer.
BeCorrect is available at https://doi.org/10.1038/s41598-020-66998-4.


### 5.1 - Convert shifted BAM to BEDGRAPH with deepTools
Script `5.1_BAM_to_BEDGRAPH.sh`.
```{bash}
#!/bin/bash 

# deepTools v3.5.1

# DEHP vs DMSO
Samples=(DEHP_1 DEHP_2 DEHP_3 DEHP_4 DEHP_5 DEHP_6 DMSO_1 DMSO_2 DMSO_3 DMSO_4 DMSO_5 DMSO_6)
fileNames=(/project_folder/ATAC-Seq/shifted_bam/DEHP*ATACshift.bam /project_folder/shifted_bam/DMSO*ATACshift.bam)

COMPOUND="DEHP"
READ_TYPE="5-prime-reads" # pe-reads, 5-prime-reads
WINDOW="window50"
MAIN_DIR="/project_folder/ATAC-Seq/DA_analysis/${COMPOUND}_${READ_TYPE}_${WINDOW}/"

mkdir -p "${MAIN_DIR}BeCorrect/"


for S in ${!Samples[@]}; do

echo "Converting sample: ${Samples[S]}"

/home/Tools/deeptools3.5.1/bin/bamCoverage --bam ${fileNames[S]} -o ${MAIN_DIR}/BeCorrect/1_raw_bedgraph/${Samples[S]}.bg --outFileFormat "bedgraph" --binSize 1 --numberOfProcessors 25

done


# DINP vs DMSO
Samples=(DINP_1 DINP_2 DINP_3 DINP_4 DINP_5 DINP_6 DMSO_1 DMSO_2 DMSO_3 DMSO_4 DMSO_5 DMSO_6)
fileNames=(/project_folder/DA_analysis/DINP*ATACshift.bam /project_folder/DA_analysis/DMSO*ATACshift.bam)

COMPOUND="DINP"
READ_TYPE="5-prime-reads" # pe-reads, 5-prime-reads
WINDOW="window50"
MAIN_DIR="/project_folder/ATAC-Seq/DA_analysis/${COMPOUND}_${READ_TYPE}_${WINDOW}/"

mkdir -p "${MAIN_DIR}BeCorrect/"


for S in ${!Samples[@]}; do

echo "Converting sample: ${Samples[S]}"

/home/Tools/deeptools3.5.1/bin/bamCoverage --bam ${fileNames[S]} -o ${MAIN_DIR}/BeCorrect/1_raw_bedgraph/${Samples[S]}.bg --outFileFormat "bedgraph" --binSize 1 --numberOfProcessors 25

done
```


### 5.2 - Retrieve raw and normalized windows counts for BeCorrect visualization
```{r}
library(tidyverse)

# load image saved on step 3
COMPOUND = 'DINP' # or 'DEHP'
READ_TYPE = '5-prime-reads' # pe-reads, 5-prime-reads
WINDOW = 'window50'
MAIN_DIR = paste0('/project_folder/ATAC-Seq/DA_analysis/', COMPOUND, '_', READ_TYPE, '_', WINDOW, '/')

setwd(MAIN_DIR)

load(paste0(COMPOUND, '.RData'))

# retrieve counts table with chr start end as columns. Remove any chromosome that doesn't start with 'chr'
normalized_counts = do.call('cbind', list( 
				chr = as.data.frame(winResults_Data@rowRanges@seqnames),
				start = as.data.frame(winResults_Data@rowRanges@ranges)[, 1],
				end = as.data.frame(winResults_Data@rowRanges@ranges)[, 2],
				set_TMM@assayData$normalizedCounts)) %>% 
	dplyr::rename(., chr = 'value') %>% # if column isn't automatically named `chr`, set it manually
	dplyr::filter(str_detect(chr, 'chr'))

# normalized count file 
write.table(normalized_counts, paste0(COMPOUND, '_normalized_counts.txt'), quote = F, row.names = F, col.names = T, sep = '\t')

# raw count file
raw_count = do.call('cbind', list( 
				chr = as.data.frame(winResults_Data@rowRanges@seqnames),
				start = as.data.frame(winResults_Data@rowRanges@ranges)[, 1],
				end = as.data.frame(winResults_Data@rowRanges@ranges)[, 2],
				set.TMM@assayData$counts)) %>% 
	dplyr::rename(., chr = 'value') %>% # if column isn't automatically named `chr`, set it manually
	dplyr::filter(str_detect(chr, 'chr'))


write.table(raw_count, paste0(COMPOUND, '_raw_counts.txt'), quote = F, row.names = F, col.names = T, sep = '\t')
```


### 5.3 - Make chromosomes names file 
TXT file with names of chromosomes. Retrieve from the file used for PEPATAC (`/project_folder/genomes/ATAC-Seq_GRCh38/GRCh38_chr_size.txt`). If needed, modify file to only retain autosomes and sex chromosomes.


### 5.4 - Retain only regions that fall within autosomes and sex chromosomes from the BEDGRAPH files exported in Step 3
```{r}
READ_TYPE = '5-prime-reads' # pe-reads, 5-prime-reads
WINDOW = 'window50'

COMPOUND = c('DEHP', 'DINP')

# For both DEHP and DINP
for (x in 1:2) {

	MAIN_DIR = paste0('/project_folder/ATAC-Seq/DA_analysis/', COMPOUND[x], '_', READ_TYPE, '_', WINDOW, '/')
	
	bg_list = list.files(path = paste0(MAIN_DIR, 'BeCorrect/1_raw_bedgraph'), pattern = '.bg', full.names = T)

	for (i in seq(1, length(bg_list))) {
	
		read.table(bg_list[i], header = F) %>%
			dplyr::filter(V4 != 0) %>% # remove lines which are 0 to decrease file size
			dplyr::filter(str_detect(V1, 'chr')) %>%
			write.table(., file = bg_list[i], quote = F, row.names = F, col.names = F, sep = '\t')
	}
}
```


### 5.5 - Download BeCorrect script and compile it on the terminal
```{bash}
wget -P /project_folder/DA_analysis/ https://raw.githubusercontent.com/Zhang-lab/BeCorrect/master/BeCorrect.cpp 
g++ /project_folder/DA_analysis/BeCorrect.cpp -o /project_folder/DA_analysis/BeCorrect
```


### 5.6 - Run BeCorrect
Script `5.6_BeCorrect.sh`.
```{bash}
READ_TYPE="5-prime-reads" # pe-reads, 5-prime-reads
WINDOW="window50"

# DEHP
COMPOUND="DEHP"
MAIN_DIR="/project_folder/DA_analysis/${COMPOUND}_${READ_TYPE}_${WINDOW}/"

mkdir -p ${MAIN_DIR}BeCorrect/2_scaled_bedgraph/
cd -p ${MAIN_DIR}BeCorrect/2_scaled_bedgraph/

/project_folder/DA_analysis/BeCorrect \
/project_folder/DA_analysis/${MAIN_DIR}/DEHP_raw_counts.txt \
/project_folder/DA_analysis/${MAIN_DIR}/DEHP_normalized_counts.txt \
/project_folder/DA_analysis/Chr_names_GRCh38.txt \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DEHP_1.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DEHP_2.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DEHP_3.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DEHP_4.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DEHP_5.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DEHP_6.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DMSO_1.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DMSO_2.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DMSO_3.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DMSO_4.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DMSO_5.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DMSO_6.bg


# DINP
COMPOUND="DINP"
MAIN_DIR="/project_folder/DA_analysis/${COMPOUND}_${READ_TYPE}_${WINDOW}/"

mkdir -p ${MAIN_DIR}BeCorrect/2_scaled_bedgraph/
cd -p ${MAIN_DIR}BeCorrect/2_scaled_bedgraph/

/project_folder/DA_analysis/BeCorrect \
/project_folder/DA_analysis/${MAIN_DIR}/DINP_raw_counts.txt \
/project_folder/DA_analysis/${MAIN_DIR}/DINP_normalized_counts.txt \
/project_folder/DA_analysis/Chr_names_GRCh38.txt \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DINP_1.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DINP_2.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DINP_3.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DINP_4.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DINP_5.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DINP_6.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DMSO_1.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DMSO_2.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DMSO_3.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DMSO_4.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DMSO_5.bg \
/project_folder/DA_analysis/${MAIN_DIR}/BeCorrect/1_raw_bedgraph/DMSO_6.bg
```


### 5.7 - If the BEDGRAPH files are too big and slow down IGV, you can convert to BIGWIG (the binary equivalent)
Script `5.7_BEDGRAPH_to_BIGWIG`.
Use UCSC utility http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
```{bash}

READ_TYPE="5-prime-reads" # pe-reads, 5-prime-reads
WINDOW="window50"

## DEHP vs DMSO
COMPOUND="DEHP"
MAIN_DIR="/project_folder/DA_analysis/${COMPOUND}_${READ_TYPE}_${WINDOW}/"

cd -p ${MAIN_DIR}BeCorrect/2_scaled_bedgraph/

# rename files to change extension from .bg to .bedGraph (the official one)
for file in *.bg; do mv -- '$file' '${file%.bg}.bedGraph'; done

# convert to BIGWIG
Samples=(DEHP_1 DEHP_2 DEHP_3 DEHP_4 DEHP_5 DEHP_6 DMSO_1 DMSO_2 DMSO_3 DMSO_4 DMSO_5 DMSO_6)
fileNames=( ${MAIN_DIR}BeCorrect/2_scaled_bedgraph/* )

for S in ${!Samples[@]}; do

echo 'Converting from BEDGRAPH to BIGWIG sample: ${Samples[S]}'

# sort BEDGRAPH
sort -k1,1 -k2,2n ${fileNames[S]} > ${MAIN_DIR}BeCorrect/3_sorted_scaled_bedgraph/${Samples[S]}_sorted.bedGraph

/home/bedGraphToBigWig \
${MAIN_DIR}BeCorrect/3_sorted_scaled_bedgraph/${Samples[S]}_sorted.bedGraph \
/project_folder/genomes/ATAC-Seq_GRCh38/GRCh38_chr_size.txt \
${MAIN_DIR}BeCorrect/4_scaled_bigwig/${Samples[S]}.bw

done


## DINP vs DMSO
COMPOUND="DINP"
MAIN_DIR="/project_folder/DA_analysis/${COMPOUND}_${READ_TYPE}_${WINDOW}/"

cd -p ${MAIN_DIR}BeCorrect/2_scaled_bedgraph/

# rename files to change extension from .bg to .bedGraph (the official one)
for file in *.bg; do mv -- '$file' '${file%.bg}.bedGraph'; done

# convert to BIGWIG
Samples=(DINP_1 DINP_2 DINP_3 DINP_4 DINP_5 DINP_6 DMSO_1 DMSO_2 DMSO_3 DMSO_4 DMSO_5 DMSO_6)
fileNames=( ${MAIN_DIR}BeCorrect/2_scaled_bedgraph/* )

for S in ${!Samples[@]}; do

echo 'Converting from BEDGRAPH to BIGWIG sample: ${Samples[S]}'

# sort BEDGRAPH
sort -k1,1 -k2,2n ${fileNames[S]} > ${MAIN_DIR}BeCorrect/3_sorted_scaled_bedgraph/${Samples[S]}_sorted.bedGraph

/home/bedGraphToBigWig \
${MAIN_DIR}BeCorrect/3_sorted_scaled_bedgraph/${Samples[S]}_sorted.bedGraph \
/project_folder/genomes/ATAC-Seq_GRCh38/GRCh38_chr_size.txt \
${MAIN_DIR}BeCorrect/4_scaled_bigwig/${Samples[S]}.bw

done
```


### 5.8 - Make IGV batch script
Two scripts are created, PART1 and PART2. The usage is: run PART1, group all tracks and autoscale, then run PART2. I have not yet found a successful way to group and autoscale from within the script!
```{r}
# add -1000 and +1000 to BED file coordinates 
padding = 1000

READ_TYPE = '5-prime-reads' # pe-reads, 5-prime-reads
WINDOW = 'window50'
DIRECTION = 'All' # All, Up, Down
COMPOUND = c('DEHP', 'DINP')

# For both DEHP and DINP
for (x in 1:2) {

	MAIN_DIR = paste0('/project_folder/ATAC-Seq/DA_analysis/', COMPOUND[x], '_', READ_TYPE, '_', WINDOW, '/')

	bigwig_compound = list.files(paste0(MAIN_DIR, 'BeCorrect/4_scaled_bigwig'), pattern = compound, full.names = T)
	bigwig_dmso = list.files(paste0(MAIN_DIR, 'BeCorrect/4_scaled_bigwig'), pattern = 'DMSO', full.names = T)

	# read BED file with regions to screenshot
	mybed = read.table(paste0(MAIN_DIR, compound, '_merged_sigAll_0.01.bed'), header = F)

	## PART 1: load and color tracks
	# add LOAD GENOME and ANNOTATIONS lines
	part1 = do.call('rbind', list('genome /project_folder/genomes/ATAC-Seq_GRCh38/GRCh38.primary_assembly.genome.fa',
					'load /project_folder/genomes/ATAC-Seq_GRCh38/gencode.v38.primary_assembly.annotation.sorted.gtf',
					'load /project_folder/genomes/ATAC-Seq_GRCh38/GRCh38-ccREs.bed',
					paste0('load ', MAIN_DIR, '/', compound, '_merged_sigAll_0.01.bed')))

	# load and color compound tracks. Different colors for DEHP and DINP
	if (compound == 'DEHP') {
		for (i in seq(1, 6)) {
			part1 = do.call('rbind', list(part1, 
							paste0('load ', bigwig_compound[i]),
							paste0('setColor #85CC6F ', bigwig_compound[i]),
							paste0('setDataRange auto ', bigwig_compound[i])))
					
	}} else if (compound == 'DINP') {
		for (i in seq(1, 6)) {
			part1 = do.call('rbind', list(part1, 
							paste0('load ', bigwig_compound[i]),
							paste0('setColor #7C80D1 ', bigwig_compound[i]),
							paste0('setDataRange auto ', bigwig_compound[i])))
	}}	
	
	# load and color DMSO tracks
	for (i in seq(1, 6)) {
		part1 = do.call('rbind', list(part1, 
						paste0('load ', bigwig_dmso[i]),
						paste0('setColor #7F7F7F ', bigwig_dmso[i]),
						paste0('setDataRange auto ', bigwig_dmso[i])))
	}


	## PART 2: add take screenshots lines. The filename is structure geneName_compound_regionName.png
	part2 = paste0('snapshotDirectory ', MAIN_DIR, '/IGV_snapshots/')
			  
	for (i in seq(1, nrow(mybed))) {
		
		# gene name comes from the HOMER annotation file created in step 4
		homer = read.table(paste0(MAIN_DIR, 'HOMER/', COMPOUND, '_sig', DIRECTION, '.txt'), header = T, sep = '\t')
		
	part2 = do.call('rbind', list(part2, 
					paste0('goto ', mybed$V1[i], ':', mybed$V2[i]-padding, '-', mybed$V3[i]+padding),
					'maxPanelHeight 3000',
					paste0('snapshot ', homer[homer$PeakID == mybed$V4[i], ]$Gene_Name, '_', compound, '_', mybed$V4[i], '.png')))
					
	}

	write.table(part1, paste0(MAIN_DIR, 'PART1_script_IGV_batch_screenshot.txt'), quote = F, row.names = F, col.names = F)
	write.table(part2, paste0(MAIN_DIR, 'PART2_script_IGV_batch_screenshot.txt'), quote = F, row.names = F, col.names = F)
}
```


### 5.9 - Run script on IGV 
How to (DEHP or DINP):

1. Open IGV
2. Run part 1 script PART1_script_IGV_batch_screenshot.txt 
3. Group and autoscale manually
4. Run part 2 script PART2_script_IGV_batch_screenshot.txt 


### 5.10 - Move UP and DOWN screenshots to their respective folder 
Script `5.10_move_screenshots.sh`.
```{bash}
READ_TYPE="5-prime-reads" # pe-reads, 5-prime-reads
WINDOW="window50"

# DEHP
COMPOUND="DEHP"
MAIN_DIR="/project_folder/ATAC-Seq/DA_analysis/${COMPOUND}_${READ_TYPE}_${WINDOW}/"

cd -p ${MAIN_DIR}IGV_snapshots/
mkdir -p ./up/
mkdir -p ./down/
mv *up.png ./up/
mv *down.png ./down/


# DINP
COMPOUND="DINP"
MAIN_DIR="/project_folder/ATAC-Seq/DA_analysis/${COMPOUND}_${READ_TYPE}_${WINDOW}/"

cd -p ${MAIN_DIR}IGV_snapshots/
mkdir -p ./up/
mkdir -p ./down/
mv *up.png ./up/
mv *down.png ./down/
```
