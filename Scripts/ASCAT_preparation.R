setwd('~/Documents/MSKCC/CPNA_analysis/TCGA/')


prepare_ASCCAT = function(countmatrix.path){
  library(vroom)
  pileup <- vroom:vroom(file = countmatrix.path, delim = ',')
  col_types = c(Chromosome = "c",
                Position = 'i',
                Ref = 'c',
                Alt = 'c',
                File1R = 'd',
                File1A = 'd',
                File1E = 'd',
                File1D = 'd',
                File2R = 'd',
                File2A = 'd',
                File2E = 'd',
                File2D = 'd'),
  progress = FALSE)
  
  ## keep only finite entries in data input:
  ii = which(pileup$File1E <= Inf & pileup$File1D <= Inf &
               pileup$File2E <= Inf & pileup$File2D <= Inf)
  
  count.matrix = pileup[ii, 1:2]
  count.matrix$all_normal = pileup$File1R[ii] + pileup$File1A[ii]
  count.matrix$ref_normal = pileup$File1R[ii]
  count.matrix$all_tumor = pileup$File2R[ii] + pileup$File2A[ii]
  count.matrix$ref_tumor = pileup$File2R[ii]
  rm(ii)
  # subset chromosomes
  count.matrix$Chromosome = gsub('chr', '', count.matrix$Chromosome)
  chromosome_levels = c(1:22, 'X')
  count.matrix = count.matrix[which(count.matrix$Chromosome %in% chromosome_levels), ]
  rm(chromosome_levels)
  
  # some thresholds; review those thresholds; especially with high-coverage sequencing 
  seq_depth = (count.matrix$all_normal >= 10) & (count.matrix$all_normal <= 1000)
  count.matrix = count.matrix[seq_depth, ]
  
  # prepare output for ASCAT analysis:
  out_ascat = list()
  
  # normal BAF and LogR
  normal_BAF = data.frame(Chr = count.matrix$Chromosome,
                          Position = count.matrix$Position,
                          N1 = 1 - (count.matrix$ref_normal / count.matrix$all_normal))
  rownames(normal_BAF) = paste0('SNP', seq(1, nrow(normal_BAF), 1))
  
  normal_LogR = 1
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # tumor BAF
  cat('Calculating BAF and LogR for tumor sample')
  
  tumor_BAF = data.frame(Chr = count.matrix$Chromosome,
                         Position = count.matrix$Position,
                         S1 = 1 - (count.matrix$ref_tumor / count.matrix$all_tumor))
  
  tumor_BAF = tumor_BAF[!is.infinite(tumor_BAF$S1), ]
  rownames(tumor_BAF) = paste0('SNP', seq(1, nrow(tumor_BAF), 1))
  
  ## tumor LogR
  tumor_LogR = data.frame(Chr = count.matrix$Chromosome,
                          Position = count.matrix$Position,
                          S1 = log2((count.matrix$all_tumor / count.matrix$all_normal)) -
                            median(log2((count.matrix$all_tumor / count.matrix$all_normal))))
  
  tumor_LogR = tumor_LogR[!is.infinite(tumor_LogR$S1), ]
  rownames(tumor_LogR) = paste0('SNP', seq(1, nrow(tumor_LogR), 1))
  
  out_ascat[[1]] = normal_BAF
  out_ascat[[2]] = normal_LogR
  out_ascat[[3]] = tumor_BAF
  out_ascat[[4]] = tumor_LogR
  names(out_ascat) = c('normal_BAF', 'normal_LogR', 'tumor_BAF', 'tumor_LogR')
  
  rm(normal_BAF)
  rm(normal_LogR)
  rm(tumor_BAF)
  rm(tumor_LogR)
  rm(count.matrix)
  gc()
  return(out_ascat)
  
}



# # preparing dummy data 
# data.in = prepare_ASCCAT(countmatrix.path = '~/Desktop/tmp_TCGA/tmp_TCGA/TCGA-2A-A8VV-10A-01D-A37A-08_TCGA-2A-A8VV-01A-11D-A377-08_3emMg7IK_countsMerged___normal_tumor.dat.gz')
# data.normal = data.in$normal_BAF
# 
# set.seed(123)
# keep1 = data.normal[sample(which(data.normal$N1 < 0.3 | data.normal$N1 > 0.7), size = 100000, replace = T), ]
# keep2 = data.normal[sample(which(data.normal$N1 > 0.3 & data.normal$N1 < 0.7), size = 100000, replace = T), ]
# 
# 
# # Germline_BAF
# keep3 = rbind(keep1, keep2)
# 
# # Tumor_BAF
# tumor_baf = data.in$tumor_BAF[rownames(data.in$tumor_BAF) %in% rownames(keep3), ]
# tumor_logR = data.in$tumor_LogR[rownames(data.in$tumor_LogR) %in% rownames(keep3), ]
# 
# 
# Germline_LogR = data.frame(Chr = 1,
#                            Position = 1,
#                            N1 = 1)
# rownames(Germline_LogR) = 1
# 
# 
# write.table(keep3, 
#             file = 'tmp/Germline_BAF.txt', 
#             row.names = T, 
#             sep = '\t')
# 
# write.table(tumor_baf, 
#             file = 'tmp/Tumor_BAF.txt',
#             row.names = T,
#             sep = '\t')
# 
# write.table(tumor_logR, 
#             file = 'tmp/Tumor_LogR.txt',
#             row.names = T,
#             sep = '\t')
# 
# write.table(Germline_LogR, 
#             file = 'tmp/Germline_LogR.txt',
#             row.names = T,
#             sep = '\t')
# 
# write.table(Germline_LogR, 
#             file = 'tmp/Germline_LogR.txt',
#             row.names = T,
#             sep = '\t')
