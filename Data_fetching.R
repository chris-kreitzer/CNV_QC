# Investigate TCGA count_metrices WES with FACETS
# snp-pileups are fetched from /ifs/res/taylorlab/tcga_facets/snp-pileup
# based on NGS data
# compare Affymetrix SNP data with Facets reported CNA
# starting with TCGA-PRAD samples

library(parallel)
library(TCGAbiolinks)

# TCGAbiolinks = API for GDC data fetch
TCGA.PRAD.CODES = TCGAbiolinks::getSampleFilesSummary("TCGA-PRAD")$.id  # n = 500 TCGA-PRAD
TCGA.PRAD.CODES = paste0(TCGA.PRAD.CODES, '-10A')
snp_pileups = list.files(path = '~/Desktop/mnt/taylorlab/tcga_facets/snp-pileup/')

PRAD_Files = mclapply(unique(TCGA.PRAD.CODES), function(x){
  grep(x, snp_pileups, value = T)
})

PRAD_Files = unlist(PRAD_Files)

