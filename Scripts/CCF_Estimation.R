## The influnce of varying purity Estimates:
## First, I will look into CCF estimations:

# First I will start and see how purity estimates influence CCF estimations:
# For this we will use FacetsSuite: ccf_Annotate_MAF


maf <- vroom::vroom('../../CPNA_analysis/TCGA/RawData/mc3.v0.2.8.PUBLIC.maf', delim = '\t')
maf = as.data.frame(maf[, c('Hugo_Symbol', 'Tumor_Sample_Barcode', 
                            'Chromosome', 'Start_Position', 
                            'End_Position', 't_ref_count', 't_alt_count', 
                            'n_depth', 'n_ref_count', 'n_alt_count')])

files = list.files('~/Desktop/tmp_TCGA/tmp_TCGA/', full.names = T)
ccf_all_out = data.frame()
for(i in 1:length(files)){
  print(which(files[i]))
  
  try({
    sample = substr(files[i], start = 78, nchar(files)[i]-56)
    
    data.in = facets::readSnpMatrix(filename = files[i])
    data.processed = facetsSuite::run_facets(data.in, 
                                        cval = 100,
                                        genome = 'hg19')
    data.segs = data.processed$segs
    data.purity = data.processed$purity
    data.purity = ifelse(is.na(data.purity), 0, data.purity)
    
    # maf
    maf_ccf = maf[grep(sample, maf$Tumor_Sample_Barcode), ]
    
    CCF_tmp = facetsSuite::ccf_annotate_maf(maf = maf_ccf,
                                            segs = data.segs,
                                            purity = data.purity)
    
    ccf_all_out = rbind(ccf_all_out, CCF_tmp)
    
    rm(data.in)
    rm(data.processed)
    rm(data.segs)
    rm(data.purity)
    rm(CCF_tmp)
    rm(maf_ccf)
    rm(sample)
    
  })
}

# write.table(ccf_all_out, file = '~/Documents/MSKCC/CPNA_analysis/TCGA/tmp/ccf_out.txt', sep = '\t', row.names = F, quote = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# investigate the CCF output on n = 203 TCGA-PRAD samples;

# how many patient are monoclonal:
x = read.csv('~/Documents/MSKCC/CPNA_analysis/TCGA/tmp/ccf_out.txt', sep = '\t')

patient.out = data.frame()
for(i in unique(x$Tumor_Sample_Barcode)){
  sample.patient = x[which(x$Tumor_Sample_Barcode == i), ]
  table.out = as.data.frame(table(sample.patient$clonality))
  summary.data = data.frame(sample = i,
                            n.clonal = ifelse(identical(table.out$Freq[which(table.out$Var1 == 'CLONAL')],
                                                        integer(0)), 0, table.out$Freq[which(table.out$Var1 == 'CLONAL')]),
                            n.subclonal = ifelse(identical(table.out$Freq[which(table.out$Var1 == 'SUBCLONAL')],
                                                           integer(0)), 0, table.out$Freq[which(table.out$Var1 == 'SUBCLONAL')]),
                            n.inter = ifelse(identical(table.out$Freq[which(table.out$Var1 == 'INDETERMINATE')],
                                                       integer(0)), 0, table.out$Freq[which(table.out$Var1 == 'INDETERMINATE')]))
  patient.out = rbind(patient.out, summary.data)
  rm(sample.patient)
  rm(summary.data)
}

# assign monoclonal VS polyclonal
patient.out$diversity = ifelse(patient.out$n.clonal / 
                                 (patient.out$n.clonal + patient.out$n.subclonal + patient.out$n.inter) >= 0.80,
                               'monoclonal', 'polyclonal')




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# look into ABSOLUTE calls and estimate 
Absolute.calls = read.csv('~/Documents/MSKCC/CPNA_analysis/TCGA/RawData/TCGA_ABSOLUTE_CPN.txt', sep = '\t')
Absolute.calls = Absolute.calls[Absolute.calls$Sample %in% substr(x$Tumor_Sample_Barcode, start = 1, stop = 15), ]
TCGA.purity = TCGAbiolinks::Tumor.purity
TCGA.purity = TCGA.purity[substr(TCGA.purity$Sample.ID, start = 1, stop = 15) %in% maf.absolute$Sample, ]
TCGA.purity$Sample.ID = substr(TCGA.purity$Sample.ID, start = 1, stop = 15)
TCGA.purity$CPE = as.numeric(as.character(gsub(pattern = ',', '.', TCGA.purity$CPE)))
TCGA.purity = TCGA.purity[!is.na(TCGA.purity$CPE), ]

maf <- read.csv('../../CPNA_analysis/TCGA/RawData/mc3.v0.2.8.PUBLIC.maf', sep = '\t')
maf = as.data.frame(maf[, c('Hugo_Symbol', 'Tumor_Sample_Barcode', 
                            'Chromosome', 'Start_Position', 
                            'End_Position', 't_ref_count', 't_alt_count', 
                            'n_depth', 'n_ref_count', 'n_alt_count')])

maf.absolute = maf[maf$Tumor_Sample_Barcode %in% x$Tumor_Sample_Barcode, ]
maf.absolute$Tumor_Sample_Barcode = substr(maf.absolute$Tumor_Sample_Barcode, start = 1, stop = 15)
colnames(maf.absolute)[2] = 'Sample'

tcga.samples = Reduce(intersect, list(Absolute.calls$Sample, TCGA.purity$Sample.ID, maf.absolute$Sample))
tcga.samples = intersect(tcga.samples, Absolute.calls$Sample)
tcga.samples = unique(tcga.samples)

ccf_absolute_out = data.frame()
for(i in 1:length(tcga.samples)){
  try({
    print(i)
    maf.selected = maf.absolute[which(maf.absolute$Sample == tcga.samples[i]), ]
    maf.selected$Chromosome = as.character(maf.selected$Chromosome)
    data.table::setDT(maf.selected, key = c('Chromosome', 'Start_Position', 'End_Position'))
    
    cpn.selected = Absolute.calls[which(Absolute.calls$Sample == tcga.samples[i]), ]
    cpn.selected$Chromosome = as.character(cpn.selected$Chromosome)
    data.table::setDT(cpn.selected, key = c('Chromosome', 'Start', 'End'))
    
    # data merged:
    data.merged = data.table::foverlaps(maf.selected,
                                        cpn.selected,
                                        by.x = c('Chromosome', 'Start_Position', 'End_Position'),
                                        by.y = c('Chromosome', 'Start', 'End'),
                                        type = 'within',
                                        mult = 'first',
                                        nomatch = NA)
    
    purity = TCGA.purity[which(TCGA.purity$Sample.ID == tcga.samples[i]), 'CPE']
    
    # CCF estimation:
    #data.selected = data.merged[,c(1, 2, 9, 21, 23, 24, 25, 26, 27, 28, 29)]
    data.selected = as.data.frame(data.merged)
    
    for(mutation in 1:nrow(data.selected)){
      
      mutant_copies = expected_mutant_copies(t_var_freq = data.selected$t_alt_count[mutation] / (data.selected$t_ref_count[mutation] + data.selected$t_alt_count[mutation]), 
                                             total_copies = data.selected$Modal_Total_CN[mutation],
                                             purity = purity)
      
      out = estimate_ccf(purity = purity,
                         total_copies = data.selected$Modal_Total_CN[mutation],
                         mutant_copies = mutant_copies,
                         t_alt_count = data.selected$t_alt_count[mutation],
                         t_depth = data.selected$t_alt_count[mutation] + data.selected$t_ref_count[mutation])
      
      ccf_mutant = out[[1]]
      ccf_upper = out[[3]]
      clonality = ifelse(ccf_mutant >= 0.8 | (ccf_mutant >= 0.7 & ccf_upper >= 0.9), 'Clonal',
                         ifelse(ccf_mutant < 0.8 & ccf_mutant > 0.001 | (ccf_mutant < 0.7 | ccf_mutant > 0.001 & ccf_upper < 0.9 | ccf_upper > 0.001), 
                                'Subclonal', 'NA'))
        
      report.out = data.frame(sample = tcga.samples[i],
                              Hugo_Symbol = data.selected$Hugo_Symbol[mutation],
                              chrom = data.selected$Chromosome[mutation],
                              Start = data.selected$Start_Position[mutation],
                              purity = purity,
                              totalCPN = data.selected$Modal_Total_CN[mutation],
                              mutant_copies = mutant_copies,
                              CCF = ccf_mutant,
                              cnA1 = data.selected$Cancer_cell_frac_a1[mutation],
                              cnA2 = data.selected$Cancer_cell_frac_a2[mutation],
                              clonality = clonality)
      
      
      ccf_absolute_out = rbind(ccf_absolute_out, report.out)
      
    }
    
  })
  
  rm(maf.selected)
  rm(cpn.selected)
  rm(data.merged)
  rm(data.selected)
  rm(out)
  rm(ccf_mutant)
  rm(ccf_upper)
  rm(report.out)
  rm(clonality)
  rm(mutant_copies)
  rm(purity)
  
}

## diversity estimate on ABSOLUTE:
ccf_out = data.frame()
for(i in unique(ccf_absolute_out$sample)){
  sample.patient = ccf_absolute_out[which(ccf_absolute_out$sample == i), ]
  table.out = as.data.frame(table(sample.patient$clonality))
  summary.data = data.frame(sample = i,
                            n.Clonal = ifelse(identical(table.out$Freq[which(table.out$Var1 == 'Clonal')],
                                                        integer(0)), 0, table.out$Freq[which(table.out$Var1 == 'Clonal')]),
                            n.Subclonal = ifelse(identical(table.out$Freq[which(table.out$Var1 == 'Subclonal')],
                                                           integer(0)), 0, table.out$Freq[which(table.out$Var1 == 'Subclonal')]),
                            n.NA = length(sample.patient$clonality[which(is.na(sample.patient$clonality))]))
  ccf_out = rbind(ccf_out, summary.data)
  rm(sample.patient)
  rm(summary.data)
}

ccf_out$diversity = ifelse(ccf_out$n.Clonal / 
                                 (ccf_out$n.Clonal + ccf_out$n.Subclonal + ccf_out$n.NA) >= 0.80,
                               'monoclonal', 'polyclonal')


# make a plot to compare the CCF density of certain genes;
# now I compare CCF obtained from ABSOLUTE/TP purity with Facets data
plot(1,
     type = 'n',
     xlim = c(0, 1.3),
     ylim = c(0, 5.8),
     axes = F,
     xlab = 'CCF estimation',
     ylab = '',
     pch = 16,
     main = '',
     xaxt = 'n',
     xaxs = 'i',
     yaxs = 'i')

axis(side = 2, at = c(0, 6), labels = FALSE, tck = 0.0, lwd = 2)
mtext(side = 2, text = "Density (a.u.)", line = 1)

axis(side = 1,
     at = c(0, 0.5, 1),
     labels = c('0', '0.5', '1'),
     las = 1,
     lwd = 2)

# add Facets density of TP53
lines(density(x$ccf_expected_copies[which(x$Hugo_Symbol == 'SPOP')], na.rm = T), lwd = 3, col = '#392C8D')
lines(density(ccf_absolute_out$CCF[which(ccf_absolute_out$Hugo_Symbol == 'SPOP')], na.rm = T), lwd = 3, col = '#EAB9B8')
title(main = 'FOXA1 (n=11)', line = -1, adj = 0.1, cex.main = 1.5)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# helper functions to estimate CCF

estimate_ccf = function(purity,
                        total_copies,
                        mutant_copies,
                        t_alt_count,
                        t_depth) {
  
  
  ccfs = seq(0.001, 1, 0.001)
  expected_vaf  = function(ccf, purity, total_copies) {
    purity * ccf * mutant_copies / (2 * (1 - purity) + purity * total_copies)
  }
  
  probs = sapply(ccfs, function(c) {
    stats::dbinom(t_alt_count, t_depth, expected_vaf(c, purity, total_copies))
  })
  probs = probs / sum(probs)
  
  ccf_max = which.max(probs)
  if (identical(ccf_max, integer(0))) ccf_max = NA
  ccf_half_max = which(probs > max(probs) / 2)
  ccf_lower = max(ccf_half_max[1] - 1, 1) # closest ccf value before half-max range (within 0-1 range)
  ccf_upper = min(ccf_half_max[length(ccf_half_max)] + 1, length(ccfs)) # closest ccf value after half-max range (within 0-1 range)
  if (is.na(purity)) ccf.upper = NA 
  ccf_max = ccf_max / length(ccfs)
  ccf_lower = ccf_lower / length(ccfs)
  ccf_upper = ccf_upper / length(ccfs)
  prob95 = sum(probs[950:1000])
  prob90 = sum(probs[900:1000])
  
  list(ccf_max, ccf_lower, ccf_upper, prob95, prob90)
}

# Estimate mutant copy number, given observed VAF, purity, and local ploidy
# Based on PMID 28270531
expected_mutant_copies = function(t_var_freq,
                                  total_copies,
                                  purity) {
  
  if (is.na(total_copies)) {
    NA_real_
  } else {
    if (total_copies == 0) total_copies = 1
    mu = t_var_freq * (1 / purity) * (purity * total_copies + (1 - purity) * 2)
    alt_copies = ifelse(mu < 1, 1, abs(mu)) # mu < 1 ~ 1, mu >= 1 ~ abs(mu)
    round(alt_copies)
  }
}







