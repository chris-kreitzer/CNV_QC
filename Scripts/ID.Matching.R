# TCGA to MSK ID matching

# neccessary libraries and data
library(parallel)
library(tictoc)
library(dplyr)
library(clipr)
library(vroom)

# TCGA TSS codes; data table
setwd('~/Documents/MSKCC/CPNA_analysis/TCGA/')
TCGA.codes = vroom::vroom('tcga_code_tables/tissueSourceSite.tsv', 
                 delim = '\t')

colnames(TCGA.codes)[1] = 'TSS.Code'
colnames(TCGA.codes)[2] = 'Source.site'
colnames(TCGA.codes)[3] = 'Study.Name'

# fetch MSK TSS codes
msk.codes = TCGA.codes[grepl('.*MSK.*', TCGA.codes$Source.site), 'TSS.Code']
msk.codes = as.character(unique(msk.codes$TSS.Code))
memorial.codes = TCGA.codes[grepl('.*Memorial.*', TCGA.codes$Source.site), 'TSS.Code']
memorial.codes = as.character(unique(memorial.codes$TSS.Code))
msk.codes = as.character(unique(c(msk.codes, memorial.codes)))

# load data
# TCGA (mc3.v0.2.8.PUBLIC.maf)
tcga.maf = vroom::vroom('RawData/mc3.v0.2.8.PUBLIC.maf', delim = '\t')
tcga.maf = as.data.frame(tcga.maf[, c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification')])
tcga.maf$eval = paste(tcga.maf$Hugo_Symbol, tcga.maf$Variant_Classification, sep = ';')
tcga.maf$TSS = substr(tcga.maf$Tumor_Sample_Barcode, start = 6, stop = 7)

# IMPACT maf (fetched from 09/2019; n ~ 30,000)
impact.maf = vroom::vroom('../../00_Data/Annotated_MAF_oncokb_hugoified.txt', delim = '\t')
impact.maf = as.data.frame(impact.maf[, c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification')])
impact.maf$eval = paste(impact.maf$Hugo_Symbol, impact.maf$Variant_Classification, sep = ';')

# annotate IMPACT samples; cancer type; 
# IMPACT meta data fetched from 07/14/2020
impact.meta = vroom::vroom('~/Downloads/new.tsv', delim = '\t')
impact.maf = merge(impact.maf, 
                   impact.meta[,c('Sample', 'CancerType')],
                   by.x = 'Tumor_Sample_Barcode',
                   by.y = 'Sample',
                   all.x = T)


# data preparation:
# TCGA samples with TSS code association to MSK/Memorial
tcga.msk = tcga.maf[which(tcga.maf$TSS %in% msk.codes),, drop = F ] # n = 701 samples 
# only include IMPACT samples >= 3 mutations
freq.mut.msk = data.frame(table(impact.maf$Tumor_Sample_Barcode))
impact.maf = impact.maf[which(impact.maf$Tumor_Sample_Barcode %in% 
                                freq.mut.msk$Var1[which(freq.mut.msk$Freq >= 3)]),, drop = F]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# one example run: n = 5,000 MSK samples
# already.tested = vroom::vroom(file = '~/Desktop/matchfirst.txt')
impact.samples = unique(impact.maf$Tumor_Sample_Barcode)
impact.samples = sample(size = 5000, impact.samples)

# subset data frame to work with those samples; IMPACT
impact.maf = impact.maf[which(impact.maf$Tumor_Sample_Barcode %in% impact.samples), ]

# loop through IMPACT
match.out.list = mclapply(unique(impact.maf$Tumor_Sample_Barcode), function(x){
  s.msk.c = a[which(impact.maf$Tumor_Sample_Barcode == x), 'eval']
  
  mclapply(unique(tcga.msk$Tumor_Sample_Barcode), function(j){
    s.tcga.c = tcga.msk[which(tcga.msk$Tumor_Sample_Barcode == j), 'eval']
      
    if(sum(s.msk.c %in% s.tcga.c) > length(s.msk.c) * 0.9){ 
      out = data.frame(msk.ID = x,
                       tcga.ID = j,
                       match = 'yes',
                       n.mut.msk = length(s.msk.c),
                       n.mut.tcga = length(s.tcga.c),
                       TSS = substr(j, start = 6, stop = 7),
                       msk.type = unique(impact.maf[impact.maf$Tumor_Sample_Barcode == x, 
                                                    'CancerType']))
      }
    
    })

})
  
match.out.data = dplyr::bind_rows(match.out.list)

# merge with TCGA IDs
match.data.complete = merge(match.out.data,
                            TCGA.codes[, c('TSS.Code', 'Study.Name')],
                            by.x = 'TSS',
                            by.y = 'TSS.Code',
                            all.x = T)




# cbio = function(x){
#   # load package
#   suppressWarnings(library(clipr))
#   out = as.character(noquote(x))
#   write_clip(out)
# }


