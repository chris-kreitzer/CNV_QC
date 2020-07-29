## PURITY Estimation: Comparison between tools:
## Pan-Cancer analysis on TCGA was made here: https://doi.org/10.1038/ncomms9971

# first look at TCGA samples:
library(TCGAbiolinks)
TP.data = TCGAbiolinks::Tumor.purity
TP.data$CPE = as.numeric(as.character(gsub(pattern = ',', '.', TP.data$CPE)))
TP.data = TP.data[!is.na(TP.data$CPE), ]

for(i in unique(TP.data$Cancer.type)){
  TP.data$purity_median[TP.data$Cancer.type == i] = median(TP.data$CPE[TP.data$Cancer.type == i])
  rm(i)
}

TP.data = TP.data %>% arrange(desc(purity_median))
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

# order cancer.types based on median
a = as.data.frame(table(TP.data$Cancer.type))
b = merge(a, TP.data[,c('Cancer.type', 'purity_median')], by.x = 'Var1', by.y = 'Cancer.type', all.x = T)
b = unique(b)
b = b[order(b$Var1, b$purity_median, decreasing = T), ]

PanPurity.plot = ggplot(TP.data, aes(x = reorder(Cancer.type, purity_median), y = CPE)) +
  geom_quasirandom(dodge.width = 0.6, cex = 1, nbins = 50, alpha = 0.5) +
  stat_summary(fun = median, 
               geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..), 
               position = position_dodge(width = 0.5), 
               width = 0.8, 
               col = 'red',
               size = 1) +
  scale_x_discrete(breaks = b$Var1,
                   labels = paste0(b$Var1, '\n(n=',b$Freq, ')')) +
  theme_std(base_size = 14) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     breaks = seq(0, 1, 0.25),
                     labels = seq(0, 1, 0.25)) +
  labs(x = '', y = 'Purity estimate')

# ggsave(filename = 'PanPurity_CPE.pdf', plot = PanPurity.plot,
#        width = 14, height = 4, units = 'in', dpi = 600,
#        device = 'pdf', type = 'cairo')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Running Facets on TCGA-PRAD samples and plot differences between methods
# Facets purity estimation on default input parameter: 
# loop through every sample in directory

library(facets)
library(facetsSuite)

files = list.files(path = '~/Desktop/tmp_TCGA/tmp_TCGA/', full.names = T)
out = data.frame()
for(i in 1:length(files)){
  try({
      print(which(files[i] == files))
      data.load = readSnpMatrix(files[i])
      data.out = run_facets(data.load, cval = 100, genome = 'hg19')
      sub = data.frame(sample = substr(files[i], start = 78, nchar(files)[i]-56),
                       purity = data.out$purity,
                       ploidy = data.out$ploidy)
      
      out = rbind(out, sub)
      rm(data.load)
      rm(data.out)
      })
}

# export table
# write.table(out, file = 'tmp/purity_out.txt', sep = '\t')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merge Facets calls with TP.data (TCGA.Biolinks)
out = read.csv('../../CPNA_analysis/TCGA/tmp/purity_out.txt', sep = '\t')
purity.merged = merge(TP.data,
                      out,
                      by.x = 'Sample.ID',
                      by.y = 'sample',
                      all.y = T)

purity.merged = purity.merged[which(purity.merged$Cancer.type == 'PRAD'), ]
purity.merged$Sample.ID = as.character(as.factor(purity.merged$Sample.ID))
purity.merged$sample = substr(purity.merged$Sample.ID, start = 1, stop = nchar(purity.merged$Sample.ID) - 1)

# convert columns into numerics
for(i in c(3, 5, 6)){
  purity.merged[, i] = as.numeric(as.character(gsub(pattern = ',', '.', purity.merged[, i])))
}

purity.merged = purity.merged[!is.na(purity.merged$ESTIMATE), ]

# make plot to compare methods:
library(ggplot2)
plot.out = list()
rm(i)
for(i in c(3, 5, 6)){
  plot = ggplot(purity.merged, aes(x = purity, y = purity.merged[, i])) +
    geom_jitter(col = '#4576bb') +
    theme_std(base_size = 14) +
    scale_y_continuous(expand = c(0,0),
                       limits = c(0, 1),
                       breaks = c(0, 0.5, 1),
                       labels = c(0, 0.5, 1)) +
    scale_x_continuous(expand = c(0,0),
                       limits = c(0, 1),
                       breaks = c(0, 0.5, 1),
                       labels = c(0, 0.5, 1)) +
    theme(aspect.ratio = 1) +
    geom_abline(intercept = 0, slope = 1, col = 'grey85', size = 1.2) +
    labs(x = 'FACETS',
         y = colnames(purity.merged)[i],
         title = paste0('paired t-test; p-value = ',
                        prettyNum(t.test(purity.merged[, i], purity.merged$purity, paired = T)['p.value'])))
  print(plot)
  
  plot.out[[i]] = plot
  
}

# arrange PanPurity Estimates and individual plots
library(gridExtra)
library(egg) 




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# some in-depth analysis:
# Samples where FACETS reports on purity == NA
cbio(purity.merged$sample[which(is.na(purity.merged$purity))]) # n = 53 samples
length(purity.merged$Sample.ID[which(is.na(purity.merged$purity))])

# look at the purity density reported from FACETS
# note that FACETS shows an peak at ~ 0.3
dens = density(purity.merged$purity, na.rm = T, bw = 0.025)
par(mar = c(5,4,4,3))
plot(1,
     type = 'n',
     xlim = c(0, 1.3),
     ylim = c(0, 4.2),
     axes = F,
     xlab = 'Purity estimates',
     ylab = '',
     pch = 16,
     main = '',
     xaxt = 'n',
     xaxs = 'i',
     yaxs = 'i')

axis(side = 2, at = c(0, 4), labels = FALSE, tck = 0.0, lwd = 2)
mtext(side = 2, text = "Density (a.u.)", line = 1)

axis(side = 1,
     at = c(0, 0.5, 1),
     labels = c('0', '0.5', '1'),
     las = 1,
     lwd = 2)

# add Facets density
lines(dens, lwd = 3, col = 'grey55')
# ESTIMATE
lines(density(purity.merged$ESTIMATE, bw = 0.08), lwd = 3, col = '#83b2a2')
# LUMP
lines(density(purity.merged$LUMP, bw = 0.08, na.rm = T), lwd = 3, col = '#9FBEE6')

# how are density estimates look like on samples where FACETS == NA
# plot layout from above (run those lines of code)
lines(density(purity.merged$ESTIMATE[is.na(purity.merged$purity)], bw = 0.08, na.rm = T), lwd = 3, col = '#83b2a2')
lines(density(purity.merged$LUMP[is.na(purity.merged$purity)], bw = 0.08, na.rm = T), lwd = 3, col = '#9FBEE6')
title(main = 'ESTIMATE/LUMP purity Estimation on\n FACETS == NA samples',
      cex.main = 1.5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# looking for alternatives: TPES 
# TPES uses SNV data to estimate tumor purity: https://doi.org/10.1093/bioinformatics/btz406

# select all samples which have FACETS NA
missing.purity = purity.merged$Sample.ID[which(is.na(purity.merged$purity))]

# prepare data for TPES calculation
# IDs
# Segments
# MAF
# Ploidy

# Segments
Segmentation.data = vroom::vroom(file = '~/Documents/MSKCC/CPNA_analysis/TCGA/RawData/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg', 
                                 delim = '\t')
missing_segmentation = data.frame()
for(i in 1:length(missing.purity)){
  selected.samples = Segmentation.data[grep(missing.purity[i], Segmentation.data$Sample), ]
  missing_segmentation = rbind(missing_segmentation, selected.samples)
}
rm(Segmentation.data)
missing_segmentation$Sample = substr(missing_segmentation$Sample, start = 1, stop = 16)
colnames(missing_segmentation) = c('ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean')

# maf data
maf <- vroom::vroom('../../CPNA_analysis/TCGA/RawData/mc3.v0.2.8.PUBLIC.maf', delim = '\t')
maf$Sample = substr(maf$Tumor_Sample_Barcode, start = 1, stop = 16)
maf.TPES = maf[,c('Sample', 'Chromosome', 'Start_Position', 'End_Position', 't_ref_count', 't_alt_count')]
colnames(maf.TPES) = c('sample', 'chr', 'start', 'end', 'ref.count', 'alt.count')
maf.TPES$chr = paste0('chr', maf.TPES$chr)
maf.TPES = maf.TPES[,c(2,3,4,5,6,1)]
rm(maf)

# Ploidy data (from facets)
Ploidy = purity.merged[,c('Sample.ID', 'ploidy')]
colnames(Ploidy) = c('sample', 'ploidy')

# run function
purity.SNV = data.frame()  
for(i in 1:length(missing.purity)){
  purity.calculated = TPES::TPES_purity(ID = missing.purity[i],
                        SEGfile = missing_segmentation,
                        SNVsReadCountsFile = maf.TPES,
                        ploidy = Ploidy,
                        minSNVs = 5)
  purity.SNV = rbind(purity.SNV, purity.calculated)
}

# make density plot of TPES calculated estimates
# make plot
plot(1,
     type = 'n',
     xlim = c(0, 1.3),
     ylim = c(0, 4.2),
     axes = F,
     xlab = 'TPES purity estimates',
     ylab = '',
     pch = 16,
     main = '',
     xaxt = 'n',
     xaxs = 'i',
     yaxs = 'i')

axis(side = 2, at = c(0, 4), labels = FALSE, tck = 0.0, lwd = 2)
mtext(side = 2, text = "Density (a.u.)", line = 1)

axis(side = 1,
     at = c(0, 0.5, 1),
     labels = c('0', '0.5', '1'),
     las = 1,
     lwd = 2)

lines(density(purity.SNV$purity.max, na.rm = T, bw = 0.10), lwd = 3, col = 'grey55')





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# for(proj in 1:length(project.ids)){
#   query = GDCquery(project = project.ids[proj],
#                    data.category = "Copy Number Variation",
#                    data.type = 'Masked Copy Number Segment',
#                    platform = 'Affymetrix SNP 6.0')
#   GDCdownload(query)
#   data = GDCprepare(query)
#   Segmentation.data = rbind(Segmentation.data, data)
# }