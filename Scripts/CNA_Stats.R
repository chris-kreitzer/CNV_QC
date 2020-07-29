# Investigate TCGA-Affy data and run Facets on those samples;
# I start this example run with an example folder; TCGA-PRAD

setwd('~/Documents/MSKCC/CPNA_analysis/TCGA/')
library(parallel)
library(facets)
library(facetsSuite)

# n = 60 TCGA samples; Affymetrix and compare with FACETS;
# first look at number of segments; then look at loss vs gains than at CNA and than gene-level
Segmentation.data = vroom::vroom('RawData/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg',
                                 delim = '\t')
example.tcga = list.files('~/Documents/MSKCC/CPNA_analysis/TCGA/PTEN_counts/', 
                          full.names = T)

# fetch samples from Whitelist; segs
selected.affy = data.frame()
for(i in 2:length(example.tcga)){
  local({
    string = substr(example.tcga[i], start = 99, stop = 114)
    sub.local = Segmentation.data[grep(string, Segmentation.data$Sample), ]
    selected.affy <<- rbind(selected.affy, sub.local)
  })
}


# look at segment level; Affy data; n = 52 samples
summary.out = data.frame()
for(i in unique(selected.affy$Sample)){
  local({
    data.selected = selected.affy[which(selected.affy$Sample == i), ]
    summary.affy = data.frame(Sample = i,
                              total = length(data.selected$Segment_Mean[which(data.selected$Segment_Mean <= -0.2 | data.selected$Segment_Mean >= 0.2)]),
                              loss = length(data.selected$Segment_Mean[which(data.selected$Segment_Mean <= -0.2)]),
                              gain = length(data.selected$Segment_Mean[which(data.selected$Segment_Mean >= 0.2)]),
                              group = 'TCGA')
    
    summary.out <<- rbind(summary.out, summary.affy)
  })
}


# now let's run FACETS on those samples
facets.all.out = data.frame()
for(i in 1:length(example.tcga)){
  
  skip_to_next = FALSE
  
  tryCatch({
    
    local({
      countmatrix = facets::readSnpMatrix(filename = example.tcga[i])
      facets.out = facetsSuite::run_facets(read_counts = countmatrix,
                                           cval = 250, 
                                           genome = 'hg19')
      facets.segs = facets.out$segs
      facets.segs$seg.mean = log2(facets.segs$tcn.em / facets.out$ploidy)
      facets.summary = data.frame(Sample = substr(example.tcga[i], start = 68, stop = 83),
                                  total = length(facets.segs$seg.mean[which(facets.segs$seg.mean <= -0.2 | facets.segs$seg.mean >= 0.2)]),
                                  loss = length(facets.segs$seg.mean[which(facets.segs$seg.mean <= -0.2)]),
                                  gain = length(facets.segs$seg.mean[which(facets.segs$seg.mean >= 0.2)]),
                                  group = 'Facets')
      
      facets.all.out <<- rbind(facets.all.out, facets.summary)
    })
    
  }, error = function(e) {skip_to_next <<- TRUE})
  
  if(skip_to_next){ next }
  
}

x = facets::readSnpMatrix('~/Documents/MSKCC/CPNA_analysis/TCGA/PTEN_counts/TCGA-2A-AAYU-10A-01D-A41N-08_TCGA-2A-AAYU-01A-11D-A41K-08_ihvrFYsz_countsMerged___normal_tumor.dat.gz')
y = facetsSuite::run_facets(x)
z = facetsSuite::arm_level_changes(segs = y$segs, ploidy = y$ploidy)




