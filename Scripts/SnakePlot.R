# Snake plot; number of segments called Affy VS Facets:

# modify data for plotting
facets.all.out$Sample = substr(facets.all.out$Sample, 
                               start = 3, 
                               stop = nchar(facets.all.out$Sample)-2)
facets.all.out$group = rep('Facets', nrow(facets.all.out))
summary.out$Sample = substr(summary.out$Sample, start = 1, stop = 12)
colnames(summary.out)[2] = 'total'
summary.out$group = rep('TCGA', nrow(summary.out))

# combine data and sort according to group
data.aggregate = rbind(facets.all.out, summary.out)

group.sorting = c('TCGA', 'Facets')
data.aggregate$order = match(data.aggregate$group, group.sorting)
data.aggregate = data.aggregate[order(data.aggregate$order, data.aggregate$total), ]

# snake plot
x = c(0.0, length(unique(data.aggregate$group)) + 0.05)
ylab = '# Segments > | threshold |'
ylim = c(0, max(data.aggregate$total) + 7)

# plot layout
plot(1,
     type = 'n',
     xlim = x,
     ylim = ylim,
     axes = F,
     xlab = '',
     ylab = ylab,
     pch = 16,
     main = '',
     xaxt = 'n',
     xaxs = 'i',
     yaxs = 'i')

# axis left
axis(side = 2,
     at = seq(0, 320, 80),
     las = 2)

axis(side = 1,
     at = c(0, 1, 2),
     labels = c('', '', ''),
     las = 3)

# plot white/grey rects:
for(i in seq(1, length(unique(data.aggregate$group)), by = 2)){
  
  xleft = i
  xright = i + 1
  ybottom = -2
  ytop = max(ylim)
  
  rect(xleft, xright, ybottom, ytop, col = 'grey90', border = F)
  
}

# help lines

for(lines in seq(0, 320, 80)){
  abline(h = lines,
         lty = 'dashed',
         col = 'grey')
}


# add count data
plot.order <- unique(data.aggregate$group)
k = 1

for (i in plot.order){
  cancerTypeTable = data.aggregate[data.aggregate$group %in% i,, drop = FALSE]
  
  x.start = k - 0.75
  x.end = k - 0.25
  y.start = median(as.numeric(cancerTypeTable$total))
  y.end = y.start
  segments(x.start, 
           y.start,
           x.end,
           y.end,
           lwd = 4,
           col = 'orange')
  
  # plot points
  start = k - 0.90
  end = k - 0.1  
  
  length.out <- nrow(cancerTypeTable)
  
  xseq = seq(start, end, length.out = length.out)
  yseq = as.numeric(cancerTypeTable$total)
  
  if (i == 'TCGA') points(xseq, yseq, pch = 16, col = "black")  
  if (i == 'Facets') points(xseq, yseq, pch = 16, col = "brown2")
  
  k = k + 1
}


# # statistical test
# x = merge(facets.all.out, 
#           summary.out, 
#           by = 'Sample', 
#           all.y = T)
# 
# x = x[!is.na(x$total.x), ]
# # wilcoxon paired test
# wilcox.test(x$total.x, x$total.y, paired = T)
# write.table(facets.all.out, file = 'tmp/facets.out.txt', sep = '\t', row.names = F, quote = F)
# write.table(summary.out, file = 'tmp/Affy.out.txt', sep = '\t', row.names = F, quote = F)


