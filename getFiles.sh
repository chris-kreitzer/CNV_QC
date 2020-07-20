#! bin/bash

PRAD_Files = paste0('/ifs/res/taylorlab/tcga_facets/snp-pileup/', 
                    PRAD_Files)
write.table(PRAD_Files, 
            file = 'countdata_tcga', 
            col.names = F, 
            row.names = F, 
            quote = F)


filename='countdata_tcga'

while read p; do
    cp $p ~/tmp_TCGA
    #com="cp"
    #path=$p
    #out="~/tmp_TCGA"
    #cat "${com} ${path} ${out}"
done < $filename

# copy files to local machine

rsync -avzh xjuno:~/tmp_TCGA /tmp/myrpms