python2.7 ../bin/CrossMap.py bed  ../data/UCSC_chain/hg19ToHg38.over.chain.gz 0_hg19.bed 0_hg38.bed

python2.7 ../bin/CrossMap.py vcf ../data/UCSC_chain/hg19ToHg38.over.chain.gz 1_hg19.vcf /database/ucsc_hg38.fa 1_out.hg38.vcf

python2.7 ../bin/CrossMap.py bam ../data/UCSC_chain/hg19ToHg38.over.chain.gz 2_hg19.bam 2_hg38.bam


python2.7 ../bin/CrossMap.py bigwig  ../data/UCSC_chain/hg19ToHg38.over.chain.gz 3_hg19.bw 3_hg38.bw

python2.7 ../bin/CrossMap.py wig  ../data/UCSC_chain/hg19ToHg38.over.chain.gz 4_hg19.bgr 4_hg38.bgr

python2.7 ../bin/CrossMap.py wig  ../data/UCSC_chain/hg19ToHg38.over.chain.gz 5_hg19.fixedStep.wig 5_hg38.fixedStep.wig

python2.7 ../bin/CrossMap.py wig  ../data/UCSC_chain/hg19ToHg38.over.chain.gz 6_hg19.variableStep.wig 6_hg38.variableStep.wig 

python2.7 ../bin/CrossMap.py gff  ../data/UCSC_chain/hg19ToHg38.over.chain.gz 7_hg19.gtf 7_hg38.gtf
