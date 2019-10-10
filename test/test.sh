python3 ../bin/CrossMap.py bed  ../data/UCSC_chain/hg19ToHg38.over.chain.gz 0_hg19.bed12 0_hg38.bed

python3 ../bin/CrossMap.py vcf ../data/UCSC_chain/hg19ToHg38.over.chain.gz 1_hg19.vcf ~/Documents/database/ucsc_hg38.fa 1_out.hg38.vcf

python3 ../bin/CrossMap.py bam ../data/UCSC_chain/hg19ToHg38.over.chain.gz 2_hg19.bam 2_hg38.bam

python3 ../bin/CrossMap.py bigwig  ../data/UCSC_chain/hg19ToHg38.over.chain.gz 3_hg19.bw 3_hg38.bw

python3 ../bin/CrossMap.py wig  ../data/UCSC_chain/hg19ToHg38.over.chain.gz 4_hg19.bgr 4_hg38.bgr

python3 ../bin/CrossMap.py wig  ../data/UCSC_chain/hg19ToHg38.over.chain.gz 5_hg19.fixedStep.wig 5_hg38.fixedStep.wig

python3 ../bin/CrossMap.py wig  ../data/UCSC_chain/hg19ToHg38.over.chain.gz 6_hg19.variableStep.wig 6_hg38.variableStep.wig 

python3 ../bin/CrossMap.py gff  ../data/UCSC_chain/hg19ToHg38.over.chain.gz 7_hg19.gtf 7_hg38.gtf

python3 ../bin/CrossMap.py maf ../data/UCSC_chain/hg19ToHg38.over.chain.gz  8_hg19.maf  ~/Documents/database/ucsc_hg38.fa GRCh38 8_hg38.maf
