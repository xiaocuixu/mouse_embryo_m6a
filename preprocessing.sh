trim_galore --paired -j 7 --basename $sample  $read1  $read2
hisat2 -p 30 -x /home1/xuxc/01.ref/mm9Hisat2/mm9 --dta-cufflinks --no-mixed --no-discordant  -1 ${sample}_R1_val_1.fq.gz -2 ${sample}_R2_val_2.fq.gz -S $sample.sam
samtools sort -@ 30 -o $sample.bam $sample.sam
mv $sample.bam bam
samtools index bam/$sample.bam
bamCoverage -p max/2 --bam bam/$sample.bam -o bw/$sample.RPKM.bw --normalizeUsing RPKM
samtools flagstat bam/$sample.bam >flagstat/$sample.flagstat


macs2 callpeak --tempdir /home1/xuxc/tmp  -g mm -t bam/$ip -c bam/$input -n $name --nomodel --keep-dup all --outdir macs2


cufflinks -G /home1/xuxc/ann/mm9.genes.chr.gtf -p 6 -o cufflinks.out/$sample bam/$sample.bam
