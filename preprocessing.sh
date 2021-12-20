#trimming and mapping
trim_galore --paired -j 7 --basename $sample  $read1  $read2
hisat2 -p 30 -x /home1/xuxc/01.ref/mm9Hisat2/mm9 --dta-cufflinks --no-mixed --no-discordant  -1 ${sample}_R1_val_1.fq.gz -2 ${sample}_R2_val_2.fq.gz -S $sample.sam
samtools sort -@ 30 -o $sample.bam $sample.sam
mv $sample.bam bam
samtools index bam/$sample.bam
bamCoverage -p max/2 --bam bam/$sample.bam -o bw/$sample.RPKM.bw --normalizeUsing RPKM
samtools flagstat bam/$sample.bam >flagstat/$sample.flagstat

##calculate FPKM for genes
cufflinks -G /home1/xuxc/ann/mm9.genes.chr.gtf -p 6 -o cufflinks.out/$sample bam/$sample.bam


#peak calling
macs2 callpeak --tempdir /home1/xuxc/tmp  -g mm -t bam/$ip -c bam/$input -n $name --nomodel --keep-dup all --outdir macs2


#peak filtering and m6a modified gene idetification
for sample in ESC GV MII Zygote 2cell 4cell
do
cat macs2/${sample}-1_peaks.narrowPeak macs2/${sample}-2_peaks.narrowPeak macs2/${sample}-3_peaks.narrowPeak |sort -k1,1 -k2,2n |mergeBed -i - >${sample}.peaks.merged.bed
multiIntersectBed -i macs2/${sample}-1_peaks.narrowPeak macs2/${sample}-2_peaks.narrowPeak macs2/${sample}-3_peaks.narrowPeak |grep '1,2,3' |cut -f1,2,3 >${sample}.3sampleOverlapBase.bed
intersectBed -u -a ${sample}.peaks.merged.bed -b ${sample}.3sampleOverlapBase.bed >${sample}.3sampleOverlapPeak.bed
intersectBed -b ${sample}.3sampleOverlapPeak.bed -a /home1/xuxc/ann/mm9.chr.exon.sorted.bed |awk '{print $5"\tyes"}' |sort |uniq >$sample.genesWithExonPeak.id.txt
done


#TE RPKM calculation
awk 'BEGIN{print "chr\tstart\tend\tGV\tMII\tZygote\t2cell\t4cell\tESC"}' >TE.rpkm.average.txt 
multiBigwigSummary  BED-file  --outRawCounts rpkm.txt --BED ~/ann/mm9.rmsk.repName.sorted.bed -o out.npz -p 30 \
-b ../bw/GV-input-1.RPKM.bw ../bw/GV-input-2.RPKM.bw ../bw/GV-input-3.RPKM.bw \
../bw/MII-input-1.RPKM.bw ../bw/MII-input-2.RPKM.bw ../bw/MII-input-3.RPKM.bw \
../bw/Zygote-input-1.RPKM.bw ../bw/Zygote-input-2.RPKM.bw ../bw/Zygote-input-3.RPKM.bw \
../bw/2cell-input-1.RPKM.bw ../bw/2cell-input-2.RPKM.bw ../bw/2cell-input-3.RPKM.bw \
../bw/4cell-input-1.RPKM.bw ../bw/4cell-input-2.RPKM.bw \
../bw/ESC-input-1.RPKM.bw ../bw/ESC-input-2.RPKM.bw ../bw/ESC-input-3.RPKM.bw
rm out.npz
awk -v OFS="\t" 'NR>1{gv=($4+$5+$6)/3;m2=($7+$8+$9)/3;z=($10+$11+$12)/3;c2=($13+$14+$15)/3;c4=($16+$17)/2;es=($18+$19+$20)/3;print $1,$2,$3,gv,m2,z,c2,c4,es}' rpkm.txt  >>TE.rpkm.average.txt 


###plot profile for m6A
sample=("GV" "MII" "Zygote" "2cell" "4cell")
col1=("#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E")
col2=("#a2efd8" "#febe8e" "#bfbddc" "#f5a4ce" "#d0f0ac")
for i  in 0 1 2 3 4 
do
computeMatrix scale-regions -bs 50 -S /home/xuxc/projects/m6a/new_pipe/bw/merge/${sample[$i]}-IP.q20.bw \
/home/xuxc/projects/m6a/new_pipe/bw/merge/${sample[$i]}-input.q20.bw \
-R ~/ann/mm9.MTA.noGenic.full-length.bed  \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 3000 \
                              -o matrix.mat.gz -p max 
plotProfile -m matrix.mat.gz \
     -out ${sample[$i]}.MTA.profile.q20.pdf \
     --perGroup \
     --regionsLabel MTA  \
     --samplesLabel ${sample[$i]}-IP ${sample[$i]}-input \
     --colors ${col1[$i]} ${col2[$i]}  --yMin -100 --yMax 1000 --plotHeight 3.5 --plotWidth 8 \
     --endLabel stop_codon --startLabel TSS --plotFileFormat pdf --dpi 300
done
