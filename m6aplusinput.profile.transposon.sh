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
