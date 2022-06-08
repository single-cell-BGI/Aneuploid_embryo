#-----------------------------------------------------------------------------------------#
#                          Step0 set options & mkdir directory                            #
#-----------------------------------------------------------------------------------------#
Sample=sample_list.txt
Rawdata_path=rawdata
Out_path=out_dir
Batch=sample_batch
mkdir -p $Out_path/0_shell
mkdir -p $Out_path/1_Cleandata
mkdir -p $Out_path/2_Alignment
mkdir -p $Out_path/3_UsableBam
mkdir -p $Out_path/4_GeneCount

for i in `cat $Sample`
do echo "
#--------------------------------------------------#
#               Step1 RNA-Clean data               #
#--------------------------------------------------#
echo start at time \`date +%F'  '%H:%M\`
mkdir $Out_path/1_Cleandata/$i && /hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/bin/SOAPnuke filter -l 5 -q 0.5 -n 0.05 -Q 2 -1 $Out_path/rawdata_trimME/${i}_trimME_1.fq.gz -2 $Out_path/rawdata_trimME/${i}_trimME_2.fq.gz -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG -o $Out_path/1_Cleandata/$i -C ${i}_clean_1.fq.gz -D ${i}_clean_2.fq.gz
echo Step1 RNA-Clean data This-Work-is-Completed!
echo finish at time \`date +%F'  '%H:%M\`
#--------------------------------------------------#
#               Step2 RNA-Rm rRNA                  #
#--------------------------------------------------#
echo start at time \`date +%F'  '%H:%M\`
/hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/bin/soap_mm_gz -a $Out_path/1_Cleandata/$i/${i}_clean_1.fq.gz -b $Out_path/1_Cleandata/$i/${i}_clean_2.fq.gz -D /hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/database/Human/hg19/Rm_rRNA/index/Human_rRNA_NCBI.fa.index -m 0 -x 1000 -s 28 -l 31 -v 5 -r 1 -p 3 -o $Out_path/1_Cleandata/${i}_clean_rRNA_PEsoap.fq.gz -2 $Out_path/1_Cleandata/${i}_clean_rRNA_SEsoap.fq.gz
/share/app/perl-5.22.0/bin/perl /hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/myScript/ATAC-seq/ATAC-seq_pipeline/Upload/RNAseq/rRNAFilter.pl -fq $Out_path/1_Cleandata/$i/${i}_clean_1.fq.gz,$Out_path/1_Cleandata/$i/${i}_clean_2.fq.gz -soap $Out_path/1_Cleandata/${i}_clean_rRNA_PEsoap.fq.gz,$Out_path/1_Cleandata/${i}_clean_rRNA_SEsoap.fq.gz -output $Out_path/1_Cleandata/${i}_clean_rRNA
echo Step2 RNA-Trim data This-Work-is-Completed!
echo finish at time \`date +%F'  '%H:%M\`
#--------------------------------------------------#
#               Step3 RNA-Alignment                #
#--------------------------------------------------#
echo start at time \`date +%F'  '%H:%M\`
/hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/biosoftware/hisat2-2.1.0/hisat2 -p 8 --phred33 --sensitive --no-discordant --no-mixed -I 1 -X 1000 -x /ldfssz1/ST_MCHRI/STEMCELL/project/Abnormal_Embryo_P18Z10200N0350/Geneome/hisat2_chrALL_ERCC_index -1 $Out_path/1_Cleandata/${i}_clean_rRNA_1.fq.gz -2 $Out_path/1_Cleandata/${i}_clean_rRNA_2.fq.gz 2>$Out_path/2_Alignment/sta_Hisat2Genome_${i}.txt | samtools view -b -S - | samtools sort -@ 5 -o $Out_path/2_Alignment/${i}_sort.bam -
samtools index $Out_path/2_Alignment/${i}_sort.bam
echo Step3 RNA-Alignment This-Work-is-Completed!
echo finish at time \`date +%F'  '%H:%M\`
#--------------------------------------------------#
#               Step4 RNA-Bam process              #
#--------------------------------------------------#
echo start at time \`date +%F'  '%H:%M\`
/hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/biosoftware/Anaconda2/bin/samtools view -h -F 12 -q 30 -@ 5 -b $Out_path/2_Alignment/${i}_sort.bam > $Out_path/2_Alignment/${i}.sort.q30.bam && /hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/biosoftware/Anaconda2/bin/samtools index -@ 5 $Out_path/2_Alignment/${i}.sort.q30.bam
/hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/biosoftware/Anaconda2/bin/java -Xmx2g -jar /hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/biosoftware/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true INPUT=$Out_path/2_Alignment/${i}.sort.q30.bam OUTPUT=$Out_path/3_UsableBam/${i}.sort.q30.nodup.bam M=$Out_path/2_Alignment/sta_rmdup_${i}.txt
samtools index $Out_path/3_UsableBam/${i}.sort.q30.nodup.bam
samtools view $Out_path/3_UsableBam/${i}.sort.q30.nodup.bam | grep \"ERCC\" | awk -F \"\t\" '{print \$3}' > $Out_path/3_UsableBam/ERCC_${i}.txt
echo Step4 RNA-Bam process This-Work-is-Completed!
echo finish at time \`date +%F'  '%H:%M\`
#--------------------------------------------------#
#              Step5 RNA-Genecount                 #
#--------------------------------------------------#
echo start at time \`date +%F'  '%H:%M\`
/hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/biosoftware/Anaconda2/bin/geneBody_coverage.py -r /hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/database/Human/hg19/hg19_RefSeq.bed -i $Out_path/3_UsableBam/${i}.sort.q30.nodup.bam -o $Out_path/4_GenebodyPlot/gbPlot_${i}
/hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/biosoftware/Anaconda2/bin/Rscript /ldfssz1/ST_MCHRI/STEMCELL/project/Abnormal_Embryo_P18Z10200N0350/0_Pipeline/count_GeneReadNumber_ERCC.R /hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/database/Human/hg19/Gencode_transcriptome/gencode.v19.chr_patch_hapl_scaff.annotation.gtf $Out_path/3_UsableBam/${i}.sort.q30.nodup.bam /ldfssz1/ST_MCHRI/STEMCELL/project/Abnormal_Embryo_P18Z10200N0350/0_Pipeline/ERCC_gene_length.txt $Out_path/4_GeneCount/genecount_${i}.txt $Out_path/3_UsableBam/ERCC_${i}.txt
echo Step5 RNA-Genecount process This-Work-is-Completed!

" > $Out_path/0_shell/RNA_qsub_${i}.sh;done
sh $Out_path/0_shell/RNA_qsub_${i}.sh