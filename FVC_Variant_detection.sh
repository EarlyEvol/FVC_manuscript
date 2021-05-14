#!/bin/bash 
# 
#usage:  Variant_detection.sh  [-q map quality] [-r minimum length of reads] [-i min insert size] [-x reference seq] [-t trial] SAMPLE 
# 
#Trims reads with trimmomatic 
#Aligns with Bowtie2 
#SAM=>BAM,sort and Rmdup with samtools 
#indel realignment with GATK 
#Variant call with SNVer 



#Defaults 
REF=FVC 
TRIAL=$(date +"%Y-%m-%d-%H:%M:%S") 
MIN_READ_LENGTH=30 
MIN_INSERT_LENGTH=0 
MQ=10 

USAGE='usage: '`basename $0`' [-q map quality] [-r minimum length of reads] [-i min insert size] [-x reference seq] [-t trial] SAMPLE' 

#Setting Variables 
while [[ ${1:0:1} = '-' ]] ; do 
N=1 
L=${#1} 
while [[ $N -lt $L ]] ; do 
  case ${1:$N:1} in 
      
     'q') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
            echo $USAGE 'r' 
            exit 1 
          fi 
          MQ=${2} 
          shift ;; 

     'r') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
            echo $USAGE 'r' 
            exit 1 
          fi 
          MIN_READ_LENGTH=${2} 
          shift ;; 
    
     'i') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
            echo $USAGE 'i' 
            exit 1 
          fi 
          MIN_INSERT_LENGTH=${2} 
          shift ;; 
     
     'x') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
            echo $USAGE 
            exit 1 
          fi 
          REF=${2} 
          shift ;; 
  
     't') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
            echo $USAGE 't' 
            exit 1 
          fi 
          TRIAL=${2} 
          shift ;; 
     *) echo $USAGE 't2' 
        exit 1 ;; 
  esac 
  N=$(($N+1)) 
done 
shift 
done 
if [[ ! -n ${1} ]] ; then 
echo $USAGE 
exit 1 
fi 
SAMPLE=$1 

#reporting and exporting variables 
echo '----Running '$0'----' 
echo 'Samples:'$SAMPLE 
echo 'REF='$REF 
echo 'TRIAL'=$TRIAL 
echo 'minimum read length'=$MIN_READ_LENGTH 
echo 'minimum insert'=$MIN_INSERT_LENGTH 
echo 'Map Quality'=$MQ 
export SAMPLE 
export REF 
export TRIAL 
export MIN_READ_LENGTH 
export MIN_INSERT_LENGTH 
export MQ 

#House keeping 
wd=~/STORAGE/WORKING_DIR/
mkdir $wd
cd $wd
mkdir $SAMPLE 
mkdir $SAMPLE/$TRIAL 
mkdir $SAMPLE/$TRIAL/ALIGN 
mkdir $SAMPLE/$TRIAL/REALIGN 
mkdir $SAMPLE/$TRIAL/DNSMPL 
mkdir $SAMPLE/$TRIAL/VCF 

#Reference indexing .fai .dict and .bt2 
cd $wd/../REFS 
bowtie2-build -q $REF.fasta $REF 
samtools faidx $REF.fasta  #some flowting point exeption thing....fai not empty, proceed.  Works with GATK realigner 
rm $REF.dict 
java -jar ~/bin/jars/CreateSequenceDictionary.jar REFERENCE=$REF.fasta OUTPUT=$REF.dict 

#Trim reads Trimmomatic 
echo "Trimmomatic fastq QC" 
cd $wd 
java -classpath ~/bin/jars/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticPE -phred33 \
$wd/$SAMPLE.R1.fastq $wd/$SAMPLE.R2.fastq \
$SAMPLE/$SAMPLE.trim.R1.fastq $SAMPLE/$SAMPLE.trim_unpaired.R1.fastq $SAMPLE/$SAMPLE.trim.R2.fastq $SAMPLE/$SAMPLE.trim_unpaired.R2.fastq \
LEADING:13 TRAILING:13 SLIDINGWINDOW:4:20 MINLEN:$MIN_READ_LENGTH 

#bowtie2 alignment of '$SAMPLE' to '$REF' 
echo 
echo 'Alignment' $SAMPLE 'to' $REF 'with Bowtie2' 
cd $wd/$SAMPLE 
bowtie2 -x $wd/../REFS/$REF -1 $SAMPLE.trim.R1.fastq -2 $SAMPLE.trim.R2.fastq \
-U $SAMPLE.trim_unpaired.R1.fastq,$SAMPLE.trim_unpaired.R2.fastq \
-S $TRIAL/ALIGN/$SAMPLE.aligned.sam -p 16 --quiet --rg-id 1234567 --rg SM:$TRIAL --un-conc $TRIAL/$SAMPLE.UnAigned.fastq \
--no-discordant -I $MIN_INSERT_LENGTH -X 2000 --no-unal --local 

#sam-bam conversion, sort bam, remove PCR duplicates, index SAMPLE.aln.srt.rmd.bam 
echo 
echo "Samtools stuff and whatnot" 
cd $wd/$SAMPLE/$TRIAL/ALIGN 
samtools view -Suhb $SAMPLE.aligned.sam | samtools sort - $SAMPLE.aligned.sorted 
samtools rmdup $SAMPLE.aligned.sorted.bam $SAMPLE.aln.srt.rmd.bam 
samtools index $SAMPLE.aln.srt.rmd.bam 
#rm $SAMPLE.aligned.sam 

#indel realignment with GATK 
cd $wd/$SAMPLE/$TRIAL 
echo 
echo -----GATK_RealignerTargetCreator----- 
echo 
java -Xmx2g -jar ~/bin/jars/GenomeAnalysisTK.jar \
-I ALIGN/$SAMPLE.aln.srt.rmd.bam \
-R $wd/../REFS/$REF.fasta \
-T RealignerTargetCreator \
-o REALIGN/$SAMPLE.intervals 
echo 
echo -----GATK_IndelRealigner------ 
echo 
java -Xmx2g -jar ~/bin/jars/GenomeAnalysisTK.jar \
-I ALIGN/$SAMPLE.aln.srt.rmd.bam \
-R $wd/../REFS/$REF.fasta \
-T IndelRealigner \
-targetIntervals REALIGN/$SAMPLE.intervals \
-dcov 10000 \
-o REALIGN/$SAMPLE.realn.srt.rmd.bam 


#SNVer Pile-up 
echo 
echo "SNVer pileup" 
cd $wd/$SAMPLE/$TRIAL 
java -Xmx4g -jar ~/apps/SNVerPool.jar -i REALIGN/ -r $wd/../REFS/$REF.fasta -n 200 -p .01 -mq $MQ -o VCF/$SAMPLE.$REF.$TRIAL 

cp VCF/*.filter.vcf $wd/VCF/ 

echo __________________________________________________ 
echo ____________WE ARE FINISHED, YOU AND I____________ 

